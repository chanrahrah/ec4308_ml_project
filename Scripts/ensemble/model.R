# Ensembles of the random forest, elastic net and boosting forecasts:
#   ens_majority_vote  each member votes with its own tuned cutoff; stocks if most say stocks
#   ens_equal          mean member probability, 0.50 cutoff
#   ens_gr             Granger-Ramanathan weights fitted once on validation, 0.50 cutoff
#   ens_gr_yearly      robustness check: GR refitted each January on all outcomes known by then (test only)
ensemble_config <- function() modifyList(base_config(), list(
  model = "ensemble", members = c(random_forest = "rf", elastic_net = "elastic_net", boosting = "boosting"),
  probability_cutoff = 0.50, gr_intercept = TRUE, gr_min_rows = 24L))

MEMBER_COLS <- c("model", "h", "origin", "probability", "cutoff", "signal", "cutoff_050", "signal_050",
                 "target", "outcome_end", "classification_eligible")

# Stocks (1) when more than half the members say stocks; ties go to bills (0).
majority_vote <- function(signals) {
  signals <- as.matrix(signals)
  ensure(ncol(signals) > 0 && !anyNA(signals) && all(signals %in% c(0, 1)), "Member signals must be 0/1")
  as.integer(rowSums(signals) > ncol(signals) / 2)
}

# Check one member's forecasts for one period: complete calendar, one cutoff, consistent signals.
check_member <- function(d, member, period) {
  missing <- setdiff(MEMBER_COLS, names(d))
  ensure(!length(missing), paste(member, period, "is missing columns", paste(missing, collapse = ", ")))
  d <- as.data.frame(d)[, MEMBER_COLS]
  d$outcome_end <- as.character(d$outcome_end)
  d$outcome_end[!is.na(d$outcome_end) & !nzchar(d$outcome_end)] <- NA
  d$origin_id <- month_number(d$origin)
  ensure(!anyNA(d$origin_id) && !anyDuplicated(d[, c("h", "origin_id")]), paste(member, period, "has bad or duplicate origins"))
  ensure(!anyNA(d$probability) && all(d$probability >= 0 & d$probability <= 1), paste(member, period, "has invalid probabilities"))
  for (h in unique(d$h)) {
    x <- d[d$h == h, ]
    ensure(length(unique(x$cutoff)) == 1L && all(x$cutoff_050 == 0.5), paste(member, period, "h", h, "cutoff is not fixed"))
    ensure(all(diff(sort(x$origin_id)) == 1L), paste(member, period, "h", h, "calendar has gaps"))
  }
  ensure(identical(as.integer(d$signal), as.integer(d$probability >= d$cutoff)) &&
           identical(as.integer(d$signal_050), as.integer(d$probability >= 0.5)), paste(member, period, "signals do not match cutoffs"))
  d$member <- member
  d
}

# Merge members by (h, origin); they must agree on targets, outcome dates and eligibility.
build_vote <- function(members, period) {
  keyed <- lapply(names(members), function(m) check_member(members[[m]], m, period))
  names(keyed) <- names(members)
  base <- keyed[[1]][order(keyed[[1]]$h, keyed[[1]]$origin_id), c("h", "origin_id", "target", "outcome_end", "classification_eligible")]
  key <- paste(base$h, base$origin_id)
  S <- S050 <- P <- matrix(NA_real_, nrow(base), length(keyed), dimnames = list(NULL, names(keyed)))
  same <- function(a, b) identical(is.na(a), is.na(b)) && all(a[!is.na(a)] == b[!is.na(b)])
  for (m in names(keyed)) {
    x <- keyed[[m]]; i <- match(key, paste(x$h, x$origin_id))
    ensure(!anyNA(i) && nrow(x) == nrow(base), paste(period, m, "does not cover the same origins"))
    ensure(same(x$target[i], base$target) && same(x$outcome_end[i], base$outcome_end) &&
             same(as.logical(x$classification_eligible[i]), as.logical(base$classification_eligible)),
           paste(period, m, "disagrees on targets"))
    S[, m] <- x$signal[i]; S050[, m] <- x$signal_050[i]; P[, m] <- x$probability[i]
  }
  list(base = base, S = S, S050 = S050, P = P,
       cutoffs = do.call(rbind, lapply(keyed, function(x) unique(x[, c("member", "h", "cutoff")]))),
       vote = majority_vote(S), vote_050 = majority_vote(S050))
}

# Granger-Ramanathan: min mean((y - a - P w)^2) with w >= 0 and sum(w) = 1.
# Least squares with sum(w) = 1 on one subset of members; NULL if rank-deficient.
gr_fit_subset <- function(P, y, S, intercept) {
  k <- length(S); base <- P[, S[k]]
  X <- NULL
  if (intercept) X <- matrix(1, nrow(P), 1L)
  if (k > 1L) X <- cbind(X, P[, S[-k], drop = FALSE] - base)
  coef <- numeric(0)
  if (!is.null(X)) {
    qx <- qr(X)
    if (qx$rank < ncol(X)) return(NULL)
    coef <- qr.coef(qx, y - base)
  }
  a <- if (intercept) coef[[1]] else 0
  wz <- if (k > 1L) coef[(intercept + 1L):length(coef)] else numeric(0)
  w <- c(wz, 1 - sum(wz))
  list(S = S, w = unname(w), a = unname(a), sse = sum((y - a - as.vector(P[, S, drop = FALSE] %*% w))^2))
}

# Exact solution: the best nonnegative fit over all member subsets.
gr_weights <- function(P, y, intercept = TRUE) {
  P <- as.matrix(P); K <- ncol(P)
  ensure(!is.null(colnames(P)) && nrow(P) == length(y) && !anyNA(P) && !anyNA(y) && K <= 12L, "Invalid GR inputs")
  best <- NULL
  for (mask in seq_len(2^K - 1)) {
    fit <- gr_fit_subset(P, y, which(bitwAnd(mask, 2L^(0:(K - 1L))) > 0), intercept)
    if (is.null(fit) || any(fit$w < -1e-10)) next
    if (is.null(best) || fit$sse < best$sse - 1e-12) best <- fit
  }
  ensure(!is.null(best), "No feasible GR solution")
  w <- setNames(numeric(K), colnames(P)); w[best$S] <- pmax(best$w, 0); w <- w / sum(w)
  list(weights = w, intercept = best$a, sse = best$sse, n = length(y), status = "ok")
}

# If the solver fails, use equal weights with zero intercept and flag it.
gr_weights_safe <- function(P, y, intercept = TRUE) {
  tryCatch(gr_weights(P, y, intercept), error = function(e)
    list(weights = setNames(rep(1 / ncol(P), ncol(P)), colnames(P)), intercept = 0, n = length(y),
         status = paste("fallback_equal_weights:", conditionMessage(e))))
}

# Weighted score, clipped to [0, 1].
gr_predict <- function(P, fit) {
  P <- as.matrix(P)
  ensure(setequal(colnames(P), names(fit$weights)), "GR weights and members do not match")
  raw <- fit$intercept + as.vector(P[, names(fit$weights), drop = FALSE] %*% fit$weights)
  list(raw = raw, probability = pmin(1, pmax(0, raw)), clipped = raw < 0 | raw > 1)
}

equal_average <- function(P) rowMeans(as.matrix(P))

# GR refitted each January of year Y on every row whose outcome was known by December of Y-1.
# P, y, origin_id, outcome_id, eligible: pooled history for one horizon; pred: rows to forecast.
gr_yearly <- function(P, y, origin_id, outcome_id, eligible, pred, intercept = TRUE, min_rows = 24L) {
  P <- as.matrix(P); years <- origin_id[pred] %/% 12L
  raw <- prob <- rep(NA_real_, length(pred)); clip <- rep(NA, length(pred)); fits <- list()
  for (Y in sort(unique(years))) {
    refit <- Y * 12L - 1L
    rows <- which(eligible & !is.na(y) & !is.na(outcome_id) & outcome_id <= refit)
    fit <- if (length(rows) >= min_rows) gr_weights_safe(P[rows, , drop = FALSE], y[rows], intercept) else
      list(weights = setNames(rep(1 / ncol(P), ncol(P)), colnames(P)), intercept = 0, n = length(rows),
           status = "fallback_equal_weights: too few known labels")
    k <- which(years == Y); g <- gr_predict(P[pred[k], , drop = FALSE], fit)
    raw[k] <- g$raw; prob[k] <- g$probability; clip[k] <- g$clipped
    fits[[as.character(Y)]] <- data.frame(year = Y, member = names(fit$weights), weight = unname(fit$weights),
      intercept = fit$intercept, n_rows = fit$n,
      last_outcome = if (length(rows)) month_string(max(outcome_id[rows])) else NA_character_, status = fit$status)
  }
  list(raw = raw, probability = prob, clipped = clip, weights = do.call(rbind, fits))
}
