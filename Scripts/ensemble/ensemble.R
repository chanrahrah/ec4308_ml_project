# ensemble.R: the three group ensembles (spec of 5 Nov 2026, rules 50-58 and 66-70). Pure base-R functions.
#   Hard majority vote: each member's signal from its own frozen validation cutoff; stocks when
#     sum(signals) > M/2, ties to bills. Nothing fitted.
#   Equal-probability average and constrained Granger-Ramanathan (with intercept, clipped): defined
#     at the bottom of this file; both use a fixed 0.50 cutoff.

EXPORT_COLS <- c("model", "h", "origin", "probability", "cutoff", "signal", "cutoff_050", "signal_050",
                 "target", "outcome_end", "classification_eligible")

# "YYYY-MM" or a "YYYY-MM-DD" date -> integer month id (year*12 + month - 1).
to_month_id <- function(x) {
  x <- substr(as.character(x), 1L, 7L); out <- rep(NA_integer_, length(x)); ok <- !is.na(x) & nzchar(x)
  out[ok] <- as.integer(substr(x[ok], 1L, 4L)) * 12L + as.integer(substr(x[ok], 6L, 7L)) - 1L
  out
}
from_month_id <- function(id) sprintf("%04d-%02d", id %/% 12L, id %% 12L + 1L)

# Strict majority: stocks (1) when more than half of the M members say stocks; ties go to bills (0).
majority_vote <- function(signals) {
  signals <- as.matrix(signals)
  if (!ncol(signals)) stop("No ensemble members")
  if (anyNA(signals) || !all(signals %in% c(0, 1))) stop("Member signals must be 0/1 without NA")
  as.integer(rowSums(signals) > ncol(signals) / 2)
}

# Check one member's export (one period) and return it keyed by (h, origin).
# Rule 51: the vote uses the member's own frozen cutoff, so the exported signal must equal
# probability >= cutoff, and the cutoff must be a single value per horizon.
check_member <- function(d, member, period) {
  miss <- setdiff(EXPORT_COLS, names(d))
  if (length(miss)) stop(member, " ", period, ": missing columns ", paste(miss, collapse=", "))
  d <- as.data.frame(d)[, EXPORT_COLS]
  # Unknown outcome_end (origin + h beyond the data) is written as "" by some members and NA by others.
  d$outcome_end <- as.character(d$outcome_end); d$outcome_end[!is.na(d$outcome_end) & !nzchar(d$outcome_end)] <- NA
  d$origin_id <- to_month_id(d$origin)
  if (anyNA(d$origin_id)) stop(member, " ", period, ": unreadable origin")
  if (anyDuplicated(d[, c("h", "origin_id")])) stop(member, " ", period, ": duplicate (h, origin) keys")
  if (anyNA(d$probability) || any(d$probability < 0 | d$probability > 1)) stop(member, " ", period, ": probability outside [0, 1] or NA")
  for (h in unique(d$h)) {
    x <- d[d$h == h, ]
    if (length(unique(x$cutoff)) != 1L || length(unique(x$cutoff_050)) != 1L)
      stop(member, " ", period, " h=", h, ": cutoff is not frozen (varies across origins)")
    if (unique(x$cutoff_050) != 0.5) stop(member, " ", period, " h=", h, ": cutoff_050 is not 0.50")
    o <- sort(x$origin_id)
    if (any(diff(o) != 1L)) stop(member, " ", period, " h=", h, ": monthly calendar has gaps")
  }
  if (!identical(as.integer(d$signal), as.integer(d$probability >= d$cutoff)))
    stop(member, " ", period, ": signal differs from probability >= cutoff")
  if (!identical(as.integer(d$signal_050), as.integer(d$probability >= 0.5)))
    stop(member, " ", period, ": signal_050 differs from probability >= 0.50")
  d$member <- member
  d
}

# Merge members by (h, origin) keys, never by row order. Members must cover identical keys and agree
# on target, outcome_end and classification eligibility. Returns one row per (h, origin) with member
# signal matrices and the votes.
build_vote <- function(members, period) {
  keyed <- lapply(names(members), function(m) check_member(members[[m]], m, period))
  names(keyed) <- names(members)
  base <- keyed[[1]][order(keyed[[1]]$h, keyed[[1]]$origin_id), c("h", "origin_id", "target", "outcome_end", "classification_eligible")]
  key <- paste(base$h, base$origin_id)
  S <- S050 <- P <- matrix(NA_real_, nrow(base), length(keyed), dimnames=list(NULL, names(keyed)))
  for (m in names(keyed)) {
    x <- keyed[[m]]; i <- match(key, paste(x$h, x$origin_id))
    if (anyNA(i) || nrow(x) != nrow(base)) stop(period, ": ", m, " does not cover the same (h, origin) keys as ", names(keyed)[1])
    same <- function(a, b) identical(is.na(a), is.na(b)) && all(a[!is.na(a)] == b[!is.na(b)])
    if (!same(x$target[i], base$target)) stop(period, ": ", m, " target differs from ", names(keyed)[1])
    if (!same(as.character(x$outcome_end[i]), as.character(base$outcome_end))) stop(period, ": ", m, " outcome_end differs")
    if (!same(as.logical(x$classification_eligible[i]), as.logical(base$classification_eligible))) stop(period, ": ", m, " eligibility differs")
    S[, m] <- x$signal[i]; S050[, m] <- x$signal_050[i]; P[, m] <- x$probability[i]
  }
  list(base=base, S=S, S050=S050, P=P,
       cutoffs=do.call(rbind, lapply(names(keyed), function(m) unique(keyed[[m]][, c("member", "h", "cutoff")]))),
       vote=majority_vote(S), vote_050=majority_vote(S050))
}

# Mann-Whitney AUC (ties count one half); NA if a class is absent.
auc <- function(score, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ------------------------------------------------------------------------- probability ensembles
# Rule 52: equal-probability average. Rules 53-56: constrained Granger-Ramanathan with intercept,
#   min_(a,w) mean((y - a - P w)^2)  s.t.  w >= 0, sum(w) = 1,  a unrestricted,
# fitted by ordinary constrained least squares on eligible validation rows, then clipped to [0, 1]
# at prediction. Rule 57: both use a fixed 0.50 cutoff (no ensemble cutoff search).

# Equality-constrained least squares on the members in S (sum(w) = 1 by substitution). NULL if rank-deficient.
gr_fit_subset <- function(P, y, S, intercept) {
  k <- length(S); base <- P[, S[k]]
  X <- NULL
  if (intercept) X <- matrix(1, nrow(P), 1L)
  if (k > 1L) X <- cbind(X, P[, S[-k], drop=FALSE] - base)
  coef <- numeric(0)
  if (!is.null(X)) {
    qx <- qr(X)
    if (qx$rank < ncol(X)) return(NULL)
    coef <- qr.coef(qx, y - base)
  }
  a <- if (intercept) coef[[1]] else 0
  wz <- if (k > 1L) coef[(intercept + 1L):length(coef)] else numeric(0)
  w <- c(wz, 1 - sum(wz))
  list(S=S, w=unname(w), a=unname(a), sse=sum((y - a - as.vector(P[, S, drop=FALSE] %*% w))^2))
}

# Exact solution by enumerating supports: the optimum of this convex problem is the equality-constrained
# fit on its own support, so the best feasible subset fit is the global optimum (fine for a few members).
gr_weights <- function(P, y, intercept=TRUE) {
  P <- as.matrix(P); K <- ncol(P)
  if (is.null(colnames(P))) stop("P needs member names as column names")
  if (nrow(P) != length(y) || anyNA(P) || anyNA(y)) stop("GR needs complete, aligned P and y")
  if (K > 12L) stop("Subset enumeration is meant for a few members")
  best <- NULL
  for (mask in seq_len(2^K - 1)) {
    S <- which(bitwAnd(mask, 2L^(0:(K - 1L))) > 0)
    fit <- gr_fit_subset(P, y, S, intercept)
    if (is.null(fit) || any(fit$w < -1e-10)) next
    if (is.null(best) || fit$sse < best$sse - 1e-12) best <- fit
  }
  if (is.null(best)) stop("No feasible GR solution")
  w <- setNames(numeric(K), colnames(P)); w[best$S] <- pmax(best$w, 0); w <- w / sum(w)
  list(weights=w, intercept=best$a, sse=best$sse, n=length(y), status="ok")
}

# Rule 63 / section 9: if the solver fails, fall back to equal weights with zero intercept and flag it.
gr_weights_safe <- function(P, y, intercept=TRUE) {
  tryCatch(gr_weights(P, y, intercept), error=function(e) {
    K <- ncol(P)
    list(weights=setNames(rep(1 / K, K), colnames(P)), intercept=0, sse=NA_real_, n=length(y),
         status=paste("fallback_equal_weights:", conditionMessage(e)))
  })
}

# Raw GR score, its clipped probability and whether clipping happened (rule 56). Members matched by name.
gr_predict <- function(P, fit) {
  P <- as.matrix(P)
  if (!setequal(colnames(P), names(fit$weights))) stop("GR weights and probability columns name different members")
  raw <- fit$intercept + as.vector(P[, names(fit$weights), drop=FALSE] %*% fit$weights)
  list(raw=raw, probability=pmin(1, pmax(0, raw)), clipped=raw < 0 | raw > 1)
}

equal_average <- function(P) rowMeans(as.matrix(P))

# Expanding-window GR, refitted once a year (the professor's optional annual update; not one of the
# rule-58 ensembles; spec rules 59-65). For each calendar year Y of the prediction origins, the weights
# are fitted in January of Y using information through the preceding December: every row (from the first
# archived origin) whose outcome_end <= (Y-1)-12, and used for all origins in Y. Rows of one horizon only.
# Rule 61 check: January 2014, h = 6 -> origins 2005-01..2013-06 = 102 rows.
#   P, y, origin_id, outcome_id, eligible: the pooled history (validation + test) for one horizon
#   pred: indices of the rows to predict
gr_yearly <- function(P, y, origin_id, outcome_id, eligible, pred, intercept = TRUE, min_rows = 24L) {
  P <- as.matrix(P); years <- origin_id[pred] %/% 12L
  raw <- prob <- rep(NA_real_, length(pred)); clip <- rep(NA, length(pred)); fits <- list()
  for (Y in sort(unique(years))) {
    refit <- Y * 12L - 1L                              # December of year Y-1 (rule 60)
    rows <- which(eligible & !is.na(y) & !is.na(outcome_id) & outcome_id <= refit)
    fit <- if (length(rows) >= min_rows) gr_weights_safe(P[rows, , drop = FALSE], y[rows], intercept) else
      list(weights = setNames(rep(1 / ncol(P), ncol(P)), colnames(P)), intercept = 0, sse = NA_real_,
           n = length(rows), status = "fallback_equal_weights: too few known labels")
    k <- which(years == Y); g <- gr_predict(P[pred[k], , drop = FALSE], fit)
    raw[k] <- g$raw; prob[k] <- g$probability; clip[k] <- g$clipped
    fits[[as.character(Y)]] <- data.frame(year = Y, information_cutoff = from_month_id(refit), member = names(fit$weights),
      weight = unname(fit$weights), intercept = fit$intercept, n_rows = fit$n,
      last_known_outcome = if (length(rows)) from_month_id(max(outcome_id[rows])) else NA_character_, status = fit$status)
  }
  list(raw = raw, probability = prob, clipped = clip, weights = do.call(rbind, fits))
}
