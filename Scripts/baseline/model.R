# Baseline: probit and logit on five predictors, refitted every month on an expanding window.
baseline_config <- function() modifyList(base_config(), list(
  model = "baseline", links = c("probit", "logit"),
  features = paste0(c("bill_yld_rel36", "term_spread_change1", "log_rv", "mom6", "log_dispersion_change1"), "_lag0"),
  bootstrap_replicates = 4999L, bootstrap_mean_block = 12L, bootstrap_block_sensitivity = c(6L, 24L)))

# Standardise on the fitting rows, then fit the binary GLM.
fit_binary <- function(d, rows, link) {
  ensure(length(rows) > ncol(d$X) + 2L && length(unique(d$y[rows])) == 2L, "Insufficient observations or classes")
  x <- d$X[rows, , drop = FALSE]
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  ensure(all(is.finite(scale) & scale > 1e-12), "Constant training predictor")
  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  design <- cbind(`(Intercept)` = 1, z)
  warnings <- character()
  fit <- withCallingHandlers(
    stats::glm.fit(design, d$y[rows], family = stats::binomial(link),
                   control = stats::glm.control(epsilon = 1e-10, maxit = 100L)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  ensure(fit$converged && fit$rank == ncol(design) && all(is.finite(fit$coefficients)) && !fit$boundary,
         paste("Failed", link, "fit"))
  ensure(!length(warnings), paste("GLM warning:", paste(warnings, collapse = "; ")))
  list(link = link, coefficients = fit$coefficients, center = center, scale = scale, rows = rows)
}

predict_binary <- function(fit, x) {
  z <- sweep(sweep(x, 2L, fit$center, "-"), 2L, fit$scale, "/")
  eta <- as.vector(cbind(1, z) %*% fit$coefficients)
  if (fit$link == "probit") stats::pnorm(eta) else stats::plogis(eta)
}

# Refit at every origin of the period using only outcomes known by that origin.
expanding_forecasts <- function(d, h, link, cfg, period) {
  rows <- period_rows(d, cfg, period)
  fits <- lapply(rows, function(i) fit_binary(d, mature_rows(d, d$m[i]), link))
  p <- vapply(seq_along(rows), function(j) predict_binary(fits[[j]], d$X[rows[j], , drop = FALSE]), 0)
  train <- lapply(fits, `[[`, "rows")
  f <- forecast_table(d, rows, h, link, period, p, lengths(train),
    vapply(train, function(r) max(d$m[r]), 0L), vapply(train, function(r) max(d$end[r]), 0L),
    vapply(train, function(r) mean(d$y[r]), 0), link)
  coefficients <- bind_rows(lapply(seq_along(rows), function(j) data.frame(h = h, model = link, period = period,
    origin = d$ym[rows[j]], term = names(fits[[j]]$coefficients), coefficient_standardized = unname(fits[[j]]$coefficients))))
  list(forecasts = f, coefficients = coefficients)
}

# Circular blocks with geometric lengths, for the stationary bootstrap.
stationary_indices <- function(n, mean_block) {
  ids <- integer(n); filled <- 0L
  while (filled < n) {
    size <- min(n - filled, stats::rgeom(1L, 1 / mean_block) + 1L)
    start <- sample.int(n, 1L)
    ids[filled + seq_len(size)] <- (start - 1L + seq_len(size) - 1L) %% n + 1L
    filled <- filled + size
  }
  ids
}

# Block-bootstrap 95% intervals for AUC, conditional on the saved forecasts.
auc_intervals <- function(f, cfg) {
  f <- eligible(f)
  ensure(all(diff(month_number(f$origin)) == 1L), "Incomplete diagnostic calendar")
  out <- data.frame(period = f$period[1], h = f$h[1], model = f$model[1])
  for (block in c(cfg$bootstrap_mean_block, cfg$bootstrap_block_sensitivity)) {
    set.seed(cfg$seed + 100L * f$h[1] + block + if (f$period[1] == "test") 10000L else 0L)
    draws <- replicate(cfg$bootstrap_replicates, {
      ids <- stationary_indices(nrow(f), block)
      auc(f$target[ids], f$probability[ids])
    })
    valid <- is.finite(draws)
    ensure(sum(valid) >= .9 * length(draws), "Too few bootstrap samples contain both outcomes")
    ci <- stats::quantile(draws[valid], c(.025, .975), names = FALSE, type = 7)
    suffix <- if (block == cfg$bootstrap_mean_block) "" else paste0("_block", block)
    out[[paste0("AUC_CI_lower", suffix)]] <- ci[1]
    out[[paste0("AUC_CI_upper", suffix)]] <- ci[2]
  }
  out
}
