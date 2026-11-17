# ADF and KPSS tests of each predictor over the training period (1974-04 to 2004-12).

# ADF lag (0 to cap) with the lowest AIC, compared on a common sample.
select_adf_lag <- function(x, cap) {
  delta <- diff(x); lagmat <- embed(delta, cap + 1L)
  y <- lagmat[, 1]
  deterministic <- cbind(constant = 1, level = x[(cap + 1L):(length(x) - 1L)])
  score <- vapply(0:cap, function(p) {
    X <- if (p) cbind(deterministic, lagmat[, 2:(p + 1L), drop = FALSE]) else deterministic
    fit <- lm.fit(X, y)
    length(y) * log(sum(fit$residuals^2) / length(y)) + 2 * fit$rank
  }, numeric(1))
  which.min(score) - 1L
}

# Automatic Bartlett bandwidth for KPSS.
kpss_bandwidth <- function(residuals) {
  n <- length(residuals); q <- floor(n^(2 / 9))
  s0 <- sum(residuals^2) / n; s1 <- 0
  for (i in seq_len(q)) {
    cross <- sum(residuals[(i + 1):n] * residuals[1:(n - i)]) * 2 / n
    s0 <- s0 + cross; s1 <- s1 + i * cross
  }
  ensure(is.finite(s0) && s0 != 0, "Undefined KPSS bandwidth")
  min(n - 1L, floor(1.1447 * ((s1 / s0)^2)^(1 / 3) * n^(1 / 3)))
}

# ADF (null: unit root) and KPSS (null: stationary) at 5%, with an intercept.
stationarity_one <- function(x, m, cap = 12L) {
  assert_monthly(m); assert_numeric(x)
  observed <- which(!is.na(x))
  ensure(length(observed) > 0, "No observations")
  rows <- seq.int(min(observed), max(observed)); x <- x[rows]; m <- m[rows]
  ensure(!anyNA(x), "Missing months inside the series")
  ensure(length(x) >= 24 && diff(range(x)) > 10 * .Machine$double.eps * max(abs(x)), "Too short or constant series")
  x <- as.numeric(scale(x)); n <- length(x)
  cap <- min(cap, floor(n / 2) - 2L)
  lag <- select_adf_lag(x, cap)
  a <- urca::ur.df(x, type = "drift", lags = lag, selectlags = "Fixed")
  bw <- kpss_bandwidth(x - mean(x))
  k <- urca::ur.kpss(x, type = "mu", use.lag = bw)
  astat <- as.numeric(a@teststat[1, 1]); acv <- as.numeric(a@cval[1, "5pct"])
  kstat <- k@teststat; kcv <- as.numeric(k@cval[1, "5pct"])
  reject_adf <- astat < acv; reject_kpss <- kstat > kcv
  verdict <- if (reject_adf && !reject_kpss) "stationarity supported" else if (!reject_adf && reject_kpss)
    "nonstationarity concern" else if (reject_adf) "conflicting tests" else "inconclusive"
  lb_lag <- min(max(12, lag + 3), floor(length(a@res) / 3))
  data.table(n = n, start = month_string(min(m)), end = month_string(max(m)),
    adf_stat = astat, adf_cv_5pct = acv, adf_reject = reject_adf, adf_lag_cap = cap, adf_lags = lag,
    adf_residual_lb_p = if (lb_lag > lag) Box.test(a@res, lag = lb_lag, type = "Ljung-Box", fitdf = lag)$p.value else NA_real_,
    kpss_stat = kstat, kpss_cv_5pct = kcv, kpss_reject = reject_kpss, kpss_lags = bw, verdict = verdict)
}

# Raw and transformed predictors at each horizon; every transformed predictor must pass.
stationarity_table <- function(predictors) {
  out <- list()
  for (h in names(predictors)) for (kind in c("raw", "transformed")) {
    d <- predictors[[h]][[if (kind == "raw") "raw" else "model"]]
    d <- csv_round_trip(d[m >= START & m <= TRAIN_END])
    for (name in setdiff(names(d), "m"))
      out[[length(out) + 1L]] <- cbind(data.table(h = as.integer(h), kind = kind, predictor = name), stationarity_one(d[[name]], d$m))
  }
  out <- rbindlist(out)
  ensure(all(out[kind == "transformed"]$verdict == "stationarity supported"), "A transformed predictor fails the stationarity tests")
  out
}
