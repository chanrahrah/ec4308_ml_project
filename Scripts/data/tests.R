# Checks of the predictor, target and stationarity calculations on small known inputs.
data_checks <- function() {
  fails <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")
  eq <- function(x, y, tol = 1e-12) ensure(isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE)), "Values differ")
  m <- seq.int(month_number("2000-01"), length.out = 60)
  raw <- data.table(m = m)
  for (n in RAW_NAMES) raw[, (n) := seq_len(60) / 100]
  raw[, dispersion := seq_len(60) / 100 + .01]
  list(
    "Monthly calendar rejects gaps, duplicates and reversals" = function()
      ensure(fails(assert_monthly(m[-4])) && fails(assert_monthly(c(m, m[60]))) && fails(assert_monthly(rev(m))), "Bad calendar accepted"),
    "CRSP daily yields compound continuously" = function() {
      y <- c(0, .0001, .0002)
      eq(daily_annual(y), exp(y * 365) - 1)
      eq(daily_annual(y), (1 + (2 * (exp(y * 182.5) - 1)) / 2)^2 - 1)
      eq(daily_accrual(.0001, 33), exp(.0001 * 33) - 1)
      ensure(fails(daily_annual(-1)), "Invalid yield accepted")
    },
    "Distribution proxy uses twelve months and cannot bridge gaps" = function() {
      v <- distribution_proxy(rep(.002, 24), rep(0, 24))
      ensure(all(is.na(v[1:11])), "Proxy defined too early"); eq(v[12:24], rep(log(.024), 13))
      rx <- rep(0, 36); rx[16] <- NA; v <- distribution_proxy(rx + .002, rx)
      ensure(all(is.na(v[16:27])), "Proxy bridged a gap"); eq(v[28], log(.024))
    },
    "Stock windows need complete months and survive price splits" = function() {
      r <- rep(.01, 12); p <- 100 * 1.01^(0:11); p[7:12] <- p[7:12] / 2
      a <- stock_signals(m[1:12], r, p)
      ensure(all(is.na(a$ma[1:9])) && all(is.na(a$hi[1:11])) && a$ma[12] == 1 && a$hi[12] == 1, "Wrong stock window")
      gap <- stock_signals(m[c(1:5, 7:12)], r[-6], p[-6])
      ensure(is.na(tail(gap$ma, 1)) && is.na(tail(gap$hi, 1)), "Window bridged a gap")
      flat <- stock_signals(m[1:12], rep(0, 12), rep(10, 12)); eq(flat$hi[12] - flat$lo[12], 0); eq(flat$ma[12], 0)
    },
    "Macro values are used only after release, including the delayed 2025 release" = function() {
      months <- seq.int(month_number("2025-07"), length.out = 5)
      v <- available_macro(months, c(.041, .042, .043, NA, .046), months + 1L, setNames(month_number("2025-11"), "2025-09"))
      eq(v, c(.041, .042, .042, .043, .046))
    },
    "Transforms use differences and prior-only rolling means" = function() {
      a <- transform_predictors(raw)
      eq(a$cpi_infl_change1[40], .01)
      eq(a$bill_yld_rel36[40], raw$bill_yld_matched[40] - mean(raw$bill_yld_matched[4:39]))
      eq(a$log_dispersion_change1[40], log(raw$dispersion[40] / raw$dispersion[39]))
      b <- copy(raw); b[45, dispersion := 0]; ensure(fails(transform_predictors(b)), "Nonpositive dispersion accepted")
      b <- copy(raw); b[45, term_spread := NA_real_]; a <- transform_predictors(b)
      ensure(all(is.na(a$term_spread_change1[45:46])) && is.finite(a$term_spread_change1[47]), "Missing value not propagated")
    },
    "Transforms and lags cannot read future observations" = function() {
      a <- transform_predictors(raw); b <- copy(raw)
      b[50:60, (RAW_NAMES) := lapply(.SD, function(x) x * 10), .SDcols = RAW_NAMES]
      eq(a[1:49], transform_predictors(b)[1:49])
      l <- lag_features(a); eq(l$cpi_infl_change1_lag5[50], a$cpi_infl_change1[45])
      ensure(ncol(l) == 133, "Expected 132 lagged features")
    },
    "Stationarity tests reject gaps and constants and match OLS and Bartlett formulas" = function() {
      x <- sin(seq_along(m)); x[30] <- NA
      ensure(fails(stationarity_one(x, m)) && fails(stationarity_one(rep(1, 60), m)), "Bad series accepted")
      set.seed(4308); x <- rnorm(369); lag <- select_adf_lag(x, 12)
      a <- stationarity_one(x, seq.int(1000L, length.out = 369))
      ensure(a$adf_lags == lag && a$verdict == "stationarity supported", "Wrong ADF lag or verdict")
      z <- as.numeric(scale(x)); d <- diff(z); E <- embed(d, lag + 1); X <- cbind(1, z[(lag + 1):(length(z) - 1)])
      if (lag) X <- cbind(X, E[, 2:(lag + 1), drop = FALSE])
      fit <- lm.fit(X, E[, 1]); se <- sqrt(diag(chol2inv(qr.R(fit$qr))) * sum(fit$residuals^2) / fit$df.residual)
      eq(a$adf_stat, unname(fit$coefficients[2] / se[2]), 1e-10)
      set.seed(1); x <- rnorm(100); res <- x - mean(x)
      variance <- sum(res^2) / 100 + 2 / 100 * sum(vapply(1:4, function(j) (1 - j / 5) * sum(res[(j + 1):100] * res[1:(100 - j)]), numeric(1)))
      eq(urca::ur.kpss(x, use.lag = 4)@teststat, sum(cumsum(res)^2) / (100^2 * variance))
    },
    "Targets use only t+1 to t+h and strict outperformance" = function() {
      months <- seq.int(month_number("2000-01"), length.out = 12)
      dates <- as.Date(paste0(month_string(months + 1), "-01")) - 1
      r <- rep(.01, 12); r[1] <- .5
      z <- forward_targets(months, r, dates, rep(0, 12), 3)
      eq(z$equity_forward[1], 1.01^3 - 1)
      ensure(z$target[1] == 1 && z$outcome_end[1] == "2000-04" && all(is.na(z$target[10:12])), "Wrong target")
      ensure(forward_targets(months, rep(0, 12), dates, rep(0, 12), 1)$target[1] == 0, "Ties must be 0")
      eq(forward_targets(months, r, dates, rep(.0001, 12), 3)$bill_accrual_proxy[1], exp(as.numeric(dates[4] - dates[1]) * .0001) - 1)
    },
    "Labels whose outcome crosses a split boundary are purged" = function() {
      origin <- month_number(c("2004-09", "2004-10", "2005-01", "2012-10", "2025-10"))
      ensure(identical(model_split(origin, month_string(origin + 3), rep(1L, 5)),
                       c("train", "purged_boundary", "validation", "purged_boundary", "purged_boundary")), "Wrong split")
    })
}
