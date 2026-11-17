# Build the 22 predictors and their 132 lagged features from the raw CRSP and macro files.
DATA <- file.path(ROOT, "Data")
INPUTS <- c(stock = "crsp_monthly_stock_data/crsp_monthly_stock_data.csv",
  monthly = "crsp_monthly_index_data/crsp_monthly_index_data.csv",
  daily = "crsp_daily_index_data/crsp_daily_index_data.csv",
  bills = "crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv",
  fixed = "crsp_treasury_fixed_term_index_data/crsp_treasury_fixed_term_index_data.csv",
  cpi = "macro/CPIAUCNS.csv", unrate = "macro/UNRATE.csv", baa = "macro/BAA10YM.csv")
FLAGS <- c("PrimaryExch", "ConditionalType", "TradingStatusFlg", "USIncFlg", "IssuerType",
           "SecurityType", "SecuritySubType", "ShareType", "MthDelFlg")
NUMERIC <- c("MthRet", "MthRetx", "MthPrc", "MthVol", "ShrOut", "MthCap")
RAW_NAMES <- c("mkt_excess", "mom3", "mom6", "mom12", "log_rv", "downside_share", "drawdown12",
  "log_dp_proxy", "bill_yld_matched", "bill_yld_chg", "term_spread", "short_slope", "cpi_infl",
  "breadth", "ma_part", "hi_lo", "dispersion", "ew_vw", "tail_imb", "turnover_dt", "unrate", "baa10ym")
# Nine predictors are transformed before lagging: monthly changes or deviations from the prior 36-month mean.
SPEC <- data.table(
  raw = c("log_dp_proxy", "bill_yld_matched", "term_spread", "short_slope", "cpi_infl", "dispersion", "tail_imb", "unrate", "baa10ym"),
  model = c("log_dp_change1", "bill_yld_rel36", "term_spread_change1", "short_slope_change1", "cpi_infl_change1",
            "log_dispersion_change1", "tail_imb_change1", "unrate_rel36", "baa10ym_rel36"),
  method = c("difference", "demean", "difference", "difference", "difference", "log_difference", "difference", "demean", "demean"))
MODEL_NAMES <- RAW_NAMES
MODEL_NAMES[match(SPEC$raw, RAW_NAMES)] <- SPEC$model
FEATURE_NAMES <- unlist(lapply(MODEL_NAMES, function(n) paste0(n, "_lag", 0:5)), use.names = FALSE)
SPEC_ID <- "r_22_predictors_continuous_crsp_yields_v3"
START <- month_number("1974-04")
TRAIN_END <- month_number(base_config()$train_end)
VAL_END <- month_number(base_config()$validation_end)
TEST_END <- month_number(base_config()$test_end)

assert_monthly <- function(m, label = "series")
  ensure(length(m) && !anyNA(m) && !anyDuplicated(m) && all(diff(m) == 1L), paste(label, "needs a complete monthly calendar"))
assert_numeric <- function(x, label = "values", missing = TRUE) {
  v <- as.matrix(x)
  ensure(is.numeric(v) && !any(is.infinite(v)) && (missing || !anyNA(v)), paste(label, "has invalid numbers"))
}
lag1 <- function(x, k = 1L) shift(x, k, type = "lag")

# Rolling statistic over complete windows only.
roll <- function(x, n, fun = mean) {
  out <- rep(NA_real_, length(x))
  if (length(x) >= n) for (i in n:length(x)) {
    w <- x[(i - n + 1):i]
    if (all(is.finite(w))) out[i] <- fun(w)
  }
  out
}

read_source <- function(key, date) {
  d <- fread(file.path(DATA, INPUTS[[key]]), na.strings = c("", "NA", "."))
  d[, m := month_number(get(date))]
  setorderv(d, c("m", date), na.last = FALSE)
  d
}

# Values pass through CSV once, as in the original build, so the saved data match it exactly.
csv_round_trip <- function(d) {
  file <- tempfile(fileext = ".csv")
  on.exit(unlink(file))
  fwrite(d, file)
  fread(file)
}

# CRSP TMYLD is a continuously compounded daily yield: accrual over n days is exp(n * TMYLD) - 1.
daily_accrual <- function(x, days) {
  assert_numeric(x); assert_numeric(days)
  ensure(!any(x <= -1, na.rm = TRUE) && !any(days < 0, na.rm = TRUE), "Invalid yield or day count")
  ensure(length(days) == 1L || length(days) == length(x), "Yield and day-count lengths differ")
  expm1(days * x)
}
daily_annual <- function(x) daily_accrual(x, 365)

transform_predictors <- function(raw) {
  assert_monthly(raw$m)
  out <- copy(raw)
  for (i in seq_len(nrow(SPEC))) {
    x <- raw[[SPEC$raw[i]]]
    value <- switch(SPEC$method[i], difference = x - lag1(x), demean = x - roll(lag1(x), 36),
      log_difference = { ensure(!any(x <= 0, na.rm = TRUE), "Nonpositive dispersion"); log(x) - lag1(log(x)) })
    out[, (SPEC$raw[i]) := value]
    setnames(out, SPEC$raw[i], SPEC$model[i])
  }
  setcolorder(out, c("m", MODEL_NAMES))
  out
}

lag_features <- function(model) {
  assert_monthly(model$m)
  out <- data.table(m = model$m)
  for (n in MODEL_NAMES) for (k in 0:5) out[, (paste0(n, "_lag", k)) := lag1(model[[n]], k)]
  out
}

# Common US operating-company shares on NYSE, AMEX or NASDAQ that are actively trading.
eligible_stock <- function(s) {
  out <- s$ShareType == "NS" & s$SecurityType == "EQTY" & s$SecuritySubType == "COM" &
    s$USIncFlg == "Y" & s$IssuerType %in% c("ACOR", "CORP") & s$PrimaryExch %in% c("N", "A", "Q") &
    s$ConditionalType == "RW" & s$TradingStatusFlg == "A"
  out[is.na(out)] <- FALSE
  out
}

# Per-stock flags: price above its 10-month average, and 12-month high or low, from complete windows.
stock_signals <- function(m, rx, price) {
  ensure(!anyDuplicated(m) && all(diff(m) > 0), "Unsorted or duplicate stock months")
  full <- seq.int(min(m), max(m)); pos <- match(m, full)
  r <- p <- rep(NA_real_, length(full)); r[pos] <- rx; p[pos] <- price
  ma <- hi <- lo <- rep(NA_real_, length(full))
  for (n in c(10L, 12L)) {
    if (length(full) < n) next
    prices <- embed(p, n)[, n:1, drop = FALSE]
    returns <- embed(r, n - 1)[-1, (n - 1):1, drop = FALSE]
    valid <- rowSums(!is.finite(prices) | prices <= 0) == 0 & rowSums(!is.finite(returns) | returns <= -1) == 0
    path <- matrix(1, nrow(prices), n)
    for (j in 2:n) path[, j] <- path[, j - 1] * (1 + returns[, j - 1])
    valid <- valid & rowSums(!is.finite(path)) == 0
    end <- seq.int(n, length(full))[valid]
    path <- path[valid, , drop = FALSE]
    if (n == 10L) ma[end] <- as.numeric(path[, n] > rowMeans(path)) else {
      hi[end] <- as.numeric(path[, n] == do.call(pmax, as.data.frame(path)))
      lo[end] <- as.numeric(path[, n] == do.call(pmin, as.data.frame(path)))
    }
  }
  list(ma = ma[pos], hi = hi[pos], lo = lo[pos])
}

# Cross-sectional predictors (breadth, trend participation, highs-lows, dispersion, tails, turnover).
build_stock <- function(calendar) {
  message("Reading CRSP stock file...")
  s <- fread(file.path(DATA, INPUTS[["stock"]]), select = c("PERMNO", "MthCalDt", FLAGS, NUMERIC), na.strings = c("", "NA"))
  s[, (NUMERIC) := lapply(.SD, as.numeric), .SDcols = NUMERIC]
  s[, m := month_number(MthCalDt)]
  ensure(!anyNA(s$PERMNO) && !anyNA(s$m), "Missing stock keys")
  assert_numeric(s[, ..NUMERIC])
  ensure(!any(as.matrix(s[, .(MthRet, MthRetx)]) < -1, na.rm = TRUE), "Return below -100%")
  s[, strict := eligible_stock(.SD)]
  # Duplicate security-months may differ in classification but never in numbers; keep the eligible row.
  duplicates <- s[duplicated(s, by = c("PERMNO", "m")) | duplicated(s, by = c("PERMNO", "m"), fromLast = TRUE)]
  conflicts <- duplicates[, lapply(.SD, uniqueN, na.rm = FALSE), by = .(PERMNO, m), .SDcols = NUMERIC]
  ensure(!any(conflicts[, rowSums(as.matrix(.SD) > 1) > 0, .SDcols = NUMERIC]), "Ambiguous numeric duplicate keys")
  setorderv(s, c("PERMNO", "m", "strict"), c(1, 1, -1))
  s <- unique(s, by = c("PERMNO", "m"))
  # Keep a delisting month that directly follows an eligible month.
  s[, terminal := TradingStatusFlg == "D" & !is.na(MthDelFlg) & MthDelFlg != "N" & !is.na(MthRet) &
      (m - shift(m) == 1L) & shift(strict, fill = FALSE) & !strict, by = PERMNO]
  s[is.na(terminal), terminal := FALSE]
  s[, recommended := strict | terminal]
  message("Computing stock windows...")
  s[, c("ma", "hi", "lo") := stock_signals(m, MthRetx, MthPrc), by = PERMNO]
  e <- s[recommended & !is.na(MthRet)]
  cross <- e[, {
    bounds <- quantile(MthRet, c(.01, .99), type = 7, names = FALSE)
    movers <- sum(MthRet != 0)
    list(breadth = if (movers) sum(MthRet > 0) / movers else NA_real_, ma_part = mean(ma, na.rm = TRUE),
         hi_lo = mean(hi, na.rm = TRUE) - mean(lo, na.rm = TRUE), dispersion = sd(pmin(pmax(MthRet, bounds[1]), bounds[2])),
         tail_imb = mean(MthRet < -0.2) - mean(MthRet > 0.2))
  }, by = m]
  act <- s[strict & PrimaryExch %in% c("N", "A") & MthVol >= 0 & ShrOut > 0]
  turnover <- act[, .(volume = sum(MthVol), shares = sum(ShrOut)), by = m]
  turnover <- merge(data.table(m = calendar), turnover, by = "m", all.x = TRUE, sort = TRUE)
  turnover[, turnover_fraction := volume / (1000 * shares)]
  turnover[, turnover_dt := log(turnover_fraction / roll(lag1(turnover_fraction), 12))]
  merge(cross, turnover[, .(m, turnover_dt)], by = "m", all = TRUE, sort = TRUE)
}

# Latest macro value released by each month-end (release assumed the month after, unless overridden).
available_macro <- function(m, x, calendar, overrides = integer()) {
  assert_monthly(m, "macro source"); assert_numeric(x)
  release <- m + 1L
  for (n in names(overrides)) release[m == month_number(n)] <- overrides[[n]]
  vapply(calendar, function(t) {
    known <- which(!is.na(x) & release <= t)
    if (length(known)) x[tail(known, 1)] else NA_real_
  }, numeric(1))
}

# Log of trailing 12-month distributions relative to price (dividend-yield proxy).
distribution_proxy <- function(r, rx) {
  ensure(!any(rx <= -1, na.rm = TRUE) && !any(r < -1, na.rm = TRUE), "Invalid index returns")
  out <- rep(NA_real_, length(r))
  if (length(r) >= 12) for (i in 12:length(r)) {
    rows <- (i - 11):i
    if (any(!is.finite(c(r[rows], rx[rows])))) next
    value <- sum((r[rows] - rx[rows]) / rev(cumprod(rev(1 + rx[rows]))))
    if (value > 0) out[i] <- log(value)
  }
  out
}

# Raw, transformed and lagged predictors for each horizon (only the matched bill yield differs).
build_predictors <- function() {
  idx <- read_source("monthly", "mthcaldt"); assert_monthly(idx$m, "monthly index")
  assert_numeric(idx[, .(vwretd, vwretx, ewretd)], missing = FALSE)
  ensure(all(as.matrix(idx[, .(vwretd, vwretx, ewretd)]) > -1), "Nonpositive index growth")
  daily <- read_source("daily", "dlycaldt")
  ensure(!anyNA(daily$dlycaldt) && !anyDuplicated(daily$dlycaldt) && all(daily$vwretd > -1), "Bad daily index data")
  fixed <- read_source("fixed", "MCALDT")[KYTREASNOX == 2000007]
  assert_monthly(fixed$m, "fixed term")
  calendar <- seq.int(min(idx$m), min(max(idx$m), max(daily$m), max(fixed$m)))
  idx <- idx[match(calendar, m)]
  bills <- read_source("bills", "MCALDT")[!is.na(m) & KYTREASNOX %in% c(2000061, 2000062, 2000063)]
  ensure(!anyDuplicated(bills, by = c("KYTREASNOX", "m")), "Duplicate Treasury keys")
  by <- lapply(c(2000061, 2000062, 2000063), function(code) bills[KYTREASNOX == code][match(calendar, m)]$TMYLD)
  names(by) <- c("1", "3", "6")
  annual <- lapply(by, daily_annual)
  fy <- fixed[match(calendar, m)]$TMYTM
  quote_days <- c(NA_real_, diff(as.numeric(as.Date(idx$mthcaldt))))

  # Market predictors from monthly and daily index returns and Treasury yields.
  base <- data.table(m = calendar, mkt_excess = idx$vwretd - daily_accrual(lag1(by[["1"]]), quote_days))
  for (n in c(3, 6, 12)) base[, (paste0("mom", n)) := expm1(roll(log1p(idx$vwretd), n, sum))]
  risk <- daily[, .(n = .N, sq = sum(vwretd^2), down = sum(vwretd[vwretd < 0]^2)), by = m][match(calendar, m)]
  ensure(!anyNA(risk$n) && all(risk$n >= 15) && all(risk$sq > 0), "Incomplete daily-risk coverage")
  base[, `:=`(log_rv = .5 * log(252 * risk$sq / risk$n), downside_share = risk$down / risk$sq)]
  wealth <- exp(cumsum(log1p(idx$vwretd)))
  base[, `:=`(drawdown12 = wealth / roll(wealth, 12, max) - 1, log_dp_proxy = distribution_proxy(idx$vwretd, idx$vwretx),
    bill_yld_chg = annual[["3"]] - lag1(annual[["3"]]), term_spread = daily_annual(fy / 36500) - annual[["3"]],
    short_slope = annual[["6"]] - annual[["1"]], ew_vw = idx$ewretd - idx$vwretd)]

  # Macro predictors, using only values released by each month-end.
  cpi <- read_source("cpi", "observation_date"); assert_monthly(cpi$m, "CPI")
  ensure(!any(cpi$CPIAUCNS <= 0, na.rm = TRUE), "Invalid CPI")
  base[, cpi_infl := available_macro(cpi$m, log(cpi$CPIAUCNS / lag1(cpi$CPIAUCNS, 12)), calendar)]
  unrate <- read_source("unrate", "observation_date"); assert_monthly(unrate$m, "UNRATE")
  ensure(!any(unrate$UNRATE < 0 | unrate$UNRATE > 100, na.rm = TRUE), "Invalid UNRATE percent")
  # The September 2025 unemployment release was delayed to November.
  base[, unrate := available_macro(unrate$m, unrate$UNRATE / 100, calendar, setNames(month_number("2025-11"), "2025-09"))]
  baa <- read_source("baa", "observation_date"); assert_monthly(baa$m, "BAA10YM")
  base[, baa10ym := available_macro(baa$m, baa$BAA10YM / 100, calendar)]

  base <- merge(base, build_stock(calendar), by = "m", all.x = TRUE, sort = TRUE)
  out <- list()
  for (h in c(1L, 3L, 6L)) {
    raw <- copy(base)
    raw[, bill_yld_matched := annual[[as.character(h)]]]
    raw <- raw[, c("m", RAW_NAMES), with = FALSE]
    model <- transform_predictors(raw)
    features <- lag_features(model)[m >= START & m <= TEST_END]
    assert_numeric(features[, ..FEATURE_NAMES], missing = FALSE)
    ensure(ncol(features) == 133L, "Expected 132 features plus month")
    out[[as.character(h)]] <- list(raw = raw, model = model, features = csv_round_trip(features))
  }
  message("Built 22 predictors / 132 features at each horizon.")
  out
}
