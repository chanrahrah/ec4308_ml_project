# Forward targets, sample splits and the model-ready RDS files.

# y = 1 when compounded equity return over t+1..t+h beats the origin's matched bill accrual.
forward_targets <- function(m, returns, quote_dates, yield, h) {
  assert_monthly(m); assert_numeric(returns, missing = FALSE); assert_numeric(yield)
  ensure(h %in% c(1L, 3L, 6L) && all(returns > -1), "Invalid horizon or returns")
  date <- as.Date(quote_dates)
  ensure(!anyNA(date) && all(diff(date) > 0) && all(month_number(date) == m), "Invalid quote dates")
  out <- data.table(m = m, outcome_end = rep(NA_character_, length(m)), equity_forward = NA_real_,
                    bill_accrual_proxy = NA_real_, target = NA_integer_)
  for (i in seq_len(length(m) - h)) {
    if (!is.finite(yield[i]) || yield[i] <= -1) next
    eq <- prod(1 + returns[(i + 1):(i + h)]) - 1
    bill <- daily_accrual(yield[i], as.numeric(date[i + h] - date[i]))
    end_month <- month_string(m[i + h])
    out[i, `:=`(outcome_end = end_month, equity_forward = eq, bill_accrual_proxy = bill, target = as.integer(eq > bill))]
  }
  out
}

# A labelled row belongs to a split only if its outcome also ends inside that split.
model_split <- function(m, outcome_end, target) {
  result <- rep("unlabeled", length(m)); end <- month_number(outcome_end)
  known <- !is.na(target)
  result[known] <- "purged_boundary"
  result[known & m >= START & m <= TRAIN_END & end <= TRAIN_END] <- "train"
  result[known & m > TRAIN_END & m <= VAL_END & end <= VAL_END] <- "validation"
  result[known & m > VAL_END & m <= TEST_END & end <= TEST_END] <- "test"
  result
}

# Model inputs for horizon h: X, y, split rows and outcome metadata (never used as predictors).
model_ready <- function(features, h) {
  idx <- read_source("monthly", "mthcaldt")
  bills <- read_source("bills", "MCALDT")[!is.na(m)]
  yield <- bills[KYTREASNOX == c(`1` = 2000061, `3` = 2000062, `6` = 2000063)[[as.character(h)]]][match(idx$m, m)]$TMYLD
  targets <- forward_targets(idx$m, idx$vwretd, idx$mthcaldt, yield, h)
  full <- merge(features, targets, by = "m", all.x = TRUE, sort = TRUE)
  full[, split := model_split(m, outcome_end, target)]
  # Lagged target from origin t-h, whose outcome is known at t.
  full[, mature_target_lag := targets$target[match(m - h, targets$m)]]
  X <- as.matrix(full[, ..FEATURE_NAMES]); storage.mode(X) <- "double"
  rows <- lapply(c("train", "validation", "test"), function(s) which(full$split == s))
  names(rows) <- c("train", "validation", "test")
  list(X = X, y = full$target, ym = month_string(full$m), split = full$split, rows = rows,
    mature_target_lag = full$mature_target_lag, feature_names = FEATURE_NAMES,
    metadata = full[, .(ym = month_string(m), outcome_end, equity_forward, bill_accrual_proxy)],
    specification = SPEC_ID,
    target_definition = "1{compound equity return t+1..t+h > exp(TMYLD_origin * quote-date days)-1}; bill is a proxy")
}
