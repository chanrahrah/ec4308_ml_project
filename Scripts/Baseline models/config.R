# Freeze research choices before validation or test evaluation.
baseline_config <- function() list(
  version = "baseline_v2_supervisor_2026_11_05", horizons = c(1L, 3L, 6L), links = c("probit", "logit"),
  features = paste0(c("bill_yld_rel36", "term_spread_change1", "log_rv",
                      "mom6", "log_dispersion_change1"), "_lag0"),
  train_end = "2004-12", validation_start = "2005-01", validation_end = "2012-12",
  test_start = "2013-01", test_end = "2025-12",
  cutoffs = seq(0, 1, by = 0.01), primary_cost = 0.001,
  costs = c(0, 0.001, 0.002), risk_aversion = 3,
  validation_refit = "monthly_expanding", test_refit = "monthly_expanding",
  cutoff_objective = "sharpe_matched_ladder", comparison_cutoff = 0.50,
  cutoff_rules = c("validation_tuned", "fixed_050"),
  common_rf_security = 2000061L,
  common_rf_rule = "4-week daily log yield at prior month-end, accrued over the next monthly quote interval",
  # Initial holdings are endowed, as in the professor's initialized example.
  # No opening/liquidation fee; subsequent net equity trades incur costs.
  initial_allocation_cost = FALSE,
  lag_rule = "Five origin-t observations (dataset lag0); no additional distributed lags",
  bill_rule = "Origin-vintage matched daily log yield, accrued over actual monthly quote gaps",
  cost_rule = "Cost rate times absolute net equity trade; exact post-cost target weights"
)

month_number <- function(x) {
  x <- as.character(x)
  as.integer(substr(x, 1, 4)) * 12L + as.integer(substr(x, 6, 7)) - 1L
}
month_string <- function(x) sprintf("%04d-%02d", x %/% 12L, x %% 12L + 1L)
ensure <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
file_hash <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
