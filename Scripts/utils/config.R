# Settings shared by every model: sample splits, cutoff grid, trading costs and seed.
base_config <- function() list(
  horizons = c(1L, 3L, 6L),
  train_end = "2004-12", validation_start = "2005-01", validation_end = "2012-12",
  test_start = "2013-01", test_end = "2025-12",
  cutoffs = seq(0L, 100L, by = 5L) / 100, comparison_cutoff = 0.50,
  cutoff_objective = "sharpe_matched_ladder", primary_cost = 0.001, costs = c(0, 0.001, 0.002),
  risk_aversion = 3, common_rf_security = 2000061L, initial_allocation_cost = FALSE,
  fixed_benchmark_weight = 0.5, seed = 4308L)

RETURN_FILES <- file.path(ROOT, "Data", c("crsp_monthly_index_data/crsp_monthly_index_data.csv",
  "crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))

model_inputs <- function(cfg)
  c(file.path(ROOT, "Data", "processed", paste0("model_ready_h", cfg$horizons, ".rds")), RETURN_FILES)
