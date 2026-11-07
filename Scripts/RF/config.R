# Hyperparameter grid. simplicity_rank orders candidates from simplest to most
# flexible and is used ONLY to break exact Brier ties: larger terminal nodes
# first, then fewer split candidates, then the default bootstrap before half-subsampling.
rf_hyper_grid <- function() {
  g <- expand.grid(mtry = c(4L, 11L, 33L), min_node_size = c(10L, 30L, 60L),
                   sampling = c("boot", "sub50"), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g$replace <- g$sampling == "boot"
  g$sample_fraction <- ifelse(g$replace, 1, 0.5)
  g <- g[order(-g$min_node_size, g$mtry, !g$replace), ]
  g$config_id <- sprintf("rf_mtry%02d_node%02d_%s", g$mtry, g$min_node_size, g$sampling)
  g$simplicity_rank <- seq_len(nrow(g))
  g$is_ranger_default <- g$mtry == 11L & g$min_node_size == 10L & g$replace
  rownames(g) <- NULL
  g[, c("config_id", "mtry", "min_node_size", "replace", "sample_fraction", "simplicity_rank", "is_ranger_default")]
}

# One candidate as a list for rf_fit().
rf_theta <- function(cfg, config_id) {
  r <- cfg$hyper_grid[cfg$hyper_grid$config_id == config_id, ]
  ensure(nrow(r) == 1L, paste("Unknown config_id", config_id))
  list(config_id = r$config_id, mtry = r$mtry, min_node_size = r$min_node_size,
       replace = r$replace, sample_fraction = r$sample_fraction, num_trees = cfg$num_trees)
}

# Freeze every random-forest research choice before validation or test evaluation.
# Shared helpers (month_number, ensure, write_csv, file_hash) and the shared
# portfolio evaluator are sourced from the baseline folder by run_rf.R.
rf_config <- function() list(
  version = "rf_v2_supervisor_2026_11_05", model = "random_forest", horizons = c(1L, 3L, 6L),
  n_features = 132L,
  # Never predictive inputs (rule 3). Asserted absent from the design matrix.
  forbidden_inputs = c("ym", "outcome_end", "equity_forward", "bill_accrual_proxy",
                       "target", "split", "mature_target_lag", "row", "origin"),
  # Fixed, finite hyperparameter grid declared before evaluation (rule 24):
  # mtry x min.node.size x sampling = 3 x 3 x 2 = 18 candidates per horizon.
  # The ranger defaults (mtry 11, node 10, bootstrap) are one of the candidates.
  hyper_grid = rf_hyper_grid(),
  num_trees = 1000L,                 # fixed, not tuned: more trees only reduce Monte Carlo noise
  splitrule = "gini",
  hyperparameter_rule = "Stage 1: per horizon, lowest validation Brier across the grid (rules 25-26); exact ties (1e-12) go to the simpler candidate. Stage 2: cutoff by validation matched-ladder Sharpe with the winner's saved forecasts (rules 27-28).",
  brier_tie_tolerance = 1e-12,
  seed_base = 4308L,
  seed_rule = "seed = 4308 + 10 * month_number(origin) + h; identical across grid candidates (common random numbers)",
  train_end = "2004-12", validation_start = "2005-01", validation_end = "2012-12",
  test_start = "2013-01", test_end = "2025-12",
  # Cutoff grid for the random forest only (the baseline keeps 0.00-1.00 by 0.01):
  # 0.30, 0.35, ..., 0.90 (13 values). Built from integers so 0.50 is exact.
  cutoffs = seq(30L, 90L, by = 5L) / 100, primary_cost = 0.001,
  costs = c(0, 0.001, 0.002), risk_aversion = 3,
  validation_refit = "monthly_expanding", test_refit = "monthly_expanding",
  cutoff_objective = "sharpe_matched_ladder", comparison_cutoff = 0.50,
  cutoff_rules = c("validation_tuned", "fixed_050"),
  common_rf_security = 2000061L,
  common_rf_rule = "4-week daily log yield at prior month-end, accrued over the next monthly quote interval",
  initial_allocation_cost = FALSE,
  # Interpretation only; never used for any selection.
  importance_reps = 5L,
  importance_rule = "Grouped (all six lags of a base variable) out-of-sample permutation across origins of the same period; delta Brier on eligible labels",
  preprocessing = "None: trees are scale-invariant; no imputation (inputs complete), no feature selection",
  bill_rule = "Origin-vintage matched daily log yield, accrued over actual monthly quote gaps",
  cost_rule = "Cost rate times absolute net equity trade; exact post-cost target weights"
)
