# Ensemble configuration, frozen before any test result is inspected (spec rule 50).
# Shared settings (horizons, splits, costs, common RF, risk aversion) come from the baseline config.
ensemble_config <- function() list(
  version = "ensemble_v3_2026_11_10",
  # One member per model family. The five-input probit/logit baseline is a parsimonious benchmark,
  # not an ensemble member. Each entry: the member's rule-66 validation and test exports, plus its
  # committed performance file (used only to cross-check the shared evaluator; may be NULL).
  members = list(
    random_forest = list(validation = "Results/RF/validation_forecasts.csv",
                         test = "Results/RF/test_forecasts.csv",
                         performance = "Results/RF/performance.csv"),
    elastic_net = list(validation = "Results/elastic_net/validation_forecasts.csv",
                       test = "Results/elastic_net/test_forecasts.csv",
                       performance = "Results/elastic_net/performance.csv"),
    boosting = list(validation = "Results/boosting/boosting_validation_probabilities.csv",
                    test = "Results/boosting/boosting_test_probabilities.csv",
                    performance = "Results/boosting/boosting_performance.csv")),  # older column names are mapped in run_ensemble.R
  ensembles = c("ens_majority_vote", "ens_equal", "ens_gr"),
  probability_cutoff = 0.50,   # rule 57: fixed for the probability ensembles, no search
  gr_intercept = TRUE,         # rule 53: constrained GR with intercept, clipped at prediction (rule 56)
  gr_fit = "once per horizon on eligible validation labels (rules 53-54)",
  vote_rule = "each member's frozen validation cutoff; stocks if sum(signals) > M/2; ties to bills (rule 51)"
)
