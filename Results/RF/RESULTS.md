# Random forest results - supervisor protocol

Configuration: `rf_v2_supervisor_2026_11_05`. A ranger probability forest on all 132 predictors (22 variables x lags 0-5) is fitted separately at h = 1, 3 and 6. Every validation and test origin is refitted monthly on an expanding sample of labels whose outcomes ended by that origin. Selection follows the group standard in two stages, both on 2005-2012 validation only: (1) per horizon, the hyperparameter candidate with the lowest validation Brier among 18 (mtry x min.node.size x sampling, 1000 trees); (2) with that candidate's saved validation forecasts fixed, the cutoff with the highest validation matched-ladder Sharpe at 10 bp. The fixed 0.50 cutoff is reported alongside. All-stock and all-bill benchmarks use the same scored months.

Generated automatically by `Scripts/RF/run_rf.R`. See `INTERPRETATION.md` in this folder for the written interpretation.

## Stage 1: hyperparameter selection (validation Brier, lower is better)

| Horizon | Brier_rank | Config | mtry | min_node | Sampling | Validation_Brier | Log_loss | AUC | Prevalence_Brier | Note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 |  1 | rf_mtry11_node10_sub50 | 11 | 10 | 50% subsample | 0.2423 | 0.6772 | 0.516 | 0.2406 | **selected** |
| 1 |  2 | rf_mtry04_node60_sub50 |  4 | 60 | 50% subsample | 0.2429 | 0.6789 | 0.490 | 0.2406 |  |
| 1 |  3 | rf_mtry11_node60_boot | 11 | 60 | bootstrap | 0.2433 | 0.6795 | 0.510 | 0.2406 |  |
| 1 |  4 | rf_mtry11_node60_sub50 | 11 | 60 | 50% subsample | 0.2435 | 0.6802 | 0.490 | 0.2406 |  |
| 1 |  5 | rf_mtry04_node30_sub50 |  4 | 30 | 50% subsample | 0.2437 | 0.6805 | 0.491 | 0.2406 |  |
| 1 |  6 | rf_mtry04_node60_boot |  4 | 60 | bootstrap | 0.2438 | 0.6809 | 0.497 | 0.2406 |  |
| 1 |  7 | rf_mtry11_node10_boot | 11 | 10 | bootstrap | 0.2439 | 0.6805 | 0.509 | 0.2406 |  ranger default |
| 1 |  8 | rf_mtry11_node30_sub50 | 11 | 30 | 50% subsample | 0.2441 | 0.6816 | 0.492 | 0.2406 |  |
| 1 |  9 | rf_mtry33_node30_sub50 | 33 | 30 | 50% subsample | 0.2442 | 0.6819 | 0.502 | 0.2406 |  |
| 1 | 10 | rf_mtry04_node30_boot |  4 | 30 | bootstrap | 0.2443 | 0.6818 | 0.494 | 0.2406 |  |
| 1 | 11 | rf_mtry33_node60_sub50 | 33 | 60 | 50% subsample | 0.2444 | 0.6824 | 0.483 | 0.2406 |  |
| 1 | 12 | rf_mtry11_node30_boot | 11 | 30 | bootstrap | 0.2449 | 0.6829 | 0.500 | 0.2406 |  |
| 1 | 13 | rf_mtry04_node10_boot |  4 | 10 | bootstrap | 0.2453 | 0.6838 | 0.491 | 0.2406 |  |
| 1 | 14 | rf_mtry04_node10_sub50 |  4 | 10 | 50% subsample | 0.2455 | 0.6844 | 0.494 | 0.2406 |  |
| 1 | 15 | rf_mtry33_node60_boot | 33 | 60 | bootstrap | 0.2460 | 0.6855 | 0.492 | 0.2406 |  |
| 1 | 16 | rf_mtry33_node30_boot | 33 | 30 | bootstrap | 0.2462 | 0.6859 | 0.486 | 0.2406 |  |
| 1 | 17 | rf_mtry33_node10_sub50 | 33 | 10 | 50% subsample | 0.2465 | 0.6868 | 0.489 | 0.2406 |  |
| 1 | 18 | rf_mtry33_node10_boot | 33 | 10 | bootstrap | 0.2490 | 0.6919 | 0.463 | 0.2406 |  |
| 3 |  1 | rf_mtry04_node10_boot |  4 | 10 | bootstrap | 0.2277 | 0.6470 | 0.570 | 0.2348 | **selected** |
| 3 |  2 | rf_mtry33_node30_boot | 33 | 30 | bootstrap | 0.2280 | 0.6470 | 0.560 | 0.2348 |  |
| 3 |  3 | rf_mtry04_node10_sub50 |  4 | 10 | 50% subsample | 0.2291 | 0.6500 | 0.555 | 0.2348 |  |
| 3 |  4 | rf_mtry11_node30_boot | 11 | 30 | bootstrap | 0.2292 | 0.6502 | 0.554 | 0.2348 |  |
| 3 |  5 | rf_mtry04_node30_boot |  4 | 30 | bootstrap | 0.2298 | 0.6516 | 0.540 | 0.2348 |  |
| 3 |  6 | rf_mtry11_node10_sub50 | 11 | 10 | 50% subsample | 0.2298 | 0.6511 | 0.559 | 0.2348 |  |
| 3 |  7 | rf_mtry11_node10_boot | 11 | 10 | bootstrap | 0.2305 | 0.6528 | 0.529 | 0.2348 |  ranger default |
| 3 |  8 | rf_mtry33_node10_sub50 | 33 | 10 | 50% subsample | 0.2310 | 0.6532 | 0.551 | 0.2348 |  |
| 3 |  9 | rf_mtry04_node30_sub50 |  4 | 30 | 50% subsample | 0.2311 | 0.6543 | 0.544 | 0.2348 |  |
| 3 | 10 | rf_mtry11_node30_sub50 | 11 | 30 | 50% subsample | 0.2313 | 0.6543 | 0.547 | 0.2348 |  |
| 3 | 11 | rf_mtry04_node60_boot |  4 | 60 | bootstrap | 0.2314 | 0.6549 | 0.536 | 0.2348 |  |
| 3 | 12 | rf_mtry33_node60_boot | 33 | 60 | bootstrap | 0.2314 | 0.6543 | 0.540 | 0.2348 |  |
| 3 | 13 | rf_mtry33_node30_sub50 | 33 | 30 | 50% subsample | 0.2329 | 0.6570 | 0.527 | 0.2348 |  |
| 3 | 14 | rf_mtry33_node10_boot | 33 | 10 | bootstrap | 0.2329 | 0.6574 | 0.526 | 0.2348 |  |
| 3 | 15 | rf_mtry11_node60_boot | 11 | 60 | bootstrap | 0.2329 | 0.6578 | 0.530 | 0.2348 |  |
| 3 | 16 | rf_mtry11_node60_sub50 | 11 | 60 | 50% subsample | 0.2332 | 0.6586 | 0.528 | 0.2348 |  |
| 3 | 17 | rf_mtry04_node60_sub50 |  4 | 60 | 50% subsample | 0.2339 | 0.6602 | 0.512 | 0.2348 |  |
| 3 | 18 | rf_mtry33_node60_sub50 | 33 | 60 | 50% subsample | 0.2361 | 0.6644 | 0.508 | 0.2348 |  |
| 6 |  1 | rf_mtry04_node10_sub50 |  4 | 10 | 50% subsample | 0.2189 | 0.6274 | 0.527 | 0.2216 | **selected** |
| 6 |  2 | rf_mtry04_node30_sub50 |  4 | 30 | 50% subsample | 0.2192 | 0.6282 | 0.525 | 0.2216 |  |
| 6 |  3 | rf_mtry11_node30_sub50 | 11 | 30 | 50% subsample | 0.2198 | 0.6289 | 0.510 | 0.2216 |  |
| 6 |  4 | rf_mtry33_node60_sub50 | 33 | 60 | 50% subsample | 0.2202 | 0.6302 | 0.530 | 0.2216 |  |
| 6 |  5 | rf_mtry33_node30_sub50 | 33 | 30 | 50% subsample | 0.2204 | 0.6303 | 0.527 | 0.2216 |  |
| 6 |  6 | rf_mtry04_node10_boot |  4 | 10 | bootstrap | 0.2212 | 0.6322 | 0.499 | 0.2216 |  |
| 6 |  7 | rf_mtry11_node60_sub50 | 11 | 60 | 50% subsample | 0.2218 | 0.6338 | 0.496 | 0.2216 |  |
| 6 |  8 | rf_mtry04_node60_sub50 |  4 | 60 | 50% subsample | 0.2220 | 0.6349 | 0.505 | 0.2216 |  |
| 6 |  9 | rf_mtry33_node60_boot | 33 | 60 | bootstrap | 0.2225 | 0.6353 | 0.493 | 0.2216 |  |
| 6 | 10 | rf_mtry04_node60_boot |  4 | 60 | bootstrap | 0.2230 | 0.6369 | 0.490 | 0.2216 |  |
| 6 | 11 | rf_mtry04_node30_boot |  4 | 30 | bootstrap | 0.2230 | 0.6368 | 0.496 | 0.2216 |  |
| 6 | 12 | rf_mtry11_node60_boot | 11 | 60 | bootstrap | 0.2232 | 0.6372 | 0.485 | 0.2216 |  |
| 6 | 13 | rf_mtry33_node10_sub50 | 33 | 10 | 50% subsample | 0.2235 | 0.6369 | 0.496 | 0.2216 |  |
| 6 | 14 | rf_mtry11_node10_sub50 | 11 | 10 | 50% subsample | 0.2237 | 0.6375 | 0.483 | 0.2216 |  |
| 6 | 15 | rf_mtry11_node30_boot | 11 | 30 | bootstrap | 0.2240 | 0.6388 | 0.484 | 0.2216 |  |
| 6 | 16 | rf_mtry33_node30_boot | 33 | 30 | bootstrap | 0.2247 | 0.6399 | 0.482 | 0.2216 |  |
| 6 | 17 | rf_mtry11_node10_boot | 11 | 10 | bootstrap | 0.2279 | 0.6475 | 0.462 | 0.2216 |  ranger default |
| 6 | 18 | rf_mtry33_node10_boot | 33 | 10 | bootstrap | 0.2304 | 0.6535 | 0.449 | 0.2216 |  |

Prevalence_Brier is the Brier of forecasting the training-sample event rate known at each origin; a candidate below it has probability skill over that naive benchmark. Exact ties (1e-12) would go to larger min.node.size, then smaller mtry, then bootstrap.

## Stage 2: selected specification and cutoff (validation, 10 bp)

| Horizon | Config | mtry | min_node_size | Sampling | Validation_Brier | Default_Brier | Prevalence_Brier | Cutoff | Tuned_validation_Sharpe_matched | Fixed_050_validation_Sharpe_matched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | rf_mtry11_node10_sub50 | 11 | 10 | 50% subsample | 0.2423 | 0.2439 | 0.2406 | 0.65 | 0.689 | 0.382 |
| 3 | rf_mtry04_node10_boot |  4 | 10 | bootstrap | 0.2277 | 0.2305 | 0.2348 | 0.70 | 0.630 | 0.256 |
| 6 | rf_mtry04_node10_sub50 |  4 | 10 | 50% subsample | 0.2189 | 0.2279 | 0.2216 | 0.80 | 0.881 | 0.269 |

Cutoff grid (random forest only): 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90. The baseline searches 0.00-1.00 in 0.01 steps. These validation Sharpe ratios are selection scores, not independent performance estimates.

## Forecast quality (eligible completed labels)

| Period | Horizon | n | Event_rate | Mean_p | Brier | Brier_prevalence | Brier_skill | Log_loss | AUC | Accuracy_tuned | Accuracy_050 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| validation | 1 |  95 | 0.600 | 0.609 | 0.2423 | 0.2406 | -0.007 | 0.6772 | 0.516 | 0.537 | 0.600 |
| validation | 3 |  93 | 0.634 | 0.613 | 0.2277 | 0.2348 | 0.030 | 0.6470 | 0.570 | 0.419 | 0.634 |
| validation | 6 |  90 | 0.678 | 0.672 | 0.2189 | 0.2216 | 0.012 | 0.6274 | 0.527 | 0.367 | 0.678 |
| test | 1 | 155 | 0.677 | 0.644 | 0.2256 | 0.2242 | -0.006 | 0.6453 | 0.507 | 0.523 | 0.671 |
| test | 3 | 153 | 0.752 | 0.678 | 0.1950 | 0.2040 | 0.044 | 0.5794 | 0.546 | 0.477 | 0.739 |
| test | 6 | 150 | 0.773 | 0.741 | 0.1933 | 0.1856 | -0.041 | 0.5906 | 0.441 | 0.380 | 0.767 |

Brier_skill = 1 - Brier / Brier_prevalence; positive means better than the expanding training base rate.

## Test performance at 10 bp

| Horizon | Strategy | Cutoff_rule | Cutoff | Start | Months | CAGR | Sharpe_common_4week | Sharpe_matched | Max_drawdown | Avg_stock_weight | Annual_turnover |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | random_forest | validation_tuned | 0.65 | 2013-02 | 155 | 7.34% | 0.650 | 0.650 | 14.17% | 0.54 | 3.40 |
| 1 | random_forest | fixed_050 | 0.50 | 2013-02 | 155 | 12.31% | 0.823 | 0.823 | 21.00% | 0.94 | 0.62 |
| 1 | all_stocks | not_applicable | - | 2013-02 | 155 | 12.80% | 0.804 | 0.804 | 24.68% | 1.00 | 0.00 |
| 1 | all_bills | not_applicable | - | 2013-02 | 155 | 1.59% | N/A | N/A | 0.00% | 0.00 | 0.00 |
| 3 | random_forest | validation_tuned | 0.70 | 2013-04 | 153 | 5.80% | 0.546 | 0.537 | 14.17% | 0.43 | 0.86 |
| 3 | random_forest | fixed_050 | 0.50 | 2013-04 | 153 | 12.22% | 0.793 | 0.787 | 21.00% | 0.95 | 0.22 |
| 3 | all_stocks | not_applicable | - | 2013-04 | 153 | 12.59% | 0.786 | 0.781 | 24.68% | 1.00 | 0.00 |
| 3 | all_bills | not_applicable | - | 2013-04 | 153 | 1.69% | 1.860 | N/A | 0.00% | 0.00 | 0.00 |
| 6 | random_forest | validation_tuned | 0.80 | 2013-07 | 150 | 3.88% | 0.381 | 0.361 | 11.78% | 0.33 | 0.35 |
| 6 | random_forest | fixed_050 | 0.50 | 2013-07 | 150 | 12.50% | 0.781 | 0.772 | 23.90% | 0.98 | 0.08 |
| 6 | all_stocks | not_applicable | - | 2013-07 | 150 | 12.69% | 0.785 | 0.776 | 24.68% | 1.00 | 0.00 |
| 6 | all_bills | not_applicable | - | 2013-07 | 150 | 1.78% | 2.165 | N/A | 0.00% | 0.00 | 0.00 |

Sharpe definitions, the overlapping-cohort engine, cost convention and benchmark construction are identical to the baseline (shared `portfolio.R`). No model or cutoff rule is selected using this table.

## What the forest relies on (test period, out-of-sample grouped permutation)

| Horizon | Reliance_rank | Variable | Mean_abs_prob_change | Delta_Brier | Signal_flip_rate |
| --- | --- | --- | --- | --- | --- |
| 1 | 1 | unrate_rel36 | 0.0374 | +0.0046 | 0.223 |
| 1 | 2 | cpi_infl_change1 | 0.0183 | +0.0006 | 0.083 |
| 1 | 3 | log_dispersion_change1 | 0.0141 | +0.0002 | 0.068 |
| 1 | 4 | baa10ym_rel36 | 0.0131 | +0.0011 | 0.046 |
| 1 | 5 | mom6 | 0.0124 | +0.0001 | 0.077 |
| 3 | 1 | unrate_rel36 | 0.0494 | +0.0077 | 0.260 |
| 3 | 2 | mom12 | 0.0224 | +0.0009 | 0.060 |
| 3 | 3 | mom6 | 0.0157 | -0.0033 | 0.065 |
| 3 | 4 | cpi_infl_change1 | 0.0134 | +0.0016 | 0.065 |
| 3 | 5 | bill_yld_rel36 | 0.0133 | +0.0004 | 0.054 |
| 6 | 1 | unrate_rel36 | 0.0502 | +0.0002 | 0.182 |
| 6 | 2 | bill_yld_rel36 | 0.0172 | -0.0007 | 0.047 |
| 6 | 3 | mom12 | 0.0154 | +0.0000 | 0.046 |
| 6 | 4 | ma_part | 0.0148 | -0.0024 | 0.067 |
| 6 | 5 | log_rv | 0.0145 | -0.0001 | 0.042 |

Reliance = mean absolute change in the forecast probability when all six lags of a variable are permuted across test origins. Delta_Brier > 0 means permuting the variable worsens test Brier (it carried useful information); <= 0 means it did not help out of sample. Signal_flip_rate is the share of origins whose tuned-cutoff decision changes. These are descriptive, not causal, and correlated predictors share credit.

## Files

- performance.csv / common_window_performance.csv / monthly_portfolios.csv / benchmark_returns.csv: same schema as the baseline.
- validation_forecasts.csv / test_forecasts.csv: one probability per origin with tuned and 0.50 cutoffs/signals, fitting audit, config_id, seed and training prevalence.
- validation_cutoff_grid.csv / selected_cutoffs.csv: full cutoff search and frozen choices.
- forecast_metrics.csv / calibration.csv: probability and classification diagnostics.
- validation_hyperparameter_scores.csv / validation_candidate_forecasts.csv: every grid candidate's Brier scores and expanding validation forecasts.
- refit_audit.csv: every validation (all candidates) and test fit with training endpoint, sample size, seed and hyperparameters.
- importance_out_of_sample.csv / importance_oob_final_fits.csv / partial_dependence.csv: interpretation outputs.
- feature_names.csv, frozen_specification.rds, validation_state.rds, validation_models.rds, final_models.rds, manifest.json, input_hashes.csv, package_versions.csv, session_info.txt, tests.csv.
- Figures: test_wealth.png, validation_hyperparameters.png, validation_cutoffs.png, probabilities.png, importance.png, partial_dependence.png.

19 focused test groups passed; parallel/serial spot checks, final-origin refit reproducibility, refit-count, maturity and input-hash checks all passed.
