# Random forest results - supervisor protocol

Configuration: `rf_v2_supervisor_2026_11_05`. A ranger probability forest on all 132 predictors (22 variables x lags 0-5) is fitted separately at h = 1, 3 and 6. Every validation and test origin is refitted monthly on an expanding sample of labels whose outcomes ended by that origin. **No hyperparameters are tuned**: every fit uses ranger's defaults (mtry 11, min.node.size 10, bootstrap, 500 trees). The only validation choice is the cutoff, with the highest validation matched-ladder Sharpe at 10 bp. The fixed 0.50 cutoff is reported alongside. All-stock and all-bill benchmarks use the same scored months.

Generated automatically by `Scripts/RF/run_rf.R`. See `INTERPRETATION.md` in this folder for the written interpretation.

## Fixed specification (ranger defaults) and cutoff choice (validation, 10 bp)

| Horizon | Config | mtry | min_node_size | Sampling | Validation_Brier | Default_Brier | Prevalence_Brier | Cutoff | Tuned_validation_Sharpe_matched | Fixed_050_validation_Sharpe_matched |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | rf_mtry11_node10_boot | 11 | 10 | bootstrap | 0.2459 | 0.2459 | 0.2406 | 0.75 | 0.593 | 0.332 |
| 3 | rf_mtry11_node10_boot | 11 | 10 | bootstrap | 0.2301 | 0.2301 | 0.2348 | 0.60 | 0.470 | 0.367 |
| 6 | rf_mtry11_node10_boot | 11 | 10 | bootstrap | 0.2276 | 0.2276 | 0.2216 | 0.80 | 0.972 | 0.269 |

Cutoff grid (random forest only): 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90. The baseline searches 0.00-1.00 in 0.01 steps. These validation Sharpe ratios are selection scores, not independent performance estimates.

## Forecast quality (eligible completed labels)

| Period | Horizon | n | Event_rate | Mean_p | Brier | Brier_prevalence | Brier_skill | Log_loss | AUC | Accuracy_tuned | Accuracy_050 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| validation | 1 |  95 | 0.600 | 0.613 | 0.2459 | 0.2406 | -0.022 | 0.6848 | 0.488 | 0.432 | 0.600 |
| validation | 3 |  93 | 0.634 | 0.616 | 0.2301 | 0.2348 | 0.020 | 0.6519 | 0.528 | 0.538 | 0.667 |
| validation | 6 |  90 | 0.678 | 0.671 | 0.2276 | 0.2216 | -0.027 | 0.6468 | 0.450 | 0.400 | 0.678 |
| test | 1 | 155 | 0.677 | 0.644 | 0.2264 | 0.2242 | -0.010 | 0.6466 | 0.502 | 0.342 | 0.658 |
| test | 3 | 153 | 0.752 | 0.685 | 0.1960 | 0.2040 | 0.039 | 0.5823 | 0.558 | 0.673 | 0.752 |
| test | 6 | 150 | 0.773 | 0.743 | 0.2019 | 0.1856 | -0.087 | 0.6327 | 0.450 | 0.393 | 0.767 |

Brier_skill = 1 - Brier / Brier_prevalence; positive means better than the expanding training base rate.

## Test performance at 10 bp

| Horizon | Strategy | Cutoff_rule | Cutoff | Start | Months | CAGR | Sharpe_common_4week | Sharpe_matched | Max_drawdown | Avg_stock_weight | Annual_turnover |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | random_forest | validation_tuned | 0.75 | 2013-02 | 155 | 0.49% | -0.264 | -0.264 | 15.83% | 0.08 | 1.55 |
| 1 | random_forest | fixed_050 | 0.50 | 2013-02 | 155 | 11.81% | 0.789 | 0.789 | 21.00% | 0.94 | 0.62 |
| 1 | all_stocks | not_applicable | - | 2013-02 | 155 | 12.80% | 0.804 | 0.804 | 24.68% | 1.00 | 0.00 |
| 1 | all_bills | not_applicable | - | 2013-02 | 155 | 1.59% | N/A | N/A | 0.00% | 0.00 | 0.00 |
| 3 | random_forest | validation_tuned | 0.60 | 2013-04 | 153 | 10.93% | 0.791 | 0.784 | 21.00% | 0.76 | 0.95 |
| 3 | random_forest | fixed_050 | 0.50 | 2013-04 | 153 | 12.54% | 0.841 | 0.836 | 21.00% | 0.92 | 0.37 |
| 3 | all_stocks | not_applicable | - | 2013-04 | 153 | 12.59% | 0.786 | 0.781 | 24.68% | 1.00 | 0.00 |
| 3 | all_bills | not_applicable | - | 2013-04 | 153 | 1.69% | 1.860 | N/A | 0.00% | 0.00 | 0.00 |
| 6 | random_forest | validation_tuned | 0.80 | 2013-07 | 150 | 3.84% | 0.344 | 0.326 | 14.17% | 0.35 | 0.35 |
| 6 | random_forest | fixed_050 | 0.50 | 2013-07 | 150 | 12.43% | 0.780 | 0.772 | 22.83% | 0.97 | 0.14 |
| 6 | all_stocks | not_applicable | - | 2013-07 | 150 | 12.69% | 0.785 | 0.776 | 24.68% | 1.00 | 0.00 |
| 6 | all_bills | not_applicable | - | 2013-07 | 150 | 1.78% | 2.165 | N/A | 0.00% | 0.00 | 0.00 |

Sharpe definitions, the overlapping-cohort engine, cost convention and benchmark construction are identical to the baseline (shared `portfolio.R`). No model or cutoff rule is selected using this table.

## What the forest relies on (test period, out-of-sample grouped permutation)

| Horizon | Reliance_rank | Variable | Mean_abs_prob_change | Delta_Brier | Signal_flip_rate |
| --- | --- | --- | --- | --- | --- |
| 1 | 1 | unrate_rel36 | 0.0402 | +0.0061 | 0.077 |
| 1 | 2 | cpi_infl_change1 | 0.0207 | +0.0003 | 0.041 |
| 1 | 3 | log_dispersion_change1 | 0.0153 | -0.0002 | 0.044 |
| 1 | 4 | mom6 | 0.0142 | -0.0010 | 0.049 |
| 1 | 5 | tail_imb_change1 | 0.0141 | -0.0015 | 0.037 |
| 3 | 1 | unrate_rel36 | 0.0693 | +0.0115 | 0.162 |
| 3 | 2 | mom12 | 0.0275 | -0.0003 | 0.078 |
| 3 | 3 | mom6 | 0.0178 | -0.0031 | 0.049 |
| 3 | 4 | cpi_infl_change1 | 0.0168 | +0.0018 | 0.049 |
| 3 | 5 | term_spread_change1 | 0.0150 | +0.0011 | 0.044 |
| 6 | 1 | unrate_rel36 | 0.0780 | +0.0016 | 0.203 |
| 6 | 2 | cpi_infl_change1 | 0.0223 | +0.0016 | 0.065 |
| 6 | 3 | mom12 | 0.0200 | -0.0018 | 0.022 |
| 6 | 4 | ma_part | 0.0199 | -0.0035 | 0.056 |
| 6 | 5 | bill_yld_rel36 | 0.0171 | -0.0010 | 0.017 |

Reliance = mean absolute change in the forecast probability when all six lags of a variable are permuted across test origins. Delta_Brier > 0 means permuting the variable worsens test Brier (it carried useful information); <= 0 means it did not help out of sample. Signal_flip_rate is the share of origins whose tuned-cutoff decision changes. These are descriptive, not causal, and correlated predictors share credit.

## Files

- performance.csv / common_window_performance.csv / monthly_portfolios.csv / benchmark_returns.csv: same schema as the baseline.
- validation_forecasts.csv / test_forecasts.csv: one probability per origin with tuned and 0.50 cutoffs/signals, fitting audit, config_id, seed and training prevalence.
- validation_cutoff_grid.csv / selected_cutoffs.csv: full cutoff search and frozen choices.
- forecast_metrics.csv / calibration.csv: probability and classification diagnostics.
- validation_hyperparameter_scores.csv / validation_candidate_forecasts.csv: Brier scores and validation forecasts of the single default specification (no grid search).
- refit_audit.csv: every validation (all candidates) and test fit with training endpoint, sample size, seed and hyperparameters.
- importance_out_of_sample.csv / importance_oob_final_fits.csv / partial_dependence.csv: interpretation outputs.
- feature_names.csv, frozen_specification.rds, validation_state.rds, validation_models.rds, final_models.rds, manifest.json, input_hashes.csv, package_versions.csv, session_info.txt, tests.csv.
- Figures: test_wealth.png, validation_cutoffs.png, probabilities.png, importance.png, partial_dependence.png.

19 focused test groups passed; parallel/serial spot checks, final-origin refit reproducibility, refit-count, maturity and input-hash checks all passed.
