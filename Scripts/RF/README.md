# Random forest in R

This folder implements the random forest under the supervisor revision of 5 November 2026 ([GROUP_MODEL_STANDARDS_DRAFT.md](../Models/GROUP_MODEL_STANDARDS_DRAFT.md)). It sits beside the baseline ([`../Baseline models/`](../Baseline%20models/)) and mirrors it: same data, splits, maturity rule, monthly expanding refits, cutoff rules, Sharpe references, portfolio engine and output schema. Four things differ:

- the estimator (a ranger probability forest);
- the predictor set (all 132 columns);
- a **declared 18-candidate hyperparameter grid**, chosen by validation Brier (the baselines have no hyperparameters);
- a **coarser cutoff grid**: 0.00-1.00 in steps of 0.05, against the baseline's 0.00-1.00 in steps of 0.01.

Hyperparameters and cutoff are chosen in sequence on the validation period, as group-standard rules 24-28 describe.

## Run

From the project root:

```sh
Rscript --vanilla Scripts/RF/run_rf.R
```

Optional stages, as in the baseline:

```sh
Rscript --vanilla Scripts/RF/run_rf.R --stage=tests       # checks only
Rscript --vanilla Scripts/RF/run_rf.R --stage=validation  # choose hyperparameters + cutoff and freeze; no test forecasts
Rscript --vanilla Scripts/RF/run_rf.R --stage=test        # requires identical config, code and input hashes
```

Dependencies: `ranger`, `digest`, `jsonlite` (in `.R-library` or the user library). `EC4308_RF_WORKERS` sets the number of parallel worker processes (default: cores - 2; `1` = serial; results are identical either way). Runtime is about 10 minutes with 14 workers, almost all of it the 5,184 validation fits of the grid search. Outputs are written to `Results/RF/` and overwritten on a rerun.

## Files

| File | Responsibility |
| --- | --- |
| `config.R` | Frozen choices: hyperparameter grid (`rf_hyper_grid()`), trees, RF cutoff grid, seed rule, dates, costs |
| `models.R` | Load validated RDS inputs, mature-row rule, ranger fit/predict, per-origin task, parallel runner, Brier scoring, candidate choice (`rf_choose_candidate`), permutation interpretation |
| `tests.R` | Shared-engine identities, class ordering, the hyperparameter grid and that ranger receives each candidate's settings, Brier and cutoff tie rules, cutoff grid, feature set, maturity and future-perturbation invariance, seeded determinism |
| `report.R` | `RESULTS.md` and figures |
| `run_rf.R` | Orchestration: tests, stage 1 (hyperparameters by Brier) and stage 2 (cutoff by Sharpe) on validation, freeze, test forecasting, portfolios, interpretation, manifest |
| `../Baseline models/config.R`, `portfolio.R` | Shared helpers and the **single shared portfolio evaluator** (sourced from the sibling folder, not copied). Moving either folder breaks this link. |

## How the code runs, step by step

### 0. Setup and tests
`run_rf.R` locates the project root (walking up until it finds `Data/` and `Scripts/`), sources the shared baseline helpers/engine and the RF files, hashes all inputs and code, then runs `run_rf_tests()`. A failed test stops the run before any fitting.

### 1. Data
For each horizon `rf_load_data()` reads `Data/processed/model_ready_h{h}.rds` and uses `$X` (621 origins x 132 columns: 22 variables x lags 0-5). It asserts the columns equal the frozen `feature_names`, none are metadata/outcome columns, the calendar is complete, every label's `outcome_end` is exactly `origin + h`, and the initial train/validation index sets respect their boundaries. No scaling, imputation or feature selection is applied: trees are scale-invariant and the inputs are complete.

### 2. The model at one origin
`rf_task()` handles one forecast origin t for one hyperparameter candidate:

1. **Training rows** = every row with an observed label whose outcome ended by t (`outcome_end <= t`). Purged boundary rows re-enter once mature. At t = 2013-01 the latest training origins are 2012-12 / 2012-10 / 2012-07 for h = 1/3/6.
2. Fit `ranger(probability = TRUE, num.trees = 1000, mtry, min.node.size, replace, sample.fraction, splitrule = "gini")` with the candidate's settings. y is a factor with levels `c("0","1")`; the fit is single-threaded with seed `4308 + 10*month_number(t) + h`. Each fit asserts that ranger actually used the requested `mtry`, `min.node.size` and tree count.
3. Predict `P(Y = 1)` (column `"1"`, checked by name) for row t only.
4. Optionally predict permuted copies of row t (interpretation only; see step 6).

The seed depends only on (h, origin), so all candidates share common random numbers (their Brier differences reflect settings, not luck) and the parallel schedule cannot change any result.

### 3. Validation: stage 1 (hyperparameters), then stage 2 (cutoff)
The declared grid (`rf_hyper_grid()` in `config.R`), fixed before any evaluation:

| Argument | Values |
| --- | --- |
| `mtry` | 4, 11 (ranger default), 33 |
| `min.node.size` | 10 (ranger default), 30, 60 |
| sampling | bootstrap with replacement (`sample.fraction` 1, default) or 50% subsample without replacement |
| `num.trees` | 1000, fixed (not tuned): more trees only reduce Monte Carlo noise |

**Stage 1 (rules 25-26).** For January 2005 to December 2012, every one of the 18 candidates is refit monthly on the expanding mature sample (18 x 96 x 3 = 5,184 fits). Each candidate's 96 forecasts are scored by Brier on eligible completed labels (95/93/90 rows). Per horizon, the lowest Brier wins; exact ties (1e-12) go to the simpler candidate (larger `min.node.size`, then smaller `mtry`, then bootstrap). OOB error and random CV are never used.

**Stage 2 (rules 27-28).** With the winner's **saved** validation forecasts held fixed, each of the 21 cutoffs 0.00, 0.05, ..., 1.00 is run through the shared portfolio engine at 10 bp. The cutoff that maximizes `sharpe_matched_ladder` is chosen; undefined Sharpes are excluded, and ties go to the cutoff nearest 0.50, then the lower one. The fixed 0.50 comparator is recorded too. The coarser cutoff grid is a declared RF-specific departure from the protocol's 0.01 grid (rule 28). Searching 13 cutoffs instead of 101 on only 96 validation months limits how finely the cutoff can be fitted to individual validation episodes.

The selected hyperparameters, cutoffs, config, input hashes and code hashes are saved to `frozen_specification.rds` **before** any test origin is forecast.

### 4. Test
For January 2013 to December 2025 (156 origins), the selected candidate is refit monthly on the expanding mature sample (468 fits). Hyperparameters and both cutoffs stay frozen. After each batch, the first and last origins are refit serially in the main process and must match the parallel output bit for bit.

### 5. Portfolios and metrics
Identical to the baseline, through the same `portfolio.R`: h overlapping cohorts with equal post-cost capital, monthly rebalancing with the exact cost equation, vintage bill accrual, the first h-1 forecasts as warm-up, scored months Feb/Apr/Jul of 2005 and 2013 to December, 0/10/20 bp with signals fixed at the 10 bp choices, both Sharpe references, plus all-stock and all-bill benchmarks. `forecast_metrics.csv` adds the Brier of the expanding training base rate and a Brier skill score against it.

### 6. Interpretation (never used for any choice)
- **Out-of-sample grouped permutation** (`importance_out_of_sample.csv`): at every validation/test origin, the fitted forest also predicts copies of that origin's row with all six lags of one variable replaced by another origin's values (5 seeded permutations per period). Reports the change in Brier on eligible labels, the mean absolute change in probability (reliance) and the share of tuned-cutoff signals that would flip.
- **OOB permutation importance** of the final validation (2012-12) and final test (2025-12) fits (`importance_oob_final_fits.csv`): describes what the forest learned in-sample; overlapping labels make OOB levels optimistic, so read it as a ranking only.
- **Partial dependence** (`partial_dependence.csv`, `partial_dependence.png`) of the final test fit for the four variables the forest relied on most in the test period, with all six lags set to the same value (solid) and lag0 only (dashed).

### 7. Outputs (`Results/RF/`)
Same as the baseline (`performance.csv`, `common_window_performance.csv`, `monthly_portfolios.csv`, `benchmark_returns.csv`, `validation_forecasts.csv`, `test_forecasts.csv`, `validation_cutoff_grid.csv`, `selected_cutoffs.csv`, `forecast_metrics.csv`, `calibration.csv`, `refit_audit.csv`, `frozen_specification.rds`, `validation_models.rds`, `final_models.rds`, `manifest.json`, `input_hashes.csv`, `tests.csv`, `session_info.txt`, `package_versions.csv`, `RESULTS.md`, `test_wealth.png`, `validation_cutoffs.png`), plus the RF-specific `validation_hyperparameter_scores.csv`, `validation_candidate_forecasts.csv`, `validation_hyperparameters.png`, `feature_names.csv`, the interpretation files and figures above, and the written `INTERPRETATION.md`. The forecast files carry every column required by rule 66, plus `mtry`, `min_node_size`, `replace`, `sample_fraction`, `num_trees`, `seed` and `train_prevalence`.

## Leakage safeguards

| Risk | Safeguard |
| --- | --- |
| Metadata or outcomes as inputs | `$X` only; asserted equal to the 132 frozen names with no outcome/split columns |
| Immature labels in a fit | `outcome_end <= origin` per fit; asserted per task, per assembled forecast set and in the final audit; future-perturbation test at 2005-01 and 2013-01 |
| Test data in tuning | Validation stage forecasts only 2005-2012 origins and uses only 2005-2012 returns; choices are frozen to disk first; `--stage=test` refuses to run if config, code or inputs changed |
| Optimistic in-sample scores | Finite pre-declared grid scored by chronological expanding validation Brier; no OOB score, no random CV; validation scores reported as selection scores, not performance estimates |
| Non-reproducibility | Origin-specific seeds, single-threaded fits, serial-vs-parallel spot checks, final-origin refit reproduction, input-hash check at the end |
