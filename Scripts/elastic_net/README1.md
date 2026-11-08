# EC4308 Elastic Net — same architecture as baseline

## Folder layout

Place this directory at `PROJECT_ROOT/Scripts/elastic_net/`, alongside the existing baseline directory. It **reuses an unmodified copy** of the supplied baseline `portfolio.R` for identical cohort conventions, net trading costs and monthly performance metrics. The original CRSP CSV files and frozen RDS objects must be present in the project root as documented by the baseline.

## Run

From project root:

```sh
Rscript --vanilla Scripts/elastic_net/run_elastic_net.R --stage=tests
Rscript --vanilla Scripts/elastic_net/run_elastic_net.R --stage=validation
Rscript --vanilla Scripts/elastic_net/run_elastic_net.R --stage=test
```

Or run `--stage=all` to complete everything. R packages: `glmnet`, `digest`, `jsonlite`. Full source input data were not included with the uploaded scripts, so no backtest values are supplied here. **R was not available in the authoring environment; run and verify in your local project** before relying on any outputs.

## Hyperparameters

`alpha in {0,.25,.5,.75,1}` and `lambda_ratio in {1,.3,.1,.03,.01}` (25 candidates per horizon). `lambda_ratio` is relative to the **training-only** binomial LASSO `lambda_max` recomputed at each monthly origin. The chosen *policy* is constant, but the absolute `lambda` may differ by month. This is a proposed coarse grid, not an established optimal choice. The same five alpha settings span ridge through LASSO; `lambda_max` is anchored at alpha=1 so it remains well-defined for ridge.

No `cv.glmnet()` is used. Every candidate is evaluated by 96 monthly expanding fits; validation Brier on mature eligible labels selects one pair separately for each horizon. Then 101 cutoff options are run through the baseline's **full** monthly portfolio at 10bp, selecting the maximum matched-ladder Sharpe; the 0.50 reference is always reported. Cutoffs, selected alpha and lambda ratio are frozen before recursive test fits.

## Shared logic and data safety

- Horizon-specific frozen `model_ready_h1/h3/h6.rds` full 132-column `X` matrix.
- One model per horizon, each month's fitting set satisfies `outcome_end <= origin`.
- Training-only standardized predictors; constant columns remain zero after scaling.
- Monthly overlapping cohorts and vintage-bill returns identical to baseline `portfolio.R`.
- Report 4-week common reference and passive matched bill ladder at identical monthly frequency.
- Include all-stock/all-bill, tuned/0.50 cutoff, 0/10/20bp, common-window performance.
- Full test-origin fits stored in `test_fits_by_origin.rds` (potentially large); candidate monthly fits are refitted on demand instead of storing 25 x 96 models.

## Outputs

`Results/elastic_net/`: `validation_tuning.csv`, `selected_settings.csv`, `validation_cutoff_grid.csv`, `validation_forecasts.csv`, `test_forecasts.csv`, `refit_audit.csv`, `monthly_portfolios.csv`, `benchmark_returns.csv`, `performance.csv`, `common_window_performance.csv`, `forecast_metrics.csv`, `frozen_specification.rds`, `validation_state.rds`, `test_fits_by_origin.rds`, hashes, manifest, tests and `RESULTS.md`.

## Caveats

- If the existing dataset has nonfinite inputs, the pipeline halts; it never imputes from the future.
- The `lambda_max` formula assumes the standard mean binomial negative-log-likelihood convention. The fit's effective absolute lambda is always saved in the audit.
- Monthly refitting uses calendar origin information only. Source-vintage and idealized execution limitations from the baseline remain.
- Other models' ensembles, probability calibration and dependence-aware hypothesis tests are not included here.
- The common grid should be set before inspecting the Elastic Net test results. Do not retune from test feedback.
