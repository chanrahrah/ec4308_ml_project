# Baseline models in R

This folder implements the supervisor revision of 5 November 2026: five-variable **probit and logit**, plus **all stocks** and **all matched Treasury bills**, at **h = 1, 3 and 6 months**. Both validation and test estimation expand monthly. Report both **validation-tuned and fixed-0.50 cutoffs**, and both **common 4-week and matched-ladder Sharpe ratios**. It uses the existing validated data without rebuilding it. Primary costs are 10 bp, with 0 and 20 bp sensitivity. The current shared coding prompt is [GROUP_MODEL_STANDARDS_DRAFT.md](../../GROUP_MODEL_STANDARDS_DRAFT.md).

The old fixed-fit validation results and source are preserved in `Results/archive/baseline_v1_pre_supervisor_2026-11-05`. This is a supervisor-directed revision after the original test results were inspected, not a new untouched evaluation sample.

## Run

From the project root:

```sh
Rscript --vanilla 'Scripts/baseline models/run_baselines.R'
```

It also works from any directory with the absolute script path. Dependencies are base R, `digest` and `jsonlite`, already used by the data pipeline. If needed, run `Rscript --vanilla Scripts/00_setup.R` once.

Outputs go to `Results/baseline models/`. Start with `RESULTS.md`, `performance.csv` and `test_wealth.png`. Output files are overwritten on a rerun. No random resampling or optimization starts are used, so the numerical results are deterministic with the same R environment.

Optional stages:

```sh
Rscript --vanilla 'Scripts/baseline models/run_baselines.R' --stage=tests
Rscript --vanilla 'Scripts/baseline models/run_baselines.R' --stage=validation
Rscript --vanilla 'Scripts/baseline models/run_baselines.R' --stage=test
```

The separate test stage requires exactly matching configuration, R code hashes and input hashes from validation. `--stage=all` (the default) runs validation, freezes the choices, then runs testing. A rerun is a reproducibility check, not another opportunity to choose a specification from test results.

## Lag decision and literature

**Use one observation of each variable, measured at forecast origin t, and no additional distributed lags.** In this dataset these are `_lag0`. They precede every return in the target window t+1,...,t+h. A paper written as `Y_t ~ X_(t-1)` uses the same timing for a one-step forecast; translating that notation into `_lag1` here would impose an additional month of delay.

The baseline is `P(Y_(t,h)=1 | I_t) = F(a_h + b_h' X_t)`, with F the normal CDF for probit or logistic CDF for logit. Each horizon gets its own coefficients and cutoff. There are five slopes plus an intercept. No penalty, feature search or lag tuning is performed.

The evidence and the design judgment are distinct:

1. [Xie, Sun and Fan (2023), Return direction forecasting, equation 27](https://link.springer.com/article/10.1186/s40854-023-00489-z) present a comparison probit with the preceding observation of the predictor vector and preceding return. This supports a compact origin-information benchmark, but is not a study proving our exact five-variable specification or its optimal lag count.
2. [Hjalmarsson (2008), Predicting Global Stock Returns](https://www.federalreserve.gov/pubs/ifdp/2008/933/ifdp933.htm) studies predictive regressions using lagged financial predictors, including short interest rates and term spreads, and expanding estimation samples. It supports the predictor timing and recursive evaluation convention, rather than a six-lag rule.
3. [Campbell and Thompson (2008), Predicting Excess Stock Returns Out of Sample](https://dash.harvard.edu/bitstream/handle/1/2622619/Campbell_Predicting.pdf) evaluate predictor regressions out of sample and discuss the difficulty of estimating expected returns precisely. This motivates a restrained baseline; it does not establish five as an optimal feature count.
4. [Nyberg (2011), Forecasting the direction of the US stock market with dynamic binary probit models](https://doi.org/10.1016/j.ijforecast.2010.02.008) motivates binary directional models and explicit dynamics as an extension. Its published abstract/preview was consulted; the full article was not available in this session. Our current five-input model is **static**, unlike the earlier proposal's dynamic probit. An additional lagged target would change the agreed predictor set and is not silently added.

My design choice is therefore to hold lag0 fixed across all horizons. `mom6` already summarizes six months, the bill-yield deviation uses a trailing 36-month reference, and the differenced variables already use adjacent months. Extra lags would add highly related inputs and additional specification choices to a small benchmark. Neither the forecast horizon nor the number of overlapping portfolio cohorts determines the number of predictor lags. The literature does not supply a universal best lag count for this exact combination; no claim is made that lag0 is empirically optimal.

## Mapping the five requested variables

| Requested concept | Existing validated model column | Meaning |
| --- | --- | --- |
| Bill yield matched | `bill_yld_rel36_lag0` | Matched effective annual yield less its prior-36-month mean; 4/13/26-week yield for h=1/3/6 |
| Term spread | `term_spread_change1_lag0` | Monthly change in 10-year minus 13-week yield |
| Log realized volatility | `log_rv_lag0` | Existing log monthly realized-volatility measure |
| Six-month momentum | `mom6_lag0` | Existing compounded six-month equity return |
| Dispersion | `log_dispersion_change1_lag0` | Monthly change in log cross-sectional dispersion |

These transformed definitions preserve the earlier stationarity work, rather than returning to raw nonstationary levels. See `../VARIABLES.md` for complete definitions. Only these columns of `$X` enter fitting; future-return metadata never enters a design matrix.

## Splits, tuning and expanding refits

- Initial training boundary: April 1974-December 2004, with provided h-specific purged rows (368/366/363 labels). These are the initial-boundary counts, not a fixed sample for all validation fits.
- Validation: January 2005-December 2012 origins; coefficients and scaler refit every month from April 1974 using only labels whose outcomes end by that origin. Previously boundary-purged rows enter once mature. There are 96 refits per model/horizon (576 total); completed-label metrics use the provided 95/93/90 eligible rows.
- Test: January 2013-December 2025 origins; every month refit from April 1974 using only labels whose outcome ends by that forecast origin. This is 156 refits per model/horizon, including end-of-sample forecasts with incomplete targets. Classified test labels number 155/153/150.
- Fit centering and sample standard deviations on the exact fitting rows each time. No randomized folds or future imputation are used.
- Search cutoffs 0.00,...,1.00 in 0.01 steps using **validation monthly net Sharpe against the horizon-matched bill ladder at 10 bp**. `p >= cutoff` means stocks. Common 4-week Sharpe is reported but does not select cutoffs. Undefined candidates are excluded; if all are undefined, stop. Ties within 1e-12 choose closest to 0.50, then the lower cutoff. Brier is a loss to minimize, not maximize. These fixed baseline specifications do not have a hyperparameter search.
- Freeze all six tuned cutoffs before forecasting any test origin. Also report the fixed 0.50 cutoff using the same probabilities and portfolio engine. Report both links and both cutoff rules; do not choose a winner from test performance. Keep each rule's signals fixed across 0/10/20 bp costs.
- A forward label crossing a split boundary is excluded from that split's classification evaluation and initial fitting. It can enter a later expanding fit once it has actually matured. Monthly portfolio accounting only uses returns within its own split; late-validation signals do not bring any test returns into cutoff selection.

The supervisor revision supersedes v1's fixed-fit validation, the older markdown's single test fit and the proposal's nonoverlapping trading schedule. The underlying dates, predictors and targets are retained.

## Overlapping portfolios and costs

Reference: `Reference codes/EC4308_6step_strategy_12 months.R`. At the return month t+1, active signals are t-h+1,...,t. Each receives 1/h of **post-cost** wealth, giving stock weight equal to the mean of those h binary signals. The first h-1 origins form the initial portfolio; valuation begins with the first complete group, endowed with wealth 1 and no entry fee, as in the professor's initialized example. Initial cohort labels are available historical forecasts, not future outcomes.

Scored validation months start February/April/July 2005; test months start February/April/July 2013. All end in December of their respective split. Benchmarks use identical months. `common_window_performance.csv` additionally compares all horizons on July 2005-December 2012 and July 2013-December 2025. Historical cohort yields are valued under the accepted accrual proxy at initialization; warm-up gains/losses are outside the scored wealth path.

Stocks earn CRSP VW total returns. Each bill cohort accrues its origin's matched daily log yield over the actual difference between monthly CRSP quote dates: `expm1(origin_yield * days)`. Thus the 3/6-month all-bills comparator is a ladder of locked vintage yields, not the current monthly yield multiplied by h. Equalizing the bill slots monthly assumes frictionless transfers under this accrual proxy; interest-rate mark-to-market changes are not modeled.

At a rebalance, V is wealth, S is existing stock value, w the new target equity fraction, c the one-way cost. Solve `q = w*(V-c*abs(q))-S`; cost is `c*abs(q)`. Purchases use `(w*V-S)/(1+c*w)`; sales use `(w*V-S)/(1-c*w)`. This captures drift as well as signal changes. Opposite cohort equity trades are netted. There is no separate bill-leg fee, no bill rollover fee, and no final liquidation fee. These conventions also apply to both passive benchmarks.

Two Sharpe measures use exactly the same monthly portfolio observations:

- `sharpe_common_rf`: subtract `expm1(TMYLD_4week_at_previous_month_end * monthly_quote_days)`. The common reference is identical for a given month across models and horizons.
- `sharpe_matched_ladder`: subtract the monthly all-bills ladder return described above. This selects the individual-model cutoff. `net_Sharpe` and `bill_reference_return` are compatibility aliases for the matched statistic and matched monthly return.

Both use `sqrt(12) * mean(excess) / sd(excess)`. Neither subtracts a full 13- or 26-week holding-period return from a monthly portfolio return. No maturity conversion bug existed in v1: its matched ladder already used monthly quote gaps. `benchmark_returns.csv` exposes both references, dates, daily common yield and security IDs for verification. These are the accepted accrual proxies, not marked-to-market Treasury total returns.

All-bills excess volatility against its identical matched ladder is zero, so matched Sharpe is N/A. Against the common 4-week reference, h=3/6 bill ladders can have nonzero spread returns and a defined ratio; do not mistake low proxy spread volatility for a riskless investment opportunity. Sortino, beta and certainty-equivalent excess returns retain the matched reference, stated in `excess_metric_reference`.

Maximum drawdown and Ulcer use monthly wealth observations and are not annualized. Intramonth drawdowns are unobserved. Sharpe and volatility are conventionally annualized from monthly returns. Higher-h forward labels overlap, but no h-month cumulative label return is mistaken for one monthly portfolio return.

## Output schema and comparisons

`validation_forecasts.csv` and `test_forecasts.csv` retain one probability per model/horizon/origin. `cutoff`/`signal` identify the validation-tuned rule; `cutoff_050`/`signal_050` identify the fixed comparator. Every row includes fitting information date, latest training origin/outcome, sample size and configuration ID. The last validation-origin fits are stored in `validation_models.rds`; per-origin parameters are in the coefficient/scaler audit files.

`performance.csv`, `common_window_performance.csv` and `monthly_portfolios.csv` include `cutoff_rule` (`validation_tuned`, `fixed_050`, or `not_applicable` for passive benchmarks). Filter by this column as well as model, horizon and cost; do not accidentally stack the two probability-cutoff paths into one portfolio. Both Sharpe columns are present in performance tables. Probability calibration is exported once per model/horizon/period because it does not depend on cutoff.

This is a monthly research backtest with idealized close-to-close execution after month-end information. Next executable prices and exact historical Treasury prices are unavailable; these execution and valuation limitations remain explicit. No causality guarantee is claimed for source fields beyond the existing validated data assumptions.

## Implementation and checks

| File | Responsibility |
| --- | --- |
| `config.R` | Fixed model, cutoff, cost and date choices |
| `models.R` | Read validated inputs, fit/scale, forecast, audit maturity |
| `portfolio.R` | Cohorts, vintage yields, exact costs and metrics |
| `tests.R` | Reference parity, leakage perturbation, maturity, cash/stocks and cost identities |
| `report.R` | Results summary and figures |
| `run_baselines.R` | Reproducible orchestration, frozen selection, outputs and hashes |

The test suite checks the professor's full 12-month example, exact post-cost weights, drift costs, vintage monthly accrual, target accrual consistency, both reference definitions, same-month common-reference identity across horizons, all-stock/all-bill identities, warm-ups, future-perturbation invariance at validation and test origins, changing monthly validation fits, matched-objective cutoff selection, undefined Sharpe behavior and AUC ties. The pipeline audits all 576 validation and 936 test fits and verifies that input hashes are unchanged. Coefficients are descriptive; naive GLM significance tests are not used for overlapping targets.
