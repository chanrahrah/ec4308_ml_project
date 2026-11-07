# Random forest: interpretation of results

Numbers come from the run of 8 November 2026 (`Scripts/RF/run_rf.R`, configuration `rf_v2_supervisor_2026_11_05`). Specification:

- **Hyperparameters:** chosen by validation Brier from a fixed grid of 18 candidates, 1,000 trees.
- **Cutoff:** chosen by validation Sharpe from 0.30–0.90 in steps of 0.05.
- **Comparison rule:** the fixed 0.50 cutoff, reported alongside.

Returns and Sharpe ratios are net of 10 bp costs unless stated otherwise. "Sharpe" means the horizon-matched ladder Sharpe; the common 4-week Sharpe is almost identical (see `RESULTS.md`). If code or data change, rerun the script and update this note.

## 1. Bottom line

1. **Tuning improved the probabilities a little, but did not create skill where there was none.**
   - On validation Brier, the selected candidates beat the ranger defaults at every horizon. The margin is small at h = 1 and h = 3 (0.0016 and 0.0028) and larger at h = 6 (0.0090).
   - In the test, Brier improved at every horizon, by 0.0008 / 0.0010 / 0.0086 compared with the earlier default-forest run. The test probabilities come from the same 2013–2025 months.
   - The forest still beats the naive base-rate forecast only at h = 3 (test Brier skill +4.4%). At h = 1 it is level with the base rate (−0.6%). At h = 6 it is still worse (−4.1%, down from −8.7% with defaults).
2. **The validation-tuned cutoffs failed again in the test.** The chosen cutoffs (0.65 / 0.70 / 0.80) each won validation by holding mostly bills through the 2008 crash. In the test they produced Sharpe 0.65 / 0.54 / 0.36 and CAGR 7.3% / 5.8% / 3.9%. All stocks gave Sharpe about 0.78–0.80 and CAGR about 12.7%. Their month-to-month exposure was uncorrelated with, or slightly opposite to, the next month's market return.
3. **The fixed 0.50 cutoff is essentially buy-and-hold and remains the reliable rule.** It stays 94–98% invested and matches all stocks: Sharpe 0.82 vs 0.80 at h = 1, 0.79 vs 0.78 at h = 3, and 0.77 vs 0.78 at h = 6. Maximum drawdown is a few points smaller. This is one test period, with no dependence-robust significance test.
4. **The tuned-cutoff results are fragile.** Across three versions of this forest (defaults with a 0.01 grid, defaults with a 0.1 grid, tuned with a 0.05 grid), the tuned-cutoff test Sharpe ranged from 0.12 to 0.65 at h = 1 and from 0.54 to 0.78 at h = 3. The 0.50 rule stayed between 0.77 and 0.84.
5. **One variable still drives the forest: `unrate_rel36`** (unemployment minus its prior 36-month mean). It ranks first everywhere and is non-linear. Its post-COVID distortion pushed the forest towards bills during the 2021–2023 recovery.

## 2. Specification and selection

| Item | Value |
| --- | --- |
| Inputs | All 132 predictors (22 variables × lags 0–5); no scaling or feature selection |
| Fitting | Refit every month on all labels whose outcome had ended by that origin |
| Hyperparameter grid (fixed in advance) | `mtry` ∈ {4, 11, 33} × `min.node.size` ∈ {10, 30, 60} × sampling ∈ {bootstrap, 50% subsample without replacement}; 1,000 trees |
| Stage 1 | Per horizon, lowest validation Brier (2005–2012, 90–95 completed labels) |
| Stage 2 | With the winner's saved validation forecasts, the cutoff in {0.30, 0.35, …, 0.90} with the highest validation Sharpe at 10 bp |

### Stage 1: which settings won

| h | Selected | Validation Brier | ranger default (11/10/bootstrap) | Base-rate Brier | Runner-up gap |
| --- | --- | --- | --- | --- | --- |
| 1 | `mtry` 11, node 10, 50% subsample | 0.2423 | 0.2439 (rank 7 of 18) | 0.2406 | 0.0006 |
| 3 | `mtry` 4, node 10, bootstrap | 0.2277 | 0.2305 (rank 7) | 0.2348 | 0.0003 |
| 6 | `mtry` 4, node 10, 50% subsample | 0.2189 | 0.2279 (rank 17) | 0.2216 | 0.0003 |

How to read this (see `validation_hyperparameters.png`):

- **Fewer split candidates (`mtry` = 4) and half-subsampling help at the longer horizons.** At h = 6 the five best candidates all use 50% subsampling. At h = 1 and h = 6 the worst candidate is the most flexible one (`mtry` 33, node 10, bootstrap). With 132 mostly noisy inputs, trees that consider fewer variables per split and see different subsamples are less likely to chase noise.
- **Larger terminal nodes did not help.** Every winner uses the smallest node size (10). Smoothing by bigger leaves was not what the data rewarded.
- **The winning margins are tiny.** The gaps to the runner-up (0.0003–0.0006) are well within the noise of 90–95 overlapping labels, so the exact winner is partly luck. The broad pattern (low `mtry` and subsampling beat high `mtry` and bootstrap at h = 6) is more trustworthy than the specific choice.
- **At h = 1 nothing helps.** No candidate beats the base rate on validation. At h = 3, 17 of 18 candidates do. At h = 6, six of 18 do, including the winner.

### Stage 2: validation Sharpe across the cutoff grid (10 bp)

| Cutoff | h = 1 | h = 3 | h = 6 |
| --- | --- | --- | --- |
| 0.30–0.40 | 0.29 | 0.29 | 0.27 |
| 0.50 | 0.38 | 0.26 | 0.27 |
| 0.55 | 0.14 | 0.34 | 0.29 |
| 0.60 | 0.40 | 0.45 | 0.23 |
| 0.65 | **0.69** | 0.47 | −0.26 |
| 0.70 | 0.50 | **0.63** | 0.57 |
| 0.75 | 0.44 | 0.33 | 0.63 |
| 0.80 | 0.35 | n/a | **0.88** |
| All stocks | 0.29 | 0.29 | 0.27 |

How to read this:

- **The curves are jagged.** One step of 0.05 can move the Sharpe by 0.2 to 0.8. At h = 6, for example, it goes from −0.26 at 0.65 to 0.88 at 0.80.
- **Every winner is a low-exposure rule.** The selected cutoffs held only 33% / 6% / 5% stocks on average in validation, and 0–17% from 2008 to 2010. They score well because the validation period contains a 38% crash, not because the probabilities time the market.

## 3. Forecast quality

Eligible completed labels only. Brier skill = 1 − Brier / Brier of the expanding training base rate (positive is better). These numbers do not depend on the cutoff.

| Period | h | Event rate | Mean forecast | Brier | Base-rate Brier | Brier skill | AUC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| validation | 1 | 0.600 | 0.609 | 0.2423 | 0.2406 | −0.7% | 0.52 |
| validation | 3 | 0.634 | 0.613 | 0.2277 | 0.2348 | **+3.0%** | 0.57 |
| validation | 6 | 0.678 | 0.672 | 0.2189 | 0.2216 | +1.2% | 0.53 |
| test | 1 | 0.677 | 0.644 | 0.2256 | 0.2242 | −0.6% | 0.51 |
| test | 3 | 0.752 | 0.678 | 0.1950 | 0.2040 | **+4.4%** | 0.55 |
| test | 6 | 0.773 | 0.741 | 0.1933 | 0.1856 | −4.1% | 0.44 |

How to read this:

- **h = 3 is the only horizon with consistent skill.** It beats the base rate in both validation and test.
- **h = 6 improved most from tuning but did not carry over.** It now beats the base rate in validation (+1.2%), but not in the test (−4.1%), and its test ranking is inverted (AUC 0.44). Calibration did improve: the forecasts above 0.9 fell from 24 to 8, and their hit rate is 75% against a mean forecast of 0.92.
- **h = 1 is noise.** AUC is about 0.5 in both periods.
- **The forecasts under-predict the 2013–2025 bull market.** The mean forecast sits 3–7 points below the realized event rate at every horizon.

## 4. Portfolio results

### Test period (10 bp)

| h | Strategy | CAGR | Sharpe | Max drawdown | Avg stock weight |
| --- | --- | --- | --- | --- | --- |
| 1 | RF tuned 0.65 | 7.34% | 0.65 | 14.2% | 0.54 |
| 1 | RF 0.50 | 12.31% | **0.82** | 21.0% | 0.94 |
| 1 | All stocks | 12.80% | 0.80 | 24.7% | 1.00 |
| 3 | RF tuned 0.70 | 5.80% | 0.54 | 14.2% | 0.43 |
| 3 | RF 0.50 | 12.22% | 0.79 | 21.0% | 0.95 |
| 3 | All stocks | 12.59% | 0.78 | 24.7% | 1.00 |
| 6 | RF tuned 0.80 | 3.88% | 0.36 | 11.8% | 0.33 |
| 6 | RF 0.50 | 12.50% | 0.77 | 23.9% | 0.98 |
| 6 | All stocks | 12.69% | 0.78 | 24.7% | 1.00 |
| any | All bills | 1.6–1.8% | n/a | 0% | 0.00 |

The rankings do not change at 0 or 20 bp. For example, the h = 1 tuned Sharpe is 0.69 / 0.65 / 0.61 at 0 / 10 / 20 bp. The common July 2013–December 2025 window (`common_window_performance.csv`) gives the same picture.

### Why the tuned cutoffs fail

- **There is no timing in the exposure.** The correlation between a tuned strategy's stock weight and the next month's market return is −0.01 / −0.04 / −0.08 for h = 1 / 3 / 6. At h = 6, the market averaged +0.71% a month when the strategy was at least half in stocks and +1.28% when it was not.
- **Where the exposure came from.** The tuned strategies were fully or mostly invested in 2014–2018. They then cut stocks to 0–17% from 2021 to 2023 at every horizon, and stayed mostly out at h = 3 and h = 6 through 2025, missing most of the post-COVID rally.
- **The drawdown gain is mostly lower exposure.** Their smaller drawdowns (12–14% vs 25%) come mainly from holding fewer stocks, not from good timing.

### Three versions of the forest compared (test, 10 bp)

| h | Defaults, 0.01 grid | Defaults, 0.1 grid | **Tuned, 0.05 grid (current)** | 0.50 rule: defaults → tuned | All stocks |
| --- | --- | --- | --- | --- | --- |
| 1 | 0.68 → 0.42 | 0.70 → 0.12 | 0.65 → **0.65** | 0.79 → 0.82 | 0.80 |
| 3 | 0.71 → 0.56 | 0.60 → 0.78 | 0.70 → **0.54** | 0.84 → 0.79 | 0.78 |
| 6 | 0.80 → 0.33 | 0.80 → 0.33 | 0.80 → **0.36** | 0.77 → 0.77 | 0.78 |

Each cell shows the tuned cutoff → test Sharpe; the "defaults" columns use ranger defaults with 500 trees.

- **Tuned-cutoff results swing widely** with small specification changes: from 0.12 to 0.65 at h = 1 and from 0.54 to 0.78 at h = 3. That is the signature of a cutoff fitted to one validation episode.
- **The 0.50 rule is stable** at 0.77–0.84 in every version.
- **Report the 0.50 rule** as the more credible measure of what the forest contributes. It contributes very little beyond buy-and-hold.

## 5. What drives the forecasts

### Variable importance

`importance_out_of_sample.csv` permutes all six lags of one variable across the test origins and records three things:

- how much the probabilities move ("reliance"),
- the change in test Brier (positive means the variable carried useful information),
- how many tuned-cutoff decisions flip.

| h | Most relied-on variables (test) | Variables whose permutation worsens test Brier by more than 0.001 |
| --- | --- | --- |
| 1 | `unrate_rel36` (mean probability change 0.037), `cpi_infl_change1`, `log_dispersion_change1`, `baa10ym_rel36` | `unrate_rel36` (+0.005), `baa10ym_rel36` (+0.001) |
| 3 | `unrate_rel36` (0.049), `mom12`, `mom6`, `cpi_infl_change1` | `unrate_rel36` (+0.008), `breadth`, `cpi_infl_change1`, `bill_yld_chg`, `short_slope_change1`, `term_spread_change1` (+0.001 to +0.002 each) |
| 6 | `unrate_rel36` (0.050), `bill_yld_rel36`, `mom12`, `ma_part` | `cpi_infl_change1` (+0.002), `breadth` (+0.001); `unrate_rel36` ≈ 0 |

- **Unemployment dominates.** Permuting `unrate_rel36` alone flips 18–26% of the tuned-cutoff decisions. In-sample OOB importance of the final fits ranks it first at every horizon and at both refit dates (`importance_oob_final_fits.csv`), followed by `bill_yld_rel36`, `mom12`, `baa10ym_rel36` and `bill_yld_chg`.
- **At h = 6 it is relied on but no longer useful out of sample.** Permuting it barely changes test Brier. That fits h = 6 failing to beat the base rate in the test.
- **Most variables add noise out of sample.** For most of the 132 inputs, permutation leaves test Brier unchanged or improves it. This is consistent with low `mtry` and subsampling winning stage 1.

### Shape of the effects

From partial dependence of the December 2025 fits (`partial_dependence.png`), setting all six lags to a sustained level:

- **`unrate_rel36` is non-linear.**
  - The probability that stocks win is highest when unemployment is about 1.2 points *below* its 3-year average, a healthy expansion: 0.67 at h = 3 and 0.73 at h = 6.
  - It is lowest near or just above the average (about 0 to +0.7 pp), typically the turn into a slowdown: 0.60 at h = 3 and 0.63 at h = 6.
  - It drops again at the extreme low end (below about −1.5 pp), which historically is the late, overheated stage of a cycle.
- **Very strong past returns lower the probability**, a mean-reversion pattern. For `mom12` around +37%, the probability is 0.58 at h = 3 and 0.64 at h = 6; moderate momentum (+6% to +14%) is most favourable.
- **Rising inflation lowers the probability.** A monthly rise in annual inflation (`cpi_infl_change1`) of more than about 0.3 pp reduces it at h = 1 and h = 3.
- **Broad market strength lowers it at h = 6.** High `ma_part`, the share of stocks above their 10-month average, reads as an "overextended" signal.

Lag-0-only partial dependence (dashed lines) is much flatter. The six lags of these slow-moving variables are highly correlated, so the forest spreads the effect across them.

### The 2021–2023 episode

`unrate_rel36` subtracts the prior 36-month mean. Once the April 2020 unemployment spike entered that mean, normal 2022 unemployment looked extraordinarily *low* relative to trend: about −0.020. That is below the 5th percentile of its 1974–2012 history (−0.014) and lands in the "overheated" region where the forest predicts poorly for stocks.

Forecast probabilities fell. The h = 3 average was 0.51 in 2021 and 0.57 in 2022, against 0.77 in 2014–2018. The tuned strategies held 0–17% in stocks from 2021 to 2023.

- **2022 worked:** stocks fell 19.5%.
- **2021 and 2023 did not:** stocks rose 23% and 24%, and those gains were missed.

This is an extrapolation problem caused by a once-in-a-century shock distorting a transformed variable, not evidence that the historical pattern has stopped holding.

## 6. Comparison with the baseline probit/logit

The baseline was run from its own scripts with the same data and protocol. It uses the fine 0.01 cutoff grid and has no hyperparameters.

| h | Test Brier: RF tuned | Test Brier: probit/logit | Test AUC: RF | Test AUC: probit |
| --- | --- | --- | --- | --- |
| 1 | **0.2256** | 0.2260 | 0.51 | 0.52 |
| 3 | **0.1950** | 0.2038 | **0.55** | 0.41 |
| 6 | 0.1933 | **0.1882** | 0.44 | 0.41 |

- **h = 3:** the tuned forest's probabilities are clearly better (lower Brier, AUC above 0.5).
- **h = 1:** they are essentially tied.
- **h = 6:** the five-variable baseline is still better.
- **Portfolios:** at h = 3 and h = 6 the baseline's tuned cutoffs (0.43–0.44) keep it 100% in stocks, so its test Sharpe equals all stocks (about 0.78). The RF's 0.50 rule ends up in the same place. Neither model's validation-tuned rule beats buy-and-hold in the test.
- **The two models' tuned cutoffs come from different grids.** The fixed 0.50 rule is the like-for-like comparison.

## 7. Caveats

- **Weak evidence.** The test has 150–156 monthly observations with one major regime (a long bull market). No dependence-robust significance test was run.
- **Choices made after seeing test results.** Hyperparameter tuning, the 0.30–0.90 cutoff grid (a departure from the protocol's 0.01 grid, rule 28) and the earlier rebuild were all adopted after test results from previous versions had been seen. Present this specification as a revised analysis, not as a pre-registered one. The fixed 0.50 rule is the comparison least affected by these choices.
- **Validation scores are selection scores.** The 18-candidate search, then the 13-cutoff search on the same 96 months, makes the reported validation Brier and Sharpe optimistic (rule 30).
- **Importance is descriptive, not causal.** Correlated predictors share credit; permutation across origins breaks correlations between variables; OOB importance is optimistic because labels overlap.
- **Macro data are current-vintage proxies** (see `Scripts/VARIABLES.md`). Since `unrate_rel36` drives the forest, unemployment revisions matter more for this model than for others.
- **Monte Carlo noise.** With 1,000 trees, each probability carries roughly a 1.3-percentage-point Monte Carlo standard error. That can flip signals near a cutoff, and it is the same order as the Brier gaps between top candidates. Results are exactly reproducible under the fixed seeds.
