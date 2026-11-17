# Ensembles: majority vote, equal-probability average, Granger-Ramanathan

These are the three ensembles of the group model spec (supervisor revision, 5 November 2026, rules 50–58). The optional no-intercept GR is not implemented. The optional annual GR update is reported separately as `ens_gr_yearly`: an expanding window from 2005-01, refitted each January on the outcomes known by the preceding December (rules 60–61), test period only. Results are in `Results/ensemble/RESULTS.md`, which `Scripts/ensemble/make_results.R` generates from the output files.

## Run it

From the project root:

```sh
Rscript --vanilla Scripts/ensemble/run_ensemble.R
```

- **Packages:** `data.table`, `jsonlite`, `digest`.
- **Data:** the CRSP files in `Data/` that the shared evaluator reads, the same files the members used. `input_hashes.csv` records their SHA-256.
- **Runtime:** about 15 seconds.
- **Process:** the 13 unit-test groups in `tests.R` run first, and any failed check stops the run.

| File | Role |
|---|---|
| `config.R` | Member export paths, the ensemble list and the fixed rules |
| `ensemble.R` | Pure functions: export checks, key merge, vote, GR solver, equal average |
| `tests.R` | Unit tests (vote and GR) |
| `run_ensemble.R` | Builds the ensembles, writes rule-66 exports, and scores everything through the shared evaluator (`Scripts/Baseline models/portfolio.R`) |

## Members (rule 50)

One model per family. The members were fixed before any test result was inspected and are declared in `config.R`:

| Member | Exports |
|---|---|
| `random_forest` | `Results/RF/{validation,test}_forecasts.csv` |
| `elastic_net` | `Results/elastic_net/{validation,test}_forecasts.csv` |
| `boosting` (`boosting_xgb`) | `Results/boosting/boosting_{validation,test}_probabilities.csv` (written by `Scripts/boosting/boosting_group_model-4.R`) |

**Baseline probit/logit is excluded.** It is the parsimonious benchmark, and no results are committed for it. To add it, put one of the two (not both) into `config.R` and rerun.

**Elastic net is included as a normal member:**
- In the equal average, its probabilities enter normally.
- In GR, the validation fit sets its weight, which may be zero.
- In the vote, its thresholded signal enters. It always votes stocks in the checked sample, so the vote reduces to "random forest OR boosting".

## The three ensembles

| Ensemble | Rule | Cutoff |
|---|---|---|
| `ens_majority_vote` | Each member votes with its own frozen validation cutoff; stocks when `sum(signals) > M/2`, ties to bills. Nothing is fitted. A second vote built from the members' 0.50 signals is the `signal_050` comparator. | Members' own |
| `ens_equal` | Mean of the member probabilities (rule 52) | Fixed 0.50 |
| `ens_gr` | `min Σ(y − a − Pw)²` subject to `w ≥ 0`, `Σw = 1`, `a` free. Fitted once per horizon on eligible validation labels (95/93/90 rows), clipped to [0, 1] at prediction (rules 53–56). | Fixed 0.50 |

**Fixed 0.50 cutoff for the probability ensembles (rule 57):** tuning a cutoff on the same validation data used to fit the weights would need a further validation set.

**GR solver:**
- It is exact: it fits equality-constrained least squares on every subset of members and keeps the best feasible fit.
- If the solver fails, it falls back to equal weights with zero intercept and flags this.
- The weights are combination weights, not measures of importance (rule 63).

**Vote export mapping:** `probability` holds the vote share and `cutoff` is `(floor(M/2)+1)/M`. An evaluator that recomputes `probability >= cutoff` therefore reproduces the vote exactly. No Brier or log loss is reported for the vote.

## Checks (stop the run on failure)

- Members are merged by `(h, origin)`, never by row order.
- Members must cover identical keys and agree on `target`, `outcome_end` and `classification_eligible`.
- Each horizon has a complete monthly calendar: 96 validation and 156 test origins.
- An unknown `outcome_end` (beyond 2025-12) may be written as `NA` or as an empty string; both are treated as missing.
- `signal == (probability >= cutoff)` and `signal_050 == (probability >= 0.50)`.
- Each member's cutoff is one value per horizon and identical between the validation and test exports.
- **Evaluator cross-check:** each member's metrics re-scored here must match its own `performance.csv` (< 1e-9). All three are checked. Boosting's `boosting_performance.csv` uses older column names (`tuned`, `cagr`, `calmar`, `window`), which are mapped before comparing.

## Benchmarks

All-stocks, all-bills and the prespecified fixed 50/50 stock/bill portfolio (rule 49) are scored on the same months. The 50/50 uses the same construction as the rebuilt boosting runner and matches its figures exactly.

## Outputs (`Results/ensemble/`)

See the file table in `RESULTS.md`. Every strategy is scored at 0, 10 and 20 bp, with both Sharpe references (common 4-week bill and matched ladder), on the full period and on the common window.

## Caveats

- **Elastic net pulls the equal average toward stocks.** At 0.50 all three ensembles are near buy-and-hold.
- **GR weights are noisy.** The 90–95 overlapping validation rows contain few independent observations, and the fitted weights are corner solutions.
- **The tuned vote's validation results are optimistic,** because it inherits the members' validation-tuned cutoffs. Only test results count.
- **Choose nothing from test rankings** (rule 58), and make no significance claims without dependence-aware inference (rule 71, not implemented).
