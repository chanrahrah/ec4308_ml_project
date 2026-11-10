# Ensemble results

Built on 10 November 2026 (`ensemble_v3_2026_11_10`) with `Rscript --vanilla Scripts/ensemble/run_ensemble.R`, then `Rscript --vanilla Scripts/ensemble/make_results.R`.

**Members**, all exactly as committed on `main` @ 458758e:
- random forest (`Results/RF/`)
- elastic net (`Results/elastic_net/`)
- boosting (`Results/boosting/`, the 0.05 cutoff grid)

The ensemble uses the members' saved forecasts without modification. Boosting export hashes: validation `366cbd52b7ad4d6a…`, test `3799aff80d53d0bd…` (`input_hashes.csv`). If any member's results change, rerun both scripts.

**Evaluator cross-check:** re-scoring all three members through the shared evaluator reproduces each member's own `performance.csv` to within 2e-14 on every row (`evaluator_crosscheck.csv`).

**Member cutoffs (test, frozen from validation):**
- RF 0.65 / 0.70 / 0.80
- Elastic net 0.50 / 0.50 / 0.50
- Boosting 0.60 / 0.60 / 0.75

## GR weights (fitted once on eligible validation labels)

| h | Random forest | Elastic net | Boosting | Intercept | Rows |
|---|---:|---:|---:|---:|---:|
| 1 | 0.000 | 0.249 | 0.751 | -0.003 | 95 |
| 3 | 1.000 | 0.000 | 0.000 | 0.021 | 93 |
| 6 | 0.000 | 0.000 | 1.000 | -0.006 | 90 |

## Test period, 10 bp cost

Sorted by matched-ladder Sharpe. "Tuned" means the members' frozen validation cutoffs; "0.50" means the fixed cutoff. Equal and GR use 0.50 only (rule 57). Fixed 50/50 is the rule-49 benchmark. All-bills is omitted (rule 46).

**h = 1**

| Strategy | Sharpe (matched) | Sharpe (common) | CAGR | Max DD | Calmar | Avg stock weight |
|---|---:|---:|---:|---:|---:|---:|
| Boosting (0.50) | 0.94 | 0.94 | 14.5% | 21.0% | 0.69 | 0.97 |
| Vote (0.50) | 0.94 | 0.94 | 14.5% | 21.0% | 0.69 | 0.97 |
| Equal average (0.50) | 0.87 | 0.87 | 13.6% | 21.1% | 0.65 | 0.99 |
| GR (0.50) | 0.86 | 0.86 | 13.6% | 21.0% | 0.65 | 0.99 |
| RF (0.50) | 0.82 | 0.82 | 12.3% | 21.0% | 0.59 | 0.94 |
| All stocks | 0.80 | 0.80 | 12.8% | 24.7% | 0.52 | 1.00 |
| Elastic net (tuned) | 0.80 | 0.80 | 12.8% | 24.7% | 0.52 | 1.00 |
| Elastic net (0.50) | 0.80 | 0.80 | 12.8% | 24.7% | 0.52 | 1.00 |
| GR, yearly refit (0.50) | 0.80 | 0.80 | 12.8% | 24.7% | 0.52 | 1.00 |
| Fixed 50/50 | 0.80 | 0.80 | 7.3% | 12.6% | 0.58 | 0.50 |
| Boosting (tuned) | 0.69 | 0.69 | 10.2% | 21.5% | 0.47 | 0.83 |
| Vote (tuned) | 0.69 | 0.69 | 10.2% | 21.5% | 0.47 | 0.83 |
| RF (tuned) | 0.65 | 0.65 | 7.3% | 14.2% | 0.52 | 0.54 |

**h = 3**

| Strategy | Sharpe (matched) | Sharpe (common) | CAGR | Max DD | Calmar | Avg stock weight |
|---|---:|---:|---:|---:|---:|---:|
| RF (0.50) | 0.79 | 0.79 | 12.2% | 21.0% | 0.58 | 0.95 |
| All stocks | 0.78 | 0.79 | 12.6% | 24.7% | 0.51 | 1.00 |
| Boosting (0.50) | 0.78 | 0.79 | 12.6% | 24.7% | 0.51 | 1.00 |
| Elastic net (tuned) | 0.78 | 0.79 | 12.6% | 24.7% | 0.51 | 1.00 |
| Elastic net (0.50) | 0.78 | 0.79 | 12.6% | 24.7% | 0.51 | 1.00 |
| Equal average (0.50) | 0.78 | 0.79 | 12.6% | 24.7% | 0.51 | 1.00 |
| GR, yearly refit (0.50) | 0.78 | 0.79 | 12.6% | 24.7% | 0.51 | 1.00 |
| Vote (0.50) | 0.78 | 0.79 | 12.6% | 24.7% | 0.51 | 1.00 |
| Fixed 50/50 | 0.78 | 0.79 | 7.3% | 12.5% | 0.58 | 0.50 |
| GR (0.50) | 0.77 | 0.78 | 12.1% | 22.1% | 0.55 | 0.95 |
| Vote (tuned) | 0.75 | 0.75 | 11.4% | 21.0% | 0.54 | 0.88 |
| Boosting (tuned) | 0.71 | 0.71 | 10.6% | 21.0% | 0.50 | 0.87 |
| RF (tuned) | 0.54 | 0.55 | 5.8% | 14.2% | 0.41 | 0.43 |

**h = 6**

| Strategy | Sharpe (matched) | Sharpe (common) | CAGR | Max DD | Calmar | Avg stock weight |
|---|---:|---:|---:|---:|---:|---:|
| All stocks | 0.78 | 0.78 | 12.7% | 24.7% | 0.51 | 1.00 |
| Elastic net (tuned) | 0.78 | 0.78 | 12.7% | 24.7% | 0.51 | 1.00 |
| Elastic net (0.50) | 0.78 | 0.78 | 12.7% | 24.7% | 0.51 | 1.00 |
| Fixed 50/50 | 0.78 | 0.79 | 7.4% | 12.5% | 0.59 | 0.50 |
| Equal average (0.50) | 0.77 | 0.78 | 12.6% | 24.7% | 0.51 | 0.99 |
| Vote (0.50) | 0.77 | 0.78 | 12.5% | 23.9% | 0.52 | 0.98 |
| RF (0.50) | 0.77 | 0.78 | 12.5% | 23.9% | 0.52 | 0.98 |
| Boosting (0.50) | 0.75 | 0.76 | 12.1% | 23.9% | 0.51 | 0.95 |
| GR (0.50) | 0.75 | 0.76 | 12.1% | 23.9% | 0.51 | 0.95 |
| GR, yearly refit (0.50) | 0.74 | 0.75 | 12.0% | 24.7% | 0.49 | 0.97 |
| Boosting (tuned) | 0.66 | 0.67 | 8.1% | 18.2% | 0.45 | 0.64 |
| Vote (tuned) | 0.66 | 0.67 | 8.1% | 18.2% | 0.45 | 0.64 |
| RF (tuned) | 0.36 | 0.38 | 3.9% | 11.8% | 0.33 | 0.33 |

## Findings

1. **No ensemble beats buy-and-hold in a way that can be relied on.** The only Sharpe gain is at h = 1: the 0.50 vote reaches 0.94 against 0.80 for all-stocks. That number is boosting's own 0.50 result, because the vote reproduces boosting's signal month by month (elastic net always votes stocks, and every month boosting chose bills, RF did too). The gain rests on about five bill calls in 2020 and 2022.
2. **At 0.50, every ensemble is close to buy-and-hold** (stock weight 0.95–1.00). At h = 3 and h = 6, every ensemble is within 0.03 Sharpe of all-stocks.
3. **The tuned vote no longer offers a drawdown advantage.** With the boosting cutoffs from the 0.05 grid, its maximum drawdown is 21.5% / 21.0% / 18.2% (h = 1/3/6), against 24.7% for all-stocks. Its Sharpe is lower at every horizon (0.69 / 0.75 / 0.66).
   The fixed 50/50 benchmark beats it on both measures at every horizon: Sharpe 0.80 / 0.78 / 0.78 and maximum drawdown about 12.6%. At h = 1, the tuned vote is identical to tuned boosting.
4. **GR weights are corner solutions:** 100% RF at h = 3, 100% boosting at h = 6, and 75% boosting / 25% elastic net at h = 1. About 90 overlapping validation rows cannot pin down stable weights, and GR never beats the equal average.
5. **No member or ensemble forecasts direction.** Test AUC: vote (tuned) 0.51 / 0.51 / 0.48, equal 0.49 / 0.51 / 0.43, GR 0.48 / 0.55 / 0.45, RF 0.51 / 0.55 / 0.44, boosting 0.48 / 0.46 / 0.45. Elastic net is 0.37 / 0.33 / 0.22: it is intercept-only, so its probability is just the training base rate.
6. **No significance claims.** Over 13 years the standard error of an annualised Sharpe is about 0.28, so every gap above is within noise. Ledoit–Wolf tests are not implemented (rule 71).

## Supplementary: GR refitted every January on an expanding window (`ens_gr_yearly`)

Each January of test year Y, GR is refitted on every member forecast from 2005-01 onward whose outcome was known by December of Y−1 (rules 60–61; for example, h = 6 in January 2014 uses 102 rows). The weights are then held for year Y. This is not a rule-58 ensemble, and nothing was selected using it.

| h | Static GR | Yearly GR | All-stocks | Yearly GR stock weight | Yearly GR AUC |
|---|---:|---:|---:|---:|---:|
| 1 | 0.86 | 0.80 | 0.80 | 1.00 | 0.41 |
| 3 | 0.77 | 0.78 | 0.78 | 1.00 | 0.46 |
| 6 | 0.75 | 0.74 | 0.78 | 0.97 | 0.44 |

It collapses to all-stocks. As bull-market test labels enter the fitting window, the intercept rises and weight shifts to elastic net (the base rate), so the combined probability stays above 0.50. The yearly weights swing a lot (see `gr_yearly_weights.csv`), which is consistent with there being no stable signal to weight.

## Files

| File | Contents |
|---|---|
| `validation_forecasts.csv`, `test_forecasts.csv` | Rule-66 exports for the ensembles, with member probabilities and signals |
| `performance.csv`, `common_window_performance.csv`, `monthly_portfolios.csv` | Shared-evaluator results at 0/10/20 bp for ensembles, members, all-stocks, all-bills and fixed 50/50 |
| `forecast_metrics.csv` | Classification metrics on eligible labels |
| `gr_weights.csv`, `gr_yearly_weights.csv`, `gr_clipping.csv` | GR weights (static and yearly) and clipping counts |
| `member_cutoffs.csv`, `vote_member_diagnostics.csv` | Frozen member cutoffs; stock-signal share and pivotal share per member |
| `evaluator_crosscheck.csv`, `tests.csv`, `input_hashes.csv`, `manifest.json`, `session_info.txt` | Reproducibility record |
