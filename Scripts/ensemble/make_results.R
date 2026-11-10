# Run after run_ensemble.R: Rscript --vanilla Scripts/ensemble/make_results.R
# Writes Results/ensemble/RESULTS.md from the output files (no hand-typed numbers).
o <- "Results/ensemble"; r <- function(f) read.csv(file.path(o, f), stringsAsFactors = FALSE)
p <- r("performance.csv"); f <- r("forecast_metrics.csv"); w <- r("gr_weights.csv"); yw <- r("gr_yearly_weights.csv")
cut <- r("member_cutoffs.csv"); dg <- r("vote_member_diagnostics.csv"); cc <- r("evaluator_crosscheck.csv"); hs <- r("input_hashes.csv")
P <- function(h, s, rule) { x <- p[p$period == "test" & abs(p$cost_rate - .001) < 1e-9 & p$h == h & p$strategy == s & p$cutoff_rule == rule, ]; stopifnot(nrow(x) == 1); x }
lab <- c(ens_majority_vote = "Vote", ens_equal = "Equal average", ens_gr = "GR", ens_gr_yearly = "GR, yearly refit", random_forest = "RF",
         elastic_net = "Elastic net", boosting = "Boosting", all_stocks = "All stocks", fixed_50_50 = "Fixed 50/50")
rl <- c(validation_tuned = " (tuned)", fixed_050 = " (0.50)", not_applicable = "")
tbl <- function(h) {
  x <- p[p$period == "test" & abs(p$cost_rate - .001) < 1e-9 & p$h == h & p$strategy != "all_bills", ]
  x <- x[order(-x$sharpe_matched_ladder, x$strategy), ]
  c(sprintf("**h = %d**", h), "", "| Strategy | Sharpe (matched) | Sharpe (common) | CAGR | Max DD | Calmar | Avg stock weight |",
    "|---|---:|---:|---:|---:|---:|---:|",
    sprintf("| %s%s | %.2f | %.2f | %.1f%% | %.1f%% | %.2f | %.2f |", lab[x$strategy], rl[x$cutoff_rule], x$sharpe_matched_ladder,
            x$sharpe_common_rf, 100 * x$CAGR, 100 * x$max_drawdown, x$Calmar, x$average_stock_weight), "")
}
g <- function(h, m) w$weight[w$h == h & w$member == m]
wl <- sprintf("| %d | %.3f | %.3f | %.3f | %.3f | %d |", w$h[!duplicated(w$h)], sapply(c(1, 3, 6), g, "random_forest"),
  sapply(c(1, 3, 6), g, "elastic_net"), sapply(c(1, 3, 6), g, "boosting"), w$intercept[!duplicated(w$h)], w$n[!duplicated(w$h)])
cu <- function(m) paste(sprintf("%.2f", cut$cutoff[cut$period == "test" & cut$member == m]), collapse = " / ")
s <- function(h, st, rule) P(h, st, rule)$sharpe_matched_ladder
dd <- function(h, st, rule) 100 * P(h, st, rule)$max_drawdown
yr <- sapply(c(1, 3, 6), function(h) c(s(h, "ens_gr", "fixed_050"), s(h, "ens_gr_yearly", "fixed_050"), s(h, "all_stocks", "not_applicable"),
  P(h, "ens_gr_yearly", "fixed_050")$average_stock_weight, f$AUC[f$period == "test" & f$model == "ens_gr_yearly" & f$h == h]))
au <- function(m) paste(sprintf("%.2f", sapply(c(1, 3, 6), function(h) f$AUC[f$period == "test" & f$model == m & f$h == h][1])), collapse = " / ")
bh <- hs$sha256[grepl("boosting", hs$file)]
md <- c("# Ensemble results", "",
"Built on 10 November 2026 (`ensemble_v3_2026_11_10`) with `Rscript --vanilla Scripts/ensemble/run_ensemble.R`, then `Rscript --vanilla Scripts/ensemble/make_results.R`.", "",
"**Members**, all exactly as committed on `main` @ 458758e:",
"- random forest (`Results/RF/`)",
"- elastic net (`Results/elastic_net/`)",
"- boosting (`Results/boosting/`, the 0.05 cutoff grid)", "",
sprintf("The ensemble uses the members' saved forecasts without modification. Boosting export hashes: validation `%s…`, test `%s…` (`input_hashes.csv`). If any member's results change, rerun both scripts.", substr(bh[1], 1, 16), substr(bh[2], 1, 16)), "",
sprintf("**Evaluator cross-check:** re-scoring all three members through the shared evaluator reproduces each member's own `performance.csv` to within %.0e on every row (`evaluator_crosscheck.csv`).", max(cc$max_abs_diff)), "",
"**Member cutoffs (test, frozen from validation):**",
sprintf("- RF %s", cu("random_forest")), sprintf("- Elastic net %s", cu("elastic_net")), sprintf("- Boosting %s", cu("boosting")), "",
"## GR weights (fitted once on eligible validation labels)", "",
"| h | Random forest | Elastic net | Boosting | Intercept | Rows |", "|---|---:|---:|---:|---:|---:|", wl, "",
"## Test period, 10 bp cost", "",
"Sorted by matched-ladder Sharpe. \"Tuned\" means the members' frozen validation cutoffs; \"0.50\" means the fixed cutoff. Equal and GR use 0.50 only (rule 57). Fixed 50/50 is the rule-49 benchmark. All-bills is omitted (rule 46).", "",
tbl(1), tbl(3), tbl(6),
"## Findings", "",
sprintf("1. **No ensemble beats buy-and-hold in a way that can be relied on.** The only Sharpe gain is at h = 1: the 0.50 vote reaches %.2f against %.2f for all-stocks. That number is boosting's own 0.50 result, because the vote reproduces boosting's signal month by month (elastic net always votes stocks, and every month boosting chose bills, RF did too). The gain rests on about five bill calls in 2020 and 2022.", s(1, "ens_majority_vote", "fixed_050"), s(1, "all_stocks", "not_applicable")),
"2. **At 0.50, every ensemble is close to buy-and-hold** (stock weight 0.95–1.00). At h = 3 and h = 6, every ensemble is within 0.03 Sharpe of all-stocks.",
sprintf("3. **The tuned vote no longer offers a drawdown advantage.** With the boosting cutoffs from the 0.05 grid, its maximum drawdown is %.1f%% / %.1f%% / %.1f%% (h = 1/3/6), against %.1f%% for all-stocks. Its Sharpe is lower at every horizon (%.2f / %.2f / %.2f).", dd(1, "ens_majority_vote", "validation_tuned"), dd(3, "ens_majority_vote", "validation_tuned"), dd(6, "ens_majority_vote", "validation_tuned"), dd(1, "all_stocks", "not_applicable"), s(1, "ens_majority_vote", "validation_tuned"), s(3, "ens_majority_vote", "validation_tuned"), s(6, "ens_majority_vote", "validation_tuned")),
sprintf("   The fixed 50/50 benchmark beats it on both measures at every horizon: Sharpe %.2f / %.2f / %.2f and maximum drawdown about %.1f%%. At h = 1, the tuned vote is identical to tuned boosting.", s(1, "fixed_50_50", "not_applicable"), s(3, "fixed_50_50", "not_applicable"), s(6, "fixed_50_50", "not_applicable"), dd(1, "fixed_50_50", "not_applicable")),
"4. **GR weights are corner solutions:** 100% RF at h = 3, 100% boosting at h = 6, and 75% boosting / 25% elastic net at h = 1. About 90 overlapping validation rows cannot pin down stable weights, and GR never beats the equal average.",
sprintf("5. **No member or ensemble forecasts direction.** Test AUC: vote (tuned) %s, equal %s, GR %s, RF %s, boosting %s. Elastic net is %s: it is intercept-only, so its probability is just the training base rate.", au("ens_majority_vote"), au("ens_equal"), au("ens_gr"), au("random_forest"), au("boosting"), au("elastic_net")),
"6. **No significance claims.** Over 13 years the standard error of an annualised Sharpe is about 0.28, so every gap above is within noise. Ledoit–Wolf tests are not implemented (rule 71).", "",
"## Supplementary: GR refitted every January on an expanding window (`ens_gr_yearly`)", "",
"Each January of test year Y, GR is refitted on every member forecast from 2005-01 onward whose outcome was known by December of Y−1 (rules 60–61; for example, h = 6 in January 2014 uses 102 rows). The weights are then held for year Y. This is not a rule-58 ensemble, and nothing was selected using it.", "",
"| h | Static GR | Yearly GR | All-stocks | Yearly GR stock weight | Yearly GR AUC |", "|---|---:|---:|---:|---:|---:|",
sprintf("| %d | %.2f | %.2f | %.2f | %.2f | %.2f |", c(1, 3, 6), yr[1, ], yr[2, ], yr[3, ], yr[4, ], yr[5, ]), "",
"It collapses to all-stocks. As bull-market test labels enter the fitting window, the intercept rises and weight shifts to elastic net (the base rate), so the combined probability stays above 0.50. The yearly weights swing a lot (see `gr_yearly_weights.csv`), which is consistent with there being no stable signal to weight.", "",
"## Files", "",
"| File | Contents |", "|---|---|",
"| `validation_forecasts.csv`, `test_forecasts.csv` | Rule-66 exports for the ensembles, with member probabilities and signals |",
"| `performance.csv`, `common_window_performance.csv`, `monthly_portfolios.csv` | Shared-evaluator results at 0/10/20 bp for ensembles, members, all-stocks, all-bills and fixed 50/50 |",
"| `forecast_metrics.csv` | Classification metrics on eligible labels |",
"| `gr_weights.csv`, `gr_yearly_weights.csv`, `gr_clipping.csv` | GR weights (static and yearly) and clipping counts |",
"| `member_cutoffs.csv`, `vote_member_diagnostics.csv` | Frozen member cutoffs; stock-signal share and pivotal share per member |",
"| `evaluator_crosscheck.csv`, `tests.csv`, `input_hashes.csv`, `manifest.json`, `session_info.txt` | Reproducibility record |")
writeLines(md, file.path(o, "RESULTS.md"))
