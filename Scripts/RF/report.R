rf_markdown_table <- function(d) {
  lines <- apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"), lines)
}

write_rf_report <- function(out, cfg, selected, hyper_scores, performance, diagnostics, importance_oos, tests) {
  fmt <- function(x, d = 3) ifelse(is.na(x), "N/A", sprintf(paste0("%.", d, "f"), x))
  pct <- function(x) sprintf("%.2f%%", 100 * x)
  hs <- hyper_scores[order(hyper_scores$h, hyper_scores$Brier), ]
  hp <- data.frame(Horizon = hs$h, Brier_rank = hs$Brier_rank, Config = hs$config_id, mtry = hs$mtry,
    min_node = hs$min_node_size, Sampling = ifelse(hs$replace, "bootstrap", "50% subsample"),
    Validation_Brier = fmt(hs$Brier, 4), Log_loss = fmt(hs$log_loss, 4), AUC = fmt(hs$AUC),
    Prevalence_Brier = fmt(hs$Brier_train_prevalence, 4),
    Note = paste0(ifelse(hs$selected, "**selected**", ""), ifelse(hs$is_ranger_default, " ranger default", "")))
  choices <- data.frame(Horizon = selected$h, Config = selected$config_id, mtry = selected$mtry,
    min_node_size = selected$min_node_size, Sampling = ifelse(selected$replace, "bootstrap", "50% subsample"),
    Validation_Brier = fmt(selected$validation_Brier, 4),
    Default_Brier = fmt(selected$validation_Brier_ranger_default, 4),
    Prevalence_Brier = fmt(selected$validation_Brier_train_prevalence, 4),
    Cutoff = sprintf("%.2f", selected$cutoff),
    Tuned_validation_Sharpe_matched = fmt(selected$validation_sharpe_matched_ladder),
    Fixed_050_validation_Sharpe_matched = fmt(selected$fixed_050_validation_sharpe_matched_ladder))
  primary <- performance[performance$period == "test" & performance$cost_rate == cfg$primary_cost, ]
  table <- data.frame(Horizon = primary$h, Strategy = primary$strategy, Cutoff_rule = primary$cutoff_rule,
    Cutoff = ifelse(is.na(primary$cutoff), "-", sprintf("%.2f", primary$cutoff)),
    Start = primary$start, Months = primary$n_months, CAGR = pct(primary$CAGR),
    Sharpe_common_4week = fmt(primary$sharpe_common_rf), Sharpe_matched = fmt(primary$sharpe_matched_ladder),
    Max_drawdown = pct(primary$max_drawdown), Avg_stock_weight = sprintf("%.2f", primary$average_stock_weight),
    Annual_turnover = sprintf("%.2f", primary$annual_turnover), check.names = FALSE)
  dg <- diagnostics[diagnostics$cutoff_rule == "validation_tuned", ]
  fm <- data.frame(Period = dg$period, Horizon = dg$h, n = dg$n, Event_rate = fmt(dg$event_rate),
    Mean_p = fmt(dg$mean_probability), Brier = fmt(dg$Brier, 4),
    Brier_prevalence = fmt(dg$Brier_train_prevalence, 4),
    Brier_skill = fmt(dg$Brier_skill_vs_prevalence), Log_loss = fmt(dg$log_loss, 4), AUC = fmt(dg$AUC),
    Accuracy_tuned = fmt(dg$accuracy), Accuracy_050 = fmt(dg$accuracy_at_050))
  imp <- importance_oos[importance_oos$period == "test", ]
  top <- do.call(rbind, lapply(cfg$horizons, function(h) {
    x <- imp[imp$h == h, ]; x <- x[order(x$rank_reliance), ][1:5, ]
    data.frame(Horizon = h, Reliance_rank = x$rank_reliance, Variable = x$base_variable,
      Mean_abs_prob_change = fmt(x$mean_abs_probability_change, 4),
      Delta_Brier = sprintf("%+.4f", x$delta_Brier),
      Signal_flip_rate = fmt(x$signal_flip_rate_at_tuned_cutoff))
  }))
  lines <- c("# Random forest results - supervisor protocol", "",
    paste0("Configuration: `", cfg$version, "`. A ranger probability forest on all 132 predictors (22 variables x lags 0-5) is fitted separately at h = 1, 3 and 6. ",
           "Every validation and test origin is refitted monthly on an expanding sample of labels whose outcomes ended by that origin. ",
           "Selection follows the group standard in two stages, both on 2005-2012 validation only: (1) per horizon, the hyperparameter candidate with the lowest validation Brier among ", nrow(cfg$hyper_grid), " (mtry x min.node.size x sampling, ", cfg$num_trees, " trees); (2) with that candidate's saved validation forecasts fixed, the cutoff with the highest validation matched-ladder Sharpe at 10 bp.",
           " The fixed 0.50 cutoff is reported alongside. All-stock and all-bill benchmarks use the same scored months."), "",
    "Generated automatically by `Scripts/RF/run_rf.R`. See `INTERPRETATION.md` in this folder for the written interpretation.", "",
    "## Stage 1: hyperparameter selection (validation Brier, lower is better)", "", rf_markdown_table(hp), "",
    "Prevalence_Brier is the Brier of forecasting the training-sample event rate known at each origin; a candidate below it has probability skill over that naive benchmark. Exact ties (1e-12) would go to larger min.node.size, then smaller mtry, then bootstrap.", "",
    "## Stage 2: selected specification and cutoff (validation, 10 bp)", "", rf_markdown_table(choices), "",
    sprintf("Cutoff grid (random forest only): %s. The baseline searches 0.00-1.00 in 0.01 steps. These validation Sharpe ratios are selection scores, not independent performance estimates.", paste(sprintf("%.2f", cfg$cutoffs), collapse = ", ")), "",
    "## Forecast quality (eligible completed labels)", "", rf_markdown_table(fm), "",
    "Brier_skill = 1 - Brier / Brier_prevalence; positive means better than the expanding training base rate.", "",
    "## Test performance at 10 bp", "", rf_markdown_table(table), "",
    "Sharpe definitions, the overlapping-cohort engine, cost convention and benchmark construction are identical to the baseline (shared `portfolio.R`). No model or cutoff rule is selected using this table.", "",
    "## What the forest relies on (test period, out-of-sample grouped permutation)", "", rf_markdown_table(top), "",
    "Reliance = mean absolute change in the forecast probability when all six lags of a variable are permuted across test origins. Delta_Brier > 0 means permuting the variable worsens test Brier (it carried useful information); <= 0 means it did not help out of sample. Signal_flip_rate is the share of origins whose tuned-cutoff decision changes. These are descriptive, not causal, and correlated predictors share credit.", "",
    "## Files", "",
    "- performance.csv / common_window_performance.csv / monthly_portfolios.csv / benchmark_returns.csv: same schema as the baseline.",
    "- validation_forecasts.csv / test_forecasts.csv: one probability per origin with tuned and 0.50 cutoffs/signals, fitting audit, config_id, seed and training prevalence.",
    "- validation_cutoff_grid.csv / selected_cutoffs.csv: full cutoff search and frozen choices.",
    "- forecast_metrics.csv / calibration.csv: probability and classification diagnostics.",
    "- validation_hyperparameter_scores.csv / validation_candidate_forecasts.csv: every grid candidate's Brier scores and expanding validation forecasts.",
    "- refit_audit.csv: every validation (all candidates) and test fit with training endpoint, sample size, seed and hyperparameters.",
    "- importance_out_of_sample.csv / importance_oob_final_fits.csv / partial_dependence.csv: interpretation outputs.",
    "- feature_names.csv, frozen_specification.rds, validation_state.rds, validation_models.rds, final_models.rds, manifest.json, input_hashes.csv, package_versions.csv, session_info.txt, tests.csv.",
    "- Figures: test_wealth.png, validation_hyperparameters.png, validation_cutoffs.png, probabilities.png, importance.png, partial_dependence.png.", "",
    sprintf("%d focused test groups passed; parallel/serial spot checks, final-origin refit reproducibility, refit-count, maturity and input-hash checks all passed.", nrow(tests)))
  writeLines(lines, file.path(out, "RESULTS.md"))
}

draw_rf_figures <- function(out, paths, grid, selected, hyper_scores, forecasts, importance_oos, pd, cfg) {
  colors <- c(rf = "#2463A7", all_stocks = "#293A45", all_bills = "#78916D", accent = "#BA5A31")
  as_date <- function(x) as.Date(paste0(x, "-01"))
  panel <- function(file, width = 1800, height = 780) {
    grDevices::png(file.path(out, file), width = width, height = height, res = 150)
    graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  }
  panel("test_wealth.png")
  for (h in cfg$horizons) {
    d <- paths[paths$period == "test" & paths$h == h & paths$cost_rate == cfg$primary_cost, ]
    graphics::plot(range(as_date(d$month)), range(c(1, d$wealth)), type = "n", xlab = "Return month",
                   ylab = "Wealth (initial = 1)", main = paste(h, "month horizon"))
    graphics::grid(col = "#E7E7E7")
    line <- function(v, col, lty = 1) graphics::lines(as_date(v$month), v$wealth, col = col, lwd = 2, lty = lty)
    line(d[d$strategy == cfg$model & d$cutoff_rule == "validation_tuned", ], colors["rf"])
    line(d[d$strategy == cfg$model & d$cutoff_rule == "fixed_050", ], colors["accent"], 2)
    line(d[d$strategy == "all_stocks", ], colors["all_stocks"])
    line(d[d$strategy == "all_bills", ], colors["all_bills"])
    graphics::legend("topleft", legend = c("RF tuned", "RF 0.50", "All stocks", "All bills"),
      col = colors[c("rf", "accent", "all_stocks", "all_bills")], lty = c(1, 2, 1, 1), lwd = 2, bty = "n", cex = .75)
  }
  graphics::mtext("Random forest test portfolios | 10 bp costs", outer = TRUE, cex = 1.05)
  grDevices::dev.off()

  panel("validation_cutoffs.png")
  for (h in cfg$horizons) {
    d <- grid[grid$h == h, ]; s <- selected[selected$h == h, ]
    graphics::plot(d$cutoff, d$sharpe_matched_ladder, type = "b", pch = 19, lwd = 2, col = colors["rf"],
      xlab = "Probability cutoff", ylab = "Validation matched-bill Sharpe", main = paste(h, "month horizon"))
    graphics::grid(col = "#E7E7E7")
    graphics::abline(v = .5, lty = 3, col = "#888888")
    graphics::points(s$cutoff, s$validation_Sharpe, pch = 19, col = colors["accent"])
    graphics::text(s$cutoff, s$validation_Sharpe, sprintf("%.2f", s$cutoff), pos = 3, cex = .8)
  }
  graphics::mtext("Expanding validation | matched-bill cutoff objective | 10 bp costs", outer = TRUE, cex = 1.05)
  grDevices::dev.off()

  panel("validation_hyperparameters.png", height = 900)
  graphics::par(mar = c(7, 4, 3, 1))
  for (h in cfg$horizons) {
    d <- hyper_scores[hyper_scores$h == h, ]
    d <- d[order(d$replace, d$mtry, d$min_node_size, decreasing = c(TRUE, FALSE, FALSE), method = "radix"), ]
    rng <- range(c(d$Brier, d$Brier_train_prevalence))
    col <- ifelse(d$selected, colors["accent"], ifelse(d$replace, colors["rf"], "#7FA7D3"))
    graphics::plot(seq_len(nrow(d)), d$Brier, ylim = rng, xaxt = "n", pch = ifelse(d$is_ranger_default, 17, 19),
      col = col, cex = 1.2, xlab = "", ylab = "Validation Brier (lower is better)", main = paste(h, "month horizon"))
    graphics::axis(1, at = seq_len(nrow(d)), las = 2, cex.axis = .62,
      labels = sprintf("m%d n%d %s", d$mtry, d$min_node_size, ifelse(d$replace, "boot", "sub")))
    graphics::abline(h = d$Brier_train_prevalence[1], lty = 2, col = "#888888")
    graphics::legend("topright", legend = c("bootstrap", "50% subsample", "selected", "ranger default (triangle)", "training base rate"),
      pch = c(19, 19, 19, 17, NA), lty = c(NA, NA, NA, NA, 2), col = c(colors["rf"], "#7FA7D3", colors["accent"], "#555555", "#888888"),
      bty = "n", cex = .65)
  }
  graphics::mtext("Stage 1: validation Brier for each hyperparameter candidate (m = mtry, n = min.node.size)", outer = TRUE, cex = 1)
  grDevices::dev.off()

  grDevices::png(file.path(out, "probabilities.png"), width = 1800, height = 1500, res = 150)
  graphics::par(mfrow = c(3, 1), mar = c(3, 4, 2.5, 1))
  for (h in cfg$horizons) {
    f <- forecasts[forecasts$h == h, ]; dt <- as_date(f$origin); s <- selected[selected$h == h, ]
    graphics::plot(dt, f$probability, type = "n", ylim = c(0, 1), xlab = "", ylab = "P(stocks beat bills)",
                   main = sprintf("h = %d: forecast probability by origin (tuned cutoff %.2f)", h, s$cutoff))
    graphics::rect(as_date(cfg$test_start), -1, max(dt) + 31, 2, col = "#F4F4F4", border = NA)
    ok <- !is.na(f$target)
    graphics::points(dt[ok], ifelse(f$target[ok] == 1, .98, .02), pch = 15, cex = .4,
                     col = ifelse(f$target[ok] == 1, colors["all_bills"], colors["accent"]))
    graphics::lines(dt, f$probability, col = colors["rf"], lwd = 1.6)
    graphics::abline(h = c(s$cutoff, .5), lty = c(2, 3), col = c(colors["accent"], "#888888"))
    graphics::text(as_date(cfg$test_start), .08, "test", pos = 4, cex = .8, col = "#666666")
    graphics::text(as_date(cfg$validation_start), .08, "validation", pos = 4, cex = .8, col = "#666666")
  }
  grDevices::dev.off()

  panel("importance.png", height = 900)
  graphics::par(mar = c(4, 12, 3, 1))
  for (h in cfg$horizons) {
    d <- importance_oos[importance_oos$period == "test" & importance_oos$h == h, ]
    d <- d[order(d$mean_abs_probability_change), ]
    graphics::barplot(d$mean_abs_probability_change, names.arg = d$base_variable, horiz = TRUE, las = 1,
      cex.names = .7, col = ifelse(d$delta_Brier > 0, colors["rf"], "#C9C9C9"), border = NA,
      xlab = "Mean |change in probability| when permuted", main = paste(h, "month horizon"))
  }
  graphics::mtext("Test-period grouped permutation reliance (blue = permuting worsens test Brier)", outer = TRUE, cex = 1.05)
  grDevices::dev.off()

  grDevices::png(file.path(out, "partial_dependence.png"), width = 1800, height = 1350, res = 150)
  graphics::par(mfrow = c(3, 4), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  for (h in cfg$horizons) for (v in unique(pd$base_variable[pd$h == h])) {
    d <- pd[pd$h == h & pd$base_variable == v, ]
    rng <- range(c(d$partial_dependence, d$partial_dependence_lag0_only))
    graphics::plot(d$value, d$partial_dependence, type = "l", lwd = 2, col = colors["rf"], ylim = rng,
      xlab = paste(v, "(value held for all six lags)"), ylab = "Average P(stocks win)",
      main = sprintf("h = %d: %s", h, v), cex.main = .9, cex.lab = .8)
    graphics::lines(d$value, d$partial_dependence_lag0_only, lwd = 1.5, lty = 2, col = colors["accent"])
    graphics::rug(d$value, col = "#999999")
    if (v == unique(pd$base_variable[pd$h == h])[1])
      graphics::legend("topright", legend = c("all six lags", "lag0 only"), lty = c(1, 2), lwd = 2,
                       col = colors[c("rf", "accent")], bty = "n", cex = .7)
  }
  graphics::mtext("Partial dependence of the final test fit (2025-12), 5th-95th percentile grid", outer = TRUE, cex = .95)
  grDevices::dev.off()
}
