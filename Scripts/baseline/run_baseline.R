# Baseline probit and logit models, plus the all-stocks, all-bills and 50/50 benchmarks.
# Run: Rscript --vanilla Scripts/baseline/run_baseline.R   (or open and click Source)
script <- c(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE),
            tryCatch(sys.frame(1)$ofile, error = function(e) NULL))
source(file.path(if (length(script)) dirname(script[1]) else "Scripts/baseline", "..", "utils", "utils.R"))
for (file in c("model.R", "tests.R", "figures.R")) source(file.path(SCRIPTS, "baseline", file))

cfg <- baseline_config()
returns <- load_return_data()
checks <- run_checks(c(shared_checks(returns, cfg), baseline_checks(cfg)))

# Fit each validation origin, choose the cutoff on validation, then fit the test origins.
forecasts <- coefficients <- grids <- list()
for (h in cfg$horizons) {
  d <- load_model_data(h, cfg, cfg$features)
  for (link in cfg$links) {
    val <- expanding_forecasts(d, h, link, cfg, "validation")
    grid <- select_cutoff(val$forecasts, returns, h, cfg)
    test <- expanding_forecasts(d, h, link, cfg, "test")
    forecasts[[length(forecasts) + 1L]] <- add_signals(rbind(val$forecasts, test$forecasts), grid$cutoff[grid$selected], cfg)
    coefficients[[length(coefficients) + 1L]] <- rbind(val$coefficients, test$coefficients)
    grids[[length(grids) + 1L]] <- grid
    message(sprintf("h=%d %s: cutoff %.2f", h, link, grid$cutoff[grid$selected]))
  }
}
forecasts <- bind_rows(forecasts)

# Evaluate both cutoff rules and the benchmarks; add bootstrap AUC intervals.
results <- evaluate(by_rule(forecasts), returns, cfg, benchmarks = TRUE)
intervals <- bind_rows(lapply(split(forecasts, list(forecasts$period, forecasts$h, forecasts$model), drop = TRUE),
                              auc_intervals, cfg = cfg))
results$forecast_metrics <- merge(results$forecast_metrics, intervals, by = c("period", "h", "model"), sort = FALSE)

out <- results_dir("baseline")
write_tables(out, c(list(forecasts = forecasts, cutoff_grid = bind_rows(grids),
                         coefficients = bind_rows(coefficients)), results))
invisible(draw_wealth_figure(out, results$monthly_portfolios, cfg))
invisible(draw_probability_figure(out, forecasts, cfg))
write_run_info(out, cfg, checks)
message("Results saved in ", out)
