# Gradient boosting: tune on validation, freeze the settings and cutoff, then forecast the test period.
# Run: Rscript --vanilla Scripts/boosting/run_boosting.R   (or open and click Source)
script <- c(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE),
            tryCatch(sys.frame(1)$ofile, error = function(e) NULL))
source(file.path(if (length(script)) dirname(script[1]) else "Scripts/boosting", "..", "utils", "utils.R"))
for (file in c("model.R", "tests.R")) source(file.path(SCRIPTS, "boosting", file))
use_packages("xgboost")

cfg <- boosting_config()
grid <- boosting_grid(cfg)
returns <- load_return_data()
checks <- run_checks(c(shared_checks(returns, cfg), boosting_checks(cfg)))

forecasts <- tuning <- cutoff_grid <- list()
for (h in cfg$horizons) {
  d <- load_model_data(h, cfg)
  # Validation: every configuration; keep the lowest Brier (ties: shallower, fewer rounds, slower).
  candidates <- boosting_forecasts(d, h, "validation", grid, cfg)
  scores <- cbind(data.frame(h = h), grid, bind_rows(lapply(candidates, validation_scores)))
  best <- best_candidate(scores$validation_Brier, order(grid$max_depth, grid$nrounds, grid$eta))
  scores$selected <- seq_len(nrow(scores)) == best
  cutoffs <- select_cutoff(candidates[[best]], returns, h, cfg)
  # Test: monthly refits with the frozen configuration and cutoff.
  test <- boosting_forecasts(d, h, "test", grid[best, ], cfg)[[1]]
  forecasts[[length(forecasts) + 1L]] <- add_signals(rbind(candidates[[best]], test), cutoffs$cutoff[cutoffs$selected], cfg)
  tuning[[length(tuning) + 1L]] <- scores
  cutoff_grid[[length(cutoff_grid) + 1L]] <- cutoffs
  message(sprintf("h=%d: %s, cutoff %.2f", h, grid$config_id[best], cutoffs$cutoff[cutoffs$selected]))
}

forecasts <- bind_rows(forecasts)
results <- evaluate(by_rule(forecasts), returns, cfg)
out <- results_dir("boosting")
write_tables(out, c(list(forecasts = forecasts, tuning = bind_rows(tuning), cutoff_grid = bind_rows(cutoff_grid)), results))
write_run_info(out, cfg, checks, "xgboost")
message("Results saved in ", out)
