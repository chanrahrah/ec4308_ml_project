# Random forest: tune on validation, freeze the settings and cutoff, then forecast the test period.
# Run: Rscript --vanilla Scripts/rf/run_rf.R   (EC4308_RF_WORKERS sets the number of parallel workers)
script <- c(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE),
            tryCatch(sys.frame(1)$ofile, error = function(e) NULL))
source(file.path(if (length(script)) dirname(script[1]) else "Scripts/rf", "..", "utils", "utils.R"))
for (file in c("model.R", "tests.R")) source(file.path(SCRIPTS, "rf", file))
use_packages("ranger")

cfg <- rf_config()
grid <- rf_grid()
returns <- load_return_data()
checks <- run_checks(c(shared_checks(returns, cfg), rf_checks(cfg, grid)))
workers <- as.integer(Sys.getenv("EC4308_RF_WORKERS", max(1L, parallel::detectCores() - 2L)))
cl <- rf_make_cluster(workers)

# Parallel results must match a serial refit of the same origin.
spot_check <- function(d, tasks, f) {
  pick <- unique(c(1L, length(tasks)))
  ensure(identical(vapply(rf_run_tasks(NULL, d, tasks[pick]), `[[`, 0, "p"), f$probability[pick]),
         "Parallel and serial forecasts differ")
}

forecasts <- tuning <- cutoff_grid <- importance <- pd <- list()
for (h in cfg$horizons) {
  d <- rf_load(h, cfg)
  # Validation: every candidate refitted monthly; keep the lowest Brier.
  rows <- period_rows(d, cfg, "validation")
  donors <- rf_donors(rows, h, "validation", cfg)
  thetas <- lapply(grid$config_id, rf_theta, grid = grid, cfg = cfg)
  tasks <- lapply(thetas, function(theta) rf_make_tasks(d, rows, h, theta, cfg, donors))
  results <- split(rf_run_tasks(cl, d, unlist(tasks, recursive = FALSE)), rep(seq_along(tasks), each = length(rows)))
  candidates <- lapply(seq_along(thetas), function(k) rf_assemble(d, rows, results[[k]], h, "validation", thetas[[k]], cfg))
  scores <- cbind(data.frame(h = h), grid, bind_rows(lapply(candidates, function(a) validation_scores(a$forecasts))))
  best <- best_candidate(scores$validation_Brier, order(grid$simplicity_rank))
  scores$selected <- seq_len(nrow(scores)) == best
  val <- candidates[[best]]
  spot_check(d, tasks[[best]], val$forecasts)
  cutoffs <- select_cutoff(val$forecasts, returns, h, cfg)

  # Test: monthly refits with the frozen candidate and cutoff.
  rows_test <- period_rows(d, cfg, "test")
  tasks_test <- rf_make_tasks(d, rows_test, h, thetas[[best]], cfg, rf_donors(rows_test, h, "test", cfg))
  test <- rf_assemble(d, rows_test, rf_run_tasks(cl, d, tasks_test), h, "test", thetas[[best]], cfg)
  spot_check(d, tasks_test, test$forecasts)
  f <- add_signals(rbind(val$forecasts, test$forecasts), cutoffs$cutoff[cutoffs$selected], cfg)

  # Interpretation: permutation importance per period and partial dependence of the final fit.
  imp <- bind_rows(list(rf_importance(f[f$period == "validation", ], val$perturbed, names(d$groups)),
                        rf_importance(f[f$period == "test", ], test$perturbed, names(d$groups))))
  top <- imp[imp$period == "test", ]
  pd[[length(pd) + 1L]] <- rf_partial_dependence(d, h, thetas[[best]], cfg, f, top$base_variable[order(top$rank_reliance)][1:4])
  forecasts[[length(forecasts) + 1L]] <- f
  tuning[[length(tuning) + 1L]] <- scores
  cutoff_grid[[length(cutoff_grid) + 1L]] <- cutoffs
  importance[[length(importance) + 1L]] <- imp
  message(sprintf("h=%d: %s, cutoff %.2f", h, thetas[[best]]$config_id, cutoffs$cutoff[cutoffs$selected]))
}
if (!is.null(cl)) parallel::stopCluster(cl)

forecasts <- bind_rows(forecasts)
results <- evaluate(by_rule(forecasts), returns, cfg)
out <- results_dir("rf")
write_tables(out, c(list(forecasts = forecasts, tuning = bind_rows(tuning), cutoff_grid = bind_rows(cutoff_grid)), results,
                    list(importance = bind_rows(importance), partial_dependence = bind_rows(pd))))
write_run_info(out, cfg, checks, "ranger")
message("Results saved in ", out)
