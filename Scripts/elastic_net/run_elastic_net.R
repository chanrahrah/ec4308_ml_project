# Elastic net: tune (alpha, lambda) on validation, freeze them and the cutoff, then forecast the test period.
# Run: Rscript --vanilla Scripts/elastic_net/run_elastic_net.R   (or open and click Source)
script <- c(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE),
            tryCatch(sys.frame(1)$ofile, error = function(e) NULL))
source(file.path(if (length(script)) dirname(script[1]) else "Scripts/elastic_net", "..", "utils", "utils.R"))
for (file in c("model.R", "tests.R")) source(file.path(SCRIPTS, "elastic_net", file))
use_packages("glmnet")

cfg <- elastic_config()
returns <- load_return_data()
checks <- run_checks(c(shared_checks(returns, cfg), elastic_checks(cfg)))

forecasts <- tuning <- cutoff_grid <- list()
for (h in cfg$horizons) {
  d <- load_model_data(h, cfg)
  tuned <- tune_elastic(d, h, cfg)
  cutoffs <- select_cutoff(tuned$forecasts, returns, h, cfg)
  test <- recursive_elastic(d, h, cfg, "test", tuned$alpha, tuned$lambda)
  forecasts[[length(forecasts) + 1L]] <- add_signals(rbind(tuned$forecasts, test$forecasts), cutoffs$cutoff[cutoffs$selected], cfg)
  tuning[[length(tuning) + 1L]] <- tuned$scores
  cutoff_grid[[length(cutoff_grid) + 1L]] <- cutoffs
  message(sprintf("h=%d: alpha %.2f, lambda %g, cutoff %.2f", h, tuned$alpha, tuned$lambda, cutoffs$cutoff[cutoffs$selected]))
}

forecasts <- bind_rows(forecasts)
results <- evaluate(by_rule(forecasts), returns, cfg)
out <- results_dir("elastic_net")
write_tables(out, c(list(forecasts = forecasts, tuning = bind_rows(tuning), cutoff_grid = bind_rows(cutoff_grid)), results))
write_run_info(out, cfg, checks, "glmnet")
message("Results saved in ", out)
