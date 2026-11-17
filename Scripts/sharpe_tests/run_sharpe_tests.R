# Ledoit-Wolf tests of equal Sharpe ratios: each model's strategy against all-stocks within each horizon,
# on the test period at 10 bp, using excess returns over the matched bill ladder.
# Run after all models: Rscript --vanilla Scripts/sharpe_tests/run_sharpe_tests.R
script <- c(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE),
            tryCatch(sys.frame(1)$ofile, error = function(e) NULL))
source(file.path(if (length(script)) dirname(script[1]) else "Scripts/sharpe_tests", "..", "utils", "utils.R"))
for (file in c("model.R", "tests.R")) source(file.path(SCRIPTS, "sharpe_tests", file))

cfg <- sharpe_test_config()
checks <- run_checks(sharpe_test_checks(cfg))
files <- file.path(ROOT, "Results", cfg$members, "monthly_portfolios.csv")
ensure(all(file.exists(files)), "Run every model before the Sharpe tests")
paths <- bind_rows(lapply(files, function(f) {
  p <- utils::read.csv(f, stringsAsFactors = FALSE)
  p[p$period == "test" & p$cost_rate == cfg$primary_cost, ]
}))
paths$excess <- paths$net_return - paths$matched_ladder_return

# Compare one strategy with all-stocks over the same months.
compare <- function(h, s, rule) {
  x <- paths[paths$h == h & paths$strategy == s & paths$cutoff_rule == rule, ]
  y <- paths[paths$h == h & paths$strategy == "all_stocks", ]
  y <- y[match(x$month, y$month), ]
  ensure(nrow(x) > 0 && !anyNA(y$excess), paste("Missing months for", s, "h", h))
  row <- data.frame(h = h, strategy = s, cutoff_rule = rule, n_months = nrow(x),
    months_different = sum(abs(x$net_return - y$net_return) > 1e-12),
    sharpe = annual_sharpe(x$excess), sharpe_all_stocks = annual_sharpe(y$excess),
    sharpe_diff = 0, se = NA_real_, p_value = NA_real_, status = "identical to all-stocks")
  if (isTRUE(all.equal(x$net_return, y$net_return, tolerance = 1e-12))) return(row)
  test <- lw_test(cbind(x$excess, y$excess), cfg$block_length, cfg$replicates, cfg$seed)
  # Annualise with the sample-SD convention used in performance.csv.
  scale <- sqrt(12 * (nrow(x) - 1) / nrow(x))
  near(test$diff * scale, row$sharpe - row$sharpe_all_stocks)
  row$sharpe_diff <- test$diff * scale; row$se <- test$se * scale
  row$p_value <- test$p_value; row$status <- "tested"
  row
}

# Tuned models at both cutoff rules; the probability ensembles use 0.50 only.
strategies <- rbind(
  expand.grid(strategy = c("random_forest", "elastic_net", "boosting", "ens_majority_vote"),
              cutoff_rule = c("validation_tuned", "fixed_050"), stringsAsFactors = FALSE),
  data.frame(strategy = c("ens_equal", "ens_gr", "ens_gr_yearly"), cutoff_rule = "fixed_050"))
results <- bind_rows(lapply(cfg$horizons, function(h)
  bind_rows(lapply(seq_len(nrow(strategies)), function(k) compare(h, strategies$strategy[k], strategies$cutoff_rule[k])))))
tested <- results$status == "tested"
results$p_holm <- NA_real_
results$p_holm[tested] <- stats::p.adjust(results$p_value[tested], method = "holm")

out <- results_dir("sharpe_tests")
write_tables(out, list(sharpe_tests = results[order(results$h, results$cutoff_rule, results$strategy), ]))
write_run_info(out, cfg, checks, inputs = files)
message(sum(tested), " comparisons tested; ", sum(results$p_value < .05, na.rm = TRUE), " with p < 0.05 (",
        sum(results$p_holm < .05, na.rm = TRUE), " after Holm).")
message("Results saved in ", out)
