# Random forest under the supervisor protocol (GROUP_MODEL_STANDARDS_DRAFT.md).
# Run from the project root (or anywhere, with the absolute path):
#   Rscript --vanilla Scripts/RF/run_rf.R [--stage=all|validation|test|tests]
# Optional environment: EC4308_RF_WORKERS (default cores-2; 1 = serial, identical results),
#                       EC4308_ROOT (project root override).
args <- commandArgs(FALSE)
entry <- grep("^--file=", args, value = TRUE)
entry_path <- if (length(entry)) gsub("~+~", " ", sub("^--file=", "", entry[1]), fixed = TRUE) else ""
script_dir <- if (nzchar(entry_path)) dirname(normalizePath(entry_path, mustWork = TRUE)) else normalizePath("Scripts/RF")
find_project_root <- function(start) {
  env <- Sys.getenv("EC4308_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  dir <- start
  for (k in 1:6) {
    if (dir.exists(file.path(dir, "Data")) && dir.exists(file.path(dir, "Scripts"))) return(dir)
    dir <- dirname(dir)
  }
  stop("Cannot locate the project root; set EC4308_ROOT.")
}
root <- find_project_root(script_dir)
.libPaths(c(file.path(root, ".R-library"), .libPaths()))
ensure_packages <- c("digest", "jsonlite", "ranger")
missing <- ensure_packages[!vapply(ensure_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
  ". Install into .R-library, e.g. install.packages('ranger', lib='.R-library').")

# One shared portfolio evaluator (rule 33): reuse the baseline's helpers and engine verbatim.
# Scripts/RF and Scripts/Baseline models are sibling folders.
baseline_dir <- file.path(script_dir, "..", "Baseline models")
if (!dir.exists(baseline_dir)) stop("Shared baseline code not found at ", baseline_dir)
baseline_dir <- normalizePath(baseline_dir)
shared_files <- file.path(baseline_dir, c("config.R", "portfolio.R"))
for (file in shared_files) source(file)
rf_files <- file.path(script_dir, c("config.R", "models.R", "tests.R", "report.R"))
for (file in rf_files) source(file)
cfg <- rf_config()
stage <- sub("^--stage=", "", grep("^--stage=", commandArgs(TRUE), value = TRUE))
if (!length(stage)) stage <- "all"
ensure(length(stage) == 1L && stage %in% c("all", "validation", "test", "tests"), "Use --stage=all|validation|test|tests")
workers <- as.integer(Sys.getenv("EC4308_RF_WORKERS", max(1L, parallel::detectCores() - 2L)))
out <- file.path(root, "Results", "RF")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
started <- proc.time()[["elapsed"]]

input_files <- c(file.path(root, "Data/processed", paste0("model_ready_h", cfg$horizons, ".rds")),
                 file.path(root, "Data/crsp_monthly_index_data/crsp_monthly_index_data.csv"),
                 file.path(root, "Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))
hashes <- setNames(vapply(input_files, file_hash, character(1)), substring(input_files, nchar(root) + 2L))
code_files <- c(sort(list.files(script_dir, pattern = "[.]R$", full.names = TRUE)), shared_files)
code_hashes <- setNames(vapply(code_files, file_hash, character(1)),
                        c(basename(code_files[seq_len(length(code_files) - 2L)]), paste0("shared/", basename(shared_files))))

message("Checking information sets, determinism and the shared portfolio engine...")
tests <- run_rf_tests(root, cfg)
write_csv(tests, file.path(out, "tests.csv"))
message(nrow(tests), " focused test groups passed.")
if (stage == "tests") quit(save = "no", status = 0)

returns <- load_return_data(root)
datasets <- setNames(lapply(cfg$horizons, function(h) rf_load_data(root, h, cfg)), cfg$horizons)
group_names <- names(datasets[[1]]$groups)
bind <- function(x) do.call(rbind, x)
message("Using ", workers, " worker(s).")
cl <- rf_make_cluster(workers, c(file.path(baseline_dir, "config.R"), file.path(script_dir, "models.R")), .libPaths())

# Serial re-run of selected tasks in this process must reproduce the parallel output bit for bit.
spot_check <- function(d, tasks, forecasts) {
  pick <- unique(c(1L, length(tasks)))
  again <- vapply(rf_run_tasks(NULL, d, tasks[pick]), `[[`, 0, "p")
  ensure(identical(again, forecasts$probability[pick]), "Parallel and serial forecasts differ")
}

validation_file <- file.path(out, "validation_state.rds")
if (stage %in% c("all", "validation")) {
  message("Monthly expanding validation forecasts: ", nrow(cfg$hyper_grid), " hyperparameter candidates x ",
          length(cfg$horizons), " horizons...")
  candidate_forecasts <- scores <- forecasts <- grids <- selections <- perturbed <- list()
  for (h in cfg$horizons) {
    d <- datasets[[as.character(h)]]
    rows <- rf_period_rows(d, cfg, "validation")
    donors <- rf_donors(rows, h, "validation", cfg)
    thetas <- lapply(cfg$hyper_grid$config_id, function(id) rf_theta(cfg, id))
    task_sets <- lapply(thetas, function(theta) rf_make_tasks(d, rows, h, theta, cfg, donors))
    results <- rf_run_tasks(cl, d, unlist(task_sets, recursive = FALSE))
    chunks <- split(results, rep(seq_along(task_sets), each = length(rows)))
    assembled <- lapply(seq_along(thetas), function(k)
      rf_assemble(d, rows, chunks[[k]], h, "validation", thetas[[k]], cfg))
    # Stage 1 (rules 25-26): every candidate's 96 expanding validation forecasts, scored by Brier.
    sc <- bind(lapply(seq_along(assembled), function(k)
      cbind(data.frame(h = h), cfg$hyper_grid[k, ], rf_probability_scores(assembled[[k]]$forecasts))))
    best <- rf_choose_candidate(sc, cfg$brier_tie_tolerance)
    sc$selected <- seq_len(nrow(sc)) == best
    sc$Brier_rank <- rank(sc$Brier, ties.method = "min")
    theta <- thetas[[best]]
    spot_check(d, task_sets[[best]], assembled[[best]]$forecasts)
    f <- assembled[[best]]$forecasts
    # Stage 2 (rules 27-28): with the winner's SAVED validation forecasts fixed, choose the cutoff.
    inp <- portfolio_inputs(f, returns, h, cfg$validation_end, cfg$common_rf_security)
    grid <- bind(lapply(cfg$cutoffs, function(cutoff) {
      p <- portfolio_path(inp, as.integer(f$probability >= cutoff), cfg$primary_cost, cfg)
      cbind(data.frame(h = h, model = cfg$model, config_id = theta$config_id, cutoff = cutoff),
            portfolio_metrics(p, cfg$risk_aversion))
    }))
    winner <- grid[choose_cutoff(grid, cfg$cutoff_objective), ]
    default <- grid[abs(grid$cutoff - cfg$comparison_cutoff) < 1e-12, ]
    ensure(nrow(default) == 1L, "Comparison cutoff missing from validation grid")
    default_row <- sc[sc$is_ranger_default, ]
    selections[[as.character(h)]] <- data.frame(h = h, model = cfg$model, config_id = theta$config_id,
      mtry = theta$mtry, min_node_size = theta$min_node_size, replace = theta$replace,
      sample_fraction = theta$sample_fraction, num_trees = theta$num_trees,
      validation_Brier = sc$Brier[best], validation_Brier_ranger_default = default_row$Brier,
      validation_Brier_train_prevalence = sc$Brier_train_prevalence[best],
      cutoff = winner$cutoff, selection_objective = cfg$cutoff_objective,
      validation_Sharpe = winner$sharpe_matched_ladder,
      validation_sharpe_matched_ladder = winner$sharpe_matched_ladder,
      validation_sharpe_common_rf = winner$sharpe_common_rf,
      fixed_050_validation_sharpe_matched_ladder = default$sharpe_matched_ladder,
      fixed_050_validation_sharpe_common_rf = default$sharpe_common_rf,
      validation_CAGR = winner$CAGR, cost_rate = cfg$primary_cost,
      validation_start = winner$start, validation_end = winner$end)
    f$cutoff <- winner$cutoff; f$signal <- as.integer(f$probability >= f$cutoff)
    f$cutoff_050 <- cfg$comparison_cutoff; f$signal_050 <- as.integer(f$probability >= f$cutoff_050)
    forecasts[[as.character(h)]] <- f; grids[[as.character(h)]] <- grid; scores[[as.character(h)]] <- sc
    candidate_forecasts[[as.character(h)]] <- bind(lapply(assembled, `[[`, "forecasts"))
    perturbed[[as.character(h)]] <- assembled[[best]]$perturbed
    message(sprintf("h=%d: selected %s (validation Brier %.4f; ranger default %.4f), cutoff %.2f",
                    h, theta$config_id, sc$Brier[best], default_row$Brier, winner$cutoff))
  }
  selected <- bind(selections); grid <- bind(grids); val_forecasts <- bind(forecasts)
  hyper_scores <- bind(scores); candidate_forecasts <- bind(candidate_forecasts)
  rownames(selected) <- rownames(grid) <- rownames(val_forecasts) <- rownames(hyper_scores) <-
    rownames(candidate_forecasts) <- NULL
  # Frozen before the first test-period forecast is generated.
  frozen <- list(config = cfg, selected = selected, input_hashes = hashes, code_hashes = code_hashes,
                 feature_names = colnames(datasets[[1]]$X),
                 frozen_at = format(Sys.time(), tz = "Asia/Singapore", usetz = TRUE))
  saveRDS(frozen, file.path(out, "frozen_specification.rds"))
  state <- list(forecasts = val_forecasts, grid = grid, selected = selected, hyper_scores = hyper_scores,
                candidate_forecasts = candidate_forecasts, perturbed = perturbed, frozen = frozen)
  saveRDS(state, validation_file)
  write_csv(hyper_scores, file.path(out, "validation_hyperparameter_scores.csv"))
  write_csv(candidate_forecasts, file.path(out, "validation_candidate_forecasts.csv"))
  write_csv(grid, file.path(out, "validation_cutoff_grid.csv"))
  write_csv(selected, file.path(out, "selected_cutoffs.csv"))
  write_csv(val_forecasts, file.path(out, "validation_forecasts.csv"))
  write_csv(data.frame(feature = colnames(datasets[[1]]$X),
                       base_variable = sub("_lag[0-5]$", "", colnames(datasets[[1]]$X))),
            file.path(out, "feature_names.csv"))
  message("Validation choices frozen.")
} else {
  ensure(file.exists(validation_file), "Run --stage=validation before --stage=test")
  state <- readRDS(validation_file)
  ensure(identical(state$frozen$config, cfg) && identical(state$frozen$input_hashes, hashes) &&
           identical(state$frozen$code_hashes, code_hashes),
         "Configuration, code or inputs changed since validation. Re-run validation explicitly.")
  val_forecasts <- state$forecasts; selected <- state$selected; grid <- state$grid
  hyper_scores <- state$hyper_scores; candidate_forecasts <- state$candidate_forecasts
}
if (stage == "validation") {
  if (!is.null(cl)) parallel::stopCluster(cl)
  message("Validation complete; final-test fitting and performance have not run.")
  quit(save = "no", status = 0)
}

message("Monthly expanding test forecasts with frozen hyperparameters and cutoffs...")
test_forecasts <- test_perturbed <- list()
for (h in cfg$horizons) {
  d <- datasets[[as.character(h)]]
  s <- selected[selected$h == h, ]
  theta <- rf_theta(cfg, s$config_id)
  rows <- rf_period_rows(d, cfg, "test")
  tasks <- rf_make_tasks(d, rows, h, theta, cfg, rf_donors(rows, h, "test", cfg))
  a <- rf_assemble(d, rows, rf_run_tasks(cl, d, tasks), h, "test", theta, cfg)
  spot_check(d, tasks, a$forecasts)
  f <- a$forecasts
  f$cutoff <- s$cutoff; f$signal <- as.integer(f$probability >= f$cutoff)
  f$cutoff_050 <- cfg$comparison_cutoff; f$signal_050 <- as.integer(f$probability >= f$cutoff_050)
  test_forecasts[[as.character(h)]] <- f; test_perturbed[[as.character(h)]] <- a$perturbed
}
if (!is.null(cl)) parallel::stopCluster(cl)
test_forecasts <- bind(test_forecasts); rownames(test_forecasts) <- NULL
all_forecasts <- rbind(val_forecasts, test_forecasts)
audit_cols <- c("h", "model", "period", "origin", "config_id", "mtry", "min_node_size", "replace",
                "sample_fraction", "num_trees", "seed", "n_train", "last_train_origin",
                "last_train_outcome", "train_prevalence")
audit <- rbind(candidate_forecasts[, audit_cols], test_forecasts[, audit_cols])
names(audit)[names(audit) == "origin"] <- "forecast_origin"
audit$selected_config <- audit$period == "test" |
  paste(audit$h, audit$config_id) %in% paste(selected$h, selected$config_id)
rownames(audit) <- NULL
ensure(all(month_number(audit$last_train_outcome) <= month_number(audit$forecast_origin)), "Final timing audit failed")
ensure(all(month_number(audit$last_train_origin) <= month_number(audit$forecast_origin) - audit$h), "Immature label in fit")
ensure(sum(audit$period == "test") == 156L * length(cfg$horizons), "Missing monthly test refits")
ensure(sum(audit$period == "validation") == 96L * length(cfg$horizons) * nrow(cfg$hyper_grid), "Missing validation refits")

paths <- performance <- diagnostics <- calibration <- benchmark_returns <- list()
for (period in c("validation", "test")) for (h in cfg$horizons) {
  f <- all_forecasts[all_forecasts$period == period & all_forecasts$h == h, ]
  end <- if (period == "validation") cfg$validation_end else cfg$test_end
  inp <- portfolio_inputs(f, returns, h, end, cfg$common_rf_security)
  benchmark_returns[[paste(period, h)]] <- data.frame(period = period, h = h, month = inp$ym,
    quote_start = inp$quote_start, quote_end = inp$quote_end, quote_days = inp$days,
    common_rf_security = inp$common_rf_security, matched_bill_security = inp$matched_bill_security,
    common_rf_daily_yield = inp$common_rf_daily_yield,
    common_rf_return = inp$common_rf_return, matched_ladder_return = inp$bill_reference)
  for (strategy in c(cfg$model, "all_stocks", "all_bills")) {
    if (strategy == cfg$model) {
      variants <- list(validation_tuned = list(signals = f$signal, cutoff = unique(f$cutoff)),
                       fixed_050 = list(signals = f$signal_050, cutoff = cfg$comparison_cutoff))
      ps <- rf_probability_scores(f)
      for (rule in names(variants)) {
        diagnostics[[paste(period, h, rule)]] <- cbind(
          data.frame(period = period, h = h, model = strategy, config_id = unique(f$config_id), cutoff_rule = rule),
          forecast_metrics(f, variants[[rule]]$cutoff),
          data.frame(Brier_train_prevalence = ps$Brier_train_prevalence,
                     Brier_skill_vs_prevalence = 1 - ps$Brier / ps$Brier_train_prevalence))
      }
      eligible <- f[f$classification_eligible & !is.na(f$target), ]
      bins <- cut(eligible$probability, breaks = seq(0, 1, by = .1), include.lowest = TRUE)
      calibration[[paste(period, h)]] <- bind(lapply(levels(bins), function(bin) {
        ids <- which(bins == bin)
        data.frame(period = period, h = h, model = strategy, bin = bin, n = length(ids),
          mean_probability = if (length(ids)) mean(eligible$probability[ids]) else NA_real_,
          event_rate = if (length(ids)) mean(eligible$target[ids]) else NA_real_)
      }))
    } else variants <- list(not_applicable = list(signals = rep(as.integer(strategy == "all_stocks"), nrow(f)), cutoff = NA_real_))
    for (rule in names(variants)) for (cost in cfg$costs) {
      key <- paste(period, h, strategy, rule, cost, sep = "_")
      p <- portfolio_path(inp, variants[[rule]]$signals, cost, cfg)
      tags <- data.frame(period = period, h = h, strategy = strategy, cutoff_rule = rule,
                         cutoff = variants[[rule]]$cutoff, cost_rate = cost)
      paths[[key]] <- cbind(tags, p)
      performance[[key]] <- cbind(tags, portfolio_metrics(p, cfg$risk_aversion))
      ensure(max(month_number(p$month)) <= month_number(end), "Portfolio return crosses split boundary")
    }
  }
}
paths <- bind(paths); performance <- bind(performance); diagnostics <- bind(diagnostics)
rownames(paths) <- rownames(performance) <- rownames(diagnostics) <- NULL
common <- list()
for (period in c("validation", "test")) {
  period_start <- if (period == "validation") cfg$validation_start else cfg$test_start
  common_start <- month_string(month_number(period_start) + max(cfg$horizons))
  p <- paths[paths$period == period & paths$month >= common_start, ]
  groups <- interaction(p$h, p$strategy, p$cutoff_rule, p$cost_rate, drop = TRUE)
  for (key in unique(groups)) {
    q <- p[groups == key, ]
    common[[paste(period, key)]] <- cbind(q[1, c("period", "h", "strategy", "cutoff_rule", "cutoff", "cost_rate")],
                                          portfolio_metrics(q, cfg$risk_aversion))
  }
}

message("Interpretation: out-of-sample grouped permutation, final-fit OOB importance, partial dependence...")
perturbed_all <- list(validation = state$perturbed, test = test_perturbed)
importance_oos <- bind(lapply(c("validation", "test"), function(period) bind(lapply(cfg$horizons, function(h) {
  f <- all_forecasts[all_forecasts$period == period & all_forecasts$h == h, ]
  imp <- rf_permutation_importance(f, perturbed_all[[period]][[as.character(h)]], group_names, unique(f$cutoff))
  imp <- imp[order(-imp$delta_Brier), ]
  cbind(data.frame(period = period, h = h, config_id = unique(f$config_id)), imp,
        rank_delta_Brier = seq_len(nrow(imp)),
        rank_reliance = rank(-imp$mean_abs_probability_change, ties.method = "min"))
}))))
final_fits <- validation_fits <- list(); importance_oob <- pd <- list()
for (h in cfg$horizons) {
  d <- datasets[[as.character(h)]]
  s <- selected[selected$h == h, ]
  for (period in c("validation", "test")) {
    origin_ym <- if (period == "validation") cfg$validation_end else cfg$test_end
    i <- which(d$ym == origin_ym); rows <- rf_mature_rows(d$y, d$m, d$end, d$m[i])
    seed <- rf_seed(cfg$seed_base, h, d$m[i])
    theta <- rf_theta(cfg, s$config_id)
    fit <- rf_fit(d$X[rows, , drop = FALSE], d$y[rows], theta, seed)
    f <- all_forecasts[all_forecasts$period == period & all_forecasts$h == h & all_forecasts$origin == origin_ym, ]
    ensure(identical(rf_prob(fit, d$X[i, , drop = FALSE]), f$probability), "Final-origin refit is not reproducible")
    if (period == "validation") validation_fits[[as.character(h)]] <- fit else final_fits[[as.character(h)]] <- fit
    imp_fit <- rf_fit(d$X[rows, , drop = FALSE], d$y[rows], theta, seed, importance = "permutation")
    vi <- imp_fit$variable.importance
    agg <- tapply(vi, factor(sub("_lag[0-5]$", "", names(vi)), levels = group_names), sum)
    importance_oob[[paste(period, h)]] <- data.frame(h = h, fit_origin = origin_ym, period = period,
      n_train = length(rows), base_variable = names(agg), oob_permutation_importance = as.numeric(agg),
      rank = rank(-as.numeric(agg), ties.method = "min"))
    if (period == "test") {
      reliance <- importance_oos[importance_oos$period == "test" & importance_oos$h == h, ]
      top_vars <- reliance$base_variable[order(reliance$rank_reliance)][1:4]
      xtrain <- d$X[rows, , drop = FALSE]
      # "Sustained level" partial dependence: all six lags of the variable are set to the
      # same value (lags of these slow-moving variables are highly correlated, so moving
      # lag0 alone understates what the forest responds to). Grid = lag0 training quantiles.
      for (v in top_vars) {
        cols <- d$groups[[v]]
        qs <- seq(.05, .95, by = .05)
        values <- stats::quantile(xtrain[, paste0(v, "_lag0")], qs, names = FALSE)
        keep <- !duplicated(values)
        values <- values[keep]
        pd_value <- function(x, which_cols) { z <- xtrain; z[, which_cols] <- x; mean(rf_prob(fit, z)) }
        pd[[paste(h, v)]] <- data.frame(h = h, fit_origin = origin_ym, base_variable = v,
          feature = paste0(v, "_lag0..lag5"), value = values, quantile = qs[keep],
          partial_dependence = vapply(values, pd_value, numeric(1), which_cols = cols),
          partial_dependence_lag0_only = vapply(values, pd_value, numeric(1), which_cols = paste0(v, "_lag0")))
      }
    }
  }
}
importance_oob <- bind(importance_oob); pd <- bind(pd)
rownames(importance_oob) <- rownames(pd) <- rownames(importance_oos) <- NULL

ensure(identical(hashes, setNames(vapply(input_files, file_hash, character(1)), names(hashes))), "Input files changed during run")
write_csv(test_forecasts, file.path(out, "test_forecasts.csv"))
write_csv(paths, file.path(out, "monthly_portfolios.csv"))
write_csv(bind(benchmark_returns), file.path(out, "benchmark_returns.csv"))
write_csv(performance, file.path(out, "performance.csv"))
write_csv(bind(common), file.path(out, "common_window_performance.csv"))
write_csv(diagnostics, file.path(out, "forecast_metrics.csv"))
write_csv(bind(calibration), file.path(out, "calibration.csv"))
write_csv(audit, file.path(out, "refit_audit.csv"))
write_csv(importance_oos, file.path(out, "importance_out_of_sample.csv"))
write_csv(importance_oob, file.path(out, "importance_oob_final_fits.csv"))
write_csv(pd, file.path(out, "partial_dependence.csv"))
write_csv(data.frame(file = names(hashes), sha256 = unname(hashes)), file.path(out, "input_hashes.csv"))
saveRDS(validation_fits, file.path(out, "validation_models.rds"))
saveRDS(final_fits, file.path(out, "final_models.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
write_csv(data.frame(package = ensure_packages, version = vapply(ensure_packages, function(p) as.character(packageVersion(p)), character(1))),
          file.path(out, "package_versions.csv"))
runtime <- proc.time()[["elapsed"]] - started
manifest <- list(version = cfg$version, model = cfg$model, completed_at = format(Sys.time(), tz = "Asia/Singapore", usetz = TRUE),
  config = cfg, feature_names = colnames(datasets[[1]]$X), input_hashes = hashes, code_hashes = code_hashes,
  tests_passed = nrow(tests), workers = workers, runtime_seconds = round(runtime, 1),
  monthly_validation_refits = sum(audit$period == "validation"), monthly_test_refits = sum(audit$period == "test"),
  validation_selected = selected, test_used_for_selection = FALSE, input_hashes_unchanged = TRUE,
  parallel_serial_identical = TRUE, final_origin_refits_reproduced = TRUE,
  protocol = "Group model standards (supervisor revision of 5 November 2026); RF-specific cutoff grid 0.00-1.00 by 0.05")
jsonlite::write_json(manifest, file.path(out, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
write_rf_report(out, cfg, selected, hyper_scores, performance, diagnostics, importance_oos, tests)
draw_rf_figures(out, paths, grid, selected, hyper_scores, all_forecasts, importance_oos, pd, cfg)
message(sprintf("Complete in %.1f s. Results: %s", runtime, file.path(out, "RESULTS.md")))
