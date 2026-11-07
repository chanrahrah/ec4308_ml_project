# Random-forest model layer: read the frozen inputs, fit/predict a ranger
# probability forest on mature-label rows only, and run monthly expanding
# forecasts. Every function here is a pure function of its arguments so that
# PSOCK workers (Windows-safe) can source this file and reproduce results exactly.

# No raw predictor reconstruction: consume only the validated RDS inputs (rule 1).
rf_load_data <- function(root, h, cfg) {
  d <- readRDS(file.path(root, "Data", "processed", paste0("model_ready_h", h, ".rds")))
  ensure(ncol(d$X) == cfg$n_features, "Expected the full 132-column predictor set")
  ensure(identical(colnames(d$X), d$feature_names), "Predictor names differ from the frozen feature list")
  ensure(all(grepl("_lag[0-5]$", colnames(d$X))), "Non-lag column in predictor matrix")
  ensure(!any(colnames(d$X) %in% cfg$forbidden_inputs), "Metadata/outcome column in predictor matrix")
  d$m <- month_number(d$ym)
  d$end <- month_number(d$metadata$outcome_end)
  ensure(all(diff(d$m) == 1L) && !anyDuplicated(d$m), "Incomplete model calendar")
  ensure(all(is.finite(d$X)), "Missing/nonfinite predictor: no imputation is specified")
  known <- !is.na(d$y)
  ensure(all(d$y[known] %in% 0:1), "Invalid binary target")
  ensure(all(d$end[known] == d$m[known] + h), "Target horizon mismatch")
  ensure(identical(d$ym, as.character(d$metadata$ym)), "Metadata alignment mismatch")
  ensure(all(d$end[d$rows$train] <= month_number(cfg$train_end)), "Training boundary leak")
  ensure(all(d$end[d$rows$validation] <= month_number(cfg$validation_end)), "Validation boundary leak")
  d$groups <- split(seq_len(ncol(d$X)), factor(sub("_lag[0-5]$", "", colnames(d$X)),
                                               levels = unique(sub("_lag[0-5]$", "", colnames(d$X)))))
  ensure(length(d$groups) == 22L && all(lengths(d$groups) == 6L), "Expected 22 variables x 6 lags")
  d
}

# At the close of origin t, a label enters estimation only once its outcome has ended.
# Boundary-purged rows re-enter once mature; the split column is not used here.
rf_mature_rows <- function(y, m, end, origin) which(!is.na(y) & !is.na(end) & end <= origin & m < origin)

rf_seed <- function(cfg_seed_base, h, origin_m) as.integer(cfg_seed_base + 10L * origin_m + h)

# One deterministic probability forest for candidate theta (a list with mtry,
# min_node_size, replace, sample_fraction, num_trees). num.threads = 1 inside the
# fit; parallelism is across origins only, so serial and parallel runs are identical.
rf_fit <- function(x, y, theta, seed, importance = "none") {
  ensure(!anyNA(y) && all(y %in% 0:1), "Training labels must be 0/1 without NA")
  ensure(length(unique(y)) == 2L, "Training labels contain one class only")
  ranger::ranger(x = x, y = factor(y, levels = c("0", "1")), probability = TRUE,
                 num.trees = theta$num_trees, mtry = theta$mtry, min.node.size = theta$min_node_size,
                 replace = theta$replace, sample.fraction = theta$sample_fraction, splitrule = "gini",
                 importance = importance, seed = seed, num.threads = 1L, verbose = FALSE)
}

# P(Y = 1) i.e. equity outperformance. Class order is checked explicitly (rule 8).
rf_prob <- function(fit, x) {
  pr <- stats::predict(fit, data = x, num.threads = 1L)$predictions
  ensure(identical(colnames(pr), c("0", "1")), "Unexpected class ordering in ranger output")
  p <- unname(pr[, "1"])
  ensure(all(is.finite(p) & p >= 0 & p <= 1), "Invalid forecast probability")
  p
}

# One forecast origin for one candidate. dd holds X, y, m, end and groups; the
# task carries the origin row, seed, candidate theta and optional donor rows.
# Optional donors (row indices) produce grouped-permutation predictions used only
# for out-of-sample interpretation; they never affect the forecast itself.
rf_task <- function(dd, task) {
  i <- task$i
  rows <- rf_mature_rows(dd$y, dd$m, dd$end, dd$m[i])
  ensure(length(rows) > 0L && !(i %in% rows), "Origin row in its own training set")
  ensure(max(dd$end[rows]) <= dd$m[i], "Leakage: training outcome ends after origin")
  fit <- rf_fit(dd$X[rows, , drop = FALSE], dd$y[rows], task$theta, task$seed)
  p <- rf_prob(fit, dd$X[i, , drop = FALSE])
  perturbed <- NULL
  if (length(task$donors)) {
    G <- length(dd$groups); K <- length(task$donors)
    z <- dd$X[rep(i, G * K), , drop = FALSE]
    for (g in seq_len(G)) for (k in seq_len(K))
      z[(g - 1L) * K + k, dd$groups[[g]]] <- dd$X[task$donors[k], dd$groups[[g]]]
    perturbed <- matrix(rf_prob(fit, z), nrow = G, ncol = K, byrow = TRUE)
  }
  ensure(fit$mtry == task$theta$mtry && fit$min.node.size == task$theta$min_node_size &&
           fit$num.trees == task$theta$num_trees, "ranger did not use the requested hyperparameters")
  list(i = i, config_id = task$theta$config_id, p = p, n_train = length(rows),
       last_origin = max(dd$m[rows]), last_outcome = max(dd$end[rows]), seed = task$seed,
       train_prevalence = mean(dd$y[rows]), perturbed = perturbed)
}

rf_task_data <- function(d) list(X = d$X, y = d$y, m = d$m, end = d$end, groups = d$groups)

rf_period_rows <- function(d, cfg, period) {
  ensure(period %in% c("validation", "test"), "Unknown forecast period")
  rows <- which(d$m >= month_number(cfg[[paste0(period, "_start")]]) &
                d$m <= month_number(cfg[[paste0(period, "_end")]]))
  ensure(length(rows) > 0L, "No forecast origins")
  rows
}

# Donor rows for permutation interpretation: K seeded permutations of the same
# period's origins. Fixed before fitting and identical across grid candidates.
rf_donors <- function(rows, h, period, cfg) {
  set.seed(cfg$seed_base + 100L * h + if (period == "validation") 1L else 2L)
  t(vapply(seq_len(cfg$importance_reps), function(k) rows[sample.int(length(rows))], integer(length(rows))))
}

rf_make_tasks <- function(d, rows, h, theta, cfg, donors = NULL) {
  lapply(seq_along(rows), function(j) list(
    i = rows[j], theta = theta, seed = rf_seed(cfg$seed_base, h, d$m[rows[j]]),
    donors = if (is.null(donors)) integer() else donors[, j]))
}

rf_make_cluster <- function(workers, source_files, lib_paths) {
  if (workers <= 1L) return(NULL)
  cl <- parallel::makePSOCKcluster(workers)
  parallel::clusterCall(cl, function(libs, files) {
    .libPaths(libs); for (f in files) source(f); NULL
  }, lib_paths, source_files)
  cl
}

# Run tasks serially (cl = NULL) or over a PSOCK cluster. Each task is seeded by
# (h, origin), so scheduling cannot change results.
rf_run_tasks <- function(cl, d, tasks) {
  dd <- rf_task_data(d)
  if (is.null(cl)) return(lapply(tasks, function(task) rf_task(dd, task)))
  parallel::clusterCall(cl, function(v) { assign("RF_DATA", v, envir = globalenv()); NULL }, dd)
  run <- function(task) rf_task(RF_DATA, task)
  environment(run) <- globalenv()   # resolve RF_DATA on the worker; do not reserialize per task
  out <- parallel::parLapplyLB(cl, tasks, run)
  parallel::clusterCall(cl, function() { rm("RF_DATA", envir = globalenv()); NULL })
  out
}

# Assemble one candidate's ordered per-origin results into the shared forecast schema.
rf_assemble <- function(d, rows, results, h, period, theta, cfg) {
  ensure(identical(vapply(results, `[[`, 0L, "i"), rows), "Forecast results out of origin order")
  ensure(all(vapply(results, `[[`, "", "config_id") == theta$config_id), "Results mixed across candidates")
  last_origin <- vapply(results, `[[`, 0L, "last_origin")
  last_outcome <- vapply(results, `[[`, 0L, "last_outcome")
  ensure(all(last_outcome <= d$m[rows]) && all(last_origin <= d$m[rows] - h),
         "Expanding-window outcome leakage")
  forecasts <- data.frame(h = h, model = cfg$model, period = period, origin = d$ym[rows],
    probability = vapply(results, `[[`, 0, "p"),
    outcome_end = as.character(d$metadata$outcome_end[rows]), target = d$y[rows],
    classification_eligible = d$split[rows] == period,
    fit_information_cutoff = d$ym[rows],
    last_train_origin = month_string(last_origin), last_train_outcome = month_string(last_outcome),
    n_train = vapply(results, `[[`, 0L, "n_train"), config_id = theta$config_id,
    mtry = theta$mtry, min_node_size = theta$min_node_size, replace = theta$replace,
    sample_fraction = theta$sample_fraction, num_trees = theta$num_trees,
    seed = vapply(results, `[[`, 0L, "seed"),
    train_prevalence = vapply(results, `[[`, 0, "train_prevalence"), stringsAsFactors = FALSE)
  perturbed <- if (!is.null(results[[1]]$perturbed))
    simplify2array(lapply(results, `[[`, "perturbed")) else NULL   # G x K x n_origins
  list(forecasts = forecasts, perturbed = perturbed)
}

# Brier/log loss/AUC on eligible completed labels only (rule 22).
rf_probability_scores <- function(f) {
  f <- f[f$classification_eligible & !is.na(f$target), ]
  y <- f$target; p <- f$probability
  pc <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  auc <- if (n1 * n0 > 0) (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) else NA_real_
  data.frame(n_eligible = length(y), Brier = mean((p - y)^2),
             log_loss = -mean(y * log(pc) + (1 - y) * log1p(-pc)), AUC = auc,
             Brier_train_prevalence = mean((f$train_prevalence - y)^2),
             mean_probability = mean(p), event_rate = mean(y))
}

# Stage 1 (rules 25-26): lowest validation Brier; exact ties go to the simpler candidate.
rf_choose_candidate <- function(scores, tol = 1e-12) {
  ensure(all(is.finite(scores$Brier)), "Undefined validation Brier")
  tied <- which(scores$Brier - min(scores$Brier) <= tol)
  tied[order(scores$simplicity_rank[tied])][1L]
}

# Out-of-sample grouped permutation importance for one period (interpretation only).
rf_permutation_importance <- function(f, perturbed, group_names, cutoff) {
  ok <- f$classification_eligible & !is.na(f$target)
  base <- mean((f$probability[ok] - f$target[ok])^2)
  do.call(rbind, lapply(seq_along(group_names), function(g) {
    pk <- perturbed[g, , , drop = TRUE]          # K x n
    if (is.null(dim(pk))) pk <- matrix(pk, nrow = 1L)
    brier_k <- apply(pk[, ok, drop = FALSE], 1L, function(p) mean((p - f$target[ok])^2))
    data.frame(base_variable = group_names[g], delta_Brier = mean(brier_k) - base,
               delta_Brier_sd_across_reps = stats::sd(brier_k),
               mean_abs_probability_change = mean(abs(sweep(pk, 2L, f$probability, "-"))),
               signal_flip_rate_at_tuned_cutoff = mean(sweep(pk >= cutoff, 2L, f$probability >= cutoff, "!=")))
  }))
}
