# Random forest (ranger probability forest) on all 132 predictors.
rf_config <- function() modifyList(base_config(), list(
  model = "random_forest", num_trees = 1000L, importance_reps = 5L))

# 18 candidates, ordered simplest first (larger nodes, fewer split variables, bootstrap) for Brier ties.
rf_grid <- function() {
  g <- expand.grid(mtry = c(4L, 11L, 33L), min_node_size = c(10L, 30L, 60L),
                   sampling = c("boot", "sub50"), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g$replace <- g$sampling == "boot"
  g$sample_fraction <- ifelse(g$replace, 1, 0.5)
  g <- g[order(-g$min_node_size, g$mtry, !g$replace), ]
  g$config_id <- sprintf("rf_mtry%02d_node%02d_%s", g$mtry, g$min_node_size, g$sampling)
  g$simplicity_rank <- seq_len(nrow(g))
  g$is_ranger_default <- g$mtry == 11L & g$min_node_size == 10L & g$replace
  rownames(g) <- NULL
  g[, c("config_id", "mtry", "min_node_size", "replace", "sample_fraction", "simplicity_rank", "is_ranger_default")]
}

rf_theta <- function(config_id, grid, cfg) {
  r <- grid[grid$config_id == config_id, ]
  list(config_id = r$config_id, mtry = r$mtry, min_node_size = r$min_node_size,
       replace = r$replace, sample_fraction = r$sample_fraction, num_trees = cfg$num_trees)
}

# Load the data and group the six lags of each of the 22 variables.
rf_load <- function(h, cfg) {
  d <- load_model_data(h, cfg)
  base <- sub("_lag[0-5]$", "", colnames(d$X))
  d$groups <- split(seq_len(ncol(d$X)), factor(base, levels = unique(base)))
  ensure(length(d$groups) == 22L && all(lengths(d$groups) == 6L), "Expected 22 variables x 6 lags")
  d
}

rf_seed <- function(cfg, h, origin_m) as.integer(cfg$seed + 10L * origin_m + h)

# Single-threaded fit so that parallel and serial runs give identical forests.
rf_fit <- function(x, y, theta, seed) {
  ensure(!anyNA(y) && length(unique(y)) == 2L, "Training labels must contain both classes")
  ranger::ranger(x = x, y = factor(y, levels = c("0", "1")), probability = TRUE,
                 num.trees = theta$num_trees, mtry = theta$mtry, min.node.size = theta$min_node_size,
                 replace = theta$replace, sample.fraction = theta$sample_fraction, splitrule = "gini",
                 seed = seed, num.threads = 1L, verbose = FALSE)
}

# Probability of class 1 (stocks beat bills).
rf_prob <- function(fit, x) {
  pr <- stats::predict(fit, data = x, num.threads = 1L)$predictions
  ensure(identical(colnames(pr), c("0", "1")), "Unexpected class order")
  unname(pr[, "1"])
}

# Fit one origin. Optional donor rows give predictions with each variable group swapped
# for the donor's values, used only for permutation importance.
rf_task <- function(dd, task) {
  i <- task$i
  rows <- mature_rows(dd, dd$m[i])
  fit <- rf_fit(dd$X[rows, , drop = FALSE], dd$y[rows], task$theta, task$seed)
  perturbed <- NULL
  if (length(task$donors)) {
    G <- length(dd$groups); K <- length(task$donors)
    z <- dd$X[rep(i, G * K), , drop = FALSE]
    for (g in seq_len(G)) for (k in seq_len(K))
      z[(g - 1L) * K + k, dd$groups[[g]]] <- dd$X[task$donors[k], dd$groups[[g]]]
    perturbed <- matrix(rf_prob(fit, z), nrow = G, ncol = K, byrow = TRUE)
  }
  list(i = i, config_id = task$theta$config_id, p = rf_prob(fit, dd$X[i, , drop = FALSE]), n_train = length(rows),
       last_origin = max(dd$m[rows]), last_outcome = max(dd$end[rows]), seed = task$seed,
       train_prevalence = mean(dd$y[rows]), perturbed = perturbed)
}

rf_task_data <- function(d) list(X = d$X, y = d$y, m = d$m, end = d$end, groups = d$groups)

# Donor rows: seeded shuffles of the period's origins, shared by every candidate.
rf_donors <- function(rows, h, period, cfg) {
  set.seed(cfg$seed + 100L * h + if (period == "validation") 1L else 2L)
  t(vapply(seq_len(cfg$importance_reps), function(k) rows[sample.int(length(rows))], integer(length(rows))))
}

rf_make_tasks <- function(d, rows, h, theta, cfg, donors = NULL)
  lapply(seq_along(rows), function(j) list(i = rows[j], theta = theta, seed = rf_seed(cfg, h, d$m[rows[j]]),
    donors = if (is.null(donors)) integer() else donors[, j]))

# Worker processes load the same shared code and model functions.
rf_make_cluster <- function(workers) {
  if (workers <= 1L) return(NULL)
  cl <- parallel::makePSOCKcluster(workers)
  files <- c(file.path(UTILS_DIR, "utils.R"), file.path(SCRIPTS, "rf", "model.R"))
  parallel::clusterCall(cl, function(files) { for (f in files) source(f); NULL }, files)
  cl
}

# Run tasks serially (cl = NULL) or in parallel; each task has its own seed.
rf_run_tasks <- function(cl, d, tasks) {
  dd <- rf_task_data(d)
  if (is.null(cl)) return(lapply(tasks, function(task) rf_task(dd, task)))
  parallel::clusterCall(cl, function(v) { assign("RF_DATA", v, envir = globalenv()); NULL }, dd)
  run <- function(task) rf_task(RF_DATA, task)
  environment(run) <- globalenv()
  out <- parallel::parLapplyLB(cl, tasks, run)
  parallel::clusterCall(cl, function() { rm("RF_DATA", envir = globalenv()); NULL })
  out
}

# Collect one candidate's per-origin results into the standard forecast table.
rf_assemble <- function(d, rows, results, h, period, theta, cfg) {
  ensure(identical(vapply(results, `[[`, 0L, "i"), rows), "Results out of origin order")
  get <- function(name, type) vapply(results, `[[`, type, name)
  f <- forecast_table(d, rows, h, cfg$model, period, get("p", 0), get("n_train", 0L), get("last_origin", 0L),
                      get("last_outcome", 0L), get("train_prevalence", 0), theta$config_id)
  f$seed <- get("seed", 0L)
  perturbed <- if (!is.null(results[[1]]$perturbed)) simplify2array(lapply(results, `[[`, "perturbed"))
  list(forecasts = f, perturbed = perturbed)
}

# Out-of-sample grouped permutation importance (all six lags of a variable swapped together).
rf_importance <- function(f, perturbed, group_names) {
  ok <- f$classification_eligible & !is.na(f$target)
  base <- mean((f$probability[ok] - f$target[ok])^2)
  imp <- bind_rows(lapply(seq_along(group_names), function(g) {
    pk <- perturbed[g, , , drop = TRUE]
    if (is.null(dim(pk))) pk <- matrix(pk, nrow = 1L)
    brier_k <- apply(pk[, ok, drop = FALSE], 1L, function(p) mean((p - f$target[ok])^2))
    data.frame(base_variable = group_names[g], delta_Brier = mean(brier_k) - base,
               delta_Brier_sd_across_reps = stats::sd(brier_k),
               mean_abs_probability_change = mean(abs(sweep(pk, 2L, f$probability, "-"))),
               signal_flip_rate_at_tuned_cutoff = mean(sweep(pk >= f$cutoff[1], 2L, f$probability >= f$cutoff[1], "!=")))
  }))
  imp <- imp[order(-imp$delta_Brier), ]
  cbind(data.frame(period = f$period[1], h = f$h[1], config_id = f$config_id[1]), imp,
        rank_delta_Brier = seq_len(nrow(imp)), rank_reliance = rank(-imp$mean_abs_probability_change, ties.method = "min"))
}

# Partial dependence of the final test fit: all six lags of a variable held at each grid value.
rf_partial_dependence <- function(d, h, theta, cfg, f, variables) {
  i <- which(d$ym == cfg$test_end); rows <- mature_rows(d, d$m[i])
  fit <- rf_fit(d$X[rows, , drop = FALSE], d$y[rows], theta, rf_seed(cfg, h, d$m[i]))
  ensure(identical(rf_prob(fit, d$X[i, , drop = FALSE]), f$probability[f$period == "test" & f$origin == cfg$test_end]),
         "Final-origin refit does not reproduce its forecast")
  x <- d$X[rows, , drop = FALSE]
  bind_rows(lapply(variables, function(v) {
    qs <- seq(.05, .95, by = .05)
    values <- stats::quantile(x[, paste0(v, "_lag0")], qs, names = FALSE)
    keep <- !duplicated(values)
    pd <- function(value, cols) { z <- x; z[, cols] <- value; mean(rf_prob(fit, z)) }
    data.frame(h = h, fit_origin = cfg$test_end, base_variable = v, feature = paste0(v, "_lag0..lag5"),
      value = values[keep], quantile = qs[keep],
      partial_dependence = vapply(values[keep], pd, numeric(1), cols = d$groups[[v]]),
      partial_dependence_lag0_only = vapply(values[keep], pd, numeric(1), cols = paste0(v, "_lag0")))
  }))
}
