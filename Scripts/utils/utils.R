# Shared helpers for every model. Sourcing this file also loads the other files in utils/.

# Locate Scripts/utils from this file's path, then the project root and its package library.
UTILS_DIR <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  dirname(normalizePath(if (length(files)) files[[length(files)]] else "Scripts/utils/utils.R"))
})
SCRIPTS <- dirname(UTILS_DIR)
ROOT <- dirname(SCRIPTS)
.libPaths(unique(c(file.path(ROOT, ".R-library"), .libPaths())))

ensure <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
near <- function(a, b, tol = 1e-10)
  ensure(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)), "Numerical check failed")

# Install missing packages into the project library.
use_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    lib <- file.path(ROOT, ".R-library")
    dir.create(lib, showWarnings = FALSE)
    utils::install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
  }
  ensure(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)), "Package installation failed")
}

# Months as integers (year * 12 + month - 1) and back to "YYYY-MM".
month_number <- function(x) {
  x <- as.character(x)
  as.integer(substr(x, 1, 4)) * 12L + as.integer(substr(x, 6, 7)) - 1L
}
month_string <- function(x) sprintf("%04d-%02d", x %/% 12L, x %% 12L + 1L)

bind_rows <- function(x) {
  x <- do.call(rbind, x)
  if (is.null(x)) return(data.frame())
  rownames(x) <- NULL
  x
}
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")

results_dir <- function(name) {
  out <- file.path(ROOT, "Results", name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}
write_tables <- function(out, tables)
  for (name in names(tables)) write_csv(tables[[name]], file.path(out, paste0(name, ".csv")))

# Load the model-ready data for horizon h and check dates, targets and split boundaries.
# features = NULL keeps all 132 predictors; otherwise keep only the named columns.
load_model_data <- function(h, cfg, features = NULL) {
  d <- readRDS(file.path(ROOT, "Data", "processed", paste0("model_ready_h", h, ".rds")))
  d$X <- as.matrix(d$X)
  if (is.null(features)) {
    ensure(ncol(d$X) == 132L && identical(colnames(d$X), d$feature_names) &&
             all(grepl("_lag[0-5]$", colnames(d$X))), "Expected the 132 lagged predictors")
  } else {
    ensure(all(features %in% colnames(d$X)), "Missing predictors")
    d$X <- d$X[, features, drop = FALSE]
  }
  d$m <- month_number(d$ym)
  d$end <- month_number(d$metadata$outcome_end)
  known <- !is.na(d$y)
  ensure(all(diff(d$m) == 1L), "Incomplete model calendar")
  ensure(all(is.finite(d$X)), "Missing or nonfinite predictor")
  ensure(identical(d$ym, as.character(d$metadata$ym)), "Metadata misaligned")
  ensure(all(d$y[known] %in% 0:1) && all(d$end[known] == d$m[known] + h), "Invalid target or horizon")
  ensure(all(d$y[known] == (d$metadata$equity_forward[known] > d$metadata$bill_accrual_proxy[known])),
         "Target is not 1{equity beats bills}")
  ensure(all(d$end[d$rows$train] <= month_number(cfg$train_end)) &&
           all(d$end[d$rows$validation] <= month_number(cfg$validation_end)), "Split boundary leak")
  d
}

# Rows whose outcome is known at the close of the origin month.
mature_rows <- function(d, origin) which(!is.na(d$y) & !is.na(d$end) & d$end <= origin & d$m < origin)

period_rows <- function(d, cfg, period)
  which(d$m >= month_number(cfg[[paste0(period, "_start")]]) & d$m <= month_number(cfg[[paste0(period, "_end")]]))

# Standard forecast table: one row per origin plus the training window behind each fit.
forecast_table <- function(d, rows, h, model, period, probability, n_train, last_origin,
                           last_outcome, train_prevalence, config_id) {
  ensure(all(last_outcome <= d$m[rows]) && all(last_origin <= d$m[rows] - h), "Training outcome leaks past origin")
  ensure(all(is.finite(probability) & probability >= 0 & probability <= 1), "Invalid probability")
  data.frame(h = h, model = model, period = period, origin = d$ym[rows], probability = probability,
    train_prevalence = train_prevalence, outcome_end = as.character(d$metadata$outcome_end[rows]),
    target = d$y[rows], classification_eligible = d$split[rows] == period,
    fit_information_cutoff = d$ym[rows], last_train_origin = month_string(last_origin),
    last_train_outcome = month_string(last_outcome), n_train = n_train, config_id = config_id,
    stringsAsFactors = FALSE)
}

# Area under the ROC curve; ties count one half.
auc <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Run a named list of check functions; any failure stops the run.
run_checks <- function(checks) {
  for (name in names(checks))
    tryCatch(checks[[name]](), error = function(e) stop("Check failed: ", name, " (", conditionMessage(e), ")", call. = FALSE))
  message(length(checks), " checks passed.")
  names(checks)
}

# Record settings, package versions, input file hashes and checks for one run.
write_run_info <- function(out, cfg, checks, packages = character(), inputs = model_inputs(cfg)) {
  settings <- vapply(cfg, function(x) paste(format(x), collapse = ", "), "")
  versions <- vapply(packages, function(p) as.character(utils::packageVersion(p)), "")
  writeLines(c(paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), R.version.string,
    if (length(packages)) paste("Package", packages, versions),
    "", "Settings:", paste0("  ", names(settings), ": ", settings),
    "", "Input files (md5):", paste0("  ", substring(inputs, nchar(ROOT) + 2L), "  ", unname(tools::md5sum(inputs))),
    "", paste0("Checks passed (", length(checks), "):"), paste0("  ", checks)),
    file.path(out, "run_info.txt"))
}

for (file in c("config.R", "portfolio.R", "evaluate.R", "tests.R")) source(file.path(UTILS_DIR, file))
