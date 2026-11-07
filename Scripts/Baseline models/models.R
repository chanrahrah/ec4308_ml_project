# No raw predictor reconstruction: consume only the existing validated RDS inputs.
load_model_data <- function(root, h, cfg) {
  d <- readRDS(file.path(root, "Data", "processed", paste0("model_ready_h", h, ".rds")))
  ensure(all(cfg$features %in% colnames(d$X)), "Required baseline predictors are absent")
  d$X <- d$X[, cfg$features, drop = FALSE]
  d$m <- month_number(d$ym)
  d$end <- month_number(d$metadata$outcome_end)
  ensure(all(diff(d$m) == 1L) && !anyDuplicated(d$m), "Incomplete model calendar")
  ensure(all(is.finite(d$X)), "Missing/nonfinite baseline predictor")
  known <- !is.na(d$y)
  ensure(all(d$y[known] %in% 0:1), "Invalid binary target")
  ensure(all(d$end[known] == d$m[known] + h), "Target horizon mismatch")
  ensure(identical(d$ym, as.character(d$metadata$ym)), "Metadata alignment mismatch")
  ensure(all(d$end[d$rows$train] <= month_number(cfg$train_end)), "Training boundary leak")
  ensure(all(d$end[d$rows$validation] <= month_number(cfg$validation_end)), "Validation boundary leak")
  d
}

# At the close of t, only outcomes completed by t can enter estimation.
# Previously purged observations are restored only once their outcomes mature.
mature_rows <- function(d, origin) which(!is.na(d$y) & !is.na(d$end) & d$end <= origin & d$m < origin)

fit_binary <- function(d, rows, link) {
  ensure(length(rows) > ncol(d$X) + 2L && length(unique(d$y[rows])) == 2L,
         "Insufficient observations/classes")
  x <- d$X[rows, , drop = FALSE]
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  ensure(all(is.finite(scale) & scale > 1e-12), "Degenerate training predictor")
  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  design <- cbind(`(Intercept)` = 1, z)
  warnings <- character()
  fit <- withCallingHandlers(
    stats::glm.fit(design, d$y[rows], family = stats::binomial(link),
                   control = stats::glm.control(epsilon = 1e-10, maxit = 100L)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  ensure(fit$converged && fit$rank == ncol(design) && all(is.finite(fit$coefficients)) && !fit$boundary,
         paste("Failed/rank-deficient/boundary", link, "fit"))
  ensure(!length(warnings), paste("GLM warning:", paste(warnings, collapse = "; ")))
  list(link = link, coefficients = fit$coefficients, center = center, scale = scale,
       rows = rows, n = length(rows), last_origin = max(d$m[rows]),
       last_outcome = max(d$end[rows]), iterations = fit$iter, deviance = fit$deviance)
}

predict_binary <- function(fit, x) {
  ensure(identical(colnames(x), names(fit$center)), "Predictor order mismatch")
  z <- sweep(sweep(x, 2L, fit$center, "-"), 2L, fit$scale, "/")
  eta <- as.vector(cbind(1, z) %*% fit$coefficients)
  p <- if (fit$link == "probit") stats::pnorm(eta) else stats::plogis(eta)
  ensure(all(is.finite(p) & p >= 0 & p <= 1), "Invalid forecast probability")
  p
}

forecast_frame <- function(d, rows, p, h, link, split) data.frame(
  h = h, model = link, period = split, origin = d$ym[rows], probability = p,
  outcome_end = as.character(d$metadata$outcome_end[rows]), target = d$y[rows],
  classification_eligible = d$split[rows] == split, stringsAsFactors = FALSE)

fit_records <- function(fit, h, link, origin, period) list(
  audit = data.frame(h = h, model = link, period = period, forecast_origin = origin,
    n_train = fit$n, last_train_origin = month_string(fit$last_origin),
    last_train_outcome = month_string(fit$last_outcome), iterations = fit$iterations,
    converged = TRUE),
  coefficients = data.frame(h = h, model = link, period = period, forecast_origin = origin,
    term = names(fit$coefficients), coefficient_standardized = unname(fit$coefficients)),
  scalers = data.frame(h = h, model = link, period = period, forecast_origin = origin,
    feature = names(fit$center), center = unname(fit$center), scale = unname(fit$scale)))

expanding_forecasts <- function(d, h, link, cfg, period = "test") {
  ensure(period %in% c("validation", "test"), "Unknown forecast period")
  start <- cfg[[paste0(period, "_start")]]; end <- cfg[[paste0(period, "_end")]]
  rows <- which(d$m >= month_number(start) & d$m <= month_number(end))
  ensure(length(rows) > 0, "No forecast origins")
  p <- numeric(length(rows)); records <- vector("list", length(rows))
  for (j in seq_along(rows)) {
    i <- rows[j]
    train <- mature_rows(d, d$m[i])
    fit <- fit_binary(d, train, link)
    ensure(fit$last_outcome <= d$m[i] && fit$last_origin <= d$m[i] - h,
           "Expanding-window outcome leakage")
    p[j] <- predict_binary(fit, d$X[i, , drop = FALSE])
    records[[j]] <- fit_records(fit, h, link, d$ym[i], period)
  }
  forecasts <- forecast_frame(d, rows, p, h, link, period)
  audit <- do.call(rbind, lapply(records, `[[`, "audit"))
  forecasts$fit_information_cutoff <- forecasts$origin
  forecasts$last_train_origin <- audit$last_train_origin
  forecasts$last_train_outcome <- audit$last_train_outcome
  forecasts$n_train <- audit$n_train
  forecasts$config_id <- cfg$version
  list(forecasts = forecasts, records = records, final_fit = fit)
}
