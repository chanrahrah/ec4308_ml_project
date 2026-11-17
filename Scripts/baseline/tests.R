# Baseline checks: no look-ahead, recursive refits and the bootstrap helper.
baseline_checks <- function(cfg) {
  checks <- list()
  for (h in cfg$horizons) checks[[paste("Future data cannot change a forecast h", h)]] <- local({ h <- h; function() {
    d <- load_model_data(h, cfg, cfg$features)
    for (link in cfg$links)
      perturbation_check(d, h, function(d, i) predict_binary(fit_binary(d, mature_rows(d, d$m[i]), link), d$X[i, , drop = FALSE]))
  }})
  checks[["Validation refits expand by one month and match an independent fit"]] <- function() {
    d <- load_model_data(3L, cfg, cfg$features)
    short <- cfg; short$validation_end <- "2005-03"
    f <- expanding_forecasts(d, 3L, "probit", short, "validation")$forecasts
    ensure(nrow(f) == 3L && all(diff(f$n_train) == 1L), "Validation history does not expand")
    ensure(all(month_number(f$last_train_origin) == month_number(f$origin) - 3L), "Immature targets in fit")
    feb <- which(d$ym == "2005-02"); rows <- mature_rows(d, d$m[feb])
    near(f$probability[2], predict_binary(fit_binary(d, rows, "probit"), d$X[feb, , drop = FALSE]))
    near(f$train_prevalence[2], mean(d$y[rows]))
  }
  checks[["Bootstrap indices stay in range"]] <- function() {
    set.seed(cfg$seed); ids <- stationary_indices(48, 12)
    ensure(length(ids) == 48 && all(ids %in% 1:48), "Invalid bootstrap row indices")
  }
  checks
}
