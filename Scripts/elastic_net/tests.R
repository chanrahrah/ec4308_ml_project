# Elastic net checks: penalty conversion, the grid and no look-ahead.
elastic_checks <- function(cfg) {
  checks <- list(
    "Penalty converts exactly to glmnet's convention" = function() {
      for (a in cfg$alphas) for (lam in cfg$lambdas) {
        ag <- (1 - a) / (1 + a); lg <- lam * (1 + a)
        near(lg * ag, lam * (1 - a)); near(lg * (1 - ag) / 2, lam * a)
      }
    },
    "Grid has 5 x 8 = 40 distinct candidates" = function()
      ensure(length(unique(cfg$alphas)) == 5L && length(unique(cfg$lambdas)) == 8L && all(cfg$lambdas > 0), "Unexpected grid"))
  for (h in cfg$horizons) checks[[paste("Future data cannot change a forecast h", h)]] <- local({ h <- h; function() {
    d <- load_model_data(h, cfg)
    perturbation_check(d, h, function(d, i) predict_elastic(fit_elastic(d, mature_rows(d, d$m[i]), .5, .01), d$X[i, , drop = FALSE]))
  }})
  checks[["Validation refits expand by one month"]] <- function() {
    short <- cfg; short$validation_end <- "2005-03"
    f <- recursive_elastic(load_model_data(3L, cfg), 3L, short, "validation", .5, .01)$forecasts
    ensure(nrow(f) == 3L && all(diff(f$n_train) == 1L), "Recursive window incorrect")
  }
  checks
}
