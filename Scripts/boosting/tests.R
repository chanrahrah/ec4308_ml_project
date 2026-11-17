# Boosting checks: no look-ahead and seeded determinism.
boosting_checks <- function(cfg) {
  checks <- list()
  for (h in cfg$horizons) checks[[paste("Future data cannot change a forecast; seeds are reproducible h", h)]] <- local({ h <- h; function() {
    d <- load_model_data(h, cfg)
    predict_at <- function(d, i) fit_predict(d, mature_rows(d, d$m[i]), i, 2L, 0.1, 99L, c(50L, 100L), cfg)
    perturbation_check(d, h, predict_at, tol = 0)
    i <- which(d$ym == "2008-06"); rows <- mature_rows(d, d$m[i])
    a <- fit_predict(d, rows, i, 2L, 0.1, 99L, 50L, cfg)
    ensure(identical(a, fit_predict(d, rows, i, 2L, 0.1, 99L, 50L, cfg)), "Same seed gave a different fit")
    ensure(!identical(a, fit_predict(d, rows, i, 2L, 0.1, 100L, 50L, cfg)), "Seed is ignored")
  }})
  checks
}
