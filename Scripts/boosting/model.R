# Gradient boosting (xgboost) on all 132 predictors.
boosting_config <- function() modifyList(base_config(), list(
  model = "boosting", depths = c(1L, 2L, 3L), etas = c(0.02, 0.05, 0.10), nrounds = c(50L, 100L, 200L, 400L),
  subsample = 0.7, colsample = 0.5, min_child_weight = 5, l2_penalty = 1))

boosting_grid <- function(cfg) {
  g <- expand.grid(max_depth = cfg$depths, eta = cfg$etas, nrounds = cfg$nrounds)
  g$config_id <- sprintf("d%d_eta%.2f_n%d", g$max_depth, g$eta, g$nrounds)
  g
}

# One single-threaded fit; returns the forecast after each requested number of rounds.
fit_predict <- function(d, rows, i, depth, eta, seed, rounds, cfg) {
  ensure(length(rows) > 50, "Too few fitting rows")
  train <- xgboost::xgb.DMatrix(d$X[rows, , drop = FALSE], label = d$y[rows], missing = NA)
  params <- list(objective = "binary:logistic", max_depth = depth, eta = eta,
                 subsample = cfg$subsample, colsample_bytree = cfg$colsample,
                 min_child_weight = cfg$min_child_weight, lambda = cfg$l2_penalty, nthread = 1L, seed = seed)
  fit <- xgboost::xgb.train(params, train, nrounds = max(rounds), verbose = 0)
  x <- d$X[i, , drop = FALSE]
  vapply(rounds, function(n) as.numeric(stats::predict(fit, x, iterationrange = c(1L, as.integer(n)))), numeric(1))
}

# Monthly expanding forecasts for every configuration in grid. Configurations that share
# a depth and learning rate come from one fit read at different numbers of rounds.
boosting_forecasts <- function(d, h, period, grid, cfg) {
  rows <- period_rows(d, cfg, period)
  fits <- unique(grid[, c("max_depth", "eta")])
  P <- matrix(NA_real_, length(rows), nrow(grid))
  train <- lapply(rows, function(i) mature_rows(d, d$m[i]))
  for (j in seq_along(rows)) {
    for (k in seq_len(nrow(fits))) {
      cols <- which(grid$max_depth == fits$max_depth[k] & grid$eta == fits$eta[k])
      P[j, cols] <- fit_predict(d, train[[j]], rows[j], fits$max_depth[k], fits$eta[k],
                                cfg$seed + h * 100003L + d$m[rows[j]], grid$nrounds[cols], cfg)
    }
    if (j %% 24 == 0) message(sprintf("h=%d %s: %d/%d origins", h, period, j, length(rows)))
  }
  lapply(seq_len(nrow(grid)), function(k) forecast_table(d, rows, h, cfg$model, period, P[, k], lengths(train),
    vapply(train, function(r) max(d$m[r]), 0L), vapply(train, function(r) max(d$end[r]), 0L),
    vapply(train, function(r) mean(d$y[r]), 0), grid$config_id[k]))
}
