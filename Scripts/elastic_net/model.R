# Logistic elastic net on all 132 predictors with a fixed (alpha, lambda) at every monthly refit.
# Penalty: lambda * (alpha * L2^2 + (1 - alpha) * L1), so alpha = 0 is lasso and alpha = 1 is ridge.
elastic_config <- function() modifyList(base_config(), list(
  model = "elastic_net", alphas = c(0, 0.25, 0.5, 0.75, 1),
  lambdas = c(0.0001, 0.0005, 0.002, 0.01, 0.05, 0.1, 0.2, 0.3)))

# Standardise on the fitting rows and fit glmnet at the requested penalty.
fit_elastic <- function(d, rows, alpha, lambda) {
  ensure(length(rows) > 40L && length(unique(d$y[rows])) == 2L, "Too few fitting rows or one class")
  x <- d$X[rows, , drop = FALSE]; y <- as.numeric(d$y[rows]); n <- length(rows)
  center <- colMeans(x); scales <- apply(x, 2, stats::sd)
  zero_scale <- scales < 1e-12
  scales[zero_scale] <- 1
  z <- sweep(sweep(x, 2, center, "-"), 2, scales, "/")
  if (any(zero_scale)) z[, zero_scale] <- 0
  # Convert to glmnet's (alpha, lambda) convention.
  alpha_glmnet <- (1 - alpha) / (1 + alpha)
  lambda_glmnet <- lambda * (1 + alpha)
  # Fit along a short path ending at lambda for numerical stability; only lambda is used.
  max_lambda <- max(abs(as.numeric(crossprod(z, y - mean(y))))) / n
  high <- max(lambda_glmnet, 1.1 * max_lambda / max(alpha_glmnet, 0.05))
  path <- exp(seq(log(high), log(lambda_glmnet), length.out = 30L))
  if (max(abs(path - lambda_glmnet)) < 1e-12) path <- lambda_glmnet
  fit <- glmnet::glmnet(x = z, y = y, family = "binomial", alpha = alpha_glmnet, lambda = path,
                        standardize = FALSE, intercept = TRUE,
                        control = list(thresh = 1e-8, maxit = 100000L))
  ensure(all(fit$jerr == 0), "glmnet did not converge")
  coefficients <- as.matrix(stats::coef(fit, s = lambda_glmnet, exact = FALSE))[, 1]
  list(object = fit, lambda_glmnet = lambda_glmnet, center = center, scale = scales,
       nonzero = sum(abs(coefficients[-1]) > 1e-8))
}

predict_elastic <- function(fit, x) {
  z <- sweep(sweep(x, 2, fit$center, "-"), 2, fit$scale, "/")
  as.numeric(stats::predict(fit$object, newx = z, s = fit$lambda_glmnet, type = "response"))
}

# Refit every month of the period with alpha and lambda held fixed.
recursive_elastic <- function(d, h, cfg, period, alpha, lambda) {
  rows <- period_rows(d, cfg, period)
  train <- lapply(rows, function(i) mature_rows(d, d$m[i]))
  fits <- lapply(train, fit_elastic, d = d, alpha = alpha, lambda = lambda)
  p <- vapply(seq_along(rows), function(j) predict_elastic(fits[[j]], d$X[rows[j], , drop = FALSE]), 0)
  f <- forecast_table(d, rows, h, cfg$model, period, p, lengths(train),
    vapply(train, function(r) max(d$m[r]), 0L), vapply(train, function(r) max(d$end[r]), 0L),
    vapply(train, function(r) mean(d$y[r]), 0), sprintf("alpha%.2f_lambda%g", alpha, lambda))
  list(forecasts = f, nonzero = vapply(fits, `[[`, 0L, "nonzero"))
}

# Score every (alpha, lambda) on validation; Brier ties go to fewer active terms,
# then larger lambda, then smaller alpha.
tune_elastic <- function(d, h, cfg) {
  grid <- expand.grid(alpha = cfg$alphas, lambda = cfg$lambdas)
  candidates <- lapply(seq_len(nrow(grid)), function(k) {
    message(sprintf("h=%d candidate %d/%d", h, k, nrow(grid)))
    recursive_elastic(d, h, cfg, "validation", grid$alpha[k], grid$lambda[k])
  })
  scores <- cbind(data.frame(h = h, config_id = vapply(candidates, function(o) o$forecasts$config_id[1], "")), grid,
                  mean_nonzero = vapply(candidates, function(o) mean(o$nonzero), 0),
                  bind_rows(lapply(candidates, function(o) validation_scores(o$forecasts))))
  best <- best_candidate(scores$validation_Brier, order(scores$mean_nonzero, -scores$lambda, scores$alpha))
  scores$selected <- seq_len(nrow(scores)) == best
  list(scores = scores, forecasts = candidates[[best]]$forecasts, alpha = grid$alpha[best], lambda = grid$lambda[best])
}
