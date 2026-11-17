# Ledoit and Wolf (2008) test of equal Sharpe ratios: studentised circular block bootstrap.
sharpe_test_config <- function() modifyList(base_config(), list(
  model = "sharpe_tests", block_length = 6L, replicates = 4999L,
  members = c("baseline", "rf", "elastic_net", "boosting", "ensemble")))

# Sharpe difference of two excess-return series and its gradient with respect to
# (mean_1, mean_2, mean square_1, mean square_2).
sharpe_diff <- function(y) {
  mu <- colMeans(y); g <- colMeans(y^2); v <- g - mu^2
  list(diff = mu[[1]] / sqrt(v[1]) - mu[[2]] / sqrt(v[2]),
       grad = c(g[1] / v[1]^1.5, -g[2] / v[2]^1.5, -mu[1] / (2 * v[1]^1.5), mu[2] / (2 * v[2]^1.5)))
}

# Centred moment series behind the two means and mean squares.
moments <- function(y) {
  mu <- colMeans(y); g <- colMeans(y^2)
  cbind(y[, 1] - mu[1], y[, 2] - mu[2], y[, 1]^2 - g[1], y[, 2]^2 - g[2])
}

# Newey-West long-run covariance with a Bartlett kernel.
newey_west <- function(v, lag) {
  n <- nrow(v)
  S <- crossprod(v) / n
  for (k in seq_len(lag)) {
    G <- crossprod(v[-(1:k), , drop = FALSE], v[1:(n - k), , drop = FALSE]) / n
    S <- S + (1 - k / (lag + 1)) * (G + t(G))
  }
  S
}

# Test H0: equal Sharpe ratios. The original standard error uses Newey-West; each bootstrap
# sample joins floor(n / b) random blocks of b months (wrapping around) for both series together,
# and its standard error comes from the block means.
lw_test <- function(y, b, replicates, seed) {
  n <- nrow(y); l <- n %/% b
  est <- sharpe_diff(y)
  se <- sqrt(drop(est$grad %*% newey_west(moments(y), b) %*% est$grad) / n)
  set.seed(seed)
  stat <- replicate(replicates, {
    ids <- (as.vector(outer(0:(b - 1L), sample.int(n, l, replace = TRUE), "+")) - 1L) %% n + 1L
    ys <- y[ids, , drop = FALSE]
    bs <- sharpe_diff(ys)
    z <- rowsum(moments(ys), rep(seq_len(l), each = b)) / sqrt(b)
    se_b <- sqrt(max(drop(bs$grad %*% (crossprod(z) / l) %*% bs$grad), 0) / (l * b))
    # A resample that misses every month where the strategies differ has zero spread; it counts as extreme.
    if (se_b > 0) abs(bs$diff - est$diff) / se_b else Inf
  })
  list(diff = est$diff, se = se, p_value = (1 + sum(stat >= abs(est$diff) / se)) / (replicates + 1))
}
