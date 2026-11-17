# Checks of the Sharpe-difference test on simulated data with known answers.
sharpe_test_checks <- function(cfg) {
  set.seed(cfg$seed)
  list(
    "Gradient matches numerical derivatives" = function() {
      y <- matrix(rnorm(200, .01, .04), 100, 2)
      f <- function(m) { v <- m[3:4] - m[1:2]^2; m[1] / sqrt(v[1]) - m[2] / sqrt(v[2]) }
      m <- c(colMeans(y), colMeans(y^2))
      num <- vapply(1:4, function(k) { e <- replace(numeric(4), k, 1e-7); (f(m + e) - f(m - e)) / 2e-7 }, 0)
      near(sharpe_diff(y)$grad, num, 1e-6); near(sharpe_diff(y)$diff, f(m))
    },
    "Newey-West with lag 0 is the sample covariance" = function() {
      v <- scale(matrix(rnorm(400), 100, 4), scale = FALSE)
      near(newey_west(v, 0), crossprod(v) / 100)
    },
    "Standard error matches the closed form for independent normal returns" = function() {
      n <- 20000; rho <- .7
      e <- matrix(rnorm(2 * n), n, 2); y <- cbind(.1 + e[, 1], .05 + rho * e[, 1] + sqrt(1 - rho^2) * e[, 2])
      s <- sharpe_diff(y); se <- sqrt(drop(s$grad %*% newey_west(moments(y), 0) %*% s$grad) / n)
      # Memmel (2003): Var = (2 - 2 rho + (SR1^2 + SR2^2 - 2 SR1 SR2 rho^2) / 2) / n
      near(se, sqrt((2 - 2 * rho + (.1^2 + .05^2 - 2 * .1 * .05 * rho^2) / 2) / n), .02)
    },
    "A large true difference is detected; swapping the series flips the sign only" = function() {
      e <- matrix(rnorm(300), 150, 2); y <- cbind(.6 + e[, 1], e[, 2])
      a <- lw_test(y, 6L, 499L, 1L); b <- lw_test(y[, 2:1], 6L, 499L, 1L)
      ensure(a$p_value < .01, "Large difference not detected")
      near(a$diff, -b$diff); near(a$se, b$se); near(a$p_value, b$p_value)
    },
    "Series that differ in only a few months still give a valid p-value" = function() {
      x <- rnorm(150, .01, .04); y <- x; y[c(20, 90)] <- 0.001
      p <- lw_test(cbind(x, y), 6L, 499L, 1L)$p_value
      ensure(is.finite(p) && p > 0 && p <= 1, "Invalid p-value for nearly identical series")
    },
    "Under equal Sharpe ratios the test rejects at about the nominal rate" = function() {
      p <- vapply(1:200, function(k) lw_test(matrix(rnorm(300, .1), 150, 2), 6L, 199L, k)$p_value, 0)
      ensure(mean(p < .05) < .10, "Test rejects too often when Sharpe ratios are equal")
    })
}
