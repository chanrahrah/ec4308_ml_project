# Ensemble checks: the vote, merging member forecasts and the GR solver.
ensemble_checks <- function() {
  # Toy member forecasts: h = 1, 12 consecutive origins.
  toy <- function(model, p, cutoff) {
    o <- month_number("2005-01") + seq_along(p) - 1L
    data.frame(model = model, h = 1L, origin = month_string(o), probability = p, cutoff = cutoff,
      signal = as.integer(p >= cutoff), cutoff_050 = .5, signal_050 = as.integer(p >= .5),
      target = rep(c(1L, 0L), length.out = length(p)), outcome_end = month_string(o + 1L), classification_eligible = TRUE)
  }
  set.seed(4308)
  pa <- runif(12); pb <- runif(12); pc <- runif(12)
  mem <- list(a = toy("a", pa, .55), b = toy("b", pb, .45), c = toy("c", pc, .6))
  fails <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")
  set.seed(11)
  n <- 400
  Q <- cbind(probit = runif(n, .2, .9), enet = runif(n, .2, .9), rf = runif(n, .2, .9))
  yq <- rbinom(n, 1, .7 * Q[, "probit"] + .3 * Q[, "enet"])

  list(
    "Majority vote: stocks only with a strict majority" = function() {
      ensure(identical(majority_vote(cbind(c(1, 1, 0, 0), c(1, 0, 1, 0), c(0, 0, 1, 0))), c(1L, 0L, 1L, 0L)) &&
               identical(majority_vote(cbind(c(1, 1, 0), c(0, 1, 0))), c(0L, 1L, 0L)) &&
               fails(majority_vote(cbind(c(1, NA), c(0, 1)))), "Vote rule failed")
    },
    "Members vote with their own cutoffs, merged by key not row order" = function() {
      v <- build_vote(mem, "validation")
      ensure(identical(v$vote, as.integer((pa >= .55) + (pb >= .45) + (pc >= .6) >= 2)) &&
               identical(v$vote_050, as.integer((pa >= .5) + (pb >= .5) + (pc >= .5) >= 2)), "Wrong vote")
      shuffled <- mem; shuffled$b <- shuffled$b[c(7:12, 1:6), ]; shuffled$c <- shuffled$c[12:1, ]
      ensure(identical(build_vote(shuffled, "validation")$vote, v$vote), "Vote depends on row order")
    },
    "Inconsistent member forecasts are refused" = function() {
      bad <- function(m, f) { x <- mem; x[[m]] <- f(x[[m]]); fails(build_vote(x, "validation")) }
      ensure(bad("b", function(x) x[-5, ]) && bad("b", function(x) rbind(x, x[3, ])) &&
               bad("c", function(x) { x$target[2] <- 1L - x$target[2]; x }) &&
               bad("c", function(x) { x$classification_eligible[4] <- FALSE; x }) &&
               bad("a", function(x) { x$signal[1] <- 1L - x$signal[1]; x }) &&
               bad("a", function(x) { x$cutoff[3] <- .7; x$signal <- as.integer(x$probability >= x$cutoff); x }) &&
               bad("a", function(x) { x$probability[1] <- 1.2; x }), "Bad member forecasts were accepted")
    },
    "Vote share >= (floor(M/2)+1)/M reproduces the vote" = function() {
      for (M in 1:6) {
        S <- as.matrix(expand.grid(rep(list(0:1), M)))
        ensure(identical(as.integer(rowMeans(S) >= (floor(M / 2) + 1) / M), majority_vote(S)), "Vote share mismatch")
      }
    },
    "GR matches a brute-force search over the weight simplex" = function() {
      grid <- expand.grid(a = seq(0, 1, by = .005), b = seq(0, 1, by = .005)); grid <- grid[grid$a + grid$b <= 1 + 1e-12, ]
      W <- cbind(grid$a, grid$b, 1 - grid$a - grid$b)
      best <- apply(W, 1, function(w) sum((yq - mean(yq - Q %*% w) - Q %*% w)^2))
      g <- gr_weights(Q, yq)
      ensure(g$sse <= min(best) + 1e-9 && max(abs(g$weights - W[which.min(best), ])) < .01, "GR is not optimal")
    },
    "GR weights are nonnegative, sum to one and recover the truth" = function() {
      m <- 20000; Z <- cbind(probit = runif(m, .1, .9), enet = runif(m, .1, .9), rf = runif(m, .1, .9))
      g <- gr_weights(Z, rbinom(m, 1, .7 * Z[, "probit"] + .3 * Z[, "enet"]))
      ensure(all(g$weights >= 0) && abs(sum(g$weights) - 1) < 1e-12 && abs(g$weights[["probit"]] - .7) < .03 &&
               abs(g$weights[["enet"]] - .3) < .03 && g$weights[["rf"]] < .03, "GR weights wrong")
      Y <- 1.5 * Z[, "probit"] - .5 * Z[, "enet"]
      near(unname(gr_weights(Z[, c("probit", "enet")], Y)$weights), c(1, 0))
    },
    "GR handles constant and duplicate members; failures fall back to equal weights" = function() {
      g <- gr_weights(cbind(Q[, c("probit", "rf")], enet = rep(.8, n)), yq)
      D <- cbind(a = Q[, 1], b = Q[, 1], c = Q[, 3])
      ensure(abs(sum(g$weights) - 1) < 1e-12 && all(g$weights >= 0) &&
               isTRUE(all.equal(gr_weights(D, yq)$sse, gr_weights(D[, c("a", "c")], yq)$sse)), "Degenerate members mishandled")
      f <- gr_weights_safe(Q, yq[-1])
      ensure(grepl("^fallback_equal_weights", f$status) && f$intercept == 0, "Fallback failed")
    },
    "GR predictions are clipped and matched to members by name" = function() {
      fit <- list(weights = c(probit = .5, enet = .5, rf = 0), intercept = .6)
      g <- gr_predict(Q[, c("rf", "enet", "probit")], fit)
      raw <- .6 + .5 * Q[, "probit"] + .5 * Q[, "enet"]
      near(g$raw, raw)
      ensure(all(g$probability <= 1) && identical(g$clipped, raw > 1) && fails(gr_predict(Q[, 1:2], fit)), "Clipping failed")
    },
    "Yearly GR only uses outcomes known by the preceding December" = function() {
      set.seed(7); n <- 120; o <- 2005L * 12L + 0:(n - 1L); oe <- o + 3L
      Z <- cbind(a = runif(n, .2, .9), b = runif(n, .2, .9)); yz <- rbinom(n, 1, Z[, "a"])
      g <- gr_yearly(Z, yz, o, oe, rep(TRUE, n), 49:n)
      yy <- yz; late <- oe > 2011L * 12L - 1L; yy[late] <- 1L - yy[late]
      k <- (o[49:n] %/% 12L) <= 2011L
      near(g$raw[k], gr_yearly(Z, yy, o, oe, rep(TRUE, n), 49:n)$raw[k])
      # January 2014, h = 6, history from 2005-01: last usable origin 2013-06, 102 rows.
      g6 <- gr_yearly(Z, yz, o, o + 6L, rep(TRUE, n), 109:120)
      ensure(g6$weights$n_rows[1] == 102L && g6$weights$last_outcome[1] == "2013-12", "Wrong yearly information set")
    })
}
