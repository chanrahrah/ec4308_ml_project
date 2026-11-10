# tests_ensemble.R: tests for the three ensembles (group spec rules 51-57 and 69-70).
# Run by 01_ensemble.R before any member export is read; a failure stops the script.
run_ensemble_tests <- function(json_path=NULL) {
  passed <- character()
  check <- function(name, expr) { force(expr); passed <<- c(passed, name) }
  fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error=function(e) e), "error")
  # Toy member export: h = 1, 12 consecutive origins.
  toy <- function(model, p, cutoff=.5, origin0="2005-01") {
    n <- length(p); o <- to_month_id(origin0) + 0:(n - 1)
    data.frame(model=model, h=1L, origin=from_month_id(o), probability=p, cutoff=cutoff,
      signal=as.integer(p >= cutoff), cutoff_050=.5, signal_050=as.integer(p >= .5),
      target=rep(c(1L, 0L), length.out=n), outcome_end=from_month_id(o + 1L), classification_eligible=TRUE)
  }
  set.seed(4308)
  pa <- runif(12); pb <- runif(12); pc <- runif(12)
  mem <- list(a=toy("a", pa, .55), b=toy("b", pb, .45), c=toy("c", pc, .6))

  check("strict majority: stocks only when more than half say stocks; ties go to bills", {
    stopifnot(identical(majority_vote(cbind(c(1, 1, 0, 0), c(1, 0, 1, 0), c(0, 0, 1, 0))), c(1L, 0L, 1L, 0L)))
    stopifnot(identical(majority_vote(cbind(c(1, 1, 0), c(0, 1, 0))), c(0L, 1L, 0L)))            # M = 2: 1-1 tie -> bills
    stopifnot(identical(majority_vote(cbind(1, 1, 0, 0)), 0L), identical(majority_vote(cbind(1, 1, 1, 0)), 1L))
    stopifnot(fails(majority_vote(cbind(c(1, NA), c(0, 1)))), fails(majority_vote(cbind(.5, 1))))
  })

  check("members vote with their own frozen cutoffs, and the 0.50 vote is built separately", {
    v <- build_vote(mem, "validation")
    stopifnot(identical(v$vote, as.integer((pa >= .55) + (pb >= .45) + (pc >= .6) >= 2)))
    stopifnot(identical(v$vote_050, as.integer((pa >= .5) + (pb >= .5) + (pc >= .5) >= 2)))
    stopifnot(identical(sort(unique(v$cutoffs$cutoff)), c(.45, .55, .6)))
  })

  check("members are merged by (h, origin) keys, not row order", {
    shuffled <- mem; shuffled$b <- shuffled$b[c(7:12, 1:6), ]; shuffled$c <- shuffled$c[12:1, ]
    stopifnot(identical(build_vote(shuffled, "validation")$vote, build_vote(mem, "validation")$vote))
  })

  check("future data cannot change an earlier vote", {
    later <- mem
    for (m in names(later)) {                      # rewrite every origin after the 6th: probabilities and targets
      x <- later[[m]]; k <- 7:12; x$probability[k] <- 1 - x$probability[k]
      x$signal[k] <- as.integer(x$probability[k] >= x$cutoff[k]); x$signal_050[k] <- as.integer(x$probability[k] >= .5)
      x$target[k] <- 1L - x$target[k]; later[[m]] <- x
    }
    stopifnot(identical(build_vote(later, "validation")$vote[1:6], build_vote(mem, "validation")$vote[1:6]))
  })

  check("inconsistent or incomplete member exports are refused", {
    bad <- function(m, f) { x <- mem; x[[m]] <- f(x[[m]]); fails(build_vote(x, "validation")) }
    stopifnot(bad("b", function(x) x[-5, ]))                                             # missing origin / gap
    stopifnot(bad("b", function(x) rbind(x, x[3, ])))                                    # duplicate key
    stopifnot(bad("c", function(x) { x$target[2] <- 1L - x$target[2]; x }))             # target disagrees
    stopifnot(bad("c", function(x) { x$classification_eligible[4] <- FALSE; x }))       # eligibility disagrees
    stopifnot(bad("a", function(x) { x$signal[1] <- 1L - x$signal[1]; x }))             # signal != p >= cutoff
    stopifnot(bad("a", function(x) { x$cutoff[3] <- .7; x$signal <- as.integer(x$probability >= x$cutoff); x }))  # cutoff not frozen
    stopifnot(bad("a", function(x) { x$probability[1] <- 1.2; x }), bad("a", function(x) x[, -4]))
  })

  check("vote share >= (floor(M/2)+1)/M reproduces the vote, so evaluators can recompute it", {
    for (M in 1:6) {
      S <- as.matrix(expand.grid(rep(list(0:1), M)))
      stopifnot(identical(as.integer(rowMeans(S) >= (floor(M / 2) + 1) / M), majority_vote(S)))
    }
  })

  # ---------------- probability ensembles ----------------
  set.seed(11)
  n <- 400
  Q <- cbind(probit=runif(n, .2, .9), enet=runif(n, .2, .9), rf=runif(n, .2, .9))
  yq <- rbinom(n, 1, .7 * Q[, "probit"] + .3 * Q[, "enet"])
  sse <- function(w, a) sum((yq - a - Q %*% w)^2)

  check("GR with intercept matches a brute-force search over the simplex", {
    grid <- expand.grid(a=seq(0, 1, by=.005), b=seq(0, 1, by=.005)); grid <- grid[grid$a + grid$b <= 1 + 1e-12, ]
    W <- cbind(grid$a, grid$b, 1 - grid$a - grid$b)
    best <- apply(W, 1, function(w) sse(w, mean(yq - Q %*% w)))   # optimal free intercept for fixed weights
    g <- gr_weights(Q, yq, intercept=TRUE)
    stopifnot(g$sse <= min(best) + 1e-9, max(abs(g$weights - W[which.min(best), ])) < .01, g$status == "ok")
    stopifnot(isTRUE(all.equal(g$sse, sse(g$weights, g$intercept))))
  })

  check("GR weights are nonnegative, sum to one, recover the truth and can be zero", {
    m <- 20000; Z <- cbind(probit=runif(m, .1, .9), enet=runif(m, .1, .9), rf=runif(m, .1, .9))
    g <- gr_weights(Z, rbinom(m, 1, .7 * Z[, "probit"] + .3 * Z[, "enet"]), intercept=TRUE)
    stopifnot(all(g$weights >= 0), abs(sum(g$weights) - 1) < 1e-12, abs(g$weights[["probit"]] - .7) < .03,
      abs(g$weights[["enet"]] - .3) < .03, g$weights[["rf"]] < .03, abs(g$intercept) < .02)
    # When the unconstrained optimum needs a negative weight, the constraint binds at a corner.
    Y <- 1.5 * Z[, "probit"] - .5 * Z[, "enet"]
    stopifnot(isTRUE(all.equal(unname(gr_weights(Z[, c("probit", "enet")], Y, intercept=TRUE)$weights), c(1, 0))))
  })

  check("the intercept absorbs a common bias without changing the weights", {
    g0 <- gr_weights(Q, yq); g1 <- gr_weights(Q + .05, yq)
    stopifnot(isTRUE(all.equal(g0$weights, g1$weights, tolerance=1e-8)), abs(g1$intercept - (g0$intercept - .05)) < 1e-8)
  })

  check("a constant member (e.g. always-high elastic net) and collinear members are handled", {
    C <- cbind(Q[, c("probit", "rf")], enet=rep(.8, n))
    g <- gr_weights(C, yq); stopifnot(abs(sum(g$weights) - 1) < 1e-12, all(g$weights >= 0))
    D <- cbind(a=Q[, 1], b=Q[, 1], c=Q[, 3]); gd <- gr_weights(D, yq)
    stopifnot(abs(sum(gd$weights) - 1) < 1e-12, isTRUE(all.equal(gd$sse, gr_weights(D[, c("a", "c")], yq)$sse)))
  })

  check("solver failure falls back to equal weights with zero intercept and is flagged", {
    f <- gr_weights_safe(Q, yq[-1])                       # misaligned input forces an error
    stopifnot(grepl("^fallback_equal_weights", f$status), isTRUE(all.equal(unname(f$weights), rep(1/3, 3))), f$intercept == 0)
    stopifnot(gr_weights_safe(Q, yq)$status == "ok")
  })

  check("GR predictions are clipped to [0, 1], flagged, and matched to members by name", {
    fit <- list(weights=c(probit=.5, enet=.5, rf=0), intercept=.6)
    g <- gr_predict(Q[, c("rf", "enet", "probit")], fit)   # columns in a different order
    raw <- .6 + .5 * Q[, "probit"] + .5 * Q[, "enet"]
    stopifnot(isTRUE(all.equal(g$raw, raw)), all(g$probability <= 1), identical(g$clipped, raw > 1), any(g$clipped))
    stopifnot(fails(gr_predict(Q[, 1:2], fit)))
    stopifnot(isTRUE(all.equal(equal_average(Q), (Q[, 1] + Q[, 2] + Q[, 3]) / 3)))
  })

  check("yearly GR uses only outcomes known by the preceding December and refits once a year", {
    set.seed(7); n <- 120; o <- 2005L * 12L + 0:(n - 1L); oe <- o + 3L    # h = 3 style labels
    Z <- cbind(a = runif(n, .2, .9), b = runif(n, .2, .9)); yz <- rbinom(n, 1, Z[, "a"])
    pred <- 49:n                                            # origins from 2009-01
    g <- gr_yearly(Z, yz, o, oe, rep(TRUE, n), pred)
    yy <- yz; late <- oe > 2011L * 12L - 1L; yy[late] <- 1L - yy[late]   # flip every label not known by 2010-12
    g2 <- gr_yearly(Z, yy, o, oe, rep(TRUE, n), pred)
    k <- (o[pred] %/% 12L) <= 2011L
    stopifnot(isTRUE(all.equal(g$raw[k], g2$raw[k])), !isTRUE(all.equal(g$raw[!k], g2$raw[!k])))
    stopifnot(nrow(g$weights) == 2L * length(unique(o[pred] %/% 12L)))
    w09 <- g$weights[g$weights$year == 2009, ]
    stopifnot(w09$n_rows[1] == sum(oe <= 2009L * 12L - 1L), w09$last_known_outcome[1] == "2008-12")
    ref <- gr_weights(Z[oe <= 2009L * 12L - 1L, ], yz[oe <= 2009L * 12L - 1L])
    stopifnot(isTRUE(all.equal(unname(w09$weight), unname(ref$weights))))
    # Rule 61: January 2014, h = 6, archive from 2005-01 -> last usable origin 2013-06, 102 rows.
    o6 <- 2005L * 12L + 0:119; g6 <- gr_yearly(Z[1:120, ], yz[1:120], o6, o6 + 6L, rep(TRUE, 120), 109:120)
    stopifnot(g6$weights$n_rows[1] == 102L, g6$weights$last_known_outcome[1] == "2013-12")
  })

  if (!is.null(json_path)) {
    dir.create(dirname(json_path), recursive=TRUE, showWarnings=FALSE)
    jsonlite::write_json(list(passed=length(passed), tests=as.list(passed), status="passed"), json_path, auto_unbox=TRUE, pretty=TRUE)
  }
  message(length(passed), " ensemble test groups passed.")
  invisible(passed)
}
