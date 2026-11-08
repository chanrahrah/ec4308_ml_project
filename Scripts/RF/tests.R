# Focused information-set, determinism and evaluator checks; run before any model fitting.
run_rf_tests <- function(root, cfg) {
  passed <- character()
  check <- function(name, expr) { force(expr); passed <<- c(passed, name) }
  near <- function(a, b, tol = 1e-10) ensure(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)), "Numerical test failed")
  quick <- cfg; quick$num_trees <- 60L             # small forests keep the suite fast
  theta <- rf_theta(quick, "rf_mtry11_node10_boot")

  check("Shared evaluator: professor's 12-month six-cohort example", {
    signals <- c(1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,0,1)
    market <- c(.03,-.02,.04,-.06,.02,.05,-.08,.03,-.01,.06,-.04,.02)
    w <- vapply(1:12, function(i) mean(signals[i:(i+5)]), numeric(1))
    p <- simulate_weights(w, market, (1-w)*.003, 0)
    value <- 100; expected <- numeric(12)
    for (i in 1:12) { value <- w[i]*value*(1+market[i]) + (1-w[i])*value*1.003; expected[i] <- value }
    near(p$wealth*100, expected)
  })
  check("Shared evaluator: exact post-cost target weights", {
    for (s in c(0, 30, 100)) for (w in c(0, .5, 1)) for (c in c(0, .001, .002)) {
      t <- rebalance_trade(100, s, w, c)
      near((s + t$trade) / t$capital, w); near(t$cost, abs(t$trade) * c)
    }
  })
  check("Class ordering: probability is P(Y = 1)", {
    set.seed(1); x <- matrix(stats::rnorm(400), 200, 2, dimnames = list(NULL, c("a", "b")))
    y <- as.integer(x[, 1] > 0)
    fit <- rf_fit(x, y, list(num_trees = 100L, mtry = 1L, min_node_size = 5L, replace = TRUE, sample_fraction = 1), 11L)
    p <- rf_prob(fit, rbind(c(a = 3, b = 0), c(a = -3, b = 0)))
    ensure(p[1] > .9 && p[2] < .1, "Probability column is not the Y = 1 class")
  })
  check("Declared 18-candidate hyperparameter grid", {
    g <- cfg$hyper_grid
    ensure(nrow(g) == 18L && !anyDuplicated(g$config_id) && sum(g$is_ranger_default) == 1L, "Unexpected grid")
    ensure(identical(sort(unique(g$mtry)), c(4L, 11L, 33L)) && identical(sort(unique(g$min_node_size)), c(10L, 30L, 60L)),
           "Unexpected grid values")
    ensure(cfg$num_trees == 1000L, "Expected 1000 trees")
    set.seed(2); x <- matrix(stats::rnorm(200 * 132), 200, 132, dimnames = list(NULL, paste0("v", 1:132)))
    y <- rep(0:1, 100)
    sub <- list(config_id = "check", mtry = 33L, min_node_size = 30L, replace = FALSE, sample_fraction = .5, num_trees = 60L)
    fit <- rf_fit(x, y, sub, 1L)
    ensure(fit$mtry == 33L && fit$min.node.size == 30L && !fit$replace && fit$num.trees == 60L,
           "ranger did not receive the candidate's settings")
    # The grid's ranger-default row must reproduce ranger's own defaults exactly.
    d0 <- rf_theta(quick, cfg$hyper_grid$config_id[cfg$hyper_grid$is_ranger_default])
    implicit <- ranger::ranger(x = x, y = factor(y, levels = c("0", "1")), probability = TRUE,
                               num.trees = quick$num_trees, seed = 1L, num.threads = 1L, verbose = FALSE)
    ensure(identical(rf_prob(rf_fit(x, y, d0, 1L), x[1:5, ]), rf_prob(implicit, x[1:5, ])),
           "Explicit default settings differ from ranger's defaults")
  })
  check("Brier tie rule favours the simpler candidate", {
    sc <- data.frame(Brier = c(.20, .19, .19 + 1e-13, .25), simplicity_rank = c(1, 3, 2, 4))
    ensure(rf_choose_candidate(sc) == 3L, "Tie should go to lower simplicity rank")
    sc$Brier[2] <- .18
    ensure(rf_choose_candidate(sc) == 2L, "Strictly lower Brier must win")
  })
  check("RF cutoff grid 0.30-0.90 by 0.05 contains the 0.50 comparator", {
    ensure(length(cfg$cutoffs) == 13L && isTRUE(all.equal(range(cfg$cutoffs), c(.3, .9))) &&
             isTRUE(all.equal(diff(cfg$cutoffs), rep(.05, 12))), "Unexpected cutoff grid")
    ensure(sum(abs(cfg$cutoffs - cfg$comparison_cutoff) < 1e-12) == 1L, "0.50 missing from grid")
  })
  check("Cutoff tie rule and undefined Sharpe exclusion (shared)", {
    g <- data.frame(cutoff = c(.2, .49, .51, 1), sharpe_matched_ladder = c(.1, .2, .2, NA))
    ensure(choose_cutoff(g) == 2L, "Cutoff tie rule failed")
    g$sharpe_matched_ladder <- NA_real_
    ensure(inherits(try(choose_cutoff(g), silent = TRUE), "try-error"), "Undefined candidates must stop")
  })
  returns <- load_return_data(root)
  for (h in cfg$horizons) {
    d <- rf_load_data(root, h, cfg)
    check(paste("Exactly the 132 frozen predictors and no metadata h", h), {
      ensure(ncol(d$X) == 132L && identical(colnames(d$X), d$feature_names), "Wrong feature set")
      ensure(!any(c("outcome_end", "target", "split", "equity_forward") %in% colnames(d$X)), "Metadata leaked")
    })
    check(paste("Boundary maturity and future-perturbation invariance h", h), {
      for (origin_string in c("2005-01", "2013-01")) {
        origin <- month_number(origin_string); i <- which(d$m == origin)
        rows <- rf_mature_rows(d$y, d$m, d$end, origin)
        ensure(max(d$m[rows]) == origin - h && max(d$end[rows]) == origin, "Wrong maturity cutoff")
        task <- rf_make_tasks(d, i, h, theta, quick)[[1]]
        a <- rf_task(rf_task_data(d), task)
        changed <- d
        unavailable <- which(is.na(d$end) | d$end > origin)
        changed$y[unavailable] <- 1L - ifelse(is.na(d$y[unavailable]), 0L, d$y[unavailable])
        changed$X[d$m > origin, ] <- changed$X[d$m > origin, ] + 1000
        b <- rf_task(rf_task_data(changed), task)
        near(a$p, b$p, 0); ensure(a$n_train == length(rows), "Wrong fitting sample")
        ensure(a$last_outcome <= origin, "Immature label in fit")
      }
    })
    check(paste("Seeded determinism and identity permutation h", h), {
      i <- which(d$ym == "2008-06")
      task <- rf_make_tasks(d, i, h, theta, quick)[[1]]
      task$donors <- c(i, i)
      a <- rf_task(rf_task_data(d), task); b <- rf_task(rf_task_data(d), task)
      ensure(identical(a$p, b$p), "Same seed gave a different forest")
      near(a$perturbed, matrix(a$p, length(d$groups), 2L), 1e-14)
      task$seed <- task$seed + 1L
      ensure(!identical(rf_task(rf_task_data(d), task)$p, a$p), "Seed is ignored")
    })
    check(paste("Passive identities on the RF forecast calendar h", h), {
      ids <- rf_period_rows(d, cfg, "validation")
      f <- data.frame(origin = d$ym[ids])
      inp <- portfolio_inputs(f, returns, h, cfg$validation_end, cfg$common_rf_security)
      ensure(inp$ym[1] == month_string(month_number(cfg$validation_start) + h), "Wrong warm-up")
      pstock <- portfolio_path(inp, rep(1L, nrow(f)), .002, cfg)
      pbill <- portfolio_path(inp, rep(0L, nrow(f)), .002, cfg)
      near(pstock$net_return, inp$market); near(pbill$net_return, inp$bill_reference)
      ensure(is.na(portfolio_metrics(pbill)$sharpe_matched_ladder), "All-bills matched Sharpe must be N/A")
    })
  }
  data.frame(test = passed, passed = TRUE)
}
