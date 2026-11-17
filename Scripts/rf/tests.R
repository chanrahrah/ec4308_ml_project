# Random forest checks: class order, the grid, no look-ahead and seeded determinism.
rf_checks <- function(cfg, grid) {
  quick <- cfg; quick$num_trees <- 60L
  theta <- rf_theta("rf_mtry11_node10_boot", grid, quick)
  checks <- list(
    "Probability is for class 1" = function() {
      set.seed(1); x <- matrix(stats::rnorm(400), 200, 2, dimnames = list(NULL, c("a", "b")))
      fit <- rf_fit(x, as.integer(x[, 1] > 0), list(num_trees = 100L, mtry = 1L, min_node_size = 5L, replace = TRUE, sample_fraction = 1), 11L)
      p <- rf_prob(fit, rbind(c(a = 3, b = 0), c(a = -3, b = 0)))
      ensure(p[1] > .9 && p[2] < .1, "Probability column is not class 1")
    },
    "Grid has 18 candidates and ranger receives each setting" = function() {
      ensure(nrow(grid) == 18L && !anyDuplicated(grid$config_id) && sum(grid$is_ranger_default) == 1L, "Unexpected grid")
      set.seed(2); x <- matrix(stats::rnorm(200 * 132), 200, 132, dimnames = list(NULL, paste0("v", 1:132)))
      y <- rep(0:1, 100)
      fit <- rf_fit(x, y, list(mtry = 33L, min_node_size = 30L, replace = FALSE, sample_fraction = .5, num_trees = 60L), 1L)
      ensure(fit$mtry == 33L && fit$min.node.size == 30L && !fit$replace && fit$num.trees == 60L, "Settings not passed to ranger")
      default <- ranger::ranger(x = x, y = factor(y, levels = c("0", "1")), probability = TRUE,
                                num.trees = 60L, seed = 1L, num.threads = 1L, verbose = FALSE)
      ensure(identical(rf_prob(rf_fit(x, y, rf_theta(grid$config_id[grid$is_ranger_default], grid, quick), 1L), x[1:5, ]),
                       rf_prob(default, x[1:5, ])), "Default candidate differs from ranger's defaults")
    })
  for (h in cfg$horizons) checks[[paste("Future data cannot change a forecast; seeds are reproducible h", h)]] <- local({ h <- h; function() {
    d <- rf_load(h, cfg)
    predict_at <- function(d, i) rf_task(rf_task_data(d), rf_make_tasks(d, i, h, theta, quick)[[1]])$p
    perturbation_check(d, h, predict_at, tol = 0)
    i <- which(d$ym == "2008-06")
    task <- rf_make_tasks(d, i, h, theta, quick)[[1]]; task$donors <- c(i, i)
    a <- rf_task(rf_task_data(d), task)
    ensure(identical(a$p, rf_task(rf_task_data(d), task)$p), "Same seed gave a different forest")
    near(a$perturbed, matrix(a$p, length(d$groups), 2L), 1e-14)
    task$seed <- task$seed + 1L
    ensure(!identical(rf_task(rf_task_data(d), task)$p, a$p), "Seed is ignored")
  }})
  checks
}
