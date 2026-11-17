# Model selection, forecast scoring and portfolio evaluation shared by every model.

# Completed labels that belong to the forecast's own split.
eligible <- function(f) f[f$classification_eligible & !is.na(f$target), ]

# Validation scores used to compare hyperparameter candidates.
validation_scores <- function(f) {
  e <- eligible(f)
  data.frame(n_eligible = nrow(e), validation_Brier = mean((e$probability - e$target)^2),
             Brier_prevalence = mean((e$train_prevalence - e$target)^2), AUC = auc(e$target, e$probability))
}

# Lowest validation Brier; exact ties go to the candidate that comes first in tie_order.
best_candidate <- function(brier, tie_order = seq_along(brier)) {
  ensure(all(is.finite(brier)), "Undefined validation Brier")
  tied <- which(brier <= min(brier) + 1e-12)
  tie_order[tie_order %in% tied][1L]
}

# Trade the validation forecasts at every cutoff (10 bp cost) and flag the best one.
select_cutoff <- function(f, returns, h, cfg) {
  f <- f[order(f$origin), ]
  inp <- portfolio_inputs(f, returns, h, cfg$validation_end, cfg$common_rf_security)
  grid <- bind_rows(lapply(cfg$cutoffs, function(cutoff) {
    p <- portfolio_path(inp, as.integer(f$probability >= cutoff), cfg$primary_cost, cfg)
    cbind(data.frame(h = h, model = f$model[1], cutoff = cutoff), portfolio_metrics(p, cfg$risk_aversion))
  }))
  grid$selected <- seq_len(nrow(grid)) == choose_cutoff(grid, cfg$cutoff_objective)
  grid
}

# Add the tuned cutoff, the fixed 0.50 cutoff and the trading signal from each.
add_signals <- function(f, cutoff, cfg) {
  f$cutoff <- cutoff
  f$signal <- as.integer(f$probability >= cutoff)
  f$cutoff_050 <- cfg$comparison_cutoff
  f$signal_050 <- as.integer(f$probability >= cfg$comparison_cutoff)
  f
}

# One block of rows per cutoff rule, so each rule is evaluated as its own strategy.
by_rule <- function(f) rbind(cbind(f, cutoff_rule = "validation_tuned"),
  cbind(transform(f, cutoff = cutoff_050, signal = signal_050), cutoff_rule = "fixed_050"))

# Probability and classification scores at a given cutoff.
forecast_metrics <- function(f, cutoff) {
  f <- eligible(f)
  y <- f$target; p <- f$probability; pred <- as.integer(p >= cutoff)
  pc <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  brier <- mean((p - y)^2)
  prevalence <- if ("train_prevalence" %in% names(f)) mean((f$train_prevalence - y)^2) else NA_real_
  ratio <- function(a, b) if (b == 0) NA_real_ else a / b
  data.frame(n = length(y), cutoff = cutoff, event_rate = mean(y), Brier = brier,
    Brier_prevalence = prevalence, Brier_skill = if (isTRUE(prevalence > 0)) 1 - brier / prevalence else NA_real_,
    log_loss = -mean(y * log(pc) + (1 - y) * log1p(-pc)), AUC = auc(y, p),
    accuracy = mean(pred == y), accuracy_at_050 = mean((p >= .5) == y),
    precision = ratio(sum(pred == 1 & y == 1), sum(pred == 1)),
    recall = ratio(sum(pred == 1 & y == 1), sum(y == 1)),
    specificity = ratio(sum(pred == 0 & y == 0), sum(y == 0)),
    mean_probability = mean(p), calibration_gap = mean(p) - mean(y))
}

# Average forecast against the realised event rate within probability deciles.
calibration_table <- function(f) {
  f <- eligible(f)
  bins <- cut(f$probability, breaks = seq(0, 1, by = .1), include.lowest = TRUE)
  bind_rows(lapply(levels(bins), function(bin) {
    ids <- which(bins == bin)
    data.frame(bin = bin, n = length(ids),
      mean_probability = if (length(ids)) mean(f$probability[ids]) else NA_real_,
      event_rate = if (length(ids)) mean(f$target[ids]) else NA_real_)
  }))
}

# Trade every (model, cutoff rule) block at every cost, by period and horizon.
# benchmarks = TRUE adds all-stocks, all-bills and the fixed 50/50 mix.
evaluate <- function(scored, returns, cfg, benchmarks = FALSE) {
  paths <- metrics <- calibration <- list()
  for (period in c("validation", "test")) for (h in cfg$horizons) {
    s <- scored[scored$period == period & scored$h == h, ]
    origins <- sort(unique(s$origin))
    end <- cfg[[paste0(period, "_end")]]
    inp <- portfolio_inputs(data.frame(origin = origins), returns, h, end, cfg$common_rf_security)
    strategies <- list()
    for (key in unique(paste(s$model, s$cutoff_rule))) {
      g <- s[paste(s$model, s$cutoff_rule) == key, ]
      g <- g[order(g$origin), ]
      ensure(identical(g$origin, origins) && length(unique(g$cutoff)) == 1L,
             "Each strategy needs every origin and one cutoff")
      strategies[[key]] <- list(strategy = g$model[1], rule = g$cutoff_rule[1], cutoff = g$cutoff[1], signals = g$signal)
      metrics[[length(metrics) + 1L]] <- cbind(data.frame(period = period, h = h, model = g$model[1],
        cutoff_rule = g$cutoff_rule[1]), forecast_metrics(g, g$cutoff[1]))
    }
    first <- s[!duplicated(s[, c("model", "origin")]), ]
    for (m in unique(first$model))
      calibration[[length(calibration) + 1L]] <- cbind(data.frame(period = period, h = h, model = m),
                                                        calibration_table(first[first$model == m, ]))
    if (benchmarks) {
      n <- length(origins)
      strategies$all_stocks <- list(strategy = "all_stocks", signals = rep(1L, n))
      strategies$all_bills <- list(strategy = "all_bills", signals = rep(0L, n))
      strategies$fixed_50_50 <- list(strategy = "fixed_50_50", signals = rep(0L, n), weight = cfg$fixed_benchmark_weight)
    }
    for (st in strategies) for (cost in cfg$costs) {
      p <- portfolio_path(inp, st$signals, cost, cfg, st$weight)
      ensure(max(month_number(p$month)) <= month_number(end), "Portfolio return crosses split boundary")
      paths[[length(paths) + 1L]] <- cbind(data.frame(period = period, h = h, strategy = st$strategy,
        cutoff_rule = if (is.null(st$rule)) "not_applicable" else st$rule,
        cutoff = if (is.null(st$cutoff)) NA_real_ else st$cutoff, cost_rate = cost), p)
    }
  }
  paths <- bind_rows(paths)
  list(performance = summarise_paths(paths, cfg), monthly_portfolios = paths,
       forecast_metrics = bind_rows(metrics), calibration = bind_rows(calibration))
}

# Performance over each full split, and over a common window that starts once every horizon trades.
summarise_paths <- function(paths, cfg) {
  tags <- c("period", "h", "strategy", "cutoff_rule", "cutoff", "cost_rate")
  out <- list()
  for (window in c("full", "common")) for (period in c("validation", "test")) {
    p <- paths[paths$period == period, ]
    if (window == "common")
      p <- p[p$month >= month_string(month_number(cfg[[paste0(period, "_start")]]) + max(cfg$horizons)), ]
    groups <- interaction(p$h, p$strategy, p$cutoff_rule, p$cost_rate, drop = TRUE)
    for (g in unique(groups)) {
      q <- p[groups == g, ]
      out[[length(out) + 1L]] <- cbind(q[1, tags], window = window, portfolio_metrics(q, cfg$risk_aversion))
    }
  }
  bind_rows(out)
}
