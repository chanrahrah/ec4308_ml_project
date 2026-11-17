# Ensembles built from the saved member forecasts; no member model is refitted.
# Run after the rf, elastic_net and boosting runs: Rscript --vanilla Scripts/ensemble/run_ensemble.R
script <- c(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE),
            tryCatch(sys.frame(1)$ofile, error = function(e) NULL))
source(file.path(if (length(script)) dirname(script[1]) else "Scripts/ensemble", "..", "utils", "utils.R"))
for (file in c("model.R", "tests.R")) source(file.path(SCRIPTS, "ensemble", file))

cfg <- ensemble_config()
checks <- run_checks(ensemble_checks())
member_files <- file.path(ROOT, "Results", cfg$members, "forecasts.csv")
ensure(all(file.exists(member_files)), paste("Run the member models first. Missing:", paste(member_files[!file.exists(member_files)], collapse = ", ")))
exports <- setNames(lapply(member_files, utils::read.csv, stringsAsFactors = FALSE), names(cfg$members))
members <- names(exports); M <- length(members)
ensure(M %% 2L == 1L, "Use an odd number of members so the vote cannot tie")

# Merge the members' forecasts for each period and check they cover every origin.
V <- list()
for (period in c("validation", "test")) {
  V[[period]] <- build_vote(lapply(exports, function(x) x[x$period == period, ]), period)
  range <- month_number(c(cfg[[paste0(period, "_start")]], cfg[[paste0(period, "_end")]]))
  for (h in cfg$horizons)
    ensure(identical(sort(V[[period]]$base$origin_id[V[[period]]$base$h == h]), range[1]:range[2]),
           paste(period, "h", h, ": members do not cover every origin"))
}
ensure(identical(V$validation$cutoffs[order(V$validation$cutoffs$member, V$validation$cutoffs$h), "cutoff"],
                 V$test$cutoffs[order(V$test$cutoffs$member, V$test$cutoffs$h), "cutoff"]),
       "A member's cutoff differs between validation and test")

# Static GR weights: fitted once per horizon on eligible validation labels.
vb <- V$validation$base; gr <- list(); weights <- list()
for (h in cfg$horizons) {
  rows <- which(vb$h == h & vb$classification_eligible & !is.na(vb$target))
  fit <- gr_weights_safe(V$validation$P[rows, , drop = FALSE], vb$target[rows], cfg$gr_intercept)
  if (fit$status != "ok") warning("GR solver failed at h = ", h, "; using equal weights")
  fit$last_outcome <- as.character(vb$outcome_end[rows[which.max(vb$origin_id[rows])]])
  gr[[as.character(h)]] <- fit
  weights[[length(weights) + 1L]] <- data.frame(fit = "validation", h = h, year = NA_integer_, member = names(fit$weights),
    weight = unname(fit$weights), intercept = fit$intercept, n_rows = fit$n, last_outcome = fit$last_outcome, status = fit$status)
}

# Ensemble forecasts for each period.
vote_cut <- (floor(M / 2) + 1) / M
forecasts <- list()
for (period in names(V)) {
  v <- V[[period]]; b <- v$base
  share <- rowMeans(v$S)
  ensure(identical(as.integer(share >= vote_cut), v$vote), "Vote share does not reproduce the vote")
  p_equal <- equal_average(v$P)
  g <- list(raw = numeric(nrow(b)), probability = numeric(nrow(b)), clipped = logical(nrow(b)))
  for (h in cfg$horizons) {
    i <- b$h == h; gh <- gr_predict(v$P[i, , drop = FALSE], gr[[as.character(h)]])
    g$raw[i] <- gh$raw; g$probability[i] <- gh$probability; g$clipped[i] <- gh$clipped
  }
  member_p <- setNames(as.data.frame(v$P), paste0("p_", members))
  make <- function(model, p, cutoff, signal, cutoff_050, signal_050, config_id)
    cbind(data.frame(h = b$h, model = model, period = period, origin = month_string(b$origin_id), probability = p,
      outcome_end = b$outcome_end, target = b$target, classification_eligible = b$classification_eligible,
      config_id = config_id, cutoff = cutoff, signal = signal, cutoff_050 = cutoff_050, signal_050 = signal_050,
      vote_count = NA_integer_, vote_count_050 = NA_integer_, raw_score = NA_real_, clipped = NA,
      n_train = NA_integer_, last_train_outcome = NA_character_, stringsAsFactors = FALSE), member_p)
  vote <- make("ens_majority_vote", share, vote_cut, v$vote, vote_cut, v$vote_050, "majority_vote|ties_to_bills")
  vote$vote_count <- rowSums(v$S); vote$vote_count_050 <- rowSums(v$S050)
  half <- function(p) as.integer(p >= cfg$probability_cutoff)
  equal <- make("ens_equal", p_equal, .5, half(p_equal), .5, half(p_equal), "equal_probability|cutoff=0.50")
  static <- make("ens_gr", g$probability, .5, half(g$probability), .5, half(g$probability), "gr_intercept_clipped|validation_fit|cutoff=0.50")
  static$raw_score <- g$raw; static$clipped <- g$clipped
  static$n_train <- vapply(b$h, function(h) gr[[as.character(h)]]$n, 0L)
  static$last_train_outcome <- vapply(b$h, function(h) gr[[as.character(h)]]$last_outcome, "")
  forecasts[[period]] <- list(vote, equal, static)
}

# Robustness check: GR refitted every January on all outcomes known by the previous December (test only).
pool_b <- rbind(V$validation$base, V$test$base); pool_P <- rbind(V$validation$P, V$test$P)
pool_period <- rep(c("validation", "test"), c(nrow(V$validation$base), nrow(V$test$base)))
yearly <- forecasts$test[[3]]
yearly$model <- "ens_gr_yearly"; yearly$config_id <- "gr_intercept_clipped|refit_each_january|cutoff=0.50"
for (h in cfg$horizons) {
  j <- which(pool_b$h == h); pred <- which(pool_period[j] == "test")
  g <- gr_yearly(pool_P[j, , drop = FALSE], pool_b$target[j], pool_b$origin_id[j], month_number(pool_b$outcome_end[j]),
                 !is.na(pool_b$target[j]), pred, cfg$gr_intercept, cfg$gr_min_rows)
  ensure(all(g$weights$status == "ok"), paste("Yearly GR fell back to equal weights at h =", h))
  i <- which(yearly$h == h); i <- i[match(pool_b$origin_id[j][pred], month_number(yearly$origin[i]))]
  yearly$raw_score[i] <- g$raw; yearly$probability[i] <- g$probability; yearly$clipped[i] <- g$clipped
  first <- g$weights[!duplicated(g$weights$year), ]
  year <- month_number(yearly$origin[i]) %/% 12L
  yearly$n_train[i] <- first$n_rows[match(year, first$year)]
  yearly$last_train_outcome[i] <- first$last_outcome[match(year, first$year)]
  weights[[length(weights) + 1L]] <- cbind(fit = "yearly", h = h, g$weights)
}
yearly$signal <- yearly$signal_050 <- as.integer(yearly$probability >= cfg$probability_cutoff)
forecasts$test[[4]] <- yearly
forecasts <- bind_rows(unlist(forecasts, recursive = FALSE))

# Score the vote at both member-cutoff rules and the probability ensembles at 0.50.
vote_rows <- forecasts$model == "ens_majority_vote"
scored <- rbind(cbind(forecasts[vote_rows, ], cutoff_rule = "validation_tuned"),
                cbind(transform(forecasts[vote_rows, ], probability = vote_count_050 / M, signal = signal_050), cutoff_rule = "fixed_050"),
                cbind(forecasts[!vote_rows, ], cutoff_rule = "fixed_050"))
returns <- load_return_data()
results <- evaluate(scored, returns, cfg)
# A vote share is not a probability, so probability scores are blank for the vote.
vote_metrics <- results$forecast_metrics$model == "ens_majority_vote"
results$forecast_metrics[vote_metrics, c("Brier", "log_loss", "mean_probability", "calibration_gap", "accuracy_at_050")] <- NA
results$calibration <- results$calibration[results$calibration$model != "ens_majority_vote", ]

out <- results_dir("ensemble")
write_tables(out, c(list(forecasts = forecasts, gr_weights = bind_rows(weights)), results))
write_run_info(out, cfg, checks, inputs = c(member_files, RETURN_FILES))
message("Results saved in ", out)
