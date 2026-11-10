# Run from any directory: Rscript --vanilla Scripts/ensemble/run_ensemble.R
# Builds the three group ensembles (spec rules 50-58) from the members' saved rule-66 forecasts and
# scores them, the members and the all-stocks/all-bills benchmarks with the SHARED evaluator
# (Scripts/Baseline models/portfolio.R), using the baseline config's costs, splits and Sharpe references.
#   ens_majority_vote  members' frozen validation cutoffs; stocks if sum > M/2; ties to bills
#   ens_equal          mean member probability, fixed 0.50 cutoff
#   ens_gr             constrained Granger-Ramanathan with intercept, fitted once per horizon on
#                      eligible validation labels, clipped to [0, 1], fixed 0.50 cutoff
# Needs the two CRSP return files under Data/ (as for the baseline runner); no model is refitted.
args <- commandArgs(FALSE)
entry <- grep("^--file=", args, value = TRUE)
entry_path <- if (length(entry)) gsub("~+~", " ", sub("^--file=", "", entry[1]), fixed = TRUE) else ""
script_dir <- if (nzchar(entry_path)) dirname(normalizePath(entry_path, mustWork = TRUE)) else normalizePath("Scripts/ensemble")
root <- normalizePath(file.path(script_dir, "../.."))
.libPaths(c(file.path(root, ".R-library"), .libPaths()))
ensure_packages <- c("digest", "jsonlite")
missing <- ensure_packages[!vapply(ensure_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Run Scripts/00_setup.R first. Missing: ", paste(missing, collapse = ", "))
shared <- file.path(root, "Scripts", "Baseline models")
for (file in c("config.R", "portfolio.R")) source(file.path(shared, file))   # shared config helpers and evaluator
for (file in c("config.R", "ensemble.R", "tests.R")) source(file.path(script_dir, file))
cfg <- baseline_config(); ecfg <- ensemble_config()
out <- file.path(root, "Results", "ensemble")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
bind <- function(x) { x <- do.call(rbind, x); rownames(x) <- NULL; x }

# ---------------- Inputs and tests ----------------
member_files <- unlist(lapply(ecfg$members, function(m) file.path(root, c(m$validation, m$test))))
return_files <- file.path(root, c("Data/crsp_monthly_index_data/crsp_monthly_index_data.csv",
                                  "Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))
input_files <- c(member_files, return_files)
ensure(all(file.exists(input_files)), paste("Missing inputs:", paste(input_files[!file.exists(input_files)], collapse = ", ")))
hashes <- setNames(vapply(input_files, file_hash, character(1)), substring(input_files, nchar(root) + 2L))
code_files <- sort(list.files(script_dir, pattern = "[.]R$", full.names = TRUE))
code_hashes <- setNames(vapply(code_files, file_hash, character(1)), basename(code_files))
message("Running ensemble tests...")
tests <- run_ensemble_tests()
write_csv(data.frame(test = tests, status = "passed"), file.path(out, "tests.csv"))

members <- names(ecfg$members); M <- length(members)
ensure(M %% 2L == 1L, "Use an odd number of members so the vote has no ties (ties would go to bills)")
read_member <- function(path) utils::read.csv(file.path(root, path), stringsAsFactors = FALSE)
V <- list()
for (period in c("validation", "test")) {
  V[[period]] <- build_vote(setNames(lapply(ecfg$members, function(m) read_member(m[[period]])), members), period)
  rng <- month_number(if (period == "validation") c(cfg$validation_start, cfg$validation_end) else c(cfg$test_start, cfg$test_end))
  for (h in cfg$horizons) ensure(identical(sort(V[[period]]$base$origin_id[V[[period]]$base$h == h]), rng[1]:rng[2]),
    paste(period, "h =", h, ": members do not cover the full origin range"))
}
cut_val <- V$validation$cutoffs; cut_test <- V$test$cutoffs
ensure(identical(cut_val[order(cut_val$member, cut_val$h), "cutoff"], cut_test[order(cut_test$member, cut_test$h), "cutoff"]),
  "A member's tuned cutoff differs between its validation and test exports")

# ---------------- GR weights: once per horizon, eligible validation labels only ----------------
vb <- V$validation$base; gr <- list(); weight_rows <- list()
for (h in cfg$horizons) {
  rows <- which(vb$h == h & vb$classification_eligible & !is.na(vb$target))
  fit <- gr_weights_safe(V$validation$P[rows, , drop = FALSE], vb$target[rows], intercept = ecfg$gr_intercept)
  if (fit$status != "ok") warning("GR solver failed at h = ", h, "; equal weights with zero intercept used")
  gr[[as.character(h)]] <- fit
  last <- rows[which.max(vb$origin_id[rows])]
  weight_rows[[length(weight_rows) + 1]] <- data.frame(h = h, member = names(fit$weights), weight = unname(fit$weights),
    intercept = fit$intercept, n_rows = fit$n, sse = fit$sse, status = fit$status,
    first_origin = month_string(min(vb$origin_id[rows])), last_origin = month_string(vb$origin_id[last]),
    last_outcome = as.character(vb$outcome_end[last]))
  message(sprintf("h=%d GR weights (n=%d): %s; intercept %.3f", h, fit$n,
    paste(sprintf("%s %.3f", names(fit$weights), fit$weights), collapse = ", "), fit$intercept))
}
gr_weight_table <- bind(weight_rows)

# ---------------- Ensemble forecasts (rule-66 exports) ----------------
vote_cut <- (floor(M / 2) + 1) / M
ensemble_forecasts <- list(); clipping <- list()
for (period in names(V)) {
  v <- V[[period]]; b <- v$base
  share <- rowMeans(v$S); share050 <- rowMeans(v$S050)
  ensure(identical(as.integer(share >= vote_cut), v$vote), "Vote share does not reproduce the vote")
  p_equal <- equal_average(v$P); raw <- prob <- numeric(nrow(b)); clip <- logical(nrow(b))
  for (h in cfg$horizons) {
    i <- b$h == h; g <- gr_predict(v$P[i, , drop = FALSE], gr[[as.character(h)]])
    raw[i] <- g$raw; prob[i] <- g$probability; clip[i] <- g$clipped
    clipping[[paste(period, h)]] <- data.frame(period = period, h = h, n = sum(i), n_clipped = sum(g$clipped),
      n_below_0 = sum(g$raw < 0), n_above_1 = sum(g$raw > 1))
  }
  key <- data.frame(h = b$h, period = period, origin = month_string(b$origin_id))
  tail_cols <- data.frame(outcome_end = b$outcome_end, target = b$target, classification_eligible = b$classification_eligible,
    fit_information_cutoff = month_string(b$origin_id))
  mem <- as.data.frame(v$P); names(mem) <- paste0("p_", members)
  for (m in members) { mem[[paste0("signal_", m)]] <- v$S[, m]; mem[[paste0("signal_050_", m)]] <- v$S050[, m] }
  mk <- function(model, p, cutoff, signal, cutoff050, signal050, config, extra = NULL) {
    x <- cbind(key[, c("h")], data.frame(model = model), key[, c("period", "origin")], data.frame(probability = p), tail_cols,
      data.frame(last_train_origin = NA, last_train_outcome = NA, n_train = NA, config_id = config,
                 cutoff = cutoff, signal = signal, cutoff_050 = cutoff050, signal_050 = signal050))
    names(x)[1] <- "h"
    if (!is.null(extra)) x <- cbind(x, extra)
    cbind(x, mem)
  }
  ids <- paste(members, collapse = ",")
  ensemble_forecasts[[period]] <- list(
    ens_majority_vote = mk("ens_majority_vote", share, vote_cut, v$vote, vote_cut, v$vote_050,
      paste0("majority_vote|", ids, "|ties_to_bills|probability=vote_share"),
      data.frame(vote_count = rowSums(v$S), vote_count_050 = rowSums(v$S050))),
    ens_equal = mk("ens_equal", p_equal, ecfg$probability_cutoff, as.integer(p_equal >= ecfg$probability_cutoff),
      ecfg$probability_cutoff, as.integer(p_equal >= ecfg$probability_cutoff), paste0("equal_probability|", ids, "|cutoff=0.50")),
    ens_gr = mk("ens_gr", prob, ecfg$probability_cutoff, as.integer(prob >= ecfg$probability_cutoff),
      ecfg$probability_cutoff, as.integer(prob >= ecfg$probability_cutoff),
      paste0("gr_intercept_clipped|", ids, "|validation_fit_once|cutoff=0.50"), data.frame(raw_score = raw, clipped = clip)))
  ensemble_forecasts[[period]]$ens_gr$n_train <- vapply(b$h, function(h) gr[[as.character(h)]]$n, 0L)
  ensemble_forecasts[[period]]$ens_gr$last_train_outcome <- vapply(b$h, function(h)
    gr_weight_table$last_outcome[gr_weight_table$h == h][1], "")
}

# ---------------- Supplementary: expanding-window GR refitted every January (test period only) ----------------
# History = validation + test member forecasts from 2005-01; at January of year Y it uses every eligible
# row with a realised target whose outcome_end <= December of Y-1 (rules 60-61; validation rows whose outcome crosses into 2013
# are ineligible for validation scoring but are known history once realised). Not a rule-58 ensemble: reported alongside, chosen by nothing.
yearly_weights <- list()
{
  pool_b <- rbind(V$validation$base, V$test$base); pool_P <- rbind(V$validation$P, V$test$P)
  pool_period <- rep(c("validation", "test"), c(nrow(V$validation$base), nrow(V$test$base)))
  tb <- V$test$base; raw <- prob <- numeric(nrow(tb)); clip <- logical(nrow(tb)); ntr <- integer(nrow(tb)); lto <- character(nrow(tb))
  for (h in cfg$horizons) {
    j <- which(pool_b$h == h); pred <- which(pool_period[j] == "test")
    g <- gr_yearly(pool_P[j, , drop = FALSE], pool_b$target[j], pool_b$origin_id[j], to_month_id(pool_b$outcome_end[j]),
                   !is.na(pool_b$target[j]), pred, intercept = ecfg$gr_intercept)
    i <- which(tb$h == h); i <- i[match(pool_b$origin_id[j][pred], tb$origin_id[i])]
    raw[i] <- g$raw; prob[i] <- g$probability; clip[i] <- g$clipped
    yrs <- tb$origin_id[i] %/% 12L; wk <- g$weights[!duplicated(g$weights$year), ]
    ntr[i] <- wk$n_rows[match(yrs, wk$year)]; lto[i] <- wk$last_known_outcome[match(yrs, wk$year)]
    yearly_weights[[as.character(h)]] <- cbind(h = h, g$weights)
    ensure(all(g$weights$status == "ok"), paste("Yearly GR fell back to equal weights at h =", h))
  }
  yearly_weights <- bind(yearly_weights)
  x <- ensemble_forecasts$test$ens_gr
  x$model <- "ens_gr_yearly"; x$probability <- prob; x$raw_score <- raw; x$clipped <- clip
  x$signal <- x$signal_050 <- as.integer(prob >= ecfg$probability_cutoff)
  x$n_train <- ntr; x$last_train_outcome <- lto
  x$config_id <- paste0("gr_intercept_clipped|", paste(members, collapse = ","), "|expanding_from_2005-01_refit_each_january|cutoff=0.50")
  ensemble_forecasts$test$ens_gr_yearly <- x
  clipping[["test yearly"]] <- data.frame(period = "test_yearly", h = NA, n = nrow(x), n_clipped = sum(clip),
    n_below_0 = sum(raw < 0), n_above_1 = sum(raw > 1))
}

# ---------------- Shared evaluator: ensembles, members and benchmarks ----------------
message("Scoring ensembles, members and benchmarks with the shared evaluator...")
returns <- load_return_data(root)
paths <- performance <- diagnostics <- list()
for (period in c("validation", "test")) for (h in cfg$horizons) {
  end <- if (period == "validation") cfg$validation_end else cfg$test_end
  b <- V[[period]]$base; i <- which(b$h == h); i <- i[order(b$origin_id[i])]
  f0 <- data.frame(origin = month_string(b$origin_id[i]))
  inp <- portfolio_inputs(f0, returns, h, end, cfg$common_rf_security)
  ef <- lapply(ensemble_forecasts[[period]], function(x) { x <- x[x$h == h, ]; x[order(month_number(x$origin)), ] })
  strategies <- list()
  # Ensembles: vote at both member-cutoff rules; probability ensembles at the fixed 0.50 only (rule 57).
  # The 0.50 vote is scored on its own vote share (members' 0.50 signals), not the tuned share.
  vote050 <- ef$ens_majority_vote; vote050$probability <- vote050$vote_count_050 / M
  strategies$ens_majority_vote <- list(validation_tuned = list(f = ef$ens_majority_vote, s = ef$ens_majority_vote$signal, cutoff = vote_cut, prob = FALSE),
                                       fixed_050 = list(f = vote050, s = vote050$signal_050, cutoff = vote_cut, prob = FALSE))
  strategies$ens_equal <- list(fixed_050 = list(f = ef$ens_equal, s = ef$ens_equal$signal, cutoff = .5, prob = TRUE))
  strategies$ens_gr <- list(fixed_050 = list(f = ef$ens_gr, s = ef$ens_gr$signal, cutoff = .5, prob = TRUE))
  if (!is.null(ef$ens_gr_yearly))
    strategies$ens_gr_yearly <- list(fixed_050 = list(f = ef$ens_gr_yearly, s = ef$ens_gr_yearly$signal, cutoff = .5, prob = TRUE))
  # Members, re-scored through the same evaluator for one comparable table.
  for (m in members) {
    mf <- data.frame(origin = month_string(b$origin_id[i]), probability = V[[period]]$P[i, m], target = b$target[i],
      classification_eligible = b$classification_eligible[i])
    cm <- V[[period]]$cutoffs$cutoff[V[[period]]$cutoffs$member == m & V[[period]]$cutoffs$h == h]
    strategies[[m]] <- list(validation_tuned = list(f = mf, s = V[[period]]$S[i, m], cutoff = cm, prob = TRUE),
                            fixed_050 = list(f = mf, s = V[[period]]$S050[i, m], cutoff = .5, prob = TRUE))
  }
  strategies$all_stocks <- list(not_applicable = list(s = rep(1L, length(i)), cutoff = NA_real_))
  strategies$all_bills <- list(not_applicable = list(s = rep(0L, length(i)), cutoff = NA_real_))
  # Rule 49: prespecified 50/50 stock/bill benchmark (same construction as the rebuilt boosting runner).
  strategies$fixed_50_50 <- list(not_applicable = list(s = rep(0L, length(i)), cutoff = NA_real_, stock_weight = .5))
  for (strategy in names(strategies)) for (rule in names(strategies[[strategy]])) {
    st <- strategies[[strategy]][[rule]]
    if (!is.null(st$f)) {
      fm <- forecast_metrics(st$f, st$cutoff)
      if (!st$prob) fm[, c("Brier", "log_loss", "mean_probability", "calibration_gap", "accuracy_at_050")] <- NA_real_
      diagnostics[[paste(period, h, strategy, rule)]] <- cbind(data.frame(period = period, h = h, model = strategy, cutoff_rule = rule), fm)
    }
    for (cost in cfg$costs) {
      p <- portfolio_path(inp, as.integer(st$s), cost, cfg)
      if (!is.null(st$stock_weight)) {
        mixed <- simulate_weights(rep(st$stock_weight, length(inp$ym)), inp$market, st$stock_weight * inp$bill_reference, cost, cfg$initial_allocation_cost)
        p[, names(mixed)] <- mixed
      }
      tags <- data.frame(period = period, h = h, strategy = strategy, cutoff_rule = rule, cutoff = st$cutoff, cost_rate = cost)
      paths[[paste(period, h, strategy, rule, cost)]] <- cbind(tags, p)
      performance[[paste(period, h, strategy, rule, cost)]] <- cbind(tags, portfolio_metrics(p, cfg$risk_aversion))
      ensure(max(month_number(p$month)) <= month_number(end), "Portfolio return crosses split boundary")
    }
  }
}
paths <- bind(paths); performance <- bind(performance); diagnostics <- bind(diagnostics)

# Common cross-horizon window (rule 40), rebased on existing paths.
common <- list()
for (period in c("validation", "test")) {
  start <- month_string(month_number(if (period == "validation") cfg$validation_start else cfg$test_start) + max(cfg$horizons))
  p <- paths[paths$period == period & paths$month >= start, ]
  groups <- interaction(p$h, p$strategy, p$cutoff_rule, p$cost_rate, drop = TRUE)
  for (g in unique(groups)) { q <- p[groups == g, ]
    common[[paste(period, g)]] <- cbind(q[1, c("period", "h", "strategy", "cutoff_rule", "cutoff", "cost_rate")], portfolio_metrics(q, cfg$risk_aversion)) }
}
common <- bind(common)

# ---------------- Cross-check: members' re-scored metrics equal their committed results ----------------
crosscheck <- list()
for (m in members) {
  pf <- ecfg$members[[m]]$performance
  if (is.null(pf) || !file.exists(file.path(root, pf))) next
  ref <- utils::read.csv(file.path(root, pf), stringsAsFactors = FALSE)
  # Older member format (boosting_performance.csv): "tuned", lowercase cagr/calmar, full/common windows.
  if ("window" %in% names(ref)) ref <- ref[ref$window == "full", ]
  ref$cutoff_rule[ref$cutoff_rule == "tuned"] <- "validation_tuned"
  if ("cagr" %in% names(ref) && !"CAGR" %in% names(ref)) names(ref)[names(ref) == "cagr"] <- "CAGR"
  if ("calmar" %in% names(ref) && !"Calmar" %in% names(ref)) names(ref)[names(ref) == "calmar"] <- "Calmar"
  ref <- ref[ref$cutoff_rule %in% c("validation_tuned", "fixed_050") & (is.null(ref$strategy) | ref$strategy %in% c(m, "boosting")), ]
  mine <- performance[performance$strategy == m, ]
  k <- c("period", "h", "cutoff_rule", "cost_rate")
  x <- merge(mine, ref, by = k, suffixes = c("", ".ref"))
  for (col in c("sharpe_matched_ladder", "sharpe_common_rf", "CAGR", "max_drawdown", "Calmar", "annual_turnover"))
    if (paste0(col, ".ref") %in% names(x)) crosscheck[[paste(m, col)]] <- data.frame(member = m, metric = col, rows = nrow(x),
      max_abs_diff = max(abs(x[[col]] - x[[paste0(col, ".ref")]]), na.rm = TRUE))
}
crosscheck <- bind(crosscheck)
if (!is.null(crosscheck)) ensure(all(crosscheck$rows > 0 & crosscheck$max_abs_diff < 1e-9),
  "Members' re-scored metrics differ from their committed performance files")

ensure(identical(hashes, setNames(vapply(input_files, file_hash, character(1)), names(hashes))), "Input files changed during run")
for (period in names(ensemble_forecasts)) {
  allcols <- unique(unlist(lapply(ensemble_forecasts[[period]], names)))
  fill <- lapply(ensemble_forecasts[[period]], function(x) { for (cn in setdiff(allcols, names(x))) x[[cn]] <- NA; x[, allcols] })
  write_csv(bind(fill), file.path(out, paste0(period, "_forecasts.csv")))
}
write_csv(gr_weight_table, file.path(out, "gr_weights.csv"))
write_csv(bind(clipping), file.path(out, "gr_clipping.csv"))
write_csv(yearly_weights, file.path(out, "gr_yearly_weights.csv"))
write_csv(performance, file.path(out, "performance.csv"))
write_csv(common, file.path(out, "common_window_performance.csv"))
write_csv(paths, file.path(out, "monthly_portfolios.csv"))
write_csv(diagnostics, file.path(out, "forecast_metrics.csv"))
if (!is.null(crosscheck)) write_csv(crosscheck, file.path(out, "evaluator_crosscheck.csv"))
write_csv(rbind(cbind(period = "validation", V$validation$cutoffs), cbind(period = "test", V$test$cutoffs)), file.path(out, "member_cutoffs.csv"))
constant <- bind(lapply(names(V), function(period) bind(lapply(cfg$horizons, function(h) {
  S <- V[[period]]$S[V[[period]]$base$h == h, , drop = FALSE]
  bind(lapply(members, function(m) data.frame(period = period, h = h, member = m, share_stock_signals = mean(S[, m]),
    pivotal_share = mean(as.integer((rowSums(S) - S[, m] + (1 - S[, m])) > M / 2) != as.integer(rowSums(S) > M / 2)))))
}))))
write_csv(constant, file.path(out, "vote_member_diagnostics.csv"))
write_csv(data.frame(file = names(hashes), sha256 = unname(hashes)), file.path(out, "input_hashes.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
jsonlite::write_json(list(version = ecfg$version, completed_at = format(Sys.time(), tz = "Asia/Singapore", usetz = TRUE),
  members = members, ensembles = ecfg$ensembles, vote_rule = ecfg$vote_rule, gr_fit = ecfg$gr_fit,
  probability_cutoff = ecfg$probability_cutoff, gr_weights = gr_weight_table, shared_config = cfg$version,
  input_hashes = as.list(hashes), code_hashes = as.list(code_hashes), tests_passed = length(tests),
  evaluator_crosscheck = if (is.null(crosscheck)) list() else crosscheck),
  file.path(out, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message("Ensemble results written to ", out)
