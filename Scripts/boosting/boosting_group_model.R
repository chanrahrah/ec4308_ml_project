#!/usr/bin/env Rscript
# =============================================================================
# EC4308 - Gradient boosting (xgboost) under the group model standards
# (GROUP_MODEL_STANDARDS_DRAFT.md, supervisor revision of 5 Nov 2026)
#
# Pipeline per horizon h in {1,3,6}:
#   1. Monthly expanding VALIDATION fits (Jan 2005-Dec 2012) for every candidate
#      theta; each fit uses only labels with outcome_end <= forecast origin.   [r.17,25]
#   2. theta* = lowest validation Brier on eligible labels (tie rule below).    [r.25-26]
#   3. Cutoff c* = argmax horizon-matched net Sharpe at 10bp on the saved
#      validation probabilities of theta* (grid seq(0,1,.01)); 0.50 reported too. [r.28-29]
#   4. Monthly expanding TEST fits (Jan 2013-Dec 2025) with theta*, c* frozen.   [r.30]
#   5. Export probabilities (rule-66 schema), audits, portfolio paths, metrics.  [r.66-68]
#
# Run:  Rscript boosting_group_model.R            (full)
#       QUICK=1 H=3 Rscript boosting_group_model.R (smoke test: tiny grid)
# Place in <project>/Scripts/boosting/. Reads model_data_h*.csv, features_h*_132.csv, base_raw_h*.csv
# from Data/processed/ (found automatically); writes to Results/boosting/.
# Env:  EC4308_ROOT = project root override; EC4308_DATA = data folder override;
#       EC4308_OUT = output folder override
# =============================================================================
suppressPackageStartupMessages(library(xgboost))

# ---- project paths -----------------------------------------------------------
# Expected layout:  <root>/Data/processed/*.csv   <root>/Scripts/boosting/this_file.R   <root>/Results/
# <root> = nearest parent folder containing both Data/ and Scripts/ (override: EC4308_ROOT).
script_dir <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1])))                       # Rscript
  for (i in rev(seq_len(sys.nframe()))) { of <- sys.frame(i)$ofile            # source()
    if (!is.null(of)) return(dirname(normalizePath(of))) }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {   # RStudio / Positron editor
    pth <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(pth)) return(dirname(normalizePath(pth))) }
  getwd()
}
find_root <- function() {
  ov <- Sys.getenv("EC4308_ROOT"); if (nzchar(ov)) return(normalizePath(ov))
  for (start in unique(c(script_dir(), getwd()))) {
    d <- normalizePath(start)
    repeat { if (dir.exists(file.path(d, "Data")) && dir.exists(file.path(d, "Scripts"))) return(d)
      up <- dirname(d); if (up == d) break; d <- up }
  }
  stop("Could not find the project root (a folder containing Data/ and Scripts/). Set EC4308_ROOT or run from inside the project.")
}
PROJECT_ROOT <- find_root()
DATA_DIR <- Sys.getenv("EC4308_DATA", file.path(PROJECT_ROOT, "Data"))
find_file <- function(name) {
  first <- file.path(DATA_DIR, "processed", name)
  if (file.exists(first)) return(first)
  hits <- list.files(DATA_DIR, pattern = paste0("^", gsub("\\.", "\\\\.", name), "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits)) return(hits[1])
  stop(sprintf("Cannot find '%s' under %s. Put it in the processed/ subfolder or another subfolder of the data directory.", name, DATA_DIR))
}

QUICK <- nzchar(Sys.getenv("QUICK"))
CFG <- list(
  model        = "boosting_xgb",
  out_dir      = Sys.getenv("EC4308_OUT", file.path(PROJECT_ROOT, "Results", "boosting")),
  horizons     = if (nzchar(Sys.getenv("H"))) as.integer(strsplit(Sys.getenv("H"), ",")[[1]]) else c(1L, 3L, 6L),
  val_start    = "2005-01", val_end  = "2012-12",
  test_start   = "2013-01", test_end = "2025-12",
  common_start = "2013-07",                      # r.40 common cross-horizon window
  # ---- hyperparameter grid, fixed BEFORE testing (r.24) -------------------
  depth_grid   = if (QUICK) c(1L, 2L) else c(1L, 2L, 3L),
  eta_grid     = if (QUICK) 0.10      else c(0.02, 0.05, 0.10),
  nrounds_grid = if (QUICK) c(50L, 100L) else c(50L, 100L, 200L, 400L),
  subsample = 0.7, colsample = 0.5, min_child_weight = 5, lambda = 1,
  seed = 4308L,
  # ---- portfolio / evaluation (r.28, 32, 42) ---------------------------------
  cost_select  = 0.001,                           # 10 bp, used for cutoff selection
  cost_report  = c(0, 0.001, 0.002),              # 0 / 10 / 20 bp sensitivity, same signals
  cutoff_step  = as.numeric(Sys.getenv("CUTOFF_STEP", "0.01")),
  cutoff_grid  = NULL,                                            # set just below
  tol          = 1e-12
)
CFG$cutoff_grid <- seq(0, 1, by = CFG$cutoff_step)
NON_FEATURES <- c("ym", "outcome_end", "equity_forward", "bill_accrual_proxy",
                  "target", "split", "mature_target_lag")      # r.3, r.13
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

ym_idx <- function(s) as.integer(substr(s, 1, 4)) * 12L + as.integer(substr(s, 6, 7)) - 1L
idx_ym <- function(i) sprintf("%04d-%02d", i %/% 12L, i %% 12L + 1L)
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n", sep = " ")

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
load_h <- function(h) {
  md <- read.csv(find_file(sprintf("model_data_h%d.csv", h)), stringsAsFactors = FALSE)
  feat <- setdiff(names(md), NON_FEATURES)
  stopifnot(length(feat) == 132L)
  ft <- read.csv(find_file(sprintf("features_h%d_132.csv", h)), stringsAsFactors = FALSE)
  stopifnot(identical(names(ft)[-1], feat), isTRUE(all.equal(as.matrix(ft[, feat]), as.matrix(md[, feat]), check.attributes = FALSE)))
  list(h = h, md = md, feat = feat, X = as.matrix(md[, feat]),
       y = as.numeric(md$target),
       o = ym_idx(md$ym),
       oe = ifelse(is.na(md$outcome_end), NA_integer_, ym_idx(md$outcome_end)))
}

# Engine inputs: monthly equity return, actual quote-day counts, entry bill yields.
# bill_yld_matched = expm1(365 * TMYLD_daily)  =>  days(t,t+1) is recovered exactly
# from the h=1 file as log1p(bill_accrual_proxy) * 365 / log1p(bill_yld_matched)
# (verified to be integer to ~1e-10 for every month).
build_engine <- function() {
  h1  <- read.csv(find_file("model_data_h1.csv"), stringsAsFactors = FALSE)
  bl1 <- read.csv(find_file("base_raw_h1.csv"), stringsAsFactors = FALSE)
  m1  <- merge(h1[, c("ym", "equity_forward", "bill_accrual_proxy")], bl1[, c("ym", "bill_yld_matched")], by = "ym", sort = TRUE)
  o   <- ym_idx(m1$ym)
  days <- round(log1p(m1$bill_accrual_proxy) * 365 / log1p(m1$bill_yld_matched))
  # Months where the 4-week yield is exactly 0 (2010-01, 2011-12, 2021-04, 2021-05) make the ratio 0/0.
  # Recover them from the 13-week file using days(t,t+3) = d(t) + d(t+1) + d(t+2).
  h3  <- read.csv(find_file("model_data_h3.csv"), stringsAsFactors = FALSE)
  bl3 <- read.csv(find_file("base_raw_h3.csv"), stringsAsFactors = FALSE)
  m3  <- merge(data.frame(ym = m1$ym), merge(h3[, c("ym", "bill_accrual_proxy")], bl3[, c("ym", "bill_yld_matched")], by = "ym"), by = "ym", sort = TRUE)
  D3  <- round(log1p(m3$bill_accrual_proxy) * 365 / log1p(m3$bill_yld_matched)); n <- length(days)
  repeat {
    na <- which(is.na(days)); changed <- FALSE
    for (i in na) {
      if (i + 2 <= n && !is.na(D3[i]) && !anyNA(days[c(i + 1, i + 2)])) { days[i] <- D3[i] - days[i + 1] - days[i + 2]; changed <- TRUE }
      else if (i - 2 >= 1 && !is.na(D3[i - 2]) && !anyNA(days[c(i - 2, i - 1)])) { days[i] <- D3[i - 2] - days[i - 2] - days[i - 1]; changed <- TRUE }
    }
    if (!changed) break
  }
  stopifnot(sum(is.na(days)) <= 1L, all(days[!is.na(days)] >= 27 & days[!is.na(days)] <= 34))   # only the final origin may stay NA
  EN <- list(o0 = o[1], req = m1$equity_forward, days = days, d = list())
  for (h in c(1L, 3L, 6L)) {
    bl <- read.csv(find_file(sprintf("base_raw_h%d.csv", h)), stringsAsFactors = FALSE)
    mm <- merge(data.frame(ym = m1$ym), bl[, c("ym", "bill_yld_matched")], by = "ym", all.x = TRUE, sort = TRUE)
    EN$d[[as.character(h)]] <- log1p(mm$bill_yld_matched) / 365        # daily cc yield at entry origin
  }
  EN
}

# -----------------------------------------------------------------------------
# 2. Model: one expanding-window xgboost fit -> probabilities for several nrounds
# -----------------------------------------------------------------------------
fit_rows_at <- function(D, t) which(!is.na(D$y) & !is.na(D$oe) & D$oe <= t)   # r.17: outcome_end <= origin

fit_predict <- function(D, fit_idx, pred_row, depth, eta, seed, nr_grid) {
  stopifnot(length(fit_idx) > 50)
  dtrain <- xgb.DMatrix(D$X[fit_idx, , drop = FALSE], label = D$y[fit_idx], missing = NA)
  params <- list(objective = "binary:logistic", max_depth = depth, eta = eta,
                 subsample = CFG$subsample, colsample_bytree = CFG$colsample,
                 min_child_weight = CFG$min_child_weight, lambda = CFG$lambda,
                 nthread = 1L, seed = seed)
  bst <- xgb.train(params, dtrain, nrounds = max(nr_grid), verbose = 0)
  xp  <- D$X[pred_row, , drop = FALSE]
  vapply(nr_grid, function(n) as.numeric(predict(bst, xp, iterationrange = c(1L, as.integer(n)))), numeric(1))
}

cfg_table <- function() {
  g <- expand.grid(max_depth = CFG$depth_grid, eta = CFG$eta_grid, nrounds = CFG$nrounds_grid)
  g$config_id <- sprintf("d%d_eta%.2f_n%d", g$max_depth, g$eta, g$nrounds)
  g
}

# Expanding monthly forecasts for all candidate configs over origins `rows`.
run_expanding <- function(D, rows, cfgs, label) {
  fits <- unique(cfgs[, c("max_depth", "eta")])
  P <- matrix(NA_real_, length(rows), nrow(cfgs), dimnames = list(NULL, cfgs$config_id))
  audit <- data.frame(row = rows, origin = D$o[rows], n_train = NA_integer_,
                      last_train_origin = NA_integer_, last_train_outcome = NA_integer_)
  for (i in seq_along(rows)) {
    r <- rows[i]; t <- D$o[r]
    fi <- fit_rows_at(D, t)
    stopifnot(max(D$oe[fi]) <= t)                                   # r.70 audit
    audit$n_train[i] <- length(fi); audit$last_train_origin[i] <- max(D$o[fi]); audit$last_train_outcome[i] <- max(D$oe[fi])
    for (k in seq_len(nrow(fits))) {
      sel <- which(cfgs$max_depth == fits$max_depth[k] & cfgs$eta == fits$eta[k])
      P[i, sel] <- fit_predict(D, fi, r, fits$max_depth[k], fits$eta[k],
                               seed = CFG$seed + D$h * 100003L + t, nr_grid = cfgs$nrounds[sel])
    }
    if (i %% 24 == 0) log_msg(sprintf("  h=%d %s: %d/%d origins", D$h, label, i, length(rows)))
  }
  list(P = P, audit = audit)
}

# -----------------------------------------------------------------------------
# 3. Shared portfolio evaluator (rules 34-39, 43-44)
# -----------------------------------------------------------------------------
# sig: binary signals aligned with consecutive origins o_vec. First h-1 signals are
# warm-up history; the first complete h-cohort portfolio forms at o_vec[h] with wealth 1
# and a free initial allocation; return months run o_vec[h]+1 ... o_vec[n].
run_portfolio <- function(sig, o_vec, h, cost, EN, w_fixed = NULL) {
  n <- length(o_vec); stopifnot(all(diff(o_vec) == 1L), length(sig) == n, n > h)
  js <- h:(n - 1L); m <- length(js)
  out <- data.frame(origin = o_vec[js], stock_weight = NA_real_, net_return = NA_real_, gross_return = NA_real_,
                    common_rf_return = NA_real_, matched_ladder_return = NA_real_, equity_return = NA_real_,
                    signed_stock_trade = NA_real_, cost = NA_real_, pretrade_wealth = NA_real_, wealth = NA_real_,
                    quote_days = NA_real_)
  dh <- EN$d[[as.character(h)]]; d1 <- EN$d[["1"]]
  w <- if (is.null(w_fixed)) mean(sig[js[1]:(js[1] - h + 1L)]) else w_fixed
  V <- 1; S <- w * V                                                 # free initial allocation (r.39)
  for (a in seq_len(m)) {
    j <- js[a]; t <- o_vec[j]; pos <- (j - h + 1L):j; ti <- t - EN$o0 + 1L
    w <- if (is.null(w_fixed)) mean(sig[pos]) else w_fixed
    need <- w * V - S
    q <- if (need >= 0) need / (1 + cost * w) else need / (1 - cost * w)   # r.38 exact post-cost trade
    cst <- cost * abs(q); Vp <- V - cst; Sp <- S + q
    stopifnot(abs(Sp - w * Vp) < 1e-10)                                  # exact post-cost weight check
    days <- EN$days[ti]
    b <- expm1(dh[(o_vec[pos] - EN$o0 + 1L)] * days)                      # each cohort keeps its entry yield (r.36)
    matched <- mean(b)
    common  <- expm1(d1[ti] * days)                                       # r.43 common 4-week reference
    req <- EN$req[ti]
    bill_end <- if (is.null(w_fixed)) (Vp / h) * sum((1 + b)[sig[pos] == 0]) else (Vp - Sp) * (1 + matched)
    end <- Sp * (1 + req) + bill_end
    out[a, ] <- list(t, w, end / V - 1, end / Vp - 1, common, matched, req, q, cst, V, end, days)
    V <- end; S <- Sp * (1 + req)
  }
  out$month <- idx_ym(out$origin + 1L); out$quote_start <- idx_ym(out$origin); out$quote_end <- out$month
  out
}

sharpe_ann <- function(x, ref) {
  e <- x - ref; s <- sd(e)
  if (!is.finite(s) || s <= 1e-12) return(NA_real_)                     # undefined, not 0 (r.46)
  sqrt(12) * mean(e) / s
}

perf_metrics <- function(df) {
  r <- df$net_return; n <- length(r)
  W <- c(1, cumprod(1 + r))
  dd <- W / cummax(W) - 1
  exm <- r - df$matched_ladder_return
  cagr <- prod(1 + r)^(12 / n) - 1; mdd <- -min(dd)
  mkt_ex <- df$equity_return - df$matched_ladder_return
  c(n_months = n,
    sharpe_common_rf = sharpe_ann(r, df$common_rf_return),
    sharpe_matched_ladder = sharpe_ann(r, df$matched_ladder_return),
    cagr = cagr, ann_vol = sqrt(12) * sd(r), max_drawdown = mdd,
    sortino_matched = sqrt(12) * mean(exm) / sqrt(mean(pmin(exm, 0)^2)),
    calmar = if (mdd > 0) cagr / mdd else NA_real_,
    ulcer = sqrt(mean(dd[-1]^2)),
    beta_matched = if (sd(mkt_ex) > 0) cov(exm, mkt_ex) / var(mkt_ex) else NA_real_,
    exposure = mean(df$stock_weight),
    annual_turnover = 12 * mean(abs(df$signed_stock_trade) / df$pretrade_wealth),
    trading_months = sum(abs(df$signed_stock_trade) > 1e-10),
    total_cost = sum(df$cost), terminal_wealth = prod(1 + r))
}

select_cutoff <- function(P_val, o_vec, h, EN) {
  g <- CFG$cutoff_grid
  res <- t(vapply(g, function(cc) {
    pf <- run_portfolio(as.integer(P_val >= cc), o_vec, h, CFG$cost_select, EN)
    c(sharpe_matched_ladder = sharpe_ann(pf$net_return, pf$matched_ladder_return),
      sharpe_common_rf = sharpe_ann(pf$net_return, pf$common_rf_return),
      exposure = mean(pf$stock_weight), terminal_wealth = prod(1 + pf$net_return))
  }, numeric(4)))
  grid <- data.frame(h = h, cutoff = g, res)
  ok <- is.finite(grid$sharpe_matched_ladder)
  if (!any(ok)) stop(sprintf("h=%d: no cutoff has a defined matched Sharpe - stop and flag (r.29)", h))
  best <- max(grid$sharpe_matched_ladder[ok])
  cand <- which(ok & grid$sharpe_matched_ladder >= best - CFG$tol)
  cand <- cand[order(abs(grid$cutoff[cand] - 0.5), grid$cutoff[cand])]       # nearest 0.50, then lower
  list(cutoff = grid$cutoff[cand[1]], grid = grid, n_tied = length(cand))
}

# -----------------------------------------------------------------------------
# 4. Forecast metrics (eligible labels only, r.22/48)
# -----------------------------------------------------------------------------
auc_rank <- function(p, y) { n1 <- sum(y == 1); n0 <- sum(y == 0); if (!n1 || !n0) return(NA_real_)
  (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
fc_metrics <- function(p, y, cutoff, prev_fc) {
  pc <- pmin(pmax(p, 1e-15), 1 - 1e-15)                                   # clip only for log loss
  c(n_eligible = length(y), event_prevalence = mean(y), brier = mean((p - y)^2),
    brier_prevalence_benchmark = mean((prev_fc - y)^2),
    log_loss = -mean(y * log(pc) + (1 - y) * log(1 - pc)), auc = auc_rank(p, y),
    accuracy_tuned = mean((p >= cutoff) == (y == 1)), accuracy_050 = mean((p >= 0.5) == (y == 1)))
}
calib_table <- function(p, y, h, period) {
  b <- cut(p, c(0, .2, .4, .6, .8, 1), include.lowest = TRUE)
  do.call(rbind, lapply(levels(b), function(l) { i <- which(b == l)
    data.frame(h = h, period = period, bin = l, n = length(i), mean_p = if (length(i)) mean(p[i]) else NA, mean_y = if (length(i)) mean(y[i]) else NA) }))
}

# -----------------------------------------------------------------------------
# 5. Main
# -----------------------------------------------------------------------------
main <- function() {
  EN <- build_engine(); cfgs <- cfg_table()
  log_msg("candidate configs:", nrow(cfgs), "| horizons:", paste(CFG$horizons, collapse = ","))
  S <- list(valprob = list(), testprob = list(), cand = list(), cutgrid = list(), sel = list(), cand_prob = list(),
            monthly = list(), perf = list(), fcm = list(), calib = list(), audit = list())
  for (h in CFG$horizons) {
    D <- load_h(h)
    vs <- ym_idx(CFG$val_start); ve <- ym_idx(CFG$val_end); ts <- ym_idx(CFG$test_start); te <- ym_idx(CFG$test_end)
    vrows <- which(D$o >= vs & D$o <= ve); trows <- which(D$o >= ts & D$o <= te)
    stopifnot(length(vrows) == 96L, length(trows) == 156L)
    v_elig <- !is.na(D$y[vrows]) & D$oe[vrows] <= ve; t_elig <- !is.na(D$y[trows]) & D$oe[trows] <= te
    stopifnot(sum(v_elig) == c(`1` = 95, `3` = 93, `6` = 90)[[as.character(h)]],       # r.15 table
              sum(t_elig) == c(`1` = 155, `3` = 153, `6` = 150)[[as.character(h)]])
    # sanity: target definition (r.7): 1{equity_forward > bill_accrual_proxy}, tie = 0
    ok <- !is.na(D$y); stopifnot(all(D$y[ok] == as.numeric(D$md$equity_forward[ok] > D$md$bill_accrual_proxy[ok])))

    # ---- validation: monthly expanding fits for all candidates ----
    log_msg(sprintf("h=%d validation fits", h))
    V <- run_expanding(D, vrows, cfgs, "val")
    brier <- apply(V$P[v_elig, , drop = FALSE], 2, function(p) mean((p - D$y[vrows][v_elig])^2))
    cand <- cbind(h = h, cfgs, brier = brier, n_eval = sum(v_elig))
    best <- min(brier); tied <- which(brier <= best + CFG$tol)
    tied <- tied[order(cfgs$max_depth[tied], cfgs$nrounds[tied], cfgs$eta[tied])]          # simpler first
    star <- tied[1]; cid <- cfgs$config_id[star]
    log_msg(sprintf("h=%d theta*=%s | val Brier=%.5f", h, cid, brier[star]))

    # ---- cutoff on saved validation probabilities of theta* ----
    p_val <- V$P[, star]
    cs <- select_cutoff(p_val, D$o[vrows], h, EN); cstar <- cs$cutoff
    log_msg(sprintf("h=%d c*=%.2f (ties=%d)", h, cstar, cs$n_tied))

    # ---- test: monthly expanding refits with frozen theta* ----
    log_msg(sprintf("h=%d test fits", h))
    one <- cfgs[star, , drop = FALSE]
    Tt <- run_expanding(D, trows, one, "test"); p_test <- Tt$P[, 1]

    mk <- function(rows, p, aud, elig, cid, cutoff) data.frame(
      model = CFG$model, h = h, origin = D$md$ym[rows], probability = p,
      cutoff = cutoff, signal = as.integer(p >= cutoff), cutoff_050 = 0.5, signal_050 = as.integer(p >= 0.5),
      target = D$y[rows], outcome_end = D$md$outcome_end[rows], classification_eligible = elig,
      fit_information_cutoff = idx_ym(D$o[rows]), last_train_origin = idx_ym(aud$last_train_origin),
      last_train_outcome = idx_ym(aud$last_train_outcome), n_train = aud$n_train, config_id = cid)
    vdf <- mk(vrows, p_val, V$audit, v_elig, cid, cstar); tdf <- mk(trows, p_test, Tt$audit, t_elig, cid, cstar)

    # ---- portfolios ----
    run_period <- function(period, df, rows, cutoffs) {
      o_vec <- D$o[rows]; res <- list()
      add <- function(strategy, rule, cutoff, sig, cost, w_fixed = NULL) {
        pf <- run_portfolio(sig, o_vec, h, cost, EN, w_fixed)
        pf$period <- period; pf$h <- h; pf$strategy <- strategy; pf$cutoff_rule <- rule; pf$cutoff <- cutoff; pf$cost_rate <- cost
        res[[length(res) + 1L]] <<- pf
      }
      for (cst in CFG$cost_report) {
        add("boosting", "tuned", cutoffs, df$signal, cst)
        add("boosting", "fixed_050", 0.5, df$signal_050, cst)
        add("all_stock", "none", NA, rep(1L, length(o_vec)), cst)
        add("all_bills", "none", NA, rep(0L, length(o_vec)), cst)
        add("fixed_50_50", "none", NA, rep(0L, length(o_vec)), cst, w_fixed = 0.5)       # r.49 diagnostic
      }
      do.call(rbind, res)
    }
    pv <- run_period("validation", vdf, vrows, cstar); pt <- run_period("test", tdf, trows, cstar)
    pm <- rbind(pv, pt)
    keep <- c("period", "h", "strategy", "cutoff_rule", "cutoff", "cost_rate", "month", "quote_start", "quote_end", "quote_days",
              "stock_weight", "net_return", "gross_return", "common_rf_return", "matched_ladder_return", "equity_return",
              "signed_stock_trade", "cost", "pretrade_wealth", "wealth")
    pm <- pm[, keep]

    # ---- r.70 portfolio checks ----
    for (cst in CFG$cost_report) {
      a <- pm[pm$period == "test" & pm$strategy == "all_stock" & pm$cost_rate == cst, ]
      stopifnot(max(abs(a$net_return - a$equity_return)) < 1e-12)                # all-stock == market
      b <- pm[pm$period == "test" & pm$strategy == "all_bills" & pm$cost_rate == cst, ]
      stopifnot(max(abs(b$net_return - b$matched_ladder_return)) < 1e-12)        # all-bills == ladder
    }

    # ---- performance tables ----
    gsplit <- split(pm, list(pm$period, pm$strategy, pm$cutoff_rule, pm$cost_rate), drop = TRUE)
    perf <- do.call(rbind, lapply(gsplit, function(g) data.frame(
      period = g$period[1], h = h, strategy = g$strategy[1], cutoff_rule = g$cutoff_rule[1], cutoff = g$cutoff[1],
      cost_rate = g$cost_rate[1], window = "full", t(perf_metrics(g)), net_Sharpe = perf_metrics(g)[["sharpe_matched_ladder"]])))
    perf_c <- do.call(rbind, lapply(gsplit, function(g) { g <- g[g$period == "test" & g$month >= CFG$common_start, ]
      if (!nrow(g)) return(NULL)
      data.frame(period = "test", h = h, strategy = g$strategy[1], cutoff_rule = g$cutoff_rule[1], cutoff = g$cutoff[1],
                 cost_rate = g$cost_rate[1], window = "common_2013-07_2025-12", t(perf_metrics(g)),
                 net_Sharpe = perf_metrics(g)[["sharpe_matched_ladder"]]) }))

    # ---- forecast metrics (same eligible rows) ----
    prev_fc <- function(rows) vapply(rows, function(r) mean(D$y[fit_rows_at(D, D$o[r])]), numeric(1))
    fcm <- rbind(
      data.frame(period = "validation", h = h, t(fc_metrics(p_val[v_elig], D$y[vrows][v_elig], cstar, prev_fc(vrows[v_elig])))),
      data.frame(period = "test", h = h, t(fc_metrics(p_test[t_elig], D$y[trows][t_elig], cstar, prev_fc(trows[t_elig])))))
    cal <- rbind(calib_table(p_val[v_elig], D$y[vrows][v_elig], h, "validation"), calib_table(p_test[t_elig], D$y[trows][t_elig], h, "test"))

    S$valprob[[as.character(h)]] <- vdf; S$testprob[[as.character(h)]] <- tdf; S$cand[[as.character(h)]] <- cand
    S$cutgrid[[as.character(h)]] <- cs$grid; S$monthly[[as.character(h)]] <- pm
    S$perf[[as.character(h)]] <- rbind(perf, perf_c); S$fcm[[as.character(h)]] <- fcm; S$calib[[as.character(h)]] <- cal
    S$sel[[as.character(h)]] <- data.frame(h = h, config_id = cid, max_depth = cfgs$max_depth[star], eta = cfgs$eta[star],
                                           nrounds = cfgs$nrounds[star], val_brier = brier[star], cutoff = cstar, n_cutoff_ties = cs$n_tied)
    cp <- data.frame(h = h, origin = D$md$ym[vrows], V$P); S$cand_prob[[as.character(h)]] <- cp
    fmt_aud <- function(a) { a$last_train_origin <- idx_ym(a$last_train_origin); a$last_train_outcome <- idx_ym(a$last_train_outcome); a[, -(1:2)] }
    S$audit[[as.character(h)]] <- rbind(data.frame(period = "validation", h = h, origin = D$md$ym[vrows], fmt_aud(V$audit)),
                                        data.frame(period = "test", h = h, origin = D$md$ym[trows], fmt_aud(Tt$audit)))
    log_msg(sprintf("h=%d done", h))
  }
  # ---- write ----
  w <- function(x, f) write.csv(do.call(rbind, x), file.path(CFG$out_dir, f), row.names = FALSE)
  w(S$valprob, "boosting_validation_probabilities.csv"); w(S$testprob, "boosting_test_probabilities.csv")
  w(S$cand, "boosting_candidate_brier.csv"); w(S$cutgrid, "boosting_cutoff_grid.csv"); w(S$sel, "boosting_selected_settings.csv")
  w(S$cand_prob, "boosting_validation_candidate_probabilities_wide.csv"); w(S$monthly, "boosting_monthly_returns.csv")
  w(S$perf, "boosting_performance.csv"); w(S$fcm, "boosting_forecast_metrics.csv"); w(S$calib, "boosting_calibration.csv")
  w(S$audit, "boosting_refit_audit.csv")
  files <- c(sprintf("model_data_h%d.csv", c(1, 3, 6)), sprintf("features_h%d_132.csv", c(1, 3, 6)), sprintf("base_raw_h%d.csv", c(1, 3, 6)))
  hash <- tools::md5sum(vapply(files, find_file, character(1)))
  D1 <- load_h(CFG$horizons[1])
  info <- c("spec: GROUP_MODEL_STANDARDS_DRAFT.md (supervisor revision 2026-11-05)",
            paste("R:", R.version.string), paste("xgboost:", as.character(packageVersion("xgboost"))),
            paste("seed base:", CFG$seed, "(per-fit seed = base + h*100003 + origin index)"),
            paste("grid depth:", toString(CFG$depth_grid), "| eta:", toString(CFG$eta_grid), "| nrounds:", toString(CFG$nrounds_grid)),
            paste("fixed: subsample", CFG$subsample, "colsample_bytree", CFG$colsample, "min_child_weight", CFG$min_child_weight, "lambda", CFG$lambda, "nthread 1"),
            "tie rule (Brier within 1e-12): lowest max_depth, then fewest rounds, then lowest eta",
            paste("cutoff grid: seq(0,1,by=", CFG$cutoff_step, ") | tie rule: nearest 0.50, then lower; undefined Sharpe excluded; stop if none defined"),
            "", "source hashes (md5):", paste(names(hash), hash), "", "feature names (132):", D1$feat)
  writeLines(info, file.path(CFG$out_dir, "boosting_run_info.txt"))
  log_msg("written to", normalizePath(CFG$out_dir))
  invisible(S)
}

# -----------------------------------------------------------------------------
# 6. r.70 look-ahead perturbation check (runs after main when CHECK=1)
# -----------------------------------------------------------------------------
future_perturbation_check <- function(h = 3L, origins = c("2008-06", "2015-03", "2021-09")) {
  D <- load_h(h); set.seed(1)
  for (om in origins) {
    r <- which(D$md$ym == om); t <- D$o[r]; fi <- fit_rows_at(D, t)
    p0 <- fit_predict(D, fi, r, 2L, 0.1, 99L, c(50L, 100L))
    D2 <- D; fut <- which(D$o > t)                                   # scramble everything after the origin
    D2$X[fut, ] <- D$X[sample(fut), ]; D2$y[fut] <- sample(D$y[fut])
    p1 <- fit_predict(D2, fit_rows_at(D2, t), r, 2L, 0.1, 99L, c(50L, 100L))
    stopifnot(isTRUE(all.equal(p0, p1, tolerance = 0)))
    log_msg("perturbation check OK at", om)
  }
}

if (!isTRUE(getOption("ec4308.no_autorun"))) {
  main()
  if (nzchar(Sys.getenv("CHECK"))) future_perturbation_check()   # header promises this; nothing called it before
}
