# LOGIT Elastic Net; teacher mixing parameter rho is named alpha in the CSV outputs.
# Teacher convention: mean binomial negative log-likelihood + lambda * sum(alpha * beta^2 + (1-alpha) * abs(beta)).
# alpha=0: Lasso; alpha=1: Ridge. Intercept unpenalized.
# All reported alpha/lambda values use this convention; glmnet values are internal conversions.
# EC4308: complete Elastic Net pipeline for h=1,3,6; cutoff step = 0.05.
# Based on EC4308_Elastic_Net_Fixed_Lambda_Max_0.2.zip (2026-11-08).
# Model/portfolio algorithms retained: 132 penalized X, unpenalized intercept.
# This latest baseline does NOT append mature_target_lag as an extra regressor.
# 40 candidates per horizon: 5 alphas x 8 shared fixed absolute lambdas.
# Selection: validation Brier, then validation net matched-ladder Sharpe.
# Both matched-ladder and common 1-month RF Sharpe are reported.
# New outputs: <project>/Results/elastic_net_logit_teacher_fixed03_cutoff005_full_pipeline/
# Prior test results have been viewed: this is an exploratory sensitivity run.
#
# Missing packages are installed automatically. Manual alternative:
# install.packages(c("glmnet", "digest", "jsonlite"))
#
# RUN: open this file and click Source ONCE.
# Runs checks -> validation -> frozen selection -> test -> saves all results.
# No Console options or second Source needed. Old stage options are ignored.
# Set PROJECT_ROOT below to the folder containing Data/ and Results/.
# For replication, provide this script together with the five data files below.
PROJECT_ROOT <- "~/Downloads/EC4308_Project"
#
# REQUIRED FILES relative to the project root (no repeated file dialogs):
# Data/processed/model_ready_h1.rds
# Data/processed/model_ready_h3.rds
# Data/processed/model_ready_h6.rds
# Data/crsp_monthly_index_data/crsp_monthly_index_data.csv
# Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv
# Raw return columns: mthcaldt,vwretd; bill columns: MCALDT,KYTREASNOX,TMYLD.
# No monthly_returns.csv / bill_vintages.csv is required in this latest pipeline.
#
# Terminal alternative (also runs both phases by default):
# Rscript EC4308_elastic_net_logit_fixed03_cutoff005_full_pipeline.R --root=/path/to/project
# Optional advanced CLI: --stage=validation, --stage=test, or --stage=tests.

# Find this exact script for frozen-code validation (RStudio Source or Rscript).
.ec4308_source_path <- local({
  ofiles <- Filter(function(z) is.character(z) && length(z)==1L && nzchar(z),
                   lapply(sys.frames(), function(f) f$ofile))
  cli <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  p <- if(length(ofiles)) tail(ofiles,1)[[1]] else if(length(cli)) sub("^--file=","",cli[1]) else ""
  if(!nzchar(p)) stop("Open this file and use Source, or run it with Rscript.")
  normalizePath(p,mustWork=TRUE)
})





# EC4308 Elastic Net: shared parameters for h=1,3,6.
elastic_config <- function() list(
  version = "elastic_net_logit_teacher_fixed03_cutoff005_v1_2026_11_17",
  horizons = c(1L,3L,6L), model = "elastic_net_logit_teacher_fixed03_cutoff005",
  train_end = "2004-12", validation_start = "2005-01", validation_end = "2012-12",
  test_start = "2013-01", test_end = "2025-12",
  penalty_convention = "teacher: lambda * (alpha * L2_squared + (1-alpha) * L1)",
  # Fixed ABSOLUTE lambda, identical at every recursive fitting origin.
  # Non-uniform shared grid; 0.3 is a practical upper bound, not a theoretical optimum.
  # Ex-post methodological sensitivity: do not describe test as untouched holdout.
  alphas = c(0,0.25,0.5,0.75,1),
  lambdas = c(0.0001,0.0005,0.002,0.01,0.05,0.1,0.2,0.3),
  # 21 candidate cutoffs on [0,1], in steps of 0.05.
  cutoffs = seq(0,1,by=0.05), comparison_cutoff = 0.5,
  cutoff_objective = "sharpe_matched_ladder", primary_cost = 0.001,
  costs = c(0,0.001,0.002), risk_aversion = 3,
  common_rf_security = 2000061L,
  initial_allocation_cost = FALSE,
  standardize = "fitting_rows_only", selection_score = "Brier",
  grid_tie_rule = "smaller number of active terms then larger absolute lambda then smaller alpha",
  seed = 4308L,
  bill_rule = "Origin-vintage matched daily log yield; actual monthly quote days",
  cost_rule = "Exact post-cost net equity trade at 10bp; 0/20bp sensitivities"
)
month_number <- function(x) {
  x <- as.character(x)
  as.integer(substr(x,1,4))*12L + as.integer(substr(x,6,7))-1L
}
month_string <- function(x) sprintf("%04d-%02d",x %/% 12L,x %% 12L+1L)
ensure <- function(ok,message) if(!isTRUE(ok)) stop(message,call.=FALSE)
write_csv <- function(x,path) utils::write.csv(x,path,row.names=FALSE,na="")
file_hash <- function(path) digest::digest(file=path,algo="sha256",serialize=FALSE)
bind_rows <- function(xs) if(length(xs)) do.call(rbind,xs) else data.frame()


# Predictors/targets are frozen upstream. No reconstruction or future metadata input.
load_elastic_data <- function(root,h,cfg) {
  path <- file.path(root,"Data","processed",paste0("model_ready_h",h,".rds"))
  d <- readRDS(path)
  ensure(is.matrix(d$X)||is.data.frame(d$X),"X must be matrix/data.frame")
  d$X <- as.matrix(d$X)
  ensure(ncol(d$X)==132L,"Expected full ML $X (132 columns: 22 x 6 lags)")
  ensure(!anyDuplicated(colnames(d$X)) && !is.null(colnames(d$X)),"Missing/duplicate feature names")
  ensure(is.numeric(d$X) && all(is.finite(d$X)),"Nonfinite ML predictors: inspect validated inputs")
  ensure(length(d$ym)==nrow(d$X) && length(d$y)==nrow(d$X),"Data lengths do not align")
  ensure(identical(as.character(d$ym),as.character(d$metadata$ym)),"Metadata origin mismatch")
  d$m <- month_number(d$ym)
  d$end <- month_number(d$metadata$outcome_end)
  ensure(all(diff(d$m)==1L)&&!anyDuplicated(d$m),"Calendar not consecutive")
  known <- !is.na(d$y)
  ensure(all(d$y[known] %in% c(0,1)) && all(d$end[known]==d$m[known]+h),"Target/horizon mismatch")
  ensure(all(d$end[d$rows$train]<=month_number(cfg$train_end)),"Initial train boundary leakage")
  ensure(all(d$end[d$rows$validation]<=month_number(cfg$validation_end)),"Validation boundary leakage")
  d
}
# At origin t only y_s with s+h<=t are known.
mature_rows <- function(d,origin) which(!is.na(d$y)&!is.na(d$end)&d$end<=origin&d$m<origin)

# Fixed absolute-lambda model: standardized using current mature training rows ONLY.
# Training gradient diagnostics never set the candidate lambda.
fit_elastic <- function(d,rows,alpha,lambda) {
  ensure(length(rows)>40L && length(unique(d$y[rows]))==2L,"Too few fitting rows or one class")
  ensure(is.finite(lambda) && lambda>0 && alpha>=0 && alpha<=1,"Invalid fixed lambda or alpha")
  x <- d$X[rows,,drop=FALSE]; y <- as.numeric(d$y[rows]); n <- length(rows)
  center <- colMeans(x); scales <- apply(x,2,stats::sd)
  ensure(all(is.finite(scales)),"Invalid predictor scaling")
  zero_scale <- scales<1e-12
  scales[zero_scale] <- 1
  z <- sweep(sweep(x,2,center,"-"),2,scales,"/")
  if(any(zero_scale)) z[,zero_scale] <- 0
  max_lam <- max(abs(as.numeric(crossprod(z,y-mean(y)))))/n
  ensure(is.finite(max_lam) && max_lam>0,"Degenerate lambda_max diagnostic")
  alpha_glmnet <- (1-alpha)/(1+alpha)
  lambda_glmnet <- lambda*(1+alpha)
  # Warm-start path for stable glmnet fitting. Only `lambda` is evaluated,
  # all other path points are numerical intermediates, NOT hyperparameters.
  # Teacher alpha=1 is Ridge: no finite all-zero threshold.
  high <- max(lambda_glmnet, 1.1*max_lam/max(alpha_glmnet,0.05))
  lam_path <- exp(seq(log(high),log(lambda_glmnet),length.out=30L))
  if(length(lam_path)>1L && max(abs(lam_path-lambda_glmnet))<1e-12) lam_path <- lambda_glmnet
  fit <- glmnet::glmnet(x=z,y=y,family="binomial",alpha=alpha_glmnet,
                       lambda=lam_path,standardize=FALSE,intercept=TRUE,
                       thresh=1e-8,maxit=100000L)
  ensure(all(fit$jerr == 0), "glmnet returned a fitting/convergence error")
  co <- as.matrix(stats::coef(fit,s=lambda_glmnet,exact=FALSE))[,1]
  ensure(all(is.finite(co)),"Nonfinite elastic-net coefficients")
  list(object=fit,alpha=alpha,lambda=lambda,alpha_glmnet=alpha_glmnet,lambda_glmnet=lambda_glmnet,
       gradient_max=max_lam,lambda_zero_threshold=if(alpha<1) max_lam/(1-alpha) else Inf,
       center=center,scale=scales,n=n,
       last_origin=max(d$m[rows]),last_outcome=max(d$end[rows]),
       coefficients=co,nonzero=sum(abs(co[-1])>1e-8),zero_scale=sum(zero_scale))
}
predict_elastic <- function(fit,x) {
  ensure(identical(colnames(x),names(fit$center)),"Feature names/order changed")
  z <- sweep(sweep(x,2,fit$center,"-"),2,fit$scale,"/")
  p <- as.numeric(stats::predict(fit$object,newx=z,s=fit$lambda_glmnet,type="response"))
  ensure(all(is.finite(p)&p>=0&p<=1),"Invalid probabilities")
  p
}

# One fixed candidate (alpha,lambda) is held constant across all monthly refits.
# Only coefficients and training-row standardization are re-estimated.
recursive_elastic <- function(d,h,cfg,period,alpha,lambda,keep_records=TRUE) {
  ensure(period %in% c("validation","test"),"Invalid period")
  first <- month_number(cfg[[paste0(period,"_start")]])
  last <- month_number(cfg[[paste0(period,"_end")]])
  rows <- which(d$m>=first & d$m<=last)
  ensure(length(rows)>0,"No forecast origins")
  prevalence <- numeric(length(rows))
  preds <- numeric(length(rows)); audit <- vector("list",length(rows))
  fits <- if(keep_records) vector("list",length(rows)) else NULL
  for(j in seq_along(rows)) {
    i <- rows[j]; fit_rows <- mature_rows(d,d$m[i])
    prevalence[j] <- mean(d$y[fit_rows])
    fit <- fit_elastic(d,fit_rows,alpha,lambda)
    ensure(fit$last_outcome<=d$m[i] && fit$last_origin<=d$m[i]-h,"Future outcome leakage")
    preds[j] <- predict_elastic(fit,d$X[i,,drop=FALSE])
    audit[[j]] <- data.frame(model=cfg$model,h=h,period=period,forecast_origin=d$ym[i],
      fit_information_cutoff=d$ym[i],last_train_origin=month_string(fit$last_origin),
      last_train_outcome=month_string(fit$last_outcome),n_train=fit$n,
      alpha=alpha,lambda_absolute=fit$lambda,
      alpha_glmnet=fit$alpha_glmnet,lambda_glmnet=fit$lambda_glmnet,
      gradient_max=fit$gradient_max,lambda_zero_threshold=fit$lambda_zero_threshold,nonzero=fit$nonzero,zero_scale=fit$zero_scale)
    if(keep_records) fits[[j]] <- list(origin=d$ym[i],coefficients=fit$coefficients,
       center=fit$center,scale=fit$scale,alpha=fit$alpha,lambda=fit$lambda,
       alpha_glmnet=fit$alpha_glmnet,lambda_glmnet=fit$lambda_glmnet)
  }
  forecasts <- data.frame(h=h,model=cfg$model,period=period,origin=d$ym[rows],
    probability=preds,prevalence_probability=prevalence,outcome_end=as.character(d$metadata$outcome_end[rows]),
    target=d$y[rows],classification_eligible=d$split[rows]==period,
    fit_information_cutoff=d$ym[rows],
    last_train_origin=vapply(audit,`[[`,character(1),"last_train_origin"),
    last_train_outcome=vapply(audit,`[[`,character(1),"last_train_outcome"),
    n_train=vapply(audit,`[[`,integer(1),"n_train"),
    config_id=cfg$version,stringsAsFactors=FALSE)
  list(forecasts=forecasts,audit=bind_rows(audit),fits=fits)
}

# Optimize predictive hyperparameters on eligible validation labels ONLY.
# No random folds and no test outcomes. Ties prefer fewer active coefficients,
# then larger absolute lambda, then smaller alpha.
tune_elastic <- function(d,h,cfg) {
  grid <- expand.grid(alpha=cfg$alphas,lambda=cfg$lambdas)
  candidates <- vector("list",nrow(grid)); score <- vector("list",nrow(grid))
  for(k in seq_len(nrow(grid))) {
    alpha <- grid$alpha[k]; lam <- grid$lambda[k]
    message(sprintf("h=%d candidate %d/%d alpha=%.2f fixed_lambda=%.6g",h,k,nrow(grid),alpha,lam))
    o <- recursive_elastic(d,h,cfg,"validation",alpha,lam,keep_records=FALSE)
    f <- o$forecasts; eligible <- f$classification_eligible & !is.na(f$target)
    bs <- mean((f$probability[eligible]-f$target[eligible])^2)
    score[[k]] <- data.frame(h=h,model=cfg$model,alpha=alpha,lambda=lam,
      validation_Brier=bs,lambda_at_upper_bound=lam==max(cfg$lambdas),n_validation=sum(eligible),
      mean_nonzero=mean(o$audit$nonzero))
    candidates[[k]] <- o
  }
  scores <- bind_rows(score)
  best <- min(scores$validation_Brier)
  eligible <- which(scores$validation_Brier <= best+1e-12)
  rank <- eligible[order(scores$mean_nonzero[eligible],
                         -scores$lambda[eligible],scores$alpha[eligible])]
  idx <- rank[1]
  if(scores$lambda_at_upper_bound[idx]) message(sprintf("h=%d: selected lambda is at grid upper bound (0.3).",h))
  list(scores=scores,selected=scores[idx,,drop=FALSE],result=candidates[[idx]])
}


# A cohort's decision persists for h monthly return periods; h cohorts receive
# equal post-cost capital at each month-end, following the professor's example.
load_return_data <- function(root) {
  market <- utils::read.csv(file.path(root, "Data/crsp_monthly_index_data/crsp_monthly_index_data.csv"))
  market$date <- as.Date(market$mthcaldt)
  market$m <- month_number(market$date)
  market <- market[order(market$m), ]
  ensure(!anyDuplicated(market$m) && all(diff(market$m) == 1L), "Market calendar mismatch")
  bills <- utils::read.csv(file.path(root, "Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))
  bills$m <- month_number(bills$MCALDT)
  bills <- bills[!is.na(bills$m), ]
  list(market = market, bills = bills)
}

portfolio_inputs <- function(forecasts, returns, h, end, common_rf_security = 2000061L) {
  origin <- month_number(forecasts$origin)
  ensure(all(diff(origin) == 1L), "Forecasts must have a complete monthly calendar")
  b <- returns$bills[returns$bills$KYTREASNOX == c(`1` = 2000061, `3` = 2000062, `6` = 2000063)[[as.character(h)]], ]
  ensure(!anyDuplicated(b$m), "Duplicate matched bill quote")
  yield <- b$TMYLD[match(origin, b$m)]
  ensure(all(is.finite(yield)), "Missing matched bill yield")
  # Require all h active origins to belong to this split. The first h-1
  # forecasts initialize the portfolio; no pre-split signals are imported.
  j <- which(seq_along(origin) >= h & origin + 1L <= month_number(end))
  ensure(length(j) > 1L, "Insufficient complete portfolio months")
  k <- match(origin[j] + 1L, returns$market$m)
  prev <- match(origin[j], returns$market$m)
  ensure(!anyNA(c(k, prev)), "Missing monthly return quote")
  days <- as.numeric(returns$market$date[k] - returns$market$date[prev])
  common <- returns$bills[returns$bills$KYTREASNOX == common_rf_security, ]
  ensure(!anyDuplicated(common$m), "Duplicate common bill quote")
  common_yield <- common$TMYLD[match(origin[j], common$m)]
  ensure(all(is.finite(common_yield)), "Missing common bill yield")
  # Both references cover the SAME one-month interval. Bill maturity never
  # multiplies days: h controls the active vintages, not the accrual interval.
  common_return <- expm1(common_yield * days)
  active <- t(matrix(unlist(lapply(j, function(i) seq.int(i - h + 1L, i))), nrow = h))
  bill_matrix <- expm1(sweep(matrix(yield[active], nrow = length(j), ncol = h), 1L, days, "*"))
  ensure(all(is.finite(bill_matrix)) && all(days > 0), "Invalid accrual proxy")
  list(h = h, j = j, active = active, ym = month_string(origin[j] + 1L), days = days,
       market = returns$market$vwretd[k], bill_matrix = bill_matrix,
       bill_reference = rowMeans(bill_matrix), common_rf_return = common_return,
       common_rf_daily_yield = common_yield, common_rf_security = common_rf_security,
       matched_bill_security = unique(b$KYTREASNOX),
       quote_start = returns$market$date[prev], quote_end = returns$market$date[k],
       origin = forecasts$origin[j])
}

# Solve q = w * (V - c * abs(q)) - S exactly, where S is equity
# immediately before the trade, and q is the signed dollar equity trade.
rebalance_trade <- function(wealth, stock_value, target_weight, cost_rate) {
  ensure(wealth > 0 && stock_value >= -1e-12 && stock_value <= wealth + 1e-12 &&
           target_weight >= 0 && target_weight <= 1 && cost_rate >= 0 && cost_rate < 1,
         "Invalid rebalancing inputs")
  gap <- target_weight * wealth - stock_value
  trade <- if (gap >= 0) gap / (1 + target_weight * cost_rate) else gap / (1 - target_weight * cost_rate)
  cost <- cost_rate * abs(trade)
  list(trade = trade, cost = cost, capital = wealth - cost)
}

simulate_weights <- function(weights, market, bill_contribution, cost_rate, initial_cost = FALSE) {
  n <- length(weights)
  ensure(n == length(market) && n == length(bill_contribution) && n > 1L &&
           all(is.finite(c(weights, market, bill_contribution))) && all(market > -1),
         "Invalid portfolio vectors")
  wealth <- 1; stock <- 0
  out <- vector("list", n)
  for (i in seq_len(n)) {
    before <- wealth
    if (i == 1L && !initial_cost) {
      trade <- list(trade = 0, cost = 0, capital = wealth)
      stock <- weights[i] * wealth
    } else trade <- rebalance_trade(wealth, stock, weights[i], cost_rate)
    gross <- weights[i] * market[i] + bill_contribution[i]
    stock_start <- weights[i] * trade$capital
    stock <- stock_start * (1 + market[i])
    wealth <- trade$capital * (1 + gross)
    ensure(wealth > 0 && stock <= wealth + 1e-10, "Invalid portfolio value")
    out[[i]] <- data.frame(stock_weight = weights[i], pretrade_wealth = before,
      gross_return = gross, net_return = wealth / before - 1,
      signed_stock_trade = trade$trade, traded_fraction = abs(trade$trade) / before,
      cost = trade$cost, cost_fraction = trade$cost / before, wealth = wealth,
      stock_start = stock_start, stock_end = stock)
  }
  do.call(rbind, out)
}

portfolio_path <- function(inputs, signals, cost_rate, cfg, fixed_weight=NULL) {
  ensure(length(signals) >= max(inputs$active) && all(signals %in% 0:1), "Invalid trading signals")
  active_signals <- matrix(signals[inputs$active], nrow = nrow(inputs$active), ncol = inputs$h)
  weights <- rowMeans(active_signals)
  bill <- rowMeans((1 - active_signals) * inputs$bill_matrix)
  if(!is.null(fixed_weight)) {
    ensure(fixed_weight>=0 && fixed_weight<=1,"Invalid fixed weight")
    weights <- rep(fixed_weight,length(weights))
    bill <- (1-fixed_weight)*inputs$bill_reference
  }
  p <- simulate_weights(weights, inputs$market, bill, cost_rate, cfg$initial_allocation_cost)
  cbind(data.frame(month = inputs$ym, last_signal_origin = inputs$origin,
    quote_start = inputs$quote_start, quote_end = inputs$quote_end,
    quote_days = inputs$days, market_return = inputs$market,
    common_rf_return = inputs$common_rf_return,
    matched_ladder_return = inputs$bill_reference,
    bill_reference_return = inputs$bill_reference), p)
}

annual_sharpe <- function(excess) {
  ensure(length(excess) > 1L && all(is.finite(excess)), "Invalid excess-return sample")
  sigma <- stats::sd(excess)
  if (sigma > 1e-12) sqrt(12) * mean(excess) / sigma else NA_real_
}

portfolio_metrics <- function(p, gamma = 3) {
  ensure(all(is.finite(as.matrix(p[, c("net_return", "common_rf_return", "matched_ladder_return")]))),
         "Missing/nonfinite returns or references")
  r <- p$net_return; excess <- r - p$matched_ladder_return
  sharpe <- annual_sharpe(excess)
  common_sharpe <- annual_sharpe(r - p$common_rf_return)
  wealth <- cumprod(1 + r)
  dd <- wealth / cummax(c(1, wealth))[-1L] - 1
  maxdd <- max(0, -min(dd))
  cagr <- tail(wealth, 1)^(12 / length(r)) - 1
  downside <- sqrt(mean(pmin(excess, 0)^2))
  market_excess <- p$market_return - p$matched_ladder_return
  data.frame(n_months = length(r), start = p$month[1], end = tail(p$month, 1),
    total_return = tail(wealth, 1) - 1, CAGR = cagr,
    annual_volatility = stats::sd(r) * sqrt(12),
    sharpe_common_rf = common_sharpe, sharpe_matched_ladder = sharpe,
    # Backward-compatible alias; all new selection/reporting uses explicit names.
    net_Sharpe = sharpe, excess_metric_reference = "matched_ladder",
    Sortino = if (downside > 1e-12) sqrt(12) * mean(excess) / downside else NA_real_,
    max_drawdown = maxdd, Calmar = if (maxdd > 1e-12) cagr / maxdd else NA_real_,
    Ulcer_index = sqrt(mean(dd^2)),
    beta = stats::cov(excess, market_excess) / stats::var(market_excess),
    annual_CE_excess = 12 * (mean(excess) - gamma / 2 * stats::var(excess)),
    average_stock_weight = mean(p$stock_weight),
    annual_turnover = 12 * mean(p$traded_fraction), total_cost_fraction = sum(p$cost_fraction),
    rebalance_trade_months = sum(p$traded_fraction > 1e-12),
    allocation_change_months = sum(abs(diff(p$stock_weight)) > 1e-12))
}

choose_cutoff <- function(grid, objective = "sharpe_matched_ladder") {
  ensure(objective %in% names(grid), "Cutoff objective missing from grid")
  scores <- grid[[objective]]
  eligible <- which(is.finite(scores))
  ensure(length(eligible) > 0, "No cutoff has a defined validation Sharpe")
  best <- max(scores[eligible])
  tied <- eligible[abs(scores[eligible] - best) < 1e-12]
  tied[order(abs(grid$cutoff[tied] - 0.5), grid$cutoff[tied])][1L]
}

forecast_metrics <- function(f, cutoff) {
  f <- f[f$classification_eligible & !is.na(f$target), ]
  y <- f$target; p <- f$probability; pred <- as.integer(p >= cutoff)
  eps <- 1e-15; pc <- pmin(pmax(p, eps), 1 - eps)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  auc <- if (n1 * n0 > 0) (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) else NA_real_
  ratio <- function(a, b) if (b == 0) NA_real_ else a / b
  data.frame(n = length(y), cutoff = cutoff, event_rate = mean(y),
    prevalence_Brier = mean((f$prevalence_probability-y)^2),
    Brier_skill_vs_prevalence = if(mean((f$prevalence_probability-y)^2)>0) 1-mean((p-y)^2)/mean((f$prevalence_probability-y)^2) else NA_real_,
    Brier = mean((p - y)^2), log_loss = -mean(y * log(pc) + (1-y) * log1p(-pc)),
    AUC = auc, accuracy = mean(pred == y), accuracy_at_050 = mean((p >= .5) == y),
    precision = ratio(sum(pred == 1 & y == 1), sum(pred == 1)),
    recall = ratio(sum(pred == 1 & y == 1), sum(y == 1)),
    specificity = ratio(sum(pred == 0 & y == 0), sum(y == 0)),
    mean_probability = mean(p), calibration_gap = mean(p) - mean(y))
}


run_elastic_tests <- function(root,cfg) {
  passed <- character()
  check <- function(name,expr) {force(expr); passed <<- c(passed,name)}
  near <- function(x,y,tol=1e-9) ensure(isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=tol)),"Numerical mismatch")
  check("Exact cost-aware buy/sell weights",{
    for(w in c(0,0.2,0.5,1)) for(s in c(0,25,75,100)) {
      t<-rebalance_trade(100,s,w,cfg$primary_cost)
      near((s+t$trade)/t$capital,w)
    }
  })
  check("Matched and common monthly Sharpe references",{
    p<-data.frame(net_return=c(.02,-.01,.03,.015),common_rf_return=c(.001,.002,.001,.002),
      matched_ladder_return=c(.002,.003,.002,.004),market_return=c(.01,-.02,.025,.012),
      stock_weight=c(.5,.5,.5,.5),traded_fraction=rep(0,4),cost_fraction=rep(0,4),month=paste0("2020-0",1:4))
    m<-portfolio_metrics(p)
    near(m$sharpe_common_rf,sqrt(12)*mean(p$net_return-p$common_rf_return)/sd(p$net_return-p$common_rf_return))
    near(m$sharpe_matched_ladder,sqrt(12)*mean(p$net_return-p$matched_ladder_return)/sd(p$net_return-p$matched_ladder_return))
  })
  check("Teacher penalty maps exactly to glmnet", {
    for(a in cfg$alphas) for(lam in cfg$lambdas) {
      ag <- (1-a)/(1+a); lg <- lam*(1+a)
      near(lg*ag,lam*(1-a))
      near(lg*(1-ag)/2,lam*a)
    }
  })
  check("Fixed grid has 5 x 8 = 40 distinct candidates",{
    ensure(length(cfg$alphas)==5L && length(cfg$lambdas)==8L &&
           length(unique(cfg$lambdas))==8L &&
           all(cfg$lambdas>0) && max(cfg$lambdas)==0.3,"Unexpected lambda grid")
  })
  check("Cutoff grid and selection", {
    near(cfg$cutoffs, (0:20)/20)
    ensure(length(cfg$cutoffs)==21L && sum(abs(cfg$cutoffs-.5)<1e-12)==1L,
           "Expected 21 cutoffs including 0.50")
    g <- data.frame(cutoff=c(.4,.5,.6),sharpe_matched_ladder=c(1,1,1))
    ensure(choose_cutoff(g)==2L,"Cutoff ties should prefer 0.50")
    g$sharpe_matched_ladder <- c(1,NA,2)
    ensure(choose_cutoff(g)==3L,"Nonfinite cutoff scores must be excluded")
  })
  for(h in cfg$horizons) {
    d<-load_elastic_data(root,h,cfg)
    check(paste("Mature label boundary h",h),{
      for(t in c("2005-01","2013-01")) {
        m<-month_number(t); r<-mature_rows(d,m)
        ensure(max(d$m[r])==m-h && max(d$end[r])==m,"Maturity cutoff incorrect")
      }
    })
    check(paste("Future perturbation invariant h",h),{
      now<-month_number("2005-02"); i<-which(d$m==now)
      r<-mature_rows(d,now)
      original<-fit_elastic(d,r,.5,.01)
      changed<-d
      changed$y[which(!is.na(d$end) & d$end>now)] <- 1L
      changed$X[which(d$m>now),] <- changed$X[which(d$m>now),]+100
      revised<-fit_elastic(changed,mature_rows(changed,now),.5,.01)
      near(predict_elastic(original,d$X[i,,drop=FALSE]),predict_elastic(revised,changed$X[i,,drop=FALSE]))
    })
    check(paste("Recursive refit short validation h",h),{
      short<-cfg;short$validation_end<-"2005-03"
      o<-recursive_elastic(d,h,short,"validation",.5,.01,FALSE)
      ensure(nrow(o$forecasts)==3 && all(diff(o$audit$n_train)==1),"Recursive window incorrect")
      ensure(all(month_number(o$audit$last_train_outcome)<=month_number(o$audit$forecast_origin)),"Leak in audit")
    })
    check(paste("Vintage monthly bill accrual h",h),{
      ret<-load_return_data(root)
      ids<-which(d$m>=month_number(cfg$validation_start)&d$m<=month_number(cfg$validation_end))
      f<-data.frame(origin=d$ym[ids])
      inp<-portfolio_inputs(f,ret,h,cfg$validation_end,cfg$common_rf_security)
      ensure(all(inp$days>=25 & inp$days<=35),"Wrong accrual days")
      stock<-portfolio_path(inp,rep(1,length(ids)),cfg$primary_cost,cfg)
      bill<-portfolio_path(inp,rep(0,length(ids)),cfg$primary_cost,cfg)
      near(stock$net_return,inp$market)
      near(bill$net_return,inp$bill_reference)
    })
  }
  data.frame(test=passed,passed=TRUE)
}


write_elastic_report <- function(out,cfg,selected,performance,diagnostics,tests) {
  primary<-performance[performance$period=="test" & abs(performance$cost_rate-cfg$primary_cost)<1e-12,]
  choices<-selected[,c("h","alpha","lambda","validation_Brier","cutoff","validation_sharpe_matched_ladder")]
  lines<-c("# EC4308 Elastic Net — recursive out-of-sample results","",
    paste("Configuration:",cfg$version),"",
    "Penalty convention: lambda * sum(alpha * beta^2 + (1-alpha) * abs(beta)); alpha=0 Lasso, alpha=1 Ridge. Loss is mean binomial negative log-likelihood.",
    "Predictors: complete 132-column frozen X; one independent model for each horizon.",
    "Validation candidates selected by minimum expanding-window Brier; candidate is (alpha, fixed absolute lambda).",
    "Fixed absolute lambda belongs to {0.0001,0.0005,0.002,0.01,0.05,0.1,0.2,0.3} and stays unchanged across origins.",
    "The 40-candidate grid is an exploratory sensitivity analysis. The upper bound 0.3 extends the earlier range; it is not a theoretically optimal value.",
    "After predictive tuning, cutoffs selected from 0–1 in 0.05 steps by validation net matched-ladder Sharpe at 10 bp.",
    "Test coefficients are re-estimated monthly, while alpha, fixed absolute lambda and cutoff are frozen.","",
    "## Frozen validation choices", "",
    paste(capture.output(print(choices,row.names=FALSE)),collapse="\n"),"",
    "## Test performance at 10bp", "",
    paste(capture.output(print(primary[,c("h","strategy","cutoff_rule","CAGR","sharpe_common_rf","sharpe_matched_ladder","max_drawdown","average_stock_weight")],row.names=FALSE)),collapse="\n"),"",
    "## Interpretation", "",
    "Common Sharpe uses 4-week bill accrual for each calendar month. Matched Sharpe uses horizon-specific vintage bill ladder accrual. Both are monthly returns.",
    "The first h-1 signals warm up the h-cohort portfolio. Performance months start February/April/July for h=1/3/6.",
    "Test-period comparisons do not select hyperparameters, cutoffs or model variants. Prior baseline test viewing and supervisor-driven protocol revision must be disclosed.",
    "Monthly maximum drawdown omits intramonth declines. Equity exposure is important when interpreting high Sharpe during bullish periods.",
    "Bill accrual is an idealized proxy and close-to-close execution assumes timely availability of end-of-month inputs.",
    paste(nrow(tests),"focused test groups passed."))
  writeLines(lines,file.path(out,"RESULTS.md"))
}


elastic_net_main <- function(root, stage="all", script_path=.ec4308_source_path) {
  root <- normalizePath(path.expand(root),mustWork=TRUE)
  .libPaths(c(file.path(root,".R-library"),.libPaths()))
  needed <- c("glmnet","digest","jsonlite")
  missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) {
    local_lib <- file.path(root,".R-library")
    dir.create(local_lib,recursive=TRUE,showWarnings=FALSE)
    .libPaths(c(local_lib,.libPaths()))
    utils::install.packages(missing,repos="https://cloud.r-project.org",lib=local_lib)
    ensure(all(vapply(needed,requireNamespace,logical(1),quietly=TRUE)),
           "Package installation incomplete; check the installation messages")
  }


cfg<-elastic_config(); set.seed(cfg$seed)
ensure(length(stage)==1 && stage %in% c("all","validation","test","tests"),"Stage: all|validation|test|tests")
out<-file.path(root,"Results","elastic_net_logit_teacher_fixed03_cutoff005_full_pipeline")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
message("Run stage: ", stage)
message("Results directory: ", out)
inputs<-c(file.path(root,"Data/processed",paste0("model_ready_h",cfg$horizons,".rds")),
 file.path(root,"Data/crsp_monthly_index_data/crsp_monthly_index_data.csv"),
 file.path(root,"Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))
ensure(all(file.exists(inputs)),"Required input files missing: check project layout")
input_hashes<-setNames(vapply(inputs,file_hash,character(1)),substring(inputs,nchar(root)+2L))
code_hashes <- c(single_script=file_hash(script_path))
message("Running numerical and data-timing tests...")
tests<-run_elastic_tests(root,cfg); write_csv(tests,file.path(out,"tests.csv"))
if(stage=="tests") return(invisible(tests))
returns<-load_return_data(root)
data<-setNames(lapply(cfg$horizons,function(h) load_elastic_data(root,h,cfg)),cfg$horizons)
state_path<-file.path(out,"validation_state.rds")
if(stage %in% c("all","validation")) {
  score_rows<-list(); selection_rows<-list(); forecast_rows<-list(); cutoff_rows<-list(); audit_rows<-list()
  for(h in cfg$horizons) {
    d<-data[[as.character(h)]]
    tuned<-tune_elastic(d,h,cfg)
    f<-tuned$result$forecasts
    inp<-portfolio_inputs(f,returns,h,cfg$validation_end,cfg$common_rf_security)
    cutoff_grid<-bind_rows(lapply(cfg$cutoffs,function(cutoff) {
      p<-portfolio_path(inp,as.integer(f$probability>=cutoff),cfg$primary_cost,cfg)
      cbind(data.frame(h=h,model=cfg$model,cutoff=cutoff),portfolio_metrics(p,cfg$risk_aversion))
    }))
    winner<-cutoff_grid[choose_cutoff(cutoff_grid,cfg$cutoff_objective),,drop=FALSE]
    fixed<-cutoff_grid[abs(cutoff_grid$cutoff-cfg$comparison_cutoff)<1e-12,,drop=FALSE]
    ensure(nrow(fixed)==1,"Missing 0.50 benchmark")
    sel<-cbind(tuned$selected,data.frame(cutoff=winner$cutoff,
      validation_sharpe_matched_ladder=winner$sharpe_matched_ladder,
      validation_sharpe_common_rf=winner$sharpe_common_rf,
      fixed_050_validation_sharpe_matched_ladder=fixed$sharpe_matched_ladder,
      selection_cost=cfg$primary_cost))
    f$cutoff<-sel$cutoff;f$signal<-as.integer(f$probability>=sel$cutoff)
    f$cutoff_050<-cfg$comparison_cutoff;f$signal_050<-as.integer(f$probability>=cfg$comparison_cutoff)
    key<-as.character(h)
    score_rows[[key]]<-tuned$scores;selection_rows[[key]]<-sel
    forecast_rows[[key]]<-f; cutoff_rows[[key]]<-cutoff_grid
    audit_rows[[key]]<-tuned$result$audit
  }
  scores<-bind_rows(score_rows); selected<-bind_rows(selection_rows)
  val_forecasts<-bind_rows(forecast_rows);cutoff_grid<-bind_rows(cutoff_rows)
  ensure(nrow(cutoff_grid)==63L,"Expected 21 cutoffs per horizon")
  val_audit<-bind_rows(audit_rows)
  frozen<-list(config=cfg,selected=selected,input_hashes=input_hashes,code_hashes=code_hashes,
    frozen_at=format(Sys.time(),tz="Asia/Singapore",usetz=TRUE))
  state<-list(scores=scores,selected=selected,forecasts=val_forecasts,cutoff_grid=cutoff_grid,
    audit=val_audit,frozen=frozen)
  saveRDS(state,state_path)
  saveRDS(frozen,file.path(out,"frozen_specification.rds"))
  write_csv(scores,file.path(out,"validation_tuning.csv"))
  write_csv(selected,file.path(out,"selected_settings.csv"))
  write_csv(cutoff_grid,file.path(out,"validation_cutoff_grid.csv"))
  write_csv(val_forecasts,file.path(out,"validation_forecasts.csv"))
  write_csv(val_audit,file.path(out,"validation_refit_audit.csv"))
  message("Validation frozen: selected hyperparameters and cutoffs for all horizons")
} else {
  ensure(file.exists(state_path),"No validation state: run --stage=validation first")
  state<-readRDS(state_path)
  ensure(identical(cfg,state$frozen$config) && identical(input_hashes,state$frozen$input_hashes) &&
         identical(code_hashes,state$frozen$code_hashes),"Code/data/config changed since validation: stop")
  selected<-state$selected;val_forecasts<-state$forecasts
}
if(stage=="validation") return(invisible(state))
message("Computing test predictions; selection is frozen...")
test_rows<-list(); test_audits<-list(); fit_records<-list()
for(h in cfg$horizons) {
  d<-data[[as.character(h)]]; sel<-selected[selected$h==h,,drop=FALSE]
  ensure(nrow(sel)==1,"Missing unique frozen horizon settings")
  result<-recursive_elastic(d,h,cfg,"test",sel$alpha,sel$lambda,TRUE)
  f<-result$forecasts
  f$cutoff<-sel$cutoff; f$signal<-as.integer(f$probability>=sel$cutoff)
  f$cutoff_050<-cfg$comparison_cutoff;f$signal_050<-as.integer(f$probability>=cfg$comparison_cutoff)
  key<-as.character(h); test_rows[[key]]<-f;test_audits[[key]]<-result$audit
  fit_records[[key]]<-result$fits
}
test_forecasts<-bind_rows(test_rows)
all_forecasts<-rbind(val_forecasts,test_forecasts)
audit<-rbind(state$audit,bind_rows(test_audits))
ensure(all(month_number(audit$last_train_outcome)<=month_number(audit$forecast_origin)),"Audit leakage")
ensure(nrow(state$audit)==96L*length(cfg$horizons) && nrow(bind_rows(test_audits))==156L*length(cfg$horizons),"Missing monthly refits")
paths<-list();perf<-list();diagnostics<-list();refs<-list()
for(period in c("validation","test")) for(h in cfg$horizons) {
  f<-all_forecasts[all_forecasts$period==period & all_forecasts$h==h,,drop=FALSE]
  end<-cfg[[paste0(period,"_end")]]
  inp<-portfolio_inputs(f,returns,h,end,cfg$common_rf_security)
  refs[[paste(period,h)]]<-data.frame(period=period,h=h,month=inp$ym,
    quote_start=inp$quote_start,quote_end=inp$quote_end,quote_days=inp$days,
    common_rf_return=inp$common_rf_return,matched_ladder_return=inp$bill_reference,
    common_rf_security=inp$common_rf_security,matched_bill_security=inp$matched_bill_security)
  variants<-list(validation_tuned=f$signal,fixed_050=f$signal_050,
    all_stocks=rep(1L,nrow(f)),all_bills=rep(0L,nrow(f)),fixed_50_50=rep(0L,nrow(f)))
  for(rule in names(variants)) {
    strat<-if(rule %in% c("all_stocks","all_bills","fixed_50_50")) rule else cfg$model
    cutoff_rule<-if(strat==cfg$model) rule else "not_applicable"
    cutoff<-if(rule=="validation_tuned") unique(f$cutoff) else if(rule=="fixed_050") cfg$comparison_cutoff else NA_real_
    if(strat==cfg$model) diagnostics[[paste(period,h,rule)]]<-cbind(data.frame(period=period,h=h,model=cfg$model,cutoff_rule=rule),forecast_metrics(f,cutoff))
    for(cost in cfg$costs) {
      p<-portfolio_path(inp,variants[[rule]],cost,cfg,
        fixed_weight=if(rule=="fixed_50_50") 0.5 else NULL)
      tags<-data.frame(period=period,h=h,strategy=strat,cutoff_rule=cutoff_rule,cutoff=cutoff,cost_rate=cost)
      key<-paste(period,h,rule,cost)
      paths[[key]]<-cbind(tags,p)
      perf[[key]]<-cbind(tags,portfolio_metrics(p,cfg$risk_aversion))
      ensure(max(month_number(p$month))<=month_number(end),"Portfolio leaks split boundary")
    }
  }
}
paths<-bind_rows(paths);performance<-bind_rows(perf)
common<-list()
for(period in c("validation","test")) {
  start<-cfg[[paste0(period,"_start")]]
  common_start<-month_string(month_number(start)+max(cfg$horizons))
  q<-paths[paths$period==period & paths$month>=common_start,]
  groups<-interaction(q$h,q$strategy,q$cutoff_rule,q$cost_rate,drop=TRUE)
  for(k in unique(groups)) {
    z<-q[groups==k,]
    common[[paste(period,k)]]<-cbind(z[1,c("period","h","strategy","cutoff_rule","cutoff","cost_rate")],portfolio_metrics(z,cfg$risk_aversion))
  }
}
ensure(identical(input_hashes,setNames(vapply(inputs,file_hash,character(1)),names(input_hashes))),"Inputs changed during run")
write_csv(test_forecasts,file.path(out,"test_forecasts.csv"))
write_csv(audit,file.path(out,"refit_audit.csv"))
write_csv(paths,file.path(out,"monthly_portfolios.csv"))
write_csv(performance,file.path(out,"performance.csv"))
write_csv(bind_rows(common),file.path(out,"common_window_performance.csv"))
write_csv(bind_rows(refs),file.path(out,"benchmark_returns.csv"))
write_csv(bind_rows(diagnostics),file.path(out,"forecast_metrics.csv"))
write_csv(data.frame(feature=colnames(data[["1"]]$X)),file.path(out,"feature_names.csv"))
write_csv(data.frame(file=names(input_hashes),sha256=unname(input_hashes)),file.path(out,"input_hashes.csv"))
saveRDS(fit_records,file.path(out,"test_fits_by_origin.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
jsonlite::write_json(list(version=cfg$version,selected=selected,input_hashes=input_hashes,
 code_hashes=code_hashes,tests_passed=nrow(tests),test_used_for_selection=FALSE,
 note="Fixed-lambda model specified after inspection of previous elastic-net test results; exploratory sensitivity, NOT virgin holdout"),
 file.path(out,"manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
write_elastic_report(out,cfg,selected,performance,bind_rows(diagnostics),tests)
message("Validation and test complete. Results directory: ", out)
message("Summary: ",file.path(out,"RESULTS.md"))
message("Sharpe ratios and drawdowns: ",file.path(out,"performance.csv"))



  invisible(list(selected=selected, performance=performance, output_dir=out))
}

.ec4308_cli <- commandArgs(TRUE)
.ec4308_stage <- grep("^--stage=",.ec4308_cli,value=TRUE)
.ec4308_root <- grep("^--root=",.ec4308_cli,value=TRUE)
.ec4308_result <- elastic_net_main(
  root=if(length(.ec4308_root)) sub("^--root=","",.ec4308_root[1]) else
    PROJECT_ROOT,
  stage=if(length(.ec4308_stage)) sub("^--stage=","",.ec4308_stage[1]) else
    "all"
)
