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
