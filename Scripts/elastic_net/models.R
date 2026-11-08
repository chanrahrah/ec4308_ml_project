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

# Derive lambda_max using ONLY current fitting rows. This is the LASSO
# all-zero-slopes threshold (for binomial logistic with mean-response intercept),
# used as a common positive reference for all alphas including ridge.
fit_elastic <- function(d,rows,alpha,lambda_ratio) {
  ensure(length(rows)>40L && length(unique(d$y[rows]))==2L,"Too few fitting rows or one class")
  x <- d$X[rows,,drop=FALSE]; y <- as.numeric(d$y[rows]); n <- length(rows)
  center <- colMeans(x); scales <- apply(x,2,stats::sd)
  ensure(all(is.finite(scales)),"Invalid predictor scaling")
  zero_scale <- scales<1e-12
  scales[zero_scale] <- 1
  z <- sweep(sweep(x,2,center,"-"),2,scales,"/")
  if(any(zero_scale)) z[,zero_scale] <- 0
  max_lam <- max(abs(as.numeric(crossprod(z,y-mean(y)))))/n
  ensure(is.finite(max_lam) && max_lam>0,"lambda_max degenerate")
  lam <- max_lam*lambda_ratio
  ensure(lam>0 && lambda_ratio<=1 && lambda_ratio>0,"Invalid lambda_ratio")
  # Several path points improve coordinate-descent warm starts, but they are
  # NOT extra hyperparameter candidates. The last value is the chosen lambda.
  lam_path <- if(abs(lam-max_lam)<1e-12) max_lam else exp(seq(log(max_lam),log(lam),length.out=12L))
  fit <- glmnet::glmnet(x=z,y=y,family="binomial",alpha=alpha,
                       lambda=lam_path,standardize=FALSE,intercept=TRUE,
                       thresh=1e-8,maxit=100000L)
  co <- as.matrix(stats::coef(fit,s=lam,exact=FALSE))[,1]
  ensure(all(is.finite(co)),"Nonfinite elastic-net coefficients")
  list(object=fit,alpha=alpha,lambda_ratio=lambda_ratio,lambda=lam,
       lambda_max=max_lam,center=center,scale=scales,
       n=n,last_origin=max(d$m[rows]),last_outcome=max(d$end[rows]),
       coefficients=co,nonzero=sum(abs(co[-1])>1e-8),zero_scale=sum(zero_scale))
}
predict_elastic <- function(fit,x) {
  ensure(identical(colnames(x),names(fit$center)),"Feature names/order changed")
  z <- sweep(sweep(x,2,fit$center,"-"),2,fit$scale,"/")
  p <- as.numeric(stats::predict(fit$object,newx=z,s=fit$lambda,type="response"))
  ensure(all(is.finite(p)&p>=0&p<=1),"Invalid probabilities")
  p
}

# One candidate (alpha,lambda_ratio) is held constant as an update policy.
# At each origin, fit new coefficients and a training-specific absolute lambda.
recursive_elastic <- function(d,h,cfg,period,alpha,lambda_ratio,keep_records=TRUE) {
  ensure(period %in% c("validation","test"),"Invalid period")
  first <- month_number(cfg[[paste0(period,"_start")]])
  last <- month_number(cfg[[paste0(period,"_end")]])
  rows <- which(d$m>=first & d$m<=last)
  ensure(length(rows)>0,"No forecast origins")
  preds <- numeric(length(rows)); audit <- vector("list",length(rows))
  fits <- if(keep_records) vector("list",length(rows)) else NULL
  for(j in seq_along(rows)) {
    i <- rows[j]; fit_rows <- mature_rows(d,d$m[i])
    fit <- fit_elastic(d,fit_rows,alpha,lambda_ratio)
    ensure(fit$last_outcome<=d$m[i] && fit$last_origin<=d$m[i]-h,"Future outcome leakage")
    preds[j] <- predict_elastic(fit,d$X[i,,drop=FALSE])
    audit[[j]] <- data.frame(model=cfg$model,h=h,period=period,forecast_origin=d$ym[i],
      fit_information_cutoff=d$ym[i],last_train_origin=month_string(fit$last_origin),
      last_train_outcome=month_string(fit$last_outcome),n_train=fit$n,
      alpha=alpha,lambda_ratio=lambda_ratio,lambda_absolute=fit$lambda,
      lambda_max=fit$lambda_max,nonzero=fit$nonzero,zero_scale=fit$zero_scale)
    if(keep_records) fits[[j]] <- list(origin=d$ym[i],coefficients=fit$coefficients,
       center=fit$center,scale=fit$scale,lambda=fit$lambda)
  }
  forecasts <- data.frame(h=h,model=cfg$model,period=period,origin=d$ym[rows],
    probability=preds,outcome_end=as.character(d$metadata$outcome_end[rows]),
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
# then stronger lambda ratio, then smaller alpha.
tune_elastic <- function(d,h,cfg) {
  grid <- expand.grid(alpha=cfg$alphas,lambda_ratio=cfg$lambda_ratios)
  candidates <- vector("list",nrow(grid)); score <- vector("list",nrow(grid))
  for(k in seq_len(nrow(grid))) {
    alpha <- grid$alpha[k]; ratio <- grid$lambda_ratio[k]
    message(sprintf("h=%d candidate %d/%d alpha=%.2f lambda_ratio=%.3f",h,k,nrow(grid),alpha,ratio))
    o <- recursive_elastic(d,h,cfg,"validation",alpha,ratio,keep_records=FALSE)
    f <- o$forecasts; eligible <- f$classification_eligible & !is.na(f$target)
    bs <- mean((f$probability[eligible]-f$target[eligible])^2)
    score[[k]] <- data.frame(h=h,model=cfg$model,alpha=alpha,lambda_ratio=ratio,
      validation_Brier=bs,n_validation=sum(eligible),
      mean_nonzero=mean(o$audit$nonzero))
    candidates[[k]] <- o
  }
  scores <- bind_rows(score)
  best <- min(scores$validation_Brier)
  eligible <- which(scores$validation_Brier <= best+1e-12)
  rank <- eligible[order(scores$mean_nonzero[eligible],
                         -scores$lambda_ratio[eligible],scores$alpha[eligible])]
  idx <- rank[1]
  list(scores=scores,selected=scores[idx,,drop=FALSE],result=candidates[[idx]])
}
