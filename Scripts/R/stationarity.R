# stationarity.R tests the transformed series and their lags using the training period.
# ADF lag selection includes zero and compares all lags on a common sample.
# The selected regression is refit on its full usable sample using urca.
select_adf_lag <- function(x,cap,criterion="AIC",trend=FALSE) {
  delta <- diff(x); lagmat <- embed(delta,cap+1L)
  y <- lagmat[,1]; level <- x[(cap+1L):(length(x)-1L)]
  deterministic <- cbind(constant=1,level=level)
  if(trend) deterministic <- cbind(deterministic,trend=seq_along(y))
  score <- numeric(cap+1L)
  for(p in 0:cap) {
    X <- if(p) cbind(deterministic,lagmat[,2:(p+1L),drop=FALSE]) else deterministic
    fit <- lm.fit(X,y)
    penalty <- if(criterion=="AIC") 2 else log(length(y))
    score[p+1L] <- length(y)*log(sum(fit$residuals^2)/length(y)) + penalty*fit$rank
  }
  which.min(score)-1L
}
kpss_bandwidth <- function(residuals) {
  n <- length(residuals); q <- floor(n^(2/9))
  s0 <- sum(residuals^2)/n; s1 <- 0
  for(i in seq_len(q)) {
    cross <- sum(residuals[(i+1):n]*residuals[1:(n-i)])*2/n
    s0 <- s0+cross; s1 <- s1+i*cross
  }
  if(!is.finite(s0) || s0==0) stop("Undefined KPSS bandwidth")
  min(n-1L,floor(1.1447*((s1/s0)^2)^(1/3)*n^(1/3)))
}
stationarity_one <- function(x,m,cap=12L,criterion="AIC",bandwidth="auto",trend=FALSE) {
  assert_monthly(m); assert_numeric(x)
  observed <- which(!is.na(x))
  if(!length(observed)) stop("No observations for stationarity")
  rows <- seq.int(min(observed),max(observed)); x <- x[rows]; m <- m[rows]
  if(anyNA(x)) stop("Internal missing months cannot be compressed for stationarity")
  if(length(x)<24 || diff(range(x))<=10*.Machine$double.eps*max(abs(x))) stop("Too short or constant series")
  x <- as.numeric(scale(x)); n <- length(x)
  cap <- min(cap,floor(n/2)-if(trend) 3L else 2L)
  lag <- select_adf_lag(x,cap,criterion,trend)
  a <- urca::ur.df(x,type=if(trend) "trend" else "drift",lags=lag,selectlags="Fixed")
  residual <- if(trend) residuals(lm(x~seq_along(x))) else x-mean(x)
  bw <- if(identical(bandwidth,"auto")) kpss_bandwidth(residual) else as.integer(bandwidth)
  k <- urca::ur.kpss(x,type=if(trend) "tau" else "mu",use.lag=bw)
  astat <- as.numeric(a@teststat[1,1]); acv <- as.numeric(a@cval[1,"5pct"])
  kstat <- k@teststat; kcv <- as.numeric(k@cval[1,"5pct"])
  reject_adf <- astat<acv; reject_kpss <- kstat>kcv
  verdict <- if(reject_adf && !reject_kpss) "stationarity supported" else if(!reject_adf && reject_kpss) "nonstationarity concern" else if(reject_adf) "conflicting tests" else "inconclusive"
  lb_lag <- min(max(12,lag+3),floor(length(a@res)/3))
  lb <- if(lb_lag>lag) Box.test(a@res,lag=lb_lag,type="Ljung-Box",fitdf=lag)$p.value else NA_real_
  # Report official urca critical values, not fabricated ADF p-values.
  data.table(n=n,start=month_label(min(m)),end=month_label(max(m)),deterministic=if(trend) "constant+trend" else "constant",
    adf_stat=astat,adf_cv_5pct=acv,adf_reject=reject_adf,adf_lag_cap=cap,adf_lags=lag,criterion=criterion,
    adf_residual_lb_p=lb,kpss_stat=kstat,kpss_cv_5pct=kcv,kpss_reject=reject_kpss,kpss_lags=bw,
    verdict=verdict)
}
stationarity_all <- function() {
  configs <- list(primary=list(),aic24=list(cap=24L),aic36=list(cap=36L),bic24=list(cap=24L,criterion="BIC"),
    kpss12=list(bandwidth=12L),kpss24=list(bandwidth=24L),kpss36=list(bandwidth=36L),trend24=list(cap=24L,trend=TRUE))
  results <- list(); j <- 0L
  for(h in c(1,3,6)) for(kind in c("raw","transformed")) {
    d <- read_table(paste0("base_",kind,"_h",h,".csv"))[m>=START & m<=TRAIN_END]
    for(n in setdiff(names(d),"m")) for(config in names(configs)) {
      j <- j+1L
      row <- do.call(stationarity_one,c(list(x=d[[n]],m=d$m),configs[[config]]))
      # Use set() so the sample-size column 'n' cannot mask the predictor-name loop variable.
      set(row,j="h",value=h); set(row,j="kind",value=kind)
      set(row,j="predictor",value=n); set(row,j="configuration",value=config)
      results[[j]] <- row
    }
  }
  result <- rbindlist(results)
  setcolorder(result,c("h","kind","predictor","configuration",setdiff(names(result),c("h","kind","predictor","configuration"))))
  fwrite(result,file.path(AUDIT,"stationarity_sensitivity.csv"))
  primary <- result[configuration=="primary"]
  fwrite(primary,file.path(AUDIT,"stationarity_primary.csv"))
  summary <- primary[,.(n=.N),by=.(h,kind,verdict)]
  fwrite(summary,file.path(AUDIT,"stationarity_summary.csv")); print(summary)
  lagged_results <- list()
  for(horizon in c(1,3,6)) {
    features <- read_table(paste0("features_h",horizon,"_132.csv"))[m>=START & m<=TRAIN_END]
    for(feature in FEATURE_NAMES) {
      row <- stationarity_one(features[[feature]],features$m)
      set(row,j="h",value=horizon); set(row,j="feature",value=feature)
      lagged_results[[length(lagged_results)+1L]] <- row
    }
  }
  lagged_results <- rbindlist(lagged_results)
  fwrite(lagged_results,file.path(AUDIT,"stationarity_lagged_primary.csv"))
  # Preserve the evidence behind the two additional R-stage specification changes.
  # Reconstruct both prior definitions from raw R outputs; no legacy files are read.
  review <- list()
  for(horizon in c(1,3,6)) {
    raw <- read_table(paste0("base_raw_h",horizon,".csv"))
    current <- read_table(paste0("base_transformed_h",horizon,".csv"))
    keep <- raw$m>=START & raw$m<=TRAIN_END
    variants <- list(dispersion_rel36=log(raw$dispersion/roll(lag1(raw$dispersion),36)),
      log_dispersion_change1=current$log_dispersion_change1,
      tail_imb=raw$tail_imb,tail_imb_change1=current$tail_imb_change1)
    for(name in names(variants)) for(k in 0:5) {
      row <- stationarity_one(lag1(variants[[name]],k)[keep],raw$m[keep])
      set(row,j="h",value=horizon); set(row,j="feature",value=paste0(name,"_lag",k))
      set(row,j="stage",value=if(name %in% c("dispersion_rel36","tail_imb")) "prior_definition" else "current_definition")
      review[[length(review)+1L]] <- row
    }
  }
  fwrite(rbindlist(review),file.path(AUDIT,"stationarity_additional_transform_review.csv"))
  write_json(list(training_start="1974-04",training_end="2004-12",primary_rows=nrow(primary),sensitivity_rows=nrow(result),
    all_22_primary_supported=all(primary[kind=="transformed"]$verdict=="stationarity supported"),
    lagged_primary_counts=as.list(table(lagged_results$verdict)),
    package="urca",version=as.character(packageVersion("urca")),
    note="Critical-value decisions, not proof. Residual dependence and bandwidth sensitivity reported."),"stationarity_manifest.json")
  invisible(result)
}
