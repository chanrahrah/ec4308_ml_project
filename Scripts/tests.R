run_tests <- function() {
  passed <- character()
  check <- function(name,expr) { force(expr); passed <<- c(passed,name) }
  fails <- function(expr) inherits(tryCatch({force(expr); NULL},error=function(e)e),"error")
  eq <- function(x,y,tol=1e-12) stopifnot(isTRUE(all.equal(x,y,tolerance=tol,check.attributes=FALSE)))
  m <- seq.int(month_id("2000-01"),length.out=60)
  check("monthly calendar rejects gaps, duplicates and reversals", {
    stopifnot(fails(assert_monthly(m[-4])),fails(assert_monthly(c(m,m[60]))),fails(assert_monthly(rev(m))))
  })
  check("CRSP continuous daily yield and nominal-percent conversion agree", {
    y <- c(0,.0001,.0002); eq(daily_annual(y),exp(y*365)-1)
    eq(daily_annual((y*36500)/36500),daily_annual(y)); stopifnot(fails(daily_annual(-1)))
    # CRSP's published semiannual-yield identity independently fixes the convention.
    semiannual <- 2*(exp(y*182.5)-1)
    eq(daily_annual(y),(1+semiannual/2)^2-1)
    stopifnot(max(abs(daily_annual(y)-((1+y)^365-1)))>1e-6)
    eq(daily_accrual(.0001,33),exp(.0001*33)-1)
  })
  check("distribution proxy uses twelve return pairs without extra warmup", {
    v <- distribution_proxy(rep(.002,24),rep(0,24)); stopifnot(all(is.na(v[1:11]))); eq(v[12:24],rep(log(.024),13))
  })
  check("distribution proxy cannot bridge missing months", {
    rx <- rep(0,36); rx[16]<-NA; v<-distribution_proxy(rx+.002,rx)
    stopifnot(all(is.na(v[16:27]))); eq(v[28],log(.024))
  })
  check("stock windows need observed endpoints and survive price splits", {
    r <- rep(.01,12); p <- 100*1.01^(0:11); p[7:12]<-p[7:12]/2
    a <- stock_signals(m[1:12],r,p); stopifnot(all(is.na(a$ma[1:9])),all(is.na(a$hi[1:11])),a$ma[12]==1,a$hi[12]==1)
    gap <- stock_signals(m[c(1:5,7:12)],r[-6],p[-6]); stopifnot(is.na(tail(gap$ma,1)),is.na(tail(gap$hi,1)))
  })
  check("flat price is simultaneously high and low", {
    a<-stock_signals(m[1:12],rep(0,12),rep(10,12)); eq(a$hi[12]-a$lo[12],0); eq(a$ma[12],0)
  })
  check("unemployment shutdown publication timing", {
    months<-seq.int(month_id("2025-07"),length.out=5)
    info<-available_macro(months,c(.041,.042,.043,NA,.046),months+1L,setNames(month_id("2025-11"),"2025-09"))
    stopifnot(info[m==month_id("2025-10")]$source_month=="2025-08",info[m==month_id("2025-11")]$source_month=="2025-09",
      info[m==month_id("2025-12")]$source_month=="2025-11")
  })
  raw <- data.table(m=m)
  for(n in RAW_NAMES) raw[,(n):=seq_len(60)/100]
  raw[,dispersion:=seq_len(60)/100+.01]
  check("first differences, signed spreads and prior-only rolling means", {
    a<-transform_predictors(raw); eq(a$cpi_infl_change1[40],.01)
    eq(a$bill_yld_rel36[40],raw$bill_yld_matched[40]-mean(raw$bill_yld_matched[4:39]))
    eq(a$log_dispersion_change1[40],log(raw$dispersion[40]/raw$dispersion[39]))
    b<-copy(raw); b[,term_spread:=-term_spread]; eq(transform_predictors(b)$term_spread_change1[40],-.01)
  })
  check("transformations and lags cannot read future observations", {
    a<-transform_predictors(raw); b<-copy(raw); b[50:60,(RAW_NAMES):=lapply(.SD,function(x)x*10),.SDcols=RAW_NAMES]
    eq(a[1:49],transform_predictors(b)[1:49]); l<-lag_features(a); eq(l$cpi_infl_change1_lag5[50],a$cpi_infl_change1[45])
    stopifnot(ncol(l)==133)
  })
  check("first difference propagates internal missingness", {
    b<-copy(raw); b[45,term_spread:=NA_real_]; a<-transform_predictors(b)
    stopifnot(all(is.na(a$term_spread_change1[45:46])),is.finite(a$term_spread_change1[47]))
  })
  check("log dispersion changes require positive levels", {
    b<-copy(raw); b[45,dispersion:=0]; stopifnot(fails(transform_predictors(b)))
  })
  check("stationarity rejects missing-month compression and constant series", {
    x<-sin(seq_along(m)); x[30]<-NA
    stopifnot(fails(stationarity_one(x,m)),fails(stationarity_one(rep(1,60),m)),fails(stationarity_one(sin(seq_len(59)),m[-30])))
  })
  check("ADF zero-lag selection is possible and urca statistic reproduces OLS", {
    set.seed(4308); x<-rnorm(369); lag<-select_adf_lag(x,12)
    a<-stationarity_one(x,seq.int(1000L,length.out=369)); stopifnot(a$adf_lags==lag,a$verdict=="stationarity supported")
    z<-as.numeric(scale(x)); d<-diff(z); E<-embed(d,lag+1); X<-cbind(1,z[(lag+1):(length(z)-1)])
    if(lag) X<-cbind(X,E[,2:(lag+1),drop=FALSE])
    fit<-lm.fit(X,E[,1]); se<-sqrt(diag(chol2inv(qr.R(fit$qr)))*sum(fit$residuals^2)/fit$df.residual)
    eq(a$adf_stat,unname(fit$coefficients[2]/se[2]),1e-10)
  })
  check("KPSS statistic matches independent Bartlett covariance calculation", {
    set.seed(1); x<-rnorm(100); bw<-4; res<-x-mean(x)
    variance<-sum(res^2)/100+2/100*sum(vapply(1:bw,function(j)(1-j/(bw+1))*sum(res[(j+1):100]*res[1:(100-j)]),numeric(1)))
    manual<-sum(cumsum(res)^2)/(100^2*variance)
    eq(urca::ur.kpss(x,use.lag=bw)@teststat,manual)
  })
  check("targets use only t+1 through t+h and strict outperformance", {
    months<-seq.int(month_id("2000-01"),length.out=12)
    dates<-as.Date(paste0(month_label(months+1),"-01"))-1
    r<-rep(.01,12); r[1]<-.5
    z<-forward_targets(months,r,dates,rep(0,12),3)
    eq(z$equity_forward[1],1.01^3-1); stopifnot(z$target[1]==1,z$outcome_end[1]=="2000-04",all(is.na(z$target[10:12])))
    zero<-forward_targets(months,rep(0,12),dates,rep(0,12),1); stopifnot(zero$target[1]==0)
    positive<-forward_targets(months,r,dates,rep(.0001,12),3)
    eq(positive$bill_accrual_proxy[1],exp(as.numeric(dates[4]-dates[1])*.0001)-1)
  })
  check("outcome boundary purges prevent training/validation leakage", {
    origin<-month_id(c("2004-09","2004-10","2005-01","2012-10","2025-10"))
    split<-model_split(origin,month_label(origin+3),rep(1L,5))
    stopifnot(identical(split,c("train","purged_boundary","validation","purged_boundary","purged_boundary")))
  })
  write_json(list(passed=length(passed),tests=as.list(passed),status="passed"),"R_tests.json")
  message(length(passed)," R test groups passed.")
  invisible(passed)
}
