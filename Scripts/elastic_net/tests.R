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
      original<-fit_elastic(d,r,.5,.1)
      changed<-d
      changed$y[which(!is.na(d$end) & d$end>now)] <- 1L
      changed$X[which(d$m>now),] <- changed$X[which(d$m>now),]+100
      revised<-fit_elastic(changed,mature_rows(changed,now),.5,.1)
      near(predict_elastic(original,d$X[i,,drop=FALSE]),predict_elastic(revised,changed$X[i,,drop=FALSE]))
    })
    check(paste("Recursive refit short validation h",h),{
      short<-cfg;short$validation_end<-"2005-03"
      o<-recursive_elastic(d,h,short,"validation",.5,.1,FALSE)
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
