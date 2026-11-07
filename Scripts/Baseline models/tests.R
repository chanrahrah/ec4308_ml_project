# Focused numerical and information-set checks; called before any model run.
run_baseline_tests <- function(root, cfg) {
  passed <- character()
  check <- function(name, expr) { force(expr); passed <<- c(passed, name) }
  near <- function(a, b, tol = 1e-10) ensure(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)), "Numerical test failed")

  check("Professor's 12-month six-cohort example", {
    signals <- c(1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,0,1)
    market <- c(.03,-.02,.04,-.06,.02,.05,-.08,.03,-.01,.06,-.04,.02)
    w <- vapply(1:12, function(i) mean(signals[i:(i+5)]), numeric(1))
    p <- simulate_weights(w, market, (1-w)*.003, 0)
    value <- 100; expected <- numeric(12)
    for (i in 1:12) {
      stock_end <- w[i]*value*(1+market[i]); bond_end <- (1-w[i])*value*1.003
      value <- stock_end + bond_end; expected[i] <- value
    }
    near(p$wealth*100, expected)
    near(p$signed_stock_trade[-1]*100, w[-1]*expected[-12] - w[-12]*c(100,expected[1:10])*(1+market[-12]))
  })
  check("Exact post-cost target for purchases, sales and endpoints", {
    for (s in c(0, 20, 50, 100)) for (w in c(0,.2,.5,1)) for (c in c(0,.001,.002)) {
      t <- rebalance_trade(100,s,w,c)
      near((s+t$trade)/t$capital,w)
      near(t$cost,abs(t$trade)*c)
      near(t$capital,100-t$cost)
    }
  })
  check("Constant target still requires drift rebalancing", {
    p <- simulate_weights(rep(.5,3), c(.1,-.05,.02), rep(.0015,3), .001)
    ensure(p$signed_stock_trade[2] < 0 && p$cost[2] > 0, "Drift trade omitted")
  })
  returns <- load_return_data(root)
  for (h in cfg$horizons) {
    d <- load_model_data(root,h,cfg)
    check(paste("Boundary maturity and future-perturbation invariance h",h), {
      for (origin_string in c("2005-01", "2013-01")) {
      origin <- month_number(origin_string)
      rows <- mature_rows(d,origin)
      ensure(max(d$m[rows]) == origin-h && max(d$end[rows]) == origin, "Wrong maturity cutoff")
      ensure(all(d$m[d$rows$train] <= month_number(cfg$train_end)-h), "Train purge failed")
      changed <- d
      unavailable <- which(is.na(d$end) | d$end > origin)
      changed$y[unavailable] <- 1L
      changed$X[d$m > origin,] <- changed$X[d$m > origin,] + 1000
      for (link in cfg$links) {
        fit <- fit_binary(d,rows,link); fit2 <- fit_binary(changed,mature_rows(changed,origin),link)
        near(fit$coefficients,fit2$coefficients)
        i <- which(d$m == origin)
        near(predict_binary(fit,d$X[i,,drop=FALSE]),predict_binary(fit2,changed$X[i,,drop=FALSE]))
        near(fit$center,colMeans(d$X[rows,,drop=FALSE]))
      }
      }
    })
    check(paste("Recursive validation refits and exports h",h), {
      short <- cfg; short$validation_end <- "2005-03"
      result <- expanding_forecasts(d,h,"probit",short,"validation")
      f <- result$forecasts
      ensure(nrow(f)==3L && length(result$records)==3L,"Missing validation origin fits")
      ensure(all(diff(f$n_train)==1L),"Validation history does not expand")
      ensure(all(month_number(f$last_train_origin)==month_number(f$origin)-h),"Validation includes immature targets")
      ensure(all(f$fit_information_cutoff==f$origin) && all(f$period=="validation"),"Wrong validation audit")
      # The February prediction must match an independent February information-set fit.
      feb <- which(d$ym=="2005-02")
      independent <- fit_binary(d,mature_rows(d,month_number("2005-02")),"probit")
      near(f$probability[2],predict_binary(independent,d$X[feb,,drop=FALSE]))
      ensure(!isTRUE(all.equal(result$records[[1]]$coefficients$coefficient_standardized,
                               result$records[[2]]$coefficients$coefficient_standardized)),
             "Validation parameters remained fixed")
    })
    check(paste("Vintage bill accrual, passive identities and portfolio dates h",h), {
      ids <- which(d$m >= month_number(cfg$validation_start) & d$m <= month_number(cfg$validation_end))
      f <- forecast_frame(d,ids,rep(.5,length(ids)),h,"probit","validation")
      inp <- portfolio_inputs(f,returns,h,cfg$validation_end)
      ensure(nrow(inp$active) == 96-h && ncol(inp$active) == h, "Portfolio dimensions wrong")
      ensure(inp$ym[1] == month_string(month_number(cfg$validation_start)+h), "Wrong warm-up")
      code <- c(`1`=2000061,`3`=2000062,`6`=2000063)[[as.character(h)]]
      b <- returns$bills[returns$bills$KYTREASNOX == code,]
      for (k in 1:h) {
        y <- b$TMYLD[match(month_number(f$origin[k]),b$m)]
        near(inp$bill_matrix[1,k], exp(y*inp$days[1])-1)
      }
      common <- returns$bills[returns$bills$KYTREASNOX==cfg$common_rf_security,]
      common_yield <- common$TMYLD[match(month_number(inp$origin),common$m)]
      near(inp$common_rf_return,expm1(common_yield*inp$days))
      near(inp$days,as.numeric(inp$quote_end-inp$quote_start))
      ensure(all(inp$days >= 25 & inp$days <= 35),"Accrual used full maturity rather than monthly days")
      if(h==1) near(inp$bill_reference,inp$common_rf_return)
      # Independently compound one vintage's accrual over h monthly quote gaps.
      origin <- month_number("2010-01"); mi <- match(origin, returns$market$m)
      y <- b$TMYLD[match(origin,b$m)]
      days <- as.numeric(diff(returns$market$date[mi:(mi+h)]))
      near(prod(exp(y*days))-1,d$metadata$bill_accrual_proxy[d$m==origin])
      pstock <- portfolio_path(inp,rep(1,nrow(f)),.002,cfg)
      pbill <- portfolio_path(inp,rep(0,nrow(f)),.002,cfg)
      near(pstock$net_return,inp$market); near(pbill$net_return,inp$bill_reference)
      near(pstock$cost,rep(0,nrow(pstock))); near(pbill$cost,rep(0,nrow(pbill)))
      ensure(is.na(portfolio_metrics(pbill)$net_Sharpe), "All-bills Sharpe must be undefined")
      if(h==1) ensure(is.na(portfolio_metrics(pbill)$sharpe_common_rf),"Identical common cash Sharpe must be undefined")
      if(h>1) ensure(is.finite(portfolio_metrics(pbill)$sharpe_common_rf),"Different bill comparator must not be forced to NA")
      near(portfolio_metrics(pstock)$beta,1)
    })
  }
  check("Cutoff tie rule and undefined Sharpe exclusion", {
    g <- data.frame(cutoff=c(.2,.49,.51,1),sharpe_matched_ladder=c(.1,.2,.2,NA),
                    sharpe_common_rf=c(.9,.1,.2,NA))
    ensure(choose_cutoff(g)==2L,"Cutoff tie rule failed")
    ensure(choose_cutoff(g,cfg$cutoff_objective)==2L,"Cutoff selected against wrong reference")
    ensure(choose_cutoff(g,"sharpe_common_rf")==1L,"Explicit alternative objective ignored")
    bad <- g; bad$sharpe_matched_ladder <- NA_real_
    ensure(inherits(try(choose_cutoff(bad),silent=TRUE),"try-error"),"Undefined candidates must stop selection")
  })
  check("Two reference Sharpes use matched monthly samples and sample SD", {
    p <- data.frame(month=c("2000-01","2000-02","2000-03","2000-04"),
      net_return=c(.03,-.01,.02,.04),market_return=c(.04,-.02,.03,.05),
      common_rf_return=c(.001,.002,.001,.002),matched_ladder_return=c(.002,.004,.003,.005),
      traded_fraction=0,cost_fraction=0,stock_weight=.5)
    m <- portfolio_metrics(p)
    for(ref in c("common_rf","matched_ladder")) {
      e <- p$net_return-p[[paste0(ref,"_return")]]
      expected <- sqrt(12)*sum(e)/4/sqrt(sum((e-sum(e)/4)^2)/3)
      near(m[[paste0("sharpe_",ref)]],expected)
    }
    near(m$net_Sharpe,m$sharpe_matched_ladder)
    ensure(abs(m$sharpe_common_rf-m$sharpe_matched_ladder)>1e-6,"Two references were conflated")
  })
  check("Same calendar month gives one common reference at all horizons", {
    refs <- lapply(cfg$horizons,function(h) {
      d <- load_model_data(root,h,cfg)
      ids <- which(d$m >= month_number(cfg$validation_start) & d$m <= month_number(cfg$validation_end))
      f <- forecast_frame(d,ids,rep(.5,length(ids)),h,"probit","validation")
      inp <- portfolio_inputs(f,returns,h,cfg$validation_end,cfg$common_rf_security)
      setNames(inp$common_rf_return,inp$ym)
    })
    months <- Reduce(intersect,lapply(refs,names))
    for(ref in refs[-1]) near(ref[months],refs[[1]][months])
  })
  check("AUC ties and inclusive cutoff", {
    f <- data.frame(classification_eligible=TRUE,target=c(0,1,0,1),probability=rep(.5,4))
    m <- forecast_metrics(f,.5)
    near(m$AUC,.5); near(m$Brier,.25); near(m$recall,1)
  })
  data.frame(test=passed,passed=TRUE)
}
