#validate.R checks transformations, lags, labels, splits and source hashes.
validate_all <- function() {
  checks <- list(); j<-0L
  index <- read_source("monthly","mthcaldt")
  bills <- read_source("bills","MCALDT")[!is.na(m)]
  for(h in c(1L,3L,6L)) {
    raw <- read_table(paste0("base_raw_h",h,".csv"))
    model <- read_table(paste0("base_transformed_h",h,".csv"))
    x <- read_table(paste0("features_h",h,"_132.csv"))
    assert_monthly(raw$m); assert_monthly(model$m); assert_monthly(x$m)
    stopifnot(identical(names(x),c("m",FEATURE_NAMES)),nrow(x)==621,all(is.finite(as.matrix(x[,..FEATURE_NAMES]))))
    for(n in MODEL_NAMES) for(k in 0:5) {
      expected <- model[[n]][match(x$m-k,model$m)]
      actual <- x[[paste0(n,"_lag",k)]]
      stopifnot(max(abs(actual-expected))<1e-12)
    }
    # Independent scalar construction of every transformation in the modeling sample.
    for(i in seq_len(nrow(SPEC))) {
      v <- raw[[SPEC$raw[i]]]; expected <- rep(NA_real_,length(v))
      for(t in 37:length(v)) {
        if(SPEC$method[i]=="difference") expected[t]<-v[t]-v[t-1] else if(SPEC$method[i]=="log_difference") {
          expected[t]<-log(v[t]/v[t-1])
        } else {
          mean_past<-sum(v[(t-36):(t-1)])/36
          expected[t]<-v[t]-mean_past
        }
      }
      keep<-raw$m>=START-5L
      stopifnot(max(abs(expected[keep]-model[[SPEC$model[i]]][keep]),na.rm=TRUE)<1e-11)
    }
    rds <- readRDS(file.path(PROCESSED,paste0("model_ready_h",h,".rds")))
    stopifnot(identical(colnames(rds$X),FEATURE_NAMES),nrow(rds$X)==621,all(is.finite(rds$X)),
      identical(rds$specification,SPEC_ID),identical(rds$ym,month_label(x$m)),
      max(abs(rds$X-as.matrix(x[,..FEATURE_NAMES])))<1e-12)
    for(s in c("train","validation","test")) {
      idx<-rds$rows[[s]]; stopifnot(length(idx)>0,all(!is.na(rds$y[idx])))
      boundary<-switch(s,train=TRAIN_END,validation=VAL_END,test=TEST_END)
      stopifnot(identical(idx,which(rds$split==s)),all(month_id(rds$metadata$outcome_end[idx])<=boundary))
    }
    # Independently reconstruct labels using log sums and exp rather than build helpers.
    code <- c(`1`=2000061,`3`=2000062,`6`=2000063)[[as.character(h)]]
    yield <- bills[KYTREASNOX==code][match(index$m,m)]$TMYLD
    for(row in which(!is.na(rds$y))) {
      origin <- x$m[row]; pos <- match(origin,index$m)
      eq <- expm1(sum(log1p(index$vwretd[(pos+1L):(pos+h)])))
      days <- as.numeric(as.Date(index$mthcaldt[pos+h])-as.Date(index$mthcaldt[pos]))
      safe <- exp(yield[pos]*days)-1
      stopifnot(abs(eq-rds$metadata$equity_forward[row])<1e-12,
        abs(safe-rds$metadata$bill_accrual_proxy[row])<1e-12,
        rds$y[row]==as.integer(eq>safe),
        rds$metadata$outcome_end[row]==month_label(origin+h))
    }
    # The dynamic-probit lag must have matured at the feature origin, never after it.
    for(row in seq_len(nrow(x))) {
      pos <- match(x$m[row],index$m)
      past_origin <- pos-h
      past_eq <- prod(1+index$vwretd[(past_origin+1L):pos])-1
      past_days <- as.numeric(as.Date(index$mthcaldt[pos])-as.Date(index$mthcaldt[past_origin]))
      past_bill <- exp(past_days*yield[past_origin])-1
      stopifnot(rds$mature_target_lag[row]==as.integer(past_eq>past_bill))
    }
    j<-j+1L; checks[[j]]<-data.table(h=h,rows=nrow(x),features=length(FEATURE_NAMES),lag_cells_checked=nrow(x)*length(FEATURE_NAMES),missing_feature_cells=sum(is.na(x)))
  }
  fwrite(rbindlist(checks),file.path(AUDIT,"dataset_checks.csv"))
  audit_yield_convention()
  manifest <- jsonlite::read_json(file.path(AUDIT,"build_manifest.json"),simplifyVector=TRUE)
  stopifnot(identical(unname(source_hashes()),unname(unlist(manifest$input_sha256))))
  primary <- fread(file.path(AUDIT,"stationarity_primary.csv"))
  for(horizon in c(1,3,6)) {
    current <- primary[h==horizon & kind=="transformed"]
    stopifnot(nrow(current)==22,identical(current$predictor,MODEL_NAMES),
              all(current$start=="1974-04"),all(current$end=="2004-12"),all(current$n==369),
              all(current$verdict=="stationarity supported"))
    lagged <- fread(file.path(AUDIT,"stationarity_lagged_primary.csv"))[h==horizon]
    stopifnot(nrow(lagged)==132,identical(lagged$feature,FEATURE_NAMES),
      all(lagged$n==369),all(lagged$start=="1974-04"),all(lagged$end=="2004-12"),
      all(lagged$verdict=="stationarity supported"))
  }
  current_scripts<-list.files(file.path(ROOT,"Scripts"),pattern="\\.R$",full.names=TRUE,recursive=TRUE)
  write_json(list(status="R feature and modeling-input contracts passed",specification=SPEC_ID,
    datasets=checks,raw_data_unchanged=TRUE,scripts_sha256=as.list(setNames(vapply(current_scripts,sha,character(1)),substring(current_scripts,nchar(ROOT)+2)))),"validation_manifest.json")
  message("R validation passed: calendar lags, transformations, labels, split boundaries, source hashes.")
}

audit_yield_convention <- function() {
  idx <- read_source("monthly","mthcaldt")
  b <- read_source("bills","MCALDT")[!is.na(m)]
  fixed <- read_source("fixed","MCALDT")[KYTREASNOX==2000007]
  days <- c(NA_real_,diff(as.numeric(as.Date(idx$mthcaldt))))
  calendar_days <- as.numeric(as.Date(paste0(month_label(idx$m+1L),"-01"))-as.Date(paste0(month_label(idx$m),"-01")))
  yields <- lapply(c(2000061,2000062,2000063),function(code) b[KYTREASNOX==code][match(idx$m,m)]$TMYLD)
  names(yields) <- c("1","3","6")
  # Independent plain-exp conversion: CRSP TMPCYLD identity implies continuous yields.
  annual <- lapply(yields,function(y) exp(365*y)-1)
  prior_annual <- lapply(yields,function(y) (1+y)^365-1)
  long_yield <- fixed[match(idx$m,m)]$TMYTM
  correct_long <- exp(long_yield/100)-1
  prior_long <- (1+long_yield/36500)^365-1
  changes <- labels <- list()
  for(horizon in c(1,3,6)) {
    hkey <- as.character(horizon)
    raw <- read_table(paste0("base_raw_h",horizon,".csv")); pos<-match(raw$m,idx$m)
    expected <- list(mkt_excess=idx$vwretd-(exp(lag1(yields[['1']])*days)-1),
      bill_yld_matched=annual[[hkey]],bill_yld_chg=annual[['3']]-lag1(annual[['3']]),
      term_spread=correct_long-annual[['3']],short_slope=annual[['6']]-annual[['1']])
    prior <- list(mkt_excess=idx$vwretd-((1+lag1(yields[['1']]))^calendar_days-1),
      bill_yld_matched=prior_annual[[hkey]],bill_yld_chg=prior_annual[['3']]-lag1(prior_annual[['3']]),
      term_spread=prior_long-prior_annual[['3']],short_slope=prior_annual[['6']]-prior_annual[['1']])
    for(name in names(expected)) {
      correct <- expected[[name]][pos]; previous <- prior[[name]][pos]
      stopifnot(identical(is.na(raw[[name]]),is.na(correct)),max(abs(raw[[name]]-correct),na.rm=TRUE)<1e-12)
      changes[[length(changes)+1L]] <- data.table(h=horizon,predictor=name,observations=sum(is.finite(correct)),
        independent_max_error=max(abs(raw[[name]]-correct),na.rm=TRUE),maximum_absolute_correction=max(abs(correct-previous),na.rm=TRUE))
    }
    ready <- readRDS(file.path(PROCESSED,paste0("model_ready_h",horizon,".rds")))
    rows<-which(!is.na(ready$y)); origin<-match(month_id(ready$ym[rows]),idx$m)
    holding_days<-as.numeric(as.Date(idx$mthcaldt[origin+horizon])-as.Date(idx$mthcaldt[origin]))
    y<-yields[[hkey]][origin]
    correct_bill<-exp(y*holding_days)-1; old_bill<-(1+y)^holding_days-1
    stopifnot(max(abs(correct_bill-ready$metadata$bill_accrual_proxy[rows]))<1e-12)
    old_label<-as.integer(ready$metadata$equity_forward[rows]>old_bill)
    labels[[length(labels)+1L]]<-data.table(h=horizon,completed_labels=length(rows),
      changed_labels=sum(ready$y[rows]!=old_label),maximum_bill_return_correction=max(abs(correct_bill-old_bill)))
    # Quote dates coincide in the model sample; the 1968 source exception is outside it.
    quote<-b[KYTREASNOX==c(`1`=2000061,`3`=2000062,`6`=2000063)[[hkey]]][match(idx$m[origin],m)]$MCALDT
    stopifnot(all(as.Date(quote)==as.Date(idx$mthcaldt[origin])))
  }
  fwrite(rbindlist(changes),file.path(AUDIT,"yield_feature_correction.csv"))
  fwrite(rbindlist(labels),file.path(AUDIT,"yield_target_correction.csv"))
  write_json(list(specification=SPEC_ID,correct_daily_accrual="exp(TMYLD * days) - 1",
    correct_effective_annual="exp(365 * TMYLD) - 1; equivalently exp(TMYTM/100) - 1",
    CRSP_identity="TMPCYLD = 2 * (exp(182.5 * TMYLD) - 1)",
    source="https://www.crsp.org/crsp_pdf/crsp-us-treasury-database-guide-sas-ascii-excel-r/",
    source_page=13,market_excess_day_count="actual index quote-date gaps; replaces calendar-month day count",
    corrected_predictors=names(expected),label_impact=labels,
    note="Inherited Python discrete compounding was incorrect. Historical parity is retained only for unchanged definitions."),"yield_correction_manifest.json")
}
