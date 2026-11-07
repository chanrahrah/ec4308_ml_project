# build.R reads the CRSP stock, index, daily-return and Treasury data plus CPI, UNRATE and BAA10YM. It applies stock filters, 
#builds market and cross-sectional signals, accounts for macro availability, and saves raw and transformed series.
eligible_stock <- function(s) {
  out <- s$ShareType=="NS" & s$SecurityType=="EQTY" & s$SecuritySubType=="COM" &
    s$USIncFlg=="Y" & s$IssuerType %in% c("ACOR","CORP") & s$PrimaryExch %in% c("N","A","Q") &
    s$ConditionalType=="RW" & s$TradingStatusFlg=="A"
  out[is.na(out)] <- FALSE; out
}

stock_signals <- function(m,rx,price) {
  if(anyDuplicated(m) || any(diff(m)<=0)) stop("Unsorted/duplicate stock months")
  full <- seq.int(min(m),max(m)); pos <- match(m,full)
  r <- p <- rep(NA_real_,length(full)); r[pos] <- rx; p[pos] <- price
  ma <- hi <- lo <- rep(NA_real_,length(full))
  for(n in c(10L,12L)) {
    if(length(full)<n) next
    # Matrix rows are complete calendar windows, oldest endpoint first.
    prices <- embed(p,n)[,n:1,drop=FALSE]
    returns <- embed(r,n-1)[-1,(n-1):1,drop=FALSE]
    valid <- rowSums(!is.finite(prices) | prices<=0)==0 & rowSums(!is.finite(returns) | returns<=-1)==0
    path <- matrix(1,nrow(prices),n)
    for(j in 2:n) path[,j] <- path[,j-1]*(1+returns[,j-1])
    valid <- valid & rowSums(!is.finite(path))==0
    end <- seq.int(n,length(full))[valid]
    path <- path[valid,,drop=FALSE]
    if(n==10L) ma[end] <- as.numeric(path[,n]>rowMeans(path)) else {
      hi[end] <- as.numeric(path[,n]==do.call(pmax,as.data.frame(path)))
      lo[end] <- as.numeric(path[,n]==do.call(pmin,as.data.frame(path)))
    }
  }
  list(ma=ma[pos],hi=hi[pos],lo=lo[pos])
}

build_stock <- function(calendar) {
  message("Reading CRSP stock source and validating keys/columns...")
  s <- fread(file.path(DATA,INPUTS[['stock']]),select=c("PERMNO","MthCalDt",FLAGS,NUMERIC),na.strings=c("","NA"))
  # Volume can exceed 2^31; normalize integer64 to double before mixed-column audits.
  s[,(NUMERIC):=lapply(.SD,as.numeric),.SDcols=NUMERIC]
  s[,m:=month_id(MthCalDt)]
  if(anyNA(s$PERMNO) || anyNA(s$m)) stop("Missing stock keys")
  assert_numeric(s[,..NUMERIC]); if(any(as.matrix(s[,.(MthRet,MthRetx)]) < -1,na.rm=TRUE)) stop("Return below -100%")
  s[,strict:=eligible_stock(.SD)]
  duplicates <- s[duplicated(s,by=c("PERMNO","m")) | duplicated(s,by=c("PERMNO","m"),fromLast=TRUE)]
  conflicts <- duplicates[,lapply(.SD,uniqueN,na.rm=FALSE),by=.(PERMNO,m),.SDcols=c(NUMERIC,FLAGS)]
  numeric_conflict <- conflicts[,rowSums(as.matrix(.SD)>1)>0,.SDcols=NUMERIC]
  flag_conflict <- conflicts[,rowSums(as.matrix(.SD)>1)>0,.SDcols=FLAGS]
  if(any(numeric_conflict)) stop("Ambiguous numeric duplicate keys")
  if(any(flag_conflict)) save_table(duplicates[conflicts[flag_conflict,.(PERMNO,m)],on=.(PERMNO,m)],"stock_classification_conflicts.csv",AUDIT)
  counts <- list(raw_rows=nrow(s),duplicate_keys=nrow(conflicts),numeric_conflict_keys=sum(numeric_conflict),classification_conflict_keys=sum(flag_conflict))
  column_audit <- rbindlist(lapply(NUMERIC,function(n) { x<-s[[n]]; data.table(column=n,observed=sum(!is.na(x)),missing=sum(is.na(x)),negative=sum(x<0,na.rm=TRUE),zero=sum(x==0,na.rm=TRUE),minimum=min(x,na.rm=TRUE),maximum=max(x,na.rm=TRUE)) }))
  fwrite(column_audit,file.path(AUDIT,"crsp_numeric_columns.csv"))
  raw_cap <- s[MthPrc>0 & ShrOut>0 & !is.na(MthCap)]
  raw_bad <- raw_cap[abs(MthCap/(MthPrc*ShrOut)-1)>.01]
  counts$raw_cap_identity_mismatches <- nrow(raw_bad)
  setorderv(s,c("PERMNO","m","strict"),c(1,1,-1))
  s <- unique(s,by=c("PERMNO","m"))
  counts$unique_rows <- nrow(s)
  s[,terminal:=TradingStatusFlg=="D" & !is.na(MthDelFlg) & MthDelFlg!="N" & !is.na(MthRet) &
      (m-shift(m)==1L) & shift(strict,fill=FALSE) & !strict,by=PERMNO]
  s[is.na(terminal),terminal:=FALSE]; s[,recommended:=strict | terminal]
  counts$terminal_recovered <- sum(s$terminal)
  bad <- s[MthPrc>0 & ShrOut>0 & !is.na(MthCap) & abs(MthCap/(MthPrc*ShrOut)-1)>.01]
  bad[,cap_price_shares_ratio:=MthCap/(MthPrc*ShrOut)]
  save_table(bad,"stock_capitalization_mismatches.csv",AUDIT)
  counts$unique_cap_identity_mismatches <- nrow(bad)
  counts$cap_mismatches_in_turnover_1973_on <- nrow(bad[strict & PrimaryExch %in% c("N","A") & MthVol>=0 & ShrOut>0 & m>=month_id("1973-01")])
  message("Computing stock windows in R (all securities, complete monthly calendars)...")
  s[,c("ma","hi","lo"):=stock_signals(m,MthRetx,MthPrc),by=PERMNO]
  e <- s[recommended & !is.na(MthRet)]
  cross <- e[,{
    bounds <- quantile(MthRet,c(.01,.99),type=7,names=FALSE)
    movers <- sum(MthRet!=0)
    list(breadth=if(movers) sum(MthRet>0)/movers else NA_real_,ma_part=mean(ma,na.rm=TRUE),
         hi_lo=mean(hi,na.rm=TRUE)-mean(lo,na.rm=TRUE),dispersion=sd(pmin(pmax(MthRet,bounds[1]),bounds[2])),
         tail_imb=mean(MthRet < -0.2)-mean(MthRet > 0.2),n_return=.N,n_ma=sum(!is.na(ma)),n_high_low=sum(!is.na(hi)))
  },by=m]
  act <- s[strict & PrimaryExch %in% c("N","A") & MthVol>=0 & ShrOut>0]
  turnover <- act[,.(volume=sum(MthVol),shares=sum(ShrOut)),by=m]
  turnover <- merge(data.table(m=calendar),turnover,by="m",all.x=TRUE,sort=TRUE)
  turnover[,turnover_fraction:=volume/(1000*shares)]
  turnover[,turnover_dt:=log(turnover_fraction/roll(lag1(turnover_fraction),12))]
  save_table(turnover,"turnover_audit.csv",AUDIT)
  save_table(cross[,.(m,n_return,n_ma,n_high_low)],"universe_counts.csv",AUDIT)
  write_json(counts,"stock_audit.json")
  merge(cross[,.(m,breadth,ma_part,hi_lo,dispersion,tail_imb)],turnover[,.(m,turnover_dt)],by="m",all=TRUE,sort=TRUE)
}

available_macro <- function(m,x,calendar,overrides=integer()) {
  assert_monthly(m,"macro source"); assert_numeric(x)
  release <- m+1L
  if(length(overrides)) for(n in names(overrides)) release[m==month_id(n)] <- overrides[[n]]
  value <- age <- source <- rep(NA_real_,length(calendar))
  for(i in seq_along(calendar)) {
    known <- which(!is.na(x) & release<=calendar[i])
    if(length(known)) { j <- tail(known,1); value[i]<-x[j]; source[i]<-m[j]; age[i]<-calendar[i]-m[j] }
  }
  data.table(m=calendar,value=value,source_month=ifelse(is.na(source),NA_character_,month_label(source)),source_age_months=age,
             carried_forward=!is.na(value) & !is.na(lag1(source)) & source==lag1(source))
}

distribution_proxy <- function(r,rx) {
  if(any(rx<=-1,na.rm=TRUE) || any(r < -1,na.rm=TRUE)) stop("Invalid index returns")
  out <- rep(NA_real_,length(r))
  if(length(r)>=12) for(i in 12:length(r)) {
    rows <- (i-11):i
    if(any(!is.finite(c(r[rows],rx[rows])))) next
    denominator <- rev(cumprod(rev(1+rx[rows])))
    value <- sum((r[rows]-rx[rows])/denominator)
    if(value>0) out[i] <- log(value)
  }
  out
}

build_all <- function() {
  before <- source_hashes()
  idx <- read_source("monthly","mthcaldt"); assert_monthly(idx$m,"monthly index")
  assert_numeric(idx[,.(vwretd,vwretx,ewretd)],missing=FALSE)
  if(any(as.matrix(idx[,.(vwretd,vwretx,ewretd)])<=-1)) stop("Nonpositive index growth")
  daily <- read_source("daily","dlycaldt")
  if(anyNA(daily$dlycaldt) || anyDuplicated(daily$dlycaldt)) stop("Bad daily dates")
  assert_numeric(daily$vwretd,missing=FALSE)
  if(any(daily$vwretd<=-1)) stop("Invalid daily returns")
  fixed <- read_source("fixed","MCALDT")[KYTREASNOX==2000007]
  assert_monthly(fixed$m,"fixed term")
  calendar <- seq.int(min(idx$m),min(max(idx$m),max(daily$m),max(fixed$m)))
  idx <- idx[match(calendar,m)]
  bills <- read_source("bills","MCALDT")[!is.na(m) & KYTREASNOX %in% c(2000061,2000062,2000063)]
  if(anyDuplicated(bills,by=c("KYTREASNOX","m"))) stop("Duplicate Treasury keys")
  by <- lapply(c(2000061,2000062,2000063),function(code) bills[KYTREASNOX==code][match(calendar,m)]$TMYLD)
  names(by) <- c("1","3","6"); annual <- lapply(by,daily_annual)
  fy <- fixed[match(calendar,m)]$TMYTM
  quote_days <- c(NA_real_,diff(as.numeric(as.Date(idx$mthcaldt))))
  base <- data.table(m=calendar,mkt_excess=idx$vwretd-daily_accrual(lag1(by[['1']]),quote_days))
  for(n in c(3,6,12)) base[,(paste0("mom",n)):=expm1(roll(log1p(idx$vwretd),n,sum))]
  risk <- daily[,.(n=.N,sq=sum(vwretd^2),down=sum(vwretd[vwretd<0]^2)),by=m][match(calendar,m)]
  if(anyNA(risk$n) || any(risk$n<15) || any(risk$sq<=0)) stop("Incomplete daily-risk coverage")
  save_table(risk,"daily_risk_coverage.csv",AUDIT)
  base[,`:=`(log_rv=.5*log(252*risk$sq/risk$n),downside_share=risk$down/risk$sq)]
  wealth <- exp(cumsum(log1p(idx$vwretd)))
  base[,`:=`(drawdown12=wealth/roll(wealth,12,max)-1,log_dp_proxy=distribution_proxy(idx$vwretd,idx$vwretx),
    bill_yld_chg=annual[['3']]-lag1(annual[['3']]),term_spread=daily_annual(fy/36500)-annual[['3']],
    short_slope=annual[['6']]-annual[['1']],ew_vw=idx$ewretd-idx$vwretd)]
  cpi <- read_source("cpi","observation_date"); assert_monthly(cpi$m,"CPI")
  if(any(cpi$CPIAUCNS<=0,na.rm=TRUE)) stop("Invalid CPI")
  info <- available_macro(cpi$m,log(cpi$CPIAUCNS/lag1(cpi$CPIAUCNS,12)),calendar)
  base[,cpi_infl:=info$value]; save_table(info,"cpi_availability.csv",AUDIT)
  for(key in c("unrate","baa")) {
    name <- if(key=="unrate") "UNRATE" else "BAA10YM"
    d <- read_source(key,"observation_date"); assert_monthly(d$m,name)
    if(key=="unrate" && any(d[[name]]<0 | d[[name]]>100,na.rm=TRUE)) stop("Invalid UNRATE percent")
    override <- if(key=="unrate") setNames(month_id("2025-11"),"2025-09") else integer()
    info <- available_macro(d$m,d[[name]]/100,calendar,override)
    base[,(if(key=="unrate") "unrate" else "baa10ym"):=info$value]
    save_table(info,paste0(name,"_availability.csv"),AUDIT)
  }
  cross <- build_stock(calendar)
  base <- merge(base,cross,by="m",all.x=TRUE,sort=TRUE)
  for(h in c(1,3,6)) {
    raw <- copy(base); raw[,bill_yld_matched:=annual[[as.character(h)]]]; raw <- raw[,c("m",RAW_NAMES),with=FALSE]
    model <- transform_predictors(raw)
    features <- lag_features(model)[m>=START & m<=TEST_END]
    assert_numeric(features[,..FEATURE_NAMES],missing=FALSE)
    if(ncol(features)!=133L) stop("Expected 132 features plus month")
    save_table(raw[m>=month_id("1970-01")],paste0("base_raw_h",h,".csv"))
    save_table(model[m>=month_id("1970-01")],paste0("base_transformed_h",h,".csv"))
    save_table(features,paste0("features_h",h,"_132.csv"))
    saveRDS(list(ym=month_label(features$m),X=as.matrix(features[,..FEATURE_NAMES]),specification=SPEC_ID),file.path(PROCESSED,paste0("features_h",h,"_132.rds")))
  }
  stopifnot(identical(before,source_hashes()))
  write_json(list(specification=SPEC_ID,input_paths=as.list(INPUTS),input_sha256=as.list(before),R=R.version.string,
    packages=as.list(vapply(REQUIRED,function(p) as.character(packageVersion(p)),character(1))),
    transforms=SPEC,training=c("1974-04","2004-12"),n_predictors=22,n_features=132,source_data_modified=FALSE),"build_manifest.json")
  message("R construction complete: 22 predictors / 132 features at each horizon.")
}
