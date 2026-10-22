# Shared, pure-R data contracts and frozen research specification.
# common.R defines the input files, CRSP field contracts, monthly calendar helpers, transformations, and lag names.
find_root <- function() {
  override <- Sys.getenv("EC4308_ROOT", "")
  if (nzchar(override)) return(normalizePath(override, mustWork=TRUE))
  arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  file <- if(length(arg)) sub("^--file=", "", arg[1]) else ""
  candidates <- c(if(nzchar(file) && file.exists(file)) dirname(normalizePath(file)), getwd())
  for (candidate in candidates) {
    for (i in 1:5) {
      if (dir.exists(file.path(candidate,"Data")) && dir.exists(file.path(candidate,"Scripts"))) return(candidate)
      candidate <- dirname(candidate)
    }
  }
  stop("Cannot locate project; set EC4308_ROOT.")
}
ROOT <- find_root()
.libPaths(c(file.path(ROOT,".R-library"), .libPaths()))
REQUIRED <- c("data.table","jsonlite","digest","urca")
missing_packages <- REQUIRED[!vapply(REQUIRED, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_packages)) stop("Run Rscript Scripts/00_setup.R; missing: ",paste(missing_packages,collapse=", "))
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
DATA <- file.path(ROOT,"Data")
PROCESSED <- file.path(DATA,"processed")
AUDIT <- file.path(DATA,"validation")
dir.create(PROCESSED, recursive=TRUE, showWarnings=FALSE)
dir.create(AUDIT, recursive=TRUE, showWarnings=FALSE)
INPUTS <- c(stock="crsp_monthly_stock_data/crsp_monthly_stock_data.csv",
 monthly="crsp_monthly_index_data/crsp_monthly_index_data.csv",
 daily="crsp_daily_index_data/crsp_daily_index_data.csv",
 bills="crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv",
 fixed="crsp_treasury_fixed_term_index_data/crsp_treasury_fixed_term_index_data.csv",
 cpi="macro/CPIAUCNS.csv",unrate="macro/UNRATE.csv",baa="macro/BAA10YM.csv")
FLAGS <- c("PrimaryExch","ConditionalType","TradingStatusFlg","USIncFlg","IssuerType",
           "SecurityType","SecuritySubType","ShareType","MthDelFlg")
NUMERIC <- c("MthRet","MthRetx","MthPrc","MthVol","ShrOut","MthCap")
RAW_NAMES <- c("mkt_excess","mom3","mom6","mom12","log_rv","downside_share","drawdown12",
 "log_dp_proxy","bill_yld_matched","bill_yld_chg","term_spread","short_slope","cpi_infl",
 "breadth","ma_part","hi_lo","dispersion","ew_vw","tail_imb","turnover_dt","unrate","baa10ym")
SPEC <- data.table(raw=c("log_dp_proxy","bill_yld_matched","term_spread","short_slope","cpi_infl","dispersion","tail_imb","unrate","baa10ym"),
 model=c("log_dp_change1","bill_yld_rel36","term_spread_change1","short_slope_change1","cpi_infl_change1","log_dispersion_change1","tail_imb_change1","unrate_rel36","baa10ym_rel36"),
 method=c("difference","demean","difference","difference","difference","log_difference","difference","demean","demean"))
MODEL_NAMES <- RAW_NAMES
MODEL_NAMES[match(SPEC$raw, RAW_NAMES)] <- SPEC$model
FEATURE_NAMES <- unlist(lapply(MODEL_NAMES, function(n) paste0(n,"_lag",0:5)),use.names=FALSE)
SPEC_ID <- "r_22_predictors_continuous_crsp_yields_v3"
month_id <- function(x) {
  z <- as.character(x)
  suppressWarnings(as.integer(substr(z,1,4))*12L + as.integer(substr(z,6,7))-1L)
}
month_label <- function(x) sprintf("%04d-%02d", x %/% 12L, x %% 12L+1L)
START <- month_id("1974-04")
TRAIN_END <- month_id("2004-12")
VAL_END <- month_id("2012-12")
TEST_END <- month_id("2025-12")
assert_monthly <- function(m, label="series") {
  if (!length(m) || anyNA(m) || anyDuplicated(m) || any(diff(m)!=1L)) stop(label,": require a complete, sorted, unique monthly calendar")
  invisible(TRUE)
}
assert_numeric <- function(x, label="values", missing=TRUE) {
  v <- as.matrix(x)
  if (!is.numeric(v) || any(is.infinite(v)) || (!missing && anyNA(v))) stop(label,": invalid numeric values")
  invisible(TRUE)
}
lag1 <- function(x,k=1L) shift(x,k,type="lag")
roll <- function(x,n,fun=mean) {
  out <- rep(NA_real_,length(x))
  if (length(x)>=n) for (i in n:length(x)) {
    w <- x[(i-n+1):i]
    if (all(is.finite(w))) out[i] <- fun(w)
  }
  out
}
carry <- function(x) {
  good <- which(!is.na(x)); out <- x
  if(length(good)) {
    pos <- findInterval(seq_along(x),good)
    out[pos>0] <- x[good[pos[pos>0]]]
  }
  out
}
read_source <- function(key,date) {
  d <- fread(file.path(DATA,INPUTS[[key]]),na.strings=c("","NA","."))
  d[,m:=month_id(get(date))]
  setorderv(d,c("m",date),na.last=FALSE)
  d
}
save_table <- function(d,name,folder=PROCESSED) {
  out <- copy(d)
  if ("m" %in% names(out)) { out[,ym:=month_label(m)]; out[,m:=NULL]; setcolorder(out,c("ym",setdiff(names(out),"ym"))) }
  fwrite(out,file.path(folder,name),na="",quote="auto")
}
read_table <- function(name,folder=PROCESSED) {
  d <- fread(file.path(folder,name)); d[,m:=month_id(ym)]; d[,ym:=NULL]; setcolorder(d,c("m",setdiff(names(d),"m"))); d
}
write_json <- function(x,name,folder=AUDIT) jsonlite::write_json(x,file.path(folder,name),pretty=TRUE,auto_unbox=TRUE,na="null",digits=16)
sha <- function(p) digest::digest(file=p,algo="sha256",serialize=FALSE)
source_hashes <- function() setNames(vapply(file.path(DATA,INPUTS),sha,character(1)),names(INPUTS))
# CRSP TMPCYLD = 2 * (exp(182.5 * TMYLD) - 1): TMYLD is a daily log yield.
# TMYTM = 36500 * TMYLD. Do not apply log1p to either converted yield.
daily_accrual <- function(x,days) {
  assert_numeric(x); assert_numeric(days)
  if(any(x<=-1,na.rm=TRUE) || any(days<0,na.rm=TRUE)) stop("Invalid CRSP yield/sentinel or day count")
  if(length(days)!=1L && length(days)!=length(x)) stop("Yield/day-count lengths differ")
  out <- expm1(days*x); assert_numeric(out); out
}
daily_annual <- function(x) daily_accrual(x,365)
transform_predictors <- function(raw) {
  assert_monthly(raw$m); out <- copy(raw)
  for(i in seq_len(nrow(SPEC))) {
    x <- raw[[SPEC$raw[i]]]; mu <- roll(lag1(x),36)
    value <- switch(SPEC$method[i],difference=x-lag1(x),demean=x-mu,
      log_difference={if(any(x<=0,na.rm=TRUE)) stop("Nonpositive dispersion"); log(x)-lag1(log(x))})
    out[,(SPEC$raw[i]):=value]; setnames(out,SPEC$raw[i],SPEC$model[i])
  }
  setcolorder(out,c("m",MODEL_NAMES)); out
}
lag_features <- function(model) {
  assert_monthly(model$m); out <- data.table(m=model$m)
  for(n in MODEL_NAMES) for(k in 0:5) out[,(paste0(n,"_lag",k)):=lag1(model[[n]],k)]
  out
}
