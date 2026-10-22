# model_data.R constructs forward-return labels, purges split-boundary observations, and writes CSV and RDS model inputs.
# Prepare proposal-aligned targets; no models, scalers, thresholds or performance fitted.
forward_targets <- function(m,returns,quote_dates,yield,h) {
  assert_monthly(m); assert_numeric(returns,missing=FALSE); assert_numeric(yield)
  if(!h %in% c(1L,3L,6L) || any(returns<=-1)) stop("Invalid horizon or returns")
  date <- as.Date(quote_dates)
  if(anyNA(date) || any(diff(date)<=0) || any(month_id(date)!=m)) stop("Invalid quote dates")
  out <- data.table(m=m,outcome_end=rep(NA_character_,length(m)),equity_forward=NA_real_,bill_accrual_proxy=NA_real_,target=NA_integer_)
  for(i in seq_len(length(m)-h)) {
    if(!is.finite(yield[i]) || yield[i]<=-1) next
    eq <- prod(1+returns[(i+1):(i+h)])-1
    bill <- daily_accrual(yield[i],as.numeric(date[i+h]-date[i]))
    end_month <- month_label(m[i+h])
    out[i,`:=`(outcome_end=end_month,equity_forward=eq,bill_accrual_proxy=bill,target=as.integer(eq>bill))]
  }
  out
}
model_split <- function(m,outcome_end,target) {
  result <- rep("unlabeled",length(m)); end <- month_id(outcome_end)
  known <- !is.na(target)
  result[known] <- "purged_boundary"
  result[known & m>=START & m<=TRAIN_END & end<=TRAIN_END] <- "train"
  result[known & m>TRAIN_END & m<=VAL_END & end<=VAL_END] <- "validation"
  result[known & m>VAL_END & m<=TEST_END & end<=TEST_END] <- "test"
  result
}
prepare_model_data <- function() {
  idx <- read_source("monthly","mthcaldt")
  bills <- read_source("bills","MCALDT")[!is.na(m)]
  counts <- list()
  for(h in c(1L,3L,6L)) {
    code <- c(`1`=2000061,`3`=2000062,`6`=2000063)[[as.character(h)]]
    yield <- bills[KYTREASNOX==code][match(idx$m,m)]$TMYLD
    targets <- forward_targets(idx$m,idx$vwretd,idx$mthcaldt,yield,h)
    x <- read_table(paste0("features_h",h,"_132.csv"))
    full <- merge(x,targets,by="m",all.x=TRUE,sort=TRUE)
    full[,split:=model_split(m,outcome_end,target)]
    # Mature lagged target for dynamic probit: origin t-h, whose outcome ends at t.
    full[,mature_target_lag:=targets$target[match(m-h,targets$m)]]
    save_table(full,paste0("model_data_h",h,".csv"))
    X <- as.matrix(full[,..FEATURE_NAMES]); storage.mode(X) <- "double"
    rows <- lapply(c("train","validation","test"),function(s) which(full$split==s)); names(rows)<-c("train","validation","test")
    object <- list(X=X,y=full$target,ym=month_label(full$m),split=full$split,rows=rows,
      mature_target_lag=full$mature_target_lag,feature_names=FEATURE_NAMES,
      metadata=full[,.(ym=month_label(m),outcome_end,equity_forward,bill_accrual_proxy)],
      specification=SPEC_ID,target_definition="1{compound equity return t+1..t+h > exp(TMYLD_origin * quote-date days)-1}; bill is a proxy")
    saveRDS(object,file.path(PROCESSED,paste0("model_ready_h",h,".rds")))
    counts[[as.character(h)]] <- full[,.(h=h,n=.N,first=month_label(min(m)),last=month_label(max(m))),by=split]
  }
  fwrite(rbindlist(counts),file.path(AUDIT,"model_split_counts.csv"))
}
