# Optional one-time migration comparison; the production build does not read legacy outputs.
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
scripts <- if(length(arg)) dirname(normalizePath(sub("^--file=","",arg[1]))) else normalizePath("Scripts")
source(file.path(scripts,"Build data/R/common.R"))
args <- commandArgs(TRUE)
reference <- if(length(args)) normalizePath(args[1]) else file.path(DATA,"predictor_review")
macro_reference <- if(length(args)>1) normalizePath(args[2]) else file.path(DATA,"predictor_review_audit/validated/macro")
rows <- list(); j<-0L
corrected_raw <- c("mkt_excess","bill_yld_matched","bill_yld_chg","term_spread","short_slope")
corrected_model <- c("mkt_excess","bill_yld_rel36","bill_yld_chg","term_spread_change1","short_slope_change1")
compare <- function(old,new,label,columns) {
  stopifnot(identical(old$ym,new$ym))
  for(n in columns) {
    x<-old[[n]]; y<-new[[n]]
    stopifnot(identical(is.na(x),is.na(y)))
    error<-max(abs(x-y),na.rm=TRUE)
    if(!is.finite(error)) error<-0
    stopifnot(error<1e-10)
    j <<- j+1L; rows[[j]] <<- data.table(comparison=label,variable=n,observed=sum(!is.na(x)),max_abs_error=error)
  }
}
for(h in c(1,3,6)) {
  raw_current<-read_table(paste0("base_raw_h",h,".csv"))
  prior_dispersion<-log(raw_current$dispersion/roll(lag1(raw_current$dispersion),36))
  for(pair in list(c("base_predictors_raw_h","base_raw_h"),c("base_predictors_h","base_transformed_h"))) {
    old<-fread(file.path(reference,paste0(pair[1],h,".csv")))
    new<-fread(file.path(PROCESSED,paste0(pair[2],h,".csv")))
    if(pair[2]=="base_transformed_h") {
      # Two additional R-stage transformations are deliberate specification changes.
      # Recompute the prior definitions in R to verify their underlying arithmetic.
      new[,dispersion_rel36:=prior_dispersion]; new[,tail_imb:=raw_current$tail_imb]
      new<-new[ym>="1973-01"]; old<-old[ym>="1973-01"]
    }
    compare(old,new,paste0(pair[2],h),setdiff(names(old),c("ym",corrected_raw,corrected_model)))
  }
  old<-fread(file.path(reference,paste0("features_h",h,"_120_post_nasdaq_warmup.csv")))
  new<-fread(file.path(PROCESSED,paste0("features_h",h,"_132.csv")))
  for(k in 0:5) {
    pos<-match(month_id(new$ym)-k,raw_current$m)
    new[,(paste0("dispersion_rel36_lag",k)):=prior_dispersion[pos]]
    new[,(paste0("tail_imb_lag",k)):=raw_current$tail_imb[pos]]
  }
  excluded_features <- unlist(lapply(corrected_model,function(n) paste0(n,"_lag",0:5)))
  compare(old,new,paste0("feature_h",h),setdiff(names(old),c("ym",excluded_features)))
}
macro<-fread(file.path(macro_reference,"macro_candidates.csv"))
new<-fread(file.path(PROCESSED,"base_transformed_h3.csv"))
compare(macro,new,"optional_macro",c("unrate_rel36","baa10ym_rel36"))
# Match statistics, not p-values: urca reports tabulated critical values.
legacy<-fread(file.path(dirname(macro_reference),"stationarity_summary.csv"))[kind=="model" & !predictor %in% corrected_model]
current<-fread(file.path(AUDIT,"stationarity_primary.csv"))[kind=="transformed"]
source(file.path(scripts,"Build data/R/stationarity.R"))
for(horizon in c(1,3,6)) {
  d<-read_table(paste0("base_raw_h",horizon,".csv")); keep<-d$m>=START & d$m<=TRAIN_END
  for(name in c("dispersion_rel36","tail_imb")) {
    v<-if(name=="tail_imb") d$tail_imb else log(d$dispersion/roll(lag1(d$dispersion),36))
    row<-stationarity_one(v[keep],d$m[keep]); set(row,j="h",value=horizon); set(row,j="predictor",value=name)
    current<-rbind(current,row,fill=TRUE)
  }
}
stats<-merge(legacy[,.(h,predictor,old_adf=train_adf_stat,old_kpss=train_kpss_stat,old_lags=train_adf_lags,old_bw=train_kpss_lags)],
 current[,.(h,predictor,adf_stat,kpss_stat,adf_lags,kpss_lags)],by=c("h","predictor"))
stopifnot(nrow(stats)==45,all(stats$old_lags==stats$adf_lags),all(stats$old_bw==stats$kpss_lags),
          max(abs(stats$old_adf-stats$adf_stat))<1e-8,max(abs(stats$old_kpss-stats$kpss_stat))<1e-8)
fwrite(stats,file.path(AUDIT,"stationarity_migration_comparison.csv"))
result<-rbindlist(rows); fwrite(result,file.path(AUDIT,"migration_numeric_comparison.csv"))
write_json(list(status="Unchanged R definitions agree with prior outputs; corrected yield definitions checked separately",specification=SPEC_ID,comparisons=nrow(result),numeric_cells=sum(result$observed),
 maximum_absolute_error=max(result$max_abs_error),adf_max_absolute_error=max(abs(stats$old_adf-stats$adf_stat)),
 kpss_max_absolute_error=max(abs(stats$old_kpss-stats$kpss_stat)),
 deliberate_R_changes=c("dispersion_rel36 replaced by log_dispersion_change1","tail_imb replaced by tail_imb_change1","CRSP continuous-compounding correction for five yield-based predictors and bill targets","actual quote-date day counts for market excess return"),
 corrected_yield_predictors_excluded_from_legacy_parity=corrected_model,
 comparison_note="Prior dispersion and tail definitions recomputed in R for parity; current model uses changes. Incorrect prior yield definitions are excluded from parity and independently checked against raw CRSP fields in yield_correction_manifest.json.",
 legacy_files_needed_by_normal_pipeline=FALSE),"migration_manifest.json")
print(result[,.(comparisons=.N,observations=sum(observed),max_error=max(max_abs_error))])
