# Root layout: Project/Scripts/elastic_net/{this and sibling modules}
#               Project/Data/processed/model_ready_h{1,3,6}.rds
args<-commandArgs(FALSE)
entry<-grep("^--file=",args,value=TRUE)
entry_path<-if(length(entry)) gsub("~+~"," ",sub("^--file=","",entry[1]),fixed=TRUE) else ""
script_dir<-if(nzchar(entry_path)) dirname(normalizePath(entry_path,mustWork=TRUE)) else getwd()
root<-normalizePath(file.path(script_dir,"../.."),mustWork=TRUE)
.libPaths(c(file.path(root,".R-library"),.libPaths()))
needed<-c("glmnet","digest","jsonlite")
missing<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Install missing R packages: ",paste(missing,collapse=", "))
for(f in c("config.R","models.R","portfolio.R","tests.R","report.R")) source(file.path(script_dir,f))
cfg<-elastic_config(); set.seed(cfg$seed)
args<-commandArgs(TRUE)
stage<-sub("^--stage=","",grep("^--stage=",args,value=TRUE))
if(!length(stage)) stage<-"all"
ensure(length(stage)==1 && stage %in% c("all","validation","test","tests"),"Stage: all|validation|test|tests")
out<-file.path(root,"Results","elastic_net")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
inputs<-c(file.path(root,"Data/processed",paste0("model_ready_h",cfg$horizons,".rds")),
 file.path(root,"Data/crsp_monthly_index_data/crsp_monthly_index_data.csv"),
 file.path(root,"Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))
ensure(all(file.exists(inputs)),"Required input files missing: check project layout")
input_hashes<-setNames(vapply(inputs,file_hash,character(1)),substring(inputs,nchar(root)+2L))
code_files<-sort(list.files(script_dir,pattern="[.]R$",full.names=TRUE))
code_hashes<-setNames(vapply(code_files,file_hash,character(1)),basename(code_files))
message("Running numerical and data-timing tests...")
tests<-run_elastic_tests(root,cfg); write_csv(tests,file.path(out,"tests.csv"))
if(stage=="tests") quit(save="no",status=0)
returns<-load_return_data(root)
data<-setNames(lapply(cfg$horizons,function(h) load_elastic_data(root,h,cfg)),cfg$horizons)
state_path<-file.path(out,"validation_state.rds")
if(stage %in% c("all","validation")) {
  score_rows<-list(); selection_rows<-list(); forecast_rows<-list(); cutoff_rows<-list(); audit_rows<-list()
  for(h in cfg$horizons) {
    d<-data[[as.character(h)]]
    tuned<-tune_elastic(d,h,cfg)
    f<-tuned$result$forecasts
    inp<-portfolio_inputs(f,returns,h,cfg$validation_end,cfg$common_rf_security)
    cutoff_grid<-bind_rows(lapply(cfg$cutoffs,function(cutoff) {
      p<-portfolio_path(inp,as.integer(f$probability>=cutoff),cfg$primary_cost,cfg)
      cbind(data.frame(h=h,model=cfg$model,cutoff=cutoff),portfolio_metrics(p,cfg$risk_aversion))
    }))
    winner<-cutoff_grid[choose_cutoff(cutoff_grid,cfg$cutoff_objective),,drop=FALSE]
    fixed<-cutoff_grid[abs(cutoff_grid$cutoff-cfg$comparison_cutoff)<1e-12,,drop=FALSE]
    ensure(nrow(fixed)==1,"Missing 0.50 benchmark")
    sel<-cbind(tuned$selected,data.frame(cutoff=winner$cutoff,
      validation_sharpe_matched_ladder=winner$sharpe_matched_ladder,
      validation_sharpe_common_rf=winner$sharpe_common_rf,
      fixed_050_validation_sharpe_matched_ladder=fixed$sharpe_matched_ladder,
      selection_cost=cfg$primary_cost))
    f$cutoff<-sel$cutoff;f$signal<-as.integer(f$probability>=sel$cutoff)
    f$cutoff_050<-cfg$comparison_cutoff;f$signal_050<-as.integer(f$probability>=cfg$comparison_cutoff)
    key<-as.character(h)
    score_rows[[key]]<-tuned$scores;selection_rows[[key]]<-sel
    forecast_rows[[key]]<-f; cutoff_rows[[key]]<-cutoff_grid
    audit_rows[[key]]<-tuned$result$audit
  }
  scores<-bind_rows(score_rows); selected<-bind_rows(selection_rows)
  val_forecasts<-bind_rows(forecast_rows);cutoff_grid<-bind_rows(cutoff_rows)
  val_audit<-bind_rows(audit_rows)
  frozen<-list(config=cfg,selected=selected,input_hashes=input_hashes,code_hashes=code_hashes,
    frozen_at=format(Sys.time(),tz="Asia/Singapore",usetz=TRUE))
  state<-list(scores=scores,selected=selected,forecasts=val_forecasts,cutoff_grid=cutoff_grid,
    audit=val_audit,frozen=frozen)
  saveRDS(state,state_path)
  saveRDS(frozen,file.path(out,"frozen_specification.rds"))
  write_csv(scores,file.path(out,"validation_tuning.csv"))
  write_csv(selected,file.path(out,"selected_settings.csv"))
  write_csv(cutoff_grid,file.path(out,"validation_cutoff_grid.csv"))
  write_csv(val_forecasts,file.path(out,"validation_forecasts.csv"))
  message("Validation frozen: selected hyperparameters and cutoffs for all horizons")
} else {
  ensure(file.exists(state_path),"No validation state: run --stage=validation first")
  state<-readRDS(state_path)
  ensure(identical(cfg,state$frozen$config) && identical(input_hashes,state$frozen$input_hashes) &&
         identical(code_hashes,state$frozen$code_hashes),"Code/data/config changed since validation: stop")
  selected<-state$selected;val_forecasts<-state$forecasts
}
if(stage=="validation") quit(save="no",status=0)
message("Computing test predictions; selection is frozen...")
test_rows<-list(); test_audits<-list(); fit_records<-list()
for(h in cfg$horizons) {
  d<-data[[as.character(h)]]; sel<-selected[selected$h==h,,drop=FALSE]
  ensure(nrow(sel)==1,"Missing unique frozen horizon settings")
  result<-recursive_elastic(d,h,cfg,"test",sel$alpha,sel$lambda_ratio,TRUE)
  f<-result$forecasts
  f$cutoff<-sel$cutoff; f$signal<-as.integer(f$probability>=sel$cutoff)
  f$cutoff_050<-cfg$comparison_cutoff;f$signal_050<-as.integer(f$probability>=cfg$comparison_cutoff)
  key<-as.character(h); test_rows[[key]]<-f;test_audits[[key]]<-result$audit
  fit_records[[key]]<-result$fits
}
test_forecasts<-bind_rows(test_rows)
all_forecasts<-rbind(val_forecasts,test_forecasts)
audit<-rbind(state$audit,bind_rows(test_audits))
ensure(all(month_number(audit$last_train_outcome)<=month_number(audit$forecast_origin)),"Audit leakage")
ensure(nrow(state$audit)==96L*length(cfg$horizons) && nrow(bind_rows(test_audits))==156L*length(cfg$horizons),"Missing monthly refits")
paths<-list();perf<-list();diagnostics<-list();refs<-list()
for(period in c("validation","test")) for(h in cfg$horizons) {
  f<-all_forecasts[all_forecasts$period==period & all_forecasts$h==h,,drop=FALSE]
  end<-cfg[[paste0(period,"_end")]]
  inp<-portfolio_inputs(f,returns,h,end,cfg$common_rf_security)
  refs[[paste(period,h)]]<-data.frame(period=period,h=h,month=inp$ym,
    quote_start=inp$quote_start,quote_end=inp$quote_end,quote_days=inp$days,
    common_rf_return=inp$common_rf_return,matched_ladder_return=inp$bill_reference,
    common_rf_security=inp$common_rf_security,matched_bill_security=inp$matched_bill_security)
  variants<-list(validation_tuned=f$signal,fixed_050=f$signal_050,
    all_stocks=rep(1L,nrow(f)),all_bills=rep(0L,nrow(f)))
  for(rule in names(variants)) {
    strat<-if(rule %in% c("all_stocks","all_bills")) rule else cfg$model
    cutoff_rule<-if(strat==cfg$model) rule else "not_applicable"
    cutoff<-if(rule=="validation_tuned") unique(f$cutoff) else if(rule=="fixed_050") cfg$comparison_cutoff else NA_real_
    if(strat==cfg$model) diagnostics[[paste(period,h,rule)]]<-cbind(data.frame(period=period,h=h,model=cfg$model,cutoff_rule=rule),forecast_metrics(f,cutoff))
    for(cost in cfg$costs) {
      p<-portfolio_path(inp,variants[[rule]],cost,cfg)
      tags<-data.frame(period=period,h=h,strategy=strat,cutoff_rule=cutoff_rule,cutoff=cutoff,cost_rate=cost)
      key<-paste(period,h,rule,cost)
      paths[[key]]<-cbind(tags,p)
      perf[[key]]<-cbind(tags,portfolio_metrics(p,cfg$risk_aversion))
      ensure(max(month_number(p$month))<=month_number(end),"Portfolio leaks split boundary")
    }
  }
}
paths<-bind_rows(paths);performance<-bind_rows(perf)
common<-list()
for(period in c("validation","test")) {
  start<-cfg[[paste0(period,"_start")]]
  common_start<-month_string(month_number(start)+max(cfg$horizons))
  q<-paths[paths$period==period & paths$month>=common_start,]
  groups<-interaction(q$h,q$strategy,q$cutoff_rule,q$cost_rate,drop=TRUE)
  for(k in unique(groups)) {
    z<-q[groups==k,]
    common[[paste(period,k)]]<-cbind(z[1,c("period","h","strategy","cutoff_rule","cutoff","cost_rate")],portfolio_metrics(z,cfg$risk_aversion))
  }
}
ensure(identical(input_hashes,setNames(vapply(inputs,file_hash,character(1)),names(input_hashes))),"Inputs changed during run")
write_csv(test_forecasts,file.path(out,"test_forecasts.csv"))
write_csv(audit,file.path(out,"refit_audit.csv"))
write_csv(paths,file.path(out,"monthly_portfolios.csv"))
write_csv(performance,file.path(out,"performance.csv"))
write_csv(bind_rows(common),file.path(out,"common_window_performance.csv"))
write_csv(bind_rows(refs),file.path(out,"benchmark_returns.csv"))
write_csv(bind_rows(diagnostics),file.path(out,"forecast_metrics.csv"))
write_csv(data.frame(feature=colnames(data[["1"]]$X)),file.path(out,"feature_names.csv"))
write_csv(data.frame(file=names(input_hashes),sha256=unname(input_hashes)),file.path(out,"input_hashes.csv"))
saveRDS(fit_records,file.path(out,"test_fits_by_origin.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
jsonlite::write_json(list(version=cfg$version,selected=selected,input_hashes=input_hashes,
 code_hashes=code_hashes,tests_passed=nrow(tests),test_used_for_selection=FALSE,
 note="Protocol revised after earlier baseline test inspection; do not characterize as virgin holdout"),
 file.path(out,"manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
write_elastic_report(out,cfg,selected,performance,bind_rows(diagnostics),tests)
message("Complete: ",file.path(out,"RESULTS.md"))
