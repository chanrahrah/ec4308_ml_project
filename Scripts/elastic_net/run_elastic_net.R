# EC4308 logistic Elastic Net: single entry point, validation + test by default.
# Place this folder in <project>/Scripts/elastic_net_logit/ and click Source.
# Alternatively: Rscript --vanilla Scripts/elastic_net_logit/run_elastic_net.R
# Data/ is located automatically by walking up from this script's directory.
# Optional override: replace NULL with the absolute project directory below.
PROJECT_ROOT <- NULL

# Resolve this file in both RStudio Source and Rscript, without changing getwd().
.elastic_entry <- local({
  ofiles <- Filter(function(z) is.character(z) && length(z)==1L && nzchar(z),
                   lapply(sys.frames(),function(f) f$ofile))
  cli <- grep("^--file=",commandArgs(FALSE),value=TRUE)
  path <- if(length(ofiles)) tail(ofiles,1)[[1]] else if(length(cli))
    gsub("~+~"," ",sub("^--file=","",cli[1]),fixed=TRUE) else ""
  if(!nzchar(path)) stop("Open run_elastic_net.R and click Source, or use Rscript.",call.=FALSE)
  normalizePath(path,mustWork=TRUE)
})

local({
script_dir <- dirname(.elastic_entry)
module_files <- c("config.R","models.R","portfolio.R","tests.R","report.R")
for(name in module_files) {
  path <- file.path(script_dir,name)
  if(!file.exists(path)) stop("Missing module: ",path,call.=FALSE)
  sys.source(path,envir=environment())
}
find_project_root <- function(start) {
  current <- normalizePath(start,mustWork=TRUE)
  repeat {
    if(dir.exists(file.path(current,"Data","processed"))) return(current)
    parent <- dirname(current)
    if(identical(parent,current)) break
    current <- parent
  }
  stop(paste("Cannot locate Data/processed above this script.",
    "Place the folder in <project>/Scripts/elastic_net_logit/",
    "or set PROJECT_ROOT at the top of run_elastic_net.R."),call.=FALSE)
}

elastic_net_main <- function(root, stage="all", script_dir) {
  root <- normalizePath(path.expand(root),mustWork=TRUE)
  .libPaths(c(file.path(root,".R-library"),.libPaths()))
  needed <- c("glmnet","digest","jsonlite")
  missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) {
    local_lib <- file.path(root,".R-library")
    dir.create(local_lib,recursive=TRUE,showWarnings=FALSE)
    .libPaths(c(local_lib,.libPaths()))
    utils::install.packages(missing,repos="https://cloud.r-project.org",lib=local_lib)
    ensure(all(vapply(needed,requireNamespace,logical(1),quietly=TRUE)),
           "Package installation incomplete; check the installation messages")
  }


cfg<-elastic_config(); set.seed(cfg$seed)
ensure(length(stage)==1 && stage %in% c("all","validation","test","tests"),"Stage: all|validation|test|tests")
out<-file.path(root,"Results","elastic_net_logit_modular_fixed03_cutoff005")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
message("Run stage: ", stage)
message("Results directory: ", out)
inputs<-c(file.path(root,"Data/processed",paste0("model_ready_h",cfg$horizons,".rds")),
 file.path(root,"Data/crsp_monthly_index_data/crsp_monthly_index_data.csv"),
 file.path(root,"Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))
ensure(all(file.exists(inputs)),"Required input files missing: check project layout")
input_hashes<-setNames(vapply(inputs,file_hash,character(1)),substring(inputs,nchar(root)+2L))
code_files <- file.path(script_dir,c("config.R","models.R","portfolio.R",
  "tests.R","report.R","run_elastic_net.R"))
code_hashes <- setNames(vapply(code_files,file_hash,character(1)),basename(code_files))
message("Running numerical and data-timing tests...")
tests<-run_elastic_tests(root,cfg); write_csv(tests,file.path(out,"tests.csv"))
if(stage=="tests") return(invisible(tests))
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
  ensure(nrow(cutoff_grid)==63L,"Expected 21 cutoffs per horizon")
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
  write_csv(val_audit,file.path(out,"validation_refit_audit.csv"))
  message("Validation frozen: selected hyperparameters and cutoffs for all horizons")
} else {
  ensure(file.exists(state_path),"No validation state: run --stage=validation first")
  state<-readRDS(state_path)
  ensure(identical(cfg,state$frozen$config) && identical(input_hashes,state$frozen$input_hashes) &&
         identical(code_hashes,state$frozen$code_hashes),"Code/data/config changed since validation: stop")
  selected<-state$selected;val_forecasts<-state$forecasts
}
if(stage=="validation") return(invisible(state))
message("Computing test predictions; selection is frozen...")
test_rows<-list(); test_audits<-list(); fit_records<-list()
for(h in cfg$horizons) {
  d<-data[[as.character(h)]]; sel<-selected[selected$h==h,,drop=FALSE]
  ensure(nrow(sel)==1,"Missing unique frozen horizon settings")
  result<-recursive_elastic(d,h,cfg,"test",sel$alpha,sel$lambda,TRUE)
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
    all_stocks=rep(1L,nrow(f)),all_bills=rep(0L,nrow(f)),fixed_50_50=rep(0L,nrow(f)))
  for(rule in names(variants)) {
    strat<-if(rule %in% c("all_stocks","all_bills","fixed_50_50")) rule else cfg$model
    cutoff_rule<-if(strat==cfg$model) rule else "not_applicable"
    cutoff<-if(rule=="validation_tuned") unique(f$cutoff) else if(rule=="fixed_050") cfg$comparison_cutoff else NA_real_
    if(strat==cfg$model) diagnostics[[paste(period,h,rule)]]<-cbind(data.frame(period=period,h=h,model=cfg$model,cutoff_rule=rule),forecast_metrics(f,cutoff))
    for(cost in cfg$costs) {
      p<-portfolio_path(inp,variants[[rule]],cost,cfg,
        fixed_weight=if(rule=="fixed_50_50") 0.5 else NULL)
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
write_csv(data.frame(package=needed,
  version=vapply(needed,function(p) as.character(utils::packageVersion(p)),character(1))),
  file.path(out,"package_versions.csv"))
jsonlite::write_json(list(version=cfg$version,selected=selected,input_hashes=input_hashes,
 code_hashes=code_hashes,tests_passed=nrow(tests),test_used_for_selection=FALSE,
 note="Fixed-lambda model specified after inspection of previous elastic-net test results; exploratory sensitivity, NOT virgin holdout"),
 file.path(out,"manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
write_elastic_report(out,cfg,selected,performance,bind_rows(diagnostics),tests)
message("Validation and test complete. Results directory: ", out)
message("Summary: ",file.path(out,"RESULTS.md"))
message("Sharpe ratios and drawdowns: ",file.path(out,"performance.csv"))



  invisible(list(selected=selected, performance=performance, output_dir=out))
}

cli <- commandArgs(TRUE)
root_arg <- grep("^--root=",cli,value=TRUE)
stage_arg <- grep("^--stage=",cli,value=TRUE)
root <- if(length(root_arg)) sub("^--root=","",root_arg[1]) else
  if(!is.null(PROJECT_ROOT)) PROJECT_ROOT else find_project_root(script_dir)
stage <- if(length(stage_arg)) sub("^--stage=","",stage_arg[1]) else "all"
elastic_net_main(root=root,stage=stage,script_dir=script_dir)
})
