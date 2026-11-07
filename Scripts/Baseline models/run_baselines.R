# Run from any directory: Rscript --vanilla '/path/to/Scripts/baseline models/run_baselines.R'
args <- commandArgs(FALSE)
entry <- grep("^--file=", args, value = TRUE)
# Rscript encodes spaces in --file as ~+~ on this R/macOS build.
entry_path <- if (length(entry)) gsub("~+~", " ", sub("^--file=", "", entry[1]), fixed=TRUE) else ""
script_dir <- if (nzchar(entry_path)) dirname(normalizePath(entry_path, mustWork=TRUE)) else normalizePath("Scripts/baseline models")
root <- normalizePath(file.path(script_dir, "../.."))
.libPaths(c(file.path(root,".R-library"), .libPaths()))
ensure_packages <- c("digest", "jsonlite")
missing <- ensure_packages[!vapply(ensure_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Run Scripts/00_setup.R first. Missing: ", paste(missing,collapse=", "))
for (file in c("config.R","models.R","portfolio.R","tests.R","report.R")) source(file.path(script_dir,file))
cfg <- baseline_config()
stage <- sub("^--stage=", "", grep("^--stage=", commandArgs(TRUE), value=TRUE))
if (!length(stage)) stage <- "all"
ensure(length(stage)==1L && stage %in% c("all","validation","test","tests"), "Use --stage=all|validation|test|tests")
out <- file.path(root,"Results","baseline models")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
input_files <- c(file.path(root,"Data/processed",paste0("model_ready_h",cfg$horizons,".rds")),
                 file.path(root,"Data/crsp_monthly_index_data/crsp_monthly_index_data.csv"),
                 file.path(root,"Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))
hashes <- setNames(vapply(input_files,file_hash,character(1)),substring(input_files,nchar(root)+2L))
code_files <- sort(list.files(script_dir,pattern="[.]R$",full.names=TRUE))
code_hashes <- setNames(vapply(code_files,file_hash,character(1)),basename(code_files))
message("Checking model timing and portfolio arithmetic...")
tests <- run_baseline_tests(root,cfg)
write_csv(tests,file.path(out,"tests.csv"))
message(nrow(tests)," focused test groups passed.")
if (stage=="tests") quit(save="no",status=0)

returns <- load_return_data(root)
datasets <- setNames(lapply(cfg$horizons,function(h) load_model_data(root,h,cfg)),cfg$horizons)
bind <- function(x) do.call(rbind,x)
validation_file <- file.path(out,"validation_state.rds")
if (stage %in% c("all","validation")) {
  message("Generating monthly expanding validation forecasts and selecting matched-Sharpe cutoffs...")
  forecasts <- grids <- selections <- fits <- records <- list()
  for (h in cfg$horizons) {
    d <- datasets[[as.character(h)]]
    for (link in cfg$links) {
      key <- paste(h,link,sep="_")
      result <- expanding_forecasts(d,h,link,cfg,period="validation")
      f <- result$forecasts
      inp <- portfolio_inputs(f,returns,h,cfg$validation_end,cfg$common_rf_security)
      grid <- bind(lapply(cfg$cutoffs,function(cutoff) {
        p <- portfolio_path(inp,as.integer(f$probability >= cutoff),cfg$primary_cost,cfg)
        cbind(data.frame(h=h,model=link,cutoff=cutoff),portfolio_metrics(p,cfg$risk_aversion))
      }))
      winner <- grid[choose_cutoff(grid,cfg$cutoff_objective),]
      default <- grid[abs(grid$cutoff-cfg$comparison_cutoff)<1e-12,]
      ensure(nrow(default)==1L,"Comparison cutoff missing from validation grid")
      selections[[key]] <- data.frame(h=h,model=link,cutoff=winner$cutoff,
        selection_objective=cfg$cutoff_objective,
        validation_Sharpe=winner$sharpe_matched_ladder,
        validation_sharpe_matched_ladder=winner$sharpe_matched_ladder,
        validation_sharpe_common_rf=winner$sharpe_common_rf,
        fixed_050_validation_sharpe_matched_ladder=default$sharpe_matched_ladder,
        fixed_050_validation_sharpe_common_rf=default$sharpe_common_rf,
        validation_CAGR=winner$CAGR,
        cost_rate=cfg$primary_cost,validation_start=winner$start,validation_end=winner$end)
      f$cutoff <- winner$cutoff; f$signal <- as.integer(f$probability >= f$cutoff)
      f$cutoff_050 <- cfg$comparison_cutoff; f$signal_050 <- as.integer(f$probability >= f$cutoff_050)
      forecasts[[key]] <- f; grids[[key]] <- grid; fits[[key]] <- result$final_fit
      records <- c(records,result$records)
    }
  }
  selected <- bind(selections); grid <- bind(grids); val_forecasts <- bind(forecasts)
  rownames(selected) <- rownames(grid) <- rownames(val_forecasts) <- NULL
  # This file is written before the first final-test prediction is generated.
  frozen <- list(config=cfg,selected=selected,input_hashes=hashes,code_hashes=code_hashes,
                 frozen_at=format(Sys.time(),tz="Asia/Singapore",usetz=TRUE))
  saveRDS(frozen,file.path(out,"frozen_specification.rds"))
  state <- list(forecasts=val_forecasts,grid=grid,selected=selected,fits=fits,records=records,frozen=frozen)
  saveRDS(state,validation_file)
  saveRDS(fits,file.path(out,"validation_models.rds"))
  write_csv(grid,file.path(out,"validation_cutoff_grid.csv"))
  write_csv(selected,file.path(out,"selected_cutoffs.csv"))
  write_csv(val_forecasts,file.path(out,"validation_forecasts.csv"))
  message("Validation choices frozen: ",paste(paste0("h",selected$h," ",selected$model,"=",sprintf("%.2f",selected$cutoff)),collapse=", "))
} else {
  ensure(file.exists(validation_file), "Run --stage=validation before --stage=test")
  state <- readRDS(validation_file)
  ensure(identical(state$frozen$config,cfg) && identical(state$frozen$input_hashes,hashes) &&
           identical(state$frozen$code_hashes,code_hashes),
         "Configuration, code or inputs changed since validation. Re-run validation explicitly.")
  val_forecasts <- state$forecasts; selected <- state$selected; grid <- state$grid
}
if (stage=="validation") {
  message("Validation complete; final-test fitting and performance have not run.")
  quit(save="no",status=0)
}

message("Generating monthly expanding-window test forecasts with frozen cutoffs...")
test_forecasts <- final_fits <- list(); records <- state$records
for (h in cfg$horizons) for (link in cfg$links) {
  key <- paste(h,link,sep="_")
  result <- expanding_forecasts(datasets[[as.character(h)]],h,link,cfg)
  f <- result$forecasts
  f$cutoff <- selected$cutoff[selected$h==h & selected$model==link]
  f$signal <- as.integer(f$probability >= f$cutoff)
  f$cutoff_050 <- cfg$comparison_cutoff; f$signal_050 <- as.integer(f$probability >= f$cutoff_050)
  test_forecasts[[key]] <- f; final_fits[[key]] <- result$final_fit
  records <- c(records,result$records)
}
test_forecasts <- bind(test_forecasts); rownames(test_forecasts) <- NULL
all_forecasts <- rbind(val_forecasts,test_forecasts)
audit <- bind(lapply(records,`[[`,"audit")); rownames(audit) <- NULL
ensure(all(month_number(audit$last_train_outcome) <= month_number(audit$forecast_origin)),"Final timing audit failed")
test_audit <- audit[audit$period=="test",]
validation_audit <- audit[audit$period=="validation",]
ensure(nrow(test_audit)==156*length(cfg$links)*length(cfg$horizons),"Missing monthly refits")
ensure(nrow(validation_audit)==96*length(cfg$links)*length(cfg$horizons),"Missing validation refits")
ensure(all(month_number(audit$last_train_origin) <= month_number(audit$forecast_origin)-audit$h),"Immature label in fit")

paths <- performance <- diagnostics <- calibration <- benchmark_returns <- list()
for (period in c("validation","test")) for (h in cfg$horizons) {
  f0 <- all_forecasts[all_forecasts$period==period & all_forecasts$h==h & all_forecasts$model==cfg$links[1],]
  end <- if(period=="validation") cfg$validation_end else cfg$test_end
  inp <- portfolio_inputs(f0,returns,h,end,cfg$common_rf_security)
  benchmark_returns[[paste(period,h)]] <- data.frame(period=period,h=h,month=inp$ym,
    quote_start=inp$quote_start,quote_end=inp$quote_end,quote_days=inp$days,
    common_rf_security=inp$common_rf_security,matched_bill_security=inp$matched_bill_security,
    common_rf_daily_yield=inp$common_rf_daily_yield,
    common_rf_return=inp$common_rf_return,matched_ladder_return=inp$bill_reference)
  for (strategy in c(cfg$links,"all_stocks","all_bills")) {
    if (strategy %in% cfg$links) {
      f <- all_forecasts[all_forecasts$period==period & all_forecasts$h==h & all_forecasts$model==strategy,]
      key <- paste(period,h,strategy,sep="_")
      variants <- list(validation_tuned=list(signals=f$signal,cutoff=unique(f$cutoff)),
                       fixed_050=list(signals=f$signal_050,cutoff=cfg$comparison_cutoff))
      for (rule in names(variants)) {
        diagnostics[[paste(key,rule)]] <- cbind(data.frame(period=period,h=h,model=strategy,cutoff_rule=rule),
                                               forecast_metrics(f,variants[[rule]]$cutoff))
      }
      eligible <- f[f$classification_eligible & !is.na(f$target),]
      bins <- cut(eligible$probability,breaks=seq(0,1,by=.1),include.lowest=TRUE)
      calibration[[key]] <- bind(lapply(levels(bins),function(bin) {
        ids <- which(bins==bin)
        data.frame(period=period,h=h,model=strategy,bin=bin,n=length(ids),
          mean_probability=if(length(ids)) mean(eligible$probability[ids]) else NA_real_,
          event_rate=if(length(ids)) mean(eligible$target[ids]) else NA_real_)
      }))
    } else variants <- list(not_applicable=list(signals=rep(as.integer(strategy=="all_stocks"),nrow(f0)),cutoff=NA_real_))
    for (rule in names(variants)) for (cost in cfg$costs) {
      key <- paste(period,h,strategy,rule,cost,sep="_")
      p <- portfolio_path(inp,variants[[rule]]$signals,cost,cfg)
      tags <- data.frame(period=period,h=h,strategy=strategy,cutoff_rule=rule,
                         cutoff=variants[[rule]]$cutoff,cost_rate=cost)
      paths[[key]] <- cbind(tags,p)
      performance[[key]] <- cbind(tags,portfolio_metrics(p,cfg$risk_aversion))
      ensure(max(month_number(p$month)) <= month_number(end),"Portfolio return crosses split boundary")
    }
  }
}
paths <- bind(paths); performance <- bind(performance); diagnostics <- bind(diagnostics)
rownames(paths) <- rownames(performance) <- rownames(diagnostics) <- NULL
common <- list()
for(period in c("validation","test")) {
  period_start <- if(period=="validation") cfg$validation_start else cfg$test_start
  common_start <- month_string(month_number(period_start)+max(cfg$horizons))
  p <- paths[paths$period==period & paths$month >= common_start,]
  groups <- interaction(p$h,p$strategy,p$cutoff_rule,p$cost_rate,drop=TRUE)
  for(key in unique(groups)) {
    q <- p[groups==key,]
    common[[paste(period,key)]] <- cbind(q[1,c("period","h","strategy","cutoff_rule","cutoff","cost_rate")],portfolio_metrics(q,cfg$risk_aversion))
  }
}
ensure(identical(hashes,setNames(vapply(input_files,file_hash,character(1)),names(hashes))), "Input files changed during run")
write_csv(test_forecasts,file.path(out,"test_forecasts.csv"))
write_csv(paths,file.path(out,"monthly_portfolios.csv"))
write_csv(bind(benchmark_returns),file.path(out,"benchmark_returns.csv"))
write_csv(performance,file.path(out,"performance.csv"))
write_csv(bind(common),file.path(out,"common_window_performance.csv"))
write_csv(diagnostics,file.path(out,"forecast_metrics.csv"))
write_csv(bind(calibration),file.path(out,"calibration.csv"))
write_csv(audit,file.path(out,"refit_audit.csv"))
write_csv(bind(lapply(records,`[[`,"coefficients")),file.path(out,"coefficients.csv"))
write_csv(bind(lapply(records,`[[`,"scalers")),file.path(out,"scalers.csv"))
write_csv(data.frame(file=names(hashes),sha256=unname(hashes)),file.path(out,"input_hashes.csv"))
saveRDS(final_fits,file.path(out,"final_models.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))
write_csv(data.frame(package=ensure_packages,version=vapply(ensure_packages,function(p) as.character(packageVersion(p)),character(1))),file.path(out,"package_versions.csv"))
manifest <- list(version=cfg$version,completed_at=format(Sys.time(),tz="Asia/Singapore",usetz=TRUE),
  config=cfg,input_hashes=hashes,code_hashes=code_hashes,tests_passed=nrow(tests),
  monthly_validation_refits=nrow(validation_audit),monthly_test_refits=nrow(test_audit),validation_selected_cutoffs=selected,
  test_used_for_selection=FALSE,input_hashes_unchanged=TRUE,
  revision_note="Supervisor protocol revision after baseline v1 test results had already been inspected")
jsonlite::write_json(manifest,file.path(out,"manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
write_baseline_report(out,cfg,selected,performance,diagnostics,tests)
draw_baseline_figures(out,paths,grid,selected,cfg)
message("Complete. Results: ",file.path(out,"RESULTS.md"))
