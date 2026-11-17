# Locate the project and load the baseline modules.
script_dir <- local({
  sourced <- unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names=FALSE)
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
  paths <- gsub("~+~", " ", c(rev(sourced), entry), fixed=TRUE)
  paths <- paths[basename(paths)=="run_baselines.R"]
  candidates <- dirname(paths)
  here <- getwd()
  repeat {
    candidates <- c(candidates, here, file.path(here,"baseline"), file.path(here,"Scripts","baseline"))
    parent <- dirname(here)
    if (identical(parent,here)) break
    here <- parent
  }
  modules <- c("run_baselines.R","config.R","models.R","portfolio.R","tests.R","figures.R")
  valid <- vapply(candidates, function(path) all(file.exists(file.path(path,modules))), logical(1))
  if (!any(valid)) stop("Cannot locate Scripts/baseline. Open the project folder or source run_baselines.R using its full path.", call.=FALSE)
  normalizePath(candidates[which(valid)[1]], mustWork=TRUE)
})
root <- normalizePath(file.path(script_dir, "../.."))
.libPaths(c(file.path(root,".R-library"), .libPaths()))
ensure_packages <- "digest"
missing <- ensure_packages[!vapply(ensure_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install the required R package with install.packages('digest'). Missing: ", paste(missing,collapse=", "))
for (file in c("config.R","models.R","portfolio.R","tests.R","figures.R")) source(file.path(script_dir,file))
cfg <- baseline_config()
ensure(length(commandArgs(TRUE)) == 0L, "Run without stage arguments; validation and testing execute together.")
out <- file.path(root,"Results","baseline")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
input_files <- c(file.path(root,"Data/processed",paste0("model_ready_h",cfg$horizons,".rds")),
                 file.path(root,"Data/crsp_monthly_index_data/crsp_monthly_index_data.csv"),
                 file.path(root,"Data/crsp_treasury_risk_free_series_data/crsp_treasury_risk_free_series_data.csv"))
hashes <- setNames(vapply(input_files,file_hash,character(1)),substring(input_files,nchar(root)+2L))
message("Checking model timing and portfolio arithmetic...")
tests <- run_baseline_tests(root,cfg)

# Load the monthly market returns, bill yields and horizon-specific model inputs.
returns <- load_return_data(root)
datasets <- setNames(lapply(cfg$horizons,function(h) load_model_data(root,h,cfg)),cfg$horizons)
bind <- function(x) do.call(rbind,x)
# Fit each validation origin, then select one cutoff per model and horizon.
message("Generating monthly expanding validation forecasts and selecting matched-Sharpe cutoffs...")
forecasts <- grids <- selections <- records <- list()
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
    forecasts[[key]] <- f; grids[[key]] <- grid
    records <- c(records,result$records)
  }
}
selected <- bind(selections); grid <- bind(grids); val_forecasts <- bind(forecasts)
rownames(selected) <- rownames(grid) <- rownames(val_forecasts) <- NULL
# Keep the selected cutoffs fixed while fitting the test-period models.
write_csv(grid,file.path(out,"validation_cutoff_grid.csv"))
write_csv(selected,file.path(out,"selected_cutoffs.csv"))
write_csv(val_forecasts,file.path(out,"validation_forecasts.csv"))
message("Selected cutoffs: ",paste(paste0("h",selected$h," ",selected$model,"=",sprintf("%.2f",selected$cutoff)),collapse=", "))

message("Generating monthly expanding-window test forecasts with frozen cutoffs...")
test_forecasts <- list()
for (h in cfg$horizons) for (link in cfg$links) {
  key <- paste(h,link,sep="_")
  result <- expanding_forecasts(datasets[[as.character(h)]],h,link,cfg)
  f <- result$forecasts
  f$cutoff <- selected$cutoff[selected$h==h & selected$model==link]
  f$signal <- as.integer(f$probability >= f$cutoff)
  f$cutoff_050 <- cfg$comparison_cutoff; f$signal_050 <- as.integer(f$probability >= f$cutoff_050)
  test_forecasts[[key]] <- f
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

# Evaluate both cutoff rules and passive strategies over the same return months.
paths <- performance <- diagnostics <- calibration <- list()
for (period in c("validation","test")) for (h in cfg$horizons) {
  f0 <- all_forecasts[all_forecasts$period==period & all_forecasts$h==h & all_forecasts$model==cfg$links[1],]
  end <- if(period=="validation") cfg$validation_end else cfg$test_end
  inp <- portfolio_inputs(f0,returns,h,end,cfg$common_rf_security)
  for (strategy in c(cfg$links,"all_stocks","all_bills")) {
    if (strategy %in% cfg$links) {
      f <- all_forecasts[all_forecasts$period==period & all_forecasts$h==h & all_forecasts$model==strategy,]
      key <- paste(period,h,strategy,sep="_")
      variants <- list(validation_tuned=list(signals=f$signal,cutoff=unique(f$cutoff)),
                       fixed_050=list(signals=f$signal_050,cutoff=cfg$comparison_cutoff))
      probability_summary <- baseline_probability_diagnostics(f,cfg)
      for (rule in names(variants)) {
        diagnostics[[paste(key,rule)]] <- cbind(data.frame(period=period,h=h,model=strategy,cutoff_rule=rule),
                                               forecast_metrics(f,variants[[rule]]$cutoff),probability_summary)
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
# Recalculate performance on a common calendar for comparisons across horizons.
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
# Save the forecasts, trading results and diagnostics used in the research report.
write_csv(test_forecasts,file.path(out,"test_forecasts.csv"))
write_csv(paths,file.path(out,"monthly_portfolios.csv"))
write_csv(performance,file.path(out,"performance.csv"))
write_csv(bind(common),file.path(out,"common_window_performance.csv"))
write_csv(diagnostics,file.path(out,"forecast_metrics.csv"))
write_csv(bind(calibration),file.path(out,"calibration.csv"))
write_csv(bind(lapply(records,`[[`,"coefficients")),file.path(out,"coefficients.csv"))
invisible(draw_baseline_figures(out,paths,cfg))
invisible(draw_probability_figure(out,test_forecasts,cfg))
message("PASS: ",nrow(tests)," check groups; ",nrow(validation_audit)," validation fits and ",nrow(test_audit)," test fits; timing and input checks passed.")
message("Complete. Research outputs saved in: ",out)
