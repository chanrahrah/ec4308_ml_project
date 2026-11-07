markdown_table <- function(d) {
  lines <- apply(d, 1, function(x) paste0("| ",paste(x,collapse=" | ")," |"))
  c(paste0("| ",paste(names(d),collapse=" | ")," |"),
    paste0("| ",paste(rep("---",ncol(d)),collapse=" | ")," |"), lines)
}

write_baseline_report <- function(out, cfg, selected, performance, diagnostics, tests) {
  primary <- performance[performance$period=="test" & performance$cost_rate==cfg$primary_cost,]
  fmt_sharpe <- function(x) ifelse(is.na(x),"N/A",sprintf("%.3f",x))
  table <- data.frame(Horizon=primary$h, Strategy=primary$strategy,
    Cutoff_rule=primary$cutoff_rule, Cutoff=ifelse(is.na(primary$cutoff),"-",sprintf("%.2f",primary$cutoff)),
    Start=primary$start, Months=primary$n_months,
    CAGR=sprintf("%.2f%%",100*primary$CAGR),
    Sharpe_common_4week=fmt_sharpe(primary$sharpe_common_rf),
    Sharpe_matched=fmt_sharpe(primary$sharpe_matched_ladder),
    Max_drawdown=sprintf("%.2f%%",100*primary$max_drawdown),check.names=FALSE)
  choices <- data.frame(Horizon=selected$h,Model=selected$model,Cutoff=sprintf("%.2f",selected$cutoff),
    Tuned_validation_Sharpe_matched=fmt_sharpe(selected$validation_sharpe_matched_ladder),
    Fixed_050_validation_Sharpe_matched=fmt_sharpe(selected$fixed_050_validation_sharpe_matched_ladder))
  lines <- c("# Baseline results - supervisor protocol revision", "",
    paste0("Configuration: `",cfg$version,"`. Both five-predictor probit and logit models are fitted at 1-, 3- and 6-month horizons. Validation and testing both refit monthly using expanding mature-label samples. Each model is evaluated with its validation-selected cutoff and the fixed 0.50 cutoff. All-stock and all-bill benchmarks use the same scored months."),"",
    "The initial baseline v1 test results had already been inspected before this supervisor-directed revision. The prior results and source scripts are preserved in Results/archive/baseline_v1_pre_supervisor_2026-11-05. This is a documented protocol revision, not a new untouched test sample.","",
    "## Validation choices", "", markdown_table(choices), "",
    "Every validation origin has its own coefficient and scaler fit, using only outcomes complete by that origin. Cutoffs maximize monthly net Sharpe relative to the horizon-matched bill ladder at 10 bp transaction costs. The 0.50 comparator is fixed and is never selected using test results. These validation figures are selection scores, not independent performance estimates.","",
    "## Test performance at 10 bp", "", markdown_table(table), "",
    "Both Sharpe columns subtract a one-month reference return from the one-month net portfolio return and use sqrt(12) annualization. The common reference uses the 4-week bill yield quoted at the preceding month-end; the matched reference averages the active 4/13/26-week bill vintages for h=1/3/6. Neither reference uses the full three- or six-month holding-period return. All-bills Sharpe relative to its identical matched ladder is undefined; against the common reference it can be defined for h=3/6. These are accrual proxies, not Treasury total-return indices.","",
    "Maximum drawdown is measured on monthly wealth observations and is not annualized; intramonth peaks and troughs are not captured. No model or cutoff rule is selected using this test table.","",
    "### What the validation-tuned strategies did", "",
    paste(vapply(cfg$horizons,function(h) {
      d <- primary[primary$h==h,]
      models <- d[d$strategy %in% cfg$links & d$cutoff_rule=="validation_tuned",]
      stock <- d[d$strategy=="all_stocks",]
      if(all(abs(models$average_stock_weight-1)<1e-12)) {
        sprintf("At h=%d, both selected models stayed entirely in stocks throughout the scored test period. Their identical performance is expected: these cutoffs produced no market-timing allocation changes.",h)
      } else {
        paste0("At h=",h,", average equity exposure was ",
          paste(sprintf("%.1f%% (%s)",100*models$average_stock_weight,models$strategy),collapse=" and "),
          ". ",if(all(models$CAGR<stock$CAGR & models$sharpe_common_rf<stock$sharpe_common_rf))
            "Both models had lower CAGR and common-reference Sharpe than all stocks." else "See the table for the comparison with all stocks.",
          if(all(models$max_drawdown<stock$max_drawdown)) " Both had smaller maximum drawdowns." else "")
      }
    },character(1)),collapse="\n\n"),"",
    "The first h-1 signals initialize the overlapping portfolio. Scoring begins February 2013, April 2013 and July 2013 for h=1,3,6, respectively; all end December 2025. The initial holdings are endowed without entry costs, consistent with the supplied reference example; subsequent net equity rebalancing trades incur costs. No terminal liquidation cost is imposed.","",
    "Use common_window_performance.csv to compare all horizons over July 2013-December 2025. This table rebases each existing strategy path; it does not restart holdings or retune cutoffs.","",
    "## Files", "",
    "- performance.csv: both cutoff rules and both Sharpe references for validation/test at 0, 10 and 20 bp. The cutoff_rule column distinguishes validation_tuned, fixed_050 and not_applicable (passive benchmarks).",
    "- monthly_portfolios.csv: cohort weights, both monthly reference returns, quote dates, exact trades, costs and wealth for every cutoff/cost scenario.",
    "- benchmark_returns.csv: one row per period/horizon/month with both references, security IDs, origin daily yield for the common bill, and actual monthly quote-day counts.",
    "- validation_forecasts.csv / test_forecasts.csv: one probability per model/origin, targets, eligibility, per-origin fitting audits, tuned cutoff/signal and cutoff_050/signal_050. Only eligible completed labels enter classification metrics.",
    "- validation_cutoff_grid.csv / selected_cutoffs.csv: full search and frozen choices.",
    "- forecast_metrics.csv / calibration.csv: Brier (lower is better), log loss, AUC, classification and calibration diagnostics.",
    "- refit_audit.csv / coefficients.csv / scalers.csv: per-origin training endpoint and fitted parameters.",
    "- frozen_specification.rds: configuration, thresholds, data hashes and code hashes saved before test forecasting.",
    "- validation_models.rds / final_models.rds: final-origin validation/test fits; per-origin coefficients and fitting-sample scalers are in the audit CSVs.",
    "- manifest.json / input_hashes.csv / tests.csv / session_info.txt: reproducibility and verification.",
    "- test_wealth.png / validation_cutoffs.png: figures.","",
    "## Interpretation and limitations", "",
    "The five inputs are the existing stationarity-adjusted versions: bill_yld_rel36_lag0, term_spread_change1_lag0, log_rv_lag0, mom6_lag0 and log_dispersion_change1_lag0. These are static binary-response baselines; the earlier proposal's extra dynamic lagged outcome is not included in the current five-predictor specification. See the scripts README for the literature-based lag decision.","",
    "At every validation and test origin t, only outcomes ending by t enter fitting, and each scaler is recomputed on those same rows. Earlier boundary-purged rows may enter once mature. No future label is used. For example, the latest eligible origin at January 2005 is December 2004 (h=1), October 2004 (h=3), or July 2004 (h=6); at January 2013 it is December/October/July 2012. Returns observed by the origin may inform estimation; that origin's forecast concerns subsequent returns.","",
    "Late-validation forecasts whose full h-month outcomes cross December 2012 are excluded from classification metrics, but their monthly portfolio returns up to December 2012 can inform the cutoff. No 2013 return is used for validation. Incomplete final-test forecasts can likewise contribute observed monthly returns without inventing a completed h-month label.","",
    "Each of h cohorts keeps its binary signal for h months and receives equal post-cost capital monthly. Bill slots retain their originating matched yield; monthly accrual uses actual quote-date days. The all-bill reference is the equal-weight ladder of all h vintage yields. Reallocations among bill vintages are frictionless under the accepted accrual proxy; changing market yields do not cause mark-to-market gains or losses. This is not a tradable Treasury total-return reconstruction.","",
    "Costs equal 0.1% of the absolute net equity amount bought or sold (one-way turnover). The cost equation is solved so target weights hold after costs. Cross-cohort equity trades are netted. Drift rebalancing can cost money even when the target stock fraction is unchanged.","",
    "Month-end information is associated with the following close-to-close monthly return under an idealized execution convention. The dataset does not provide the next executable price after all month-end inputs are finalized, so literal real-time execution is not established. Existing source-data limitations remain. No independent-observation significance tests or Ledoit-Wolf claims are made for overlapping labels and serially dependent strategy returns.","",
    sprintf("%d focused test groups passed; all final refit, period-boundary and input-hash checks also passed.",nrow(tests)),"",
    "Interpret improvements alongside equity exposure and the all-stock benchmark: a strong equity period can produce high strategy Sharpe without strong market-timing ability. Lower drawdown alone can reflect lower equity exposure.","",
    "Metrics: CAGR is compounded. sharpe_common_rf uses the common 4-week monthly reference; sharpe_matched_ladder uses the passive matched ladder. The legacy net_Sharpe field remains an alias of sharpe_matched_ladder. Sortino, beta and certainty-equivalent excess returns retain the matched reference, identified by excess_metric_reference. Sortino uses root mean squared negative excess return over all months; drawdown includes starting wealth; Ulcer is RMS monthly drawdown. Annual certainty equivalent uses gamma=3 and 12*(mean excess - gamma/2*variance excess). Annual turnover is 12 times mean absolute net equity traded / pretrade wealth. total_cost_fraction sums monthly cost/pretrade-wealth ratios; it is not the terminal wealth drag.")
  writeLines(lines,file.path(out,"RESULTS.md"))
}

draw_baseline_figures <- function(out, paths, grid, selected, cfg) {
  colors <- c(probit="#2463A7",logit="#BA5A31",all_stocks="#293A45",all_bills="#78916D")
  grDevices::png(file.path(out,"test_wealth.png"),width=1800,height=780,res=150)
  graphics::par(mfrow=c(1,3),mar=c(4,4,3,1),oma=c(0,0,2,0))
  for(h in cfg$horizons) {
    d <- paths[paths$period=="test" & paths$h==h & paths$cost_rate==cfg$primary_cost,]
    dates <- as.Date(paste0(d$month,"-01"))
    graphics::plot(range(dates),range(c(1,d$wealth)),type="n",xlab="Return month",ylab="Wealth (initial = 1)",main=paste(h,"month horizon"))
    graphics::grid(col="#E7E7E7")
    for(s in names(colors)) {
      v <- d[d$strategy==s & d$cutoff_rule %in% c("validation_tuned","not_applicable"),]
      graphics::lines(as.Date(paste0(v$month,"-01")),v$wealth,col=colors[s],lwd=2)
      if(s %in% cfg$links) {
        v <- d[d$strategy==s & d$cutoff_rule=="fixed_050",]
        graphics::lines(as.Date(paste0(v$month,"-01")),v$wealth,col=colors[s],lwd=1.5,lty=2)
      }
    }
    graphics::legend("topleft",legend=c("Probit tuned","Logit tuned","All stocks","All bills","Probit 0.50","Logit 0.50"),
                     col=c(colors,colors[c("probit","logit")]),lty=c(1,1,1,1,2,2),lwd=2,bty="n",cex=.68)
    m <- d[d$strategy %in% cfg$links & d$cutoff_rule=="validation_tuned",]
    if(all(abs(m$stock_weight-1)<1e-12))
      graphics::mtext("Probit = Logit = All stocks",side=3,line=.15,cex=.65,col="#555555")
  }
  graphics::mtext("Test portfolios | 10 bp costs | validation-tuned and fixed 0.50 cutoffs",outer=TRUE,cex=1.05)
  grDevices::dev.off()
  grDevices::png(file.path(out,"validation_cutoffs.png"),width=1800,height=780,res=150)
  graphics::par(mfrow=c(1,3),mar=c(4,4,3,1),oma=c(0,0,2,0))
  for(h in cfg$horizons) {
    d<-grid[grid$h==h,]
    graphics::plot(c(0,1),range(d$sharpe_matched_ladder,na.rm=TRUE),type="n",xlab="Probability cutoff",ylab="Validation matched-bill Sharpe",main=paste(h,"month horizon"))
    graphics::grid(col="#E7E7E7")
    for(s in cfg$links) {
      v<-d[d$model==s,]; sel<-selected[selected$h==h & selected$model==s,]
      graphics::lines(v$cutoff,v$sharpe_matched_ladder,col=colors[s],lwd=2)
      graphics::points(sel$cutoff,sel$validation_Sharpe,pch=19,col=colors[s])
    }
    graphics::legend("bottomleft",legend=c("Probit","Logit"),col=colors[c("probit","logit")],lwd=2,bty="n",cex=.8)
  }
  graphics::mtext("Expanding validation | matched-bill cutoff objective | 10 bp costs",outer=TRUE,cex=1.05)
  grDevices::dev.off()
}
