write_elastic_report <- function(out,cfg,selected,performance,diagnostics,tests) {
  primary<-performance[performance$period=="test" & abs(performance$cost_rate-cfg$primary_cost)<1e-12,]
  choices<-selected[,c("h","alpha","lambda_ratio","validation_Brier","cutoff","validation_sharpe_matched_ladder")]
  lines<-c("# EC4308 Elastic Net — recursive out-of-sample results","",
    paste("Configuration:",cfg$version),"",
    "Predictors: complete 132-column frozen X; one independent model for each horizon.",
    "Validation candidates selected by minimum expanding-window Brier; candidate is (alpha, relative lambda).",
    "Relative lambda = fraction of fitting-only binomial LASSO lambda_max, recomputed at every estimation origin.",
    "The 25-candidate grid was chosen for methodological interpretability, not test performance.",
    "After predictive tuning, cutoffs selected from 0–1 in 0.01 steps by validation net matched-ladder Sharpe at 10 bp.",
    "Test coefficients are re-estimated monthly, while alpha, relative lambda and cutoff are frozen.","",
    "## Frozen validation choices", "",
    paste(capture.output(print(choices,row.names=FALSE)),collapse="\n"),"",
    "## Test performance at 10bp", "",
    paste(capture.output(print(primary[,c("h","strategy","cutoff_rule","CAGR","sharpe_common_rf","sharpe_matched_ladder","max_drawdown","average_stock_weight")],row.names=FALSE)),collapse="\n"),"",
    "## Interpretation", "",
    "Common Sharpe uses 4-week bill accrual for each calendar month. Matched Sharpe uses horizon-specific vintage bill ladder accrual. Both are monthly returns.",
    "The first h-1 signals warm up the h-cohort portfolio. Performance months start February/April/July for h=1/3/6.",
    "Test-period comparisons do not select hyperparameters, cutoffs or model variants. Prior baseline test viewing and supervisor-driven protocol revision must be disclosed.",
    "Monthly maximum drawdown omits intramonth declines. Equity exposure is important when interpreting high Sharpe during bullish periods.",
    "Bill accrual is an idealized proxy and close-to-close execution assumes timely availability of end-of-month inputs.",
    paste(nrow(tests),"focused test groups passed."))
  writeLines(lines,file.path(out,"RESULTS.md"))
}
