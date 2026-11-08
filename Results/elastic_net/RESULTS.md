# EC4308 Elastic Net — recursive out-of-sample results

Configuration: elastic_net_v1_recursive_2026_11_07

Predictors: complete 132-column frozen X; one independent model for each horizon.
Validation candidates selected by minimum expanding-window Brier; candidate is (alpha, relative lambda).
Relative lambda = fraction of fitting-only binomial LASSO lambda_max, recomputed at every estimation origin.
The 25-candidate grid was chosen for methodological interpretability, not test performance.
After predictive tuning, cutoffs selected from 0–1 in 0.01 steps by validation net matched-ladder Sharpe at 10 bp.
Test coefficients are re-estimated monthly, while alpha, relative lambda and cutoff are frozen.

## Frozen validation choices

 h alpha lambda_ratio validation_Brier cutoff validation_sharpe_matched_ladder
 1     1            1        0.2406476    0.5                        0.2906256
 3     1            1        0.2347984    0.5                        0.2870421
 6     1            1        0.2215805    0.5                        0.2691803

## Test performance at 10bp

 h    strategy      cutoff_rule       CAGR sharpe_common_rf
 1 elastic_net validation_tuned 0.12796661        0.8035990
 1 elastic_net        fixed_050 0.12796661        0.8035990
 1  all_stocks   not_applicable 0.12796661        0.8035990
 1   all_bills   not_applicable 0.01594181               NA
 3 elastic_net validation_tuned 0.12594561        0.7861993
 3 elastic_net        fixed_050 0.12594561        0.7861993
 3  all_stocks   not_applicable 0.12594561        0.7861993
 3   all_bills   not_applicable 0.01694944        1.8603172
 6 elastic_net validation_tuned 0.12694289        0.7846701
 6 elastic_net        fixed_050 0.12694289        0.7846701
 6  all_stocks   not_applicable 0.12694289        0.7846701
 6   all_bills   not_applicable 0.01775471        2.1646620
 sharpe_matched_ladder max_drawdown average_stock_weight
             0.8035990 2.467634e-01                    1
             0.8035990 2.467634e-01                    1
             0.8035990 2.467634e-01                    1
                    NA 2.270791e-05                    0
             0.7809807 2.467634e-01                    1
             0.7809807 2.467634e-01                    1
             0.7809807 2.467634e-01                    1
                    NA 2.043762e-05                    0
             0.7763533 2.467634e-01                    1
             0.7763533 2.467634e-01                    1
             0.7763533 2.467634e-01                    1
                    NA 0.000000e+00                    0

## Interpretation

Common Sharpe uses 4-week bill accrual for each calendar month. Matched Sharpe uses horizon-specific vintage bill ladder accrual. Both are monthly returns.
The first h-1 signals warm up the h-cohort portfolio. Performance months start February/April/July for h=1/3/6.
Test-period comparisons do not select hyperparameters, cutoffs or model variants. Prior baseline test viewing and supervisor-driven protocol revision must be disclosed.
Monthly maximum drawdown omits intramonth declines. Equity exposure is important when interpreting high Sharpe during bullish periods.
Bill accrual is an idealized proxy and close-to-close execution assumes timely availability of end-of-month inputs.
14 focused test groups passed.
