# EC4308 Elastic Net: shared parameters for h=1,3,6.
elastic_config <- function() list(
  version = "elastic_net_logit_teacher_fixed03_cutoff005_v1_2026_11_17",
  horizons = c(1L,3L,6L), model = "elastic_net_logit_teacher_fixed03_cutoff005",
  train_end = "2004-12", validation_start = "2005-01", validation_end = "2012-12",
  test_start = "2013-01", test_end = "2025-12",
  penalty_convention = "teacher: lambda * (alpha * L2_squared + (1-alpha) * L1)",
  # Fixed ABSOLUTE lambda, identical at every recursive fitting origin.
  # Non-uniform shared grid; 0.3 is a practical upper bound, not a theoretical optimum.
  # Ex-post methodological sensitivity: do not describe test as untouched holdout.
  alphas = c(0,0.25,0.5,0.75,1),
  lambdas = c(0.0001,0.0005,0.002,0.01,0.05,0.1,0.2,0.3),
  # 21 candidate cutoffs on [0,1], in steps of 0.05.
  cutoffs = seq(0,1,by=0.05), comparison_cutoff = 0.5,
  cutoff_objective = "sharpe_matched_ladder", primary_cost = 0.001,
  costs = c(0,0.001,0.002), risk_aversion = 3,
  common_rf_security = 2000061L,
  initial_allocation_cost = FALSE,
  standardize = "fitting_rows_only", selection_score = "Brier",
  grid_tie_rule = "smaller number of active terms then larger absolute lambda then smaller alpha",
  seed = 4308L,
  bill_rule = "Origin-vintage matched daily log yield; actual monthly quote days",
  cost_rule = "Exact post-cost net equity trade at 10bp; 0/20bp sensitivities"
)
month_number <- function(x) {
  x <- as.character(x)
  as.integer(substr(x,1,4))*12L + as.integer(substr(x,6,7))-1L
}
month_string <- function(x) sprintf("%04d-%02d",x %/% 12L,x %% 12L+1L)
ensure <- function(ok,message) if(!isTRUE(ok)) stop(message,call.=FALSE)
write_csv <- function(x,path) utils::write.csv(x,path,row.names=FALSE,na="")
file_hash <- function(path) digest::digest(file=path,algo="sha256",serialize=FALSE)
bind_rows <- function(xs) if(length(xs)) do.call(rbind,xs) else data.frame()
