# One R-only command from any working directory: Rscript /path/to/Scripts/run_all.R
# The entry point in run_all.R runs the tests, builds the predictors, runs stationarity diagnostics, 
# prepares targets and splits, validates the outputs, and writes a report.

arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
scripts <- if(length(arg)) dirname(normalizePath(sub("^--file=","",arg[1]))) else normalizePath("Scripts")
source(file.path(scripts,"R/common.R"))
source(file.path(scripts,"R/build.R"))
source(file.path(scripts,"R/stationarity.R"))
source(file.path(scripts,"R/model_data.R"))
source(file.path(scripts,"tests.R"))
run_tests()
build_all()
stationarity_all()
prepare_model_data()
source(file.path(scripts,"R/validate.R"))
validate_all()
source(file.path(scripts,"R/report.R"))
write_validation_report()
writeLines(capture.output(sessionInfo()),file.path(AUDIT,"R_session_info.txt"))
fwrite(data.table(package=REQUIRED,version=vapply(REQUIRED,function(p) as.character(packageVersion(p)),character(1))),file.path(scripts,"packages.lock.csv"))
message("R pipeline complete. Read Scripts/README.md and Scripts/VALIDATION.md before modeling.")
