# Build the model inputs Data/processed/model_ready_h{1,3,6}.rds from the raw files in Data/.
# Run: Rscript --vanilla Scripts/data/run_build.R   (or open and click Source)
script <- c(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE),
            tryCatch(sys.frame(1)$ofile, error = function(e) NULL))
source(file.path(if (length(script)) dirname(script[1]) else "Scripts/data", "..", "utils", "utils.R"))
use_packages(c("data.table", "urca"))
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
for (file in c("build.R", "targets.R", "stationarity.R", "tests.R")) source(file.path(SCRIPTS, "data", file))

checks <- run_checks(data_checks())
predictors <- build_predictors()
stationarity <- stationarity_table(predictors)

processed <- file.path(DATA, "processed"); validation <- file.path(DATA, "validation")
dir.create(processed, recursive = TRUE, showWarnings = FALSE); dir.create(validation, recursive = TRUE, showWarnings = FALSE)
for (h in names(predictors))
  saveRDS(model_ready(predictors[[h]]$features, as.integer(h)), file.path(processed, paste0("model_ready_h", h, ".rds")))
fwrite(stationarity, file.path(validation, "stationarity.csv"))
write_run_info(validation, list(specification = SPEC_ID, sample = "1974-04 to 2025-12", predictors = length(MODEL_NAMES),
                                features = length(FEATURE_NAMES), stationarity_sample = "1974-04 to 2004-12"),
               c(checks, "All transformed predictors pass ADF and KPSS"), c("data.table", "urca"), file.path(DATA, INPUTS))
message("Model inputs saved in ", processed)
