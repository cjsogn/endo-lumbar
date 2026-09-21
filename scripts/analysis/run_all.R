# Run from any working directory. Data and output locations must be explicit.
file_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(file_arg)!=1L) stop("Run this entry point with Rscript.")
script <- sub("^--file=","",file_arg)
# Rscript can encode spaces in its --file command argument.
if(!file.exists(script)) script <- gsub("~+~"," ",script,fixed=TRUE)
code <- dirname(normalizePath(script,mustWork=TRUE))
Sys.setenv(ENDO_CODE_DIR=code)
source(file.path(code,"00_config.R"))
required <- c("brms","posterior","haven","cmdstanr","dplyr","tidyr","tibble",
              "tableone","loo","tmle","SuperLearner","mice","glmnet","ranger","xgboost","grf")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop(paste("Install the required packages:",paste(missing,collapse=", ")))
if(!nzchar(Sys.which("python3"))) stop("Python 3 is required for scheduling.")
cmdstanr::cmdstan_path()
if(!nzchar(RAW_DATA) || !file.exists(RAW_DATA)) stop("Set ENDO_LUMBAR_RAW_DATA to the private SPSS export.")
Sys.unsetenv(c("ENDO_MODEL_WORKER","ENDO_CHAIN_CORES","ENDO_LOAD_ONLY",
              "ENDO_MANIFEST","ENDO_SCRIPT","ENDO_BATCH"))
stages <- c("01_import_registry.R","01_prepare_data.R","02_fit_primary.R",
 "02_summarize_primary.R","04_run_calendar_models.R","05_primary_derived_sensitivities.R",
 "06_primary_model_sensitivities.R","07_tmle_matched_checks.R","08_selection_model.R",
 "09_model_and_population_audit.R","11_predictive_diagnostics.R",
 "12_finalize_result_tables.R","13_descriptive_analyses.R","22_restore_exploratory.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)) {
 if(any(!args %in% stages)) stop("Unknown stage. Supply exact stage filenames listed in run_all.R.")
 stages<-stages[stages %in% args]
}
versions<-data.frame(package=required,version=vapply(required,function(p)as.character(packageVersion(p)),character(1)))
write_csv(versions,"10_logs/package_versions.csv")
for(stage in stages) {
 cat("Running",stage,"\n");flush.console()
 status<-system2(file.path(R.home("bin"),"Rscript"),shQuote(file.path(code,stage)),
                 stdout=file.path(ROOT,"10_logs",paste0(stage,".log")),
                 stderr=file.path(ROOT,"10_logs",paste0(stage,".log")))
 if(status!=0) stop(paste("Stage failed:",stage,". Read its log in ENDO_WORK_DIR/10_logs."))
}
cat("Requested stages finished. Output is in",ROOT,"\n")
