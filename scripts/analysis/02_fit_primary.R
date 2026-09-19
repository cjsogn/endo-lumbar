source(file.path(Sys.getenv("ENDO_CODE_DIR"),"02_primary_jobs.R"))
if(Sys.getenv("ENDO_MODEL_WORKER","")=="1") {
 ids<-commandArgs(trailingOnly=TRUE)
 stopifnot(length(ids)>0L,all(ids %in% primary_ids))
 for(id in ids) {
  b<-build_primary(id)
  fit<-brm(b$formula,data=b$data,family=b$family,prior=b$prior,
   chains=4,iter=2000,warmup=1000,cores=as.integer(Sys.getenv("ENDO_CHAIN_CORES")),
   seed=SEED,control=b$control,backend="cmdstanr",refresh=500,
   file=file.path(ROOT,"03_models",paste0(id,"_calendar")),file_refit="on_change")
  dg<-diagnostics(fit,id,b$control$max_treedepth)
  if(!dg$pass) stop(paste("Unresolved primary/perioperative sampler diagnostics:",id))
 }
} else {
 seen<-character()
 for(id in primary_ids) {
  b<-build_primary(id)
  fn<-cmdstanr::write_stan_file(make_stancode(b$formula,data=b$data,family=b$family,prior=b$prior))
  if(!fn %in% seen) {cmdstanr::cmdstan_model(fn,quiet=TRUE);seen<-c(seen,fn)}
 }
 write_csv(data.frame(id=primary_ids),"08_qa/primary_models_manifest.csv")
 Sys.setenv(ENDO_MANIFEST="primary_models_manifest.csv",ENDO_SCRIPT="02_fit_primary.R",ENDO_BATCH="primary")
 status<-system2("python3",shQuote(file.path(CODE_DIR,"04_model_scheduler.py")))
 if(status!=0) stop("A primary/perioperative model failed. Read the worker logs.")
}
