source(file.path(Sys.getenv("ENDO_CODE_DIR"),"02_primary_jobs.R"))
res<-diag<-list()
for(id in primary_ids) {
 b<-build_primary(id)
 fit<-readRDS(file.path(ROOT,"03_models",paste0(id,"_calendar.rds")))
 diag[[id]]<-diagnostics(fit,id,b$control$max_treedepth)
 stopifnot(diag[[id]]$pass)
 if(id=="los_postop") draws<-data.frame(delta=-as_draws_df(fit)$b_treatmentELD) else
  draws<-gcomp(fit,b$target,id!="day_surgery",if(id=="odi_3m")100 else 1)
 res[[id]]<-cbind(id=id,summarize_effect(draws$delta,b$margin))
 saveRDS(draws,file.path(ROOT,"04_results",paste0(id,"_calendar_draws.rds")))
 if(id=="odi_3m") {
  write_csv(do.call(rbind,lapply(c(7,5,3),function(m)summarize_effect(draws$delta,m))),
   "04_results/primary_stricter_margins.csv")
  saveRDS(gcomp(fit,b$target[!is.na(b$target$odi_3m),],TRUE,100),
   file.path(ROOT,"04_results/odi_3m_calendar_observed_target_draws.rds"))
 }
}
write_csv(do.call(rbind,res),"04_results/verified_calendar_primary_perioperative.csv")
write_csv(do.call(rbind,diag),"08_qa/verified_calendar_diagnostics.csv")
