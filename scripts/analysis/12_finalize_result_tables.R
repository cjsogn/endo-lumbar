source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
d <- readRDS(file.path(ROOT,"02_data/derived/cohort_3m.rds"))
new <- read.csv(file.path(ROOT,"04_results/calendar_all_bayesian_results.csv"))
los <- new[new$id=="los_postop",]
write_csv(data.frame(estimate=exp(los$mean),lower=exp(los$lower),upper=exp(los$upper),
 summary="Exponentiated posterior mean log OR, as in submitted reporting"),"04_results/los_odds_ratio.csv")

# Unadjusted endpoint summaries with actual outcome-specific denominators.
ids <- c(new$id,"perop_comp_any")
d$operating_time[d$operating_time>=500] <- NA
d12 <- readRDS(file.path(ROOT,"02_data/derived/cohort_12m_date_eligible.rds"))
desc <- do.call(rbind,lapply(ids,function(id) {
 z<-if(grepl("12m",id))d12 else d
 do.call(rbind,lapply(c("MSD","ELD"),function(a) {
  y<-z[[id]][z$treatment==a];y<-y[!is.na(y)]
  data.frame(id=id,arm=a,eligible=sum(z$treatment==a),observed=length(y),
   events=if(all(y %in% c(0,1)))sum(y)else NA_real_,mean=mean(y),sd=sd(y),
   median=median(y),q25=unname(quantile(y,.25)),q75=unname(quantile(y,.75)))
 }))
}))
write_csv(desc,"04_results/endpoint_descriptive_verified.csv")

# Conditional association checks retain the submitted likelihood-ratio method,
# adding the shared calendar spline. These do not test away unmeasured anatomy.
fals <- lapply(c("sex","bmi_z","education"),function(v) {
 full<-glm(as.formula(paste("I(treatment=='ELD') ~",paste(c(COVS,"ct1","ct2"),collapse=" + "))),data=d,family=binomial())
 reduced<-glm(as.formula(paste("I(treatment=='ELD') ~",paste(c(setdiff(COVS,v),"ct1","ct2"),collapse=" + "))),data=d,family=binomial())
 lr<-anova(reduced,full,test="LRT")
 k<-grep(paste0("^",v),names(coef(full)))
 data.frame(variable=v,lr_chi_sq=lr$Deviance[2],df=lr$Df[2],p=lr$`Pr(>Chi)`[2],
 max_abs_coefficient=max(abs(coef(full)[k])))
})
write_csv(do.call(rbind,fals),"04_results/falsification_calendar.csv")

# Consolidate all actual sampler settings and diagnostics for reproducibility.
diag_files<-list.files(file.path(ROOT,"08_qa"),"_diagnostics.csv$",full.names=TRUE)
diag_files<-diag_files[!grepl("all_calendar|calendar_psis|verified_calendar|final_mcmc",basename(diag_files))]
dg0<-read.csv(file.path(ROOT,"08_qa/verified_calendar_diagnostics.csv"))
dg<-rbind(dg0, do.call(rbind,lapply(diag_files,function(fn) read.csv(fn)[names(dg0)])))
stopifnot(!anyDuplicated(dg$id),all(dg$pass))
write_csv(dg,"08_qa/final_mcmc_diagnostics.csv")
settings<-lapply(list.files(file.path(ROOT,"03_models"),".rds$",full.names=TRUE),function(fn) {
 x<-readRDS(fn)
 if(!inherits(x,"brmsfit") || grepl("_initial.rds$",fn)) return(NULL)
 a<-x$fit@stan_args[[1]]
 data.frame(id=basename(fn),chains=length(x$fit@stan_args),iterations=a$iter,
 warmup=a$warmup,thin=a$thin,adapt_delta=a$control$adapt_delta,
 max_treedepth=a$control$max_treedepth,seed=as.character(a$seed))
})
write_csv(do.call(rbind,settings),"08_qa/final_sampler_settings.csv")
print(do.call(rbind,fals));print(range(dg$rhat_max));print(range(dg$ess_bulk_min))
