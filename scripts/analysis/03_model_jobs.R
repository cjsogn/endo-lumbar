source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
cache<-file.path(ROOT,"03_models/stan_cache"); dir.create(cache,showWarnings=FALSE)
options(cmdstanr_write_stan_file_dir=cache)
jobs <- list()
add <- function(id,family,scale=1,lower=TRUE,margin=NA_real_,baseline=NULL) {
 jobs[[id]] <<- list(id=id,family=family,scale=scale,lower=lower,margin=margin,baseline=baseline)
}
add("odi_12m","zib",100,TRUE,7,"odi_baseline_z")
for(tp in c("3m","12m")) {
 for(site in c("back","leg")) add(paste0("nrs_",site,"_",tp),"zib",10,TRUE,1,paste0("nrs_",site,"_baseline_z"))
 add(paste0("eq5d_",tp),"eq",1,FALSE,.05)
 for(y in c("rtw","analgesic","satisfied","gpe_success")) add(paste0(y,"_",tp),"binary",1,y=="analgesic",.1)
 add(paste0("eq5d_anxiety_",tp),"negative",1,TRUE)
}
add("responder_3m","binary",1,FALSE,.1)
add("operating_time","optime",1,TRUE)
build_job <- function(j) {
 d<-readRDS(file.path(ROOT,"02_data/derived",if(grepl("12m$",j$id))
  "cohort_12m_date_eligible.rds" else "cohort_3m.rds"))
 d$y<-as.numeric(haven::zap_labels(d[[j$id]]))
 if(j$family=="optime") d$y[d$y>=500]<-NA_real_
 observed<-sum(!is.na(d$y))
 if(j$family=="zib") {
  d$y<-pmin(d$y/j$scale,1-1e-6);d$baseline_zi<-d[[j$baseline]]
  formula<-bf(as.formula(paste("y ~",RHS)),zi~treatment+baseline_zi)
  family<-zero_inflated_beta(); prior<-pri_zib()
 } else if(j$family=="binary") {
  formula<-bf(as.formula(paste("y ~",RHS)));family<-bernoulli();prior<-pri_binary
 } else {
  formula<-bf(as.formula(paste("y | mi() ~",RHS)));family<-gaussian()
  prior<-switch(j$family,eq=pri_eq,
   negative=c(prior(normal(0,2),class=b,coef=treatmentELD),prior(normal(0,2),class=b),
    prior(student_t(3,0,3),class=sigma),prior(normal(2,2),class=Intercept)),
   optime=c(prior(normal(0,30),class=b,coef=treatmentELD),prior(normal(0,2),class=b),
    prior(student_t(3,0,30),class=sigma),prior(normal(60,30),class=Intercept)))
 }
 target<-d
 if(j$family %in% c("zib","binary")) d<-d[!is.na(d$y),]
 list(job=j,data=d,target=target,observed=observed,formula=formula,family=family,prior=prior)
}
