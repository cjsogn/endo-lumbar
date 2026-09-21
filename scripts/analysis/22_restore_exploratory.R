# Calendar-adjusted exploratory subgroup and forest analyses, and reconstructed
# post-hoc weighting. All saved individual predictions remain in private output.
source(file.path(Sys.getenv("ENDO_CODE_DIR"),"00_config.R"))
suppressPackageStartupMessages(library(grf))
options(cmdstanr_write_stan_file_dir=file.path(ROOT,"03_models/stan_cache"))

read_exploratory_data <- function() {
 d <- readRDS(file.path(ROOT,"02_data/derived/cohort_3m.rds"))
 cal <- readRDS(file.path(ROOT,"02_data/derived/calendar_basis_specification.rds"))
 z <- (d$calendar_time-cal$center)/cal$scale
 b <- predict(cal$basis,newx=z)
 stopifnot(nrow(d)==421L,identical(levels(d$treatment),c("MSD","ELD")),
  sum(!is.na(d$odi_3m))==328L,length(COVS)==27L,
  all(complete.cases(d[c(COVS,"ct1","ct2")])),
  max(abs(b-as.matrix(d[c("ct1","ct2")])))<1e-12,
  all(d$odi_3m_zib==pmin(d$odi_3m/100,1-1e-6),na.rm=TRUE))
 d$subgroup_age <- factor(ifelse(d$age<50,"younger","older"),levels=c("younger","older"))
 d$subgroup_severity <- factor(ifelse(d$odi_baseline<40,"lower","higher"),levels=c("lower","higher"))
 d$subgroup_symptom <- factor(ifelse(d$symptom_duration_leg<=3,"shorter","longer"),levels=c("shorter","longer"))
 d$subgroup_levels <- factor(ifelse(d$multilevel==0,"single","multiple"),levels=c("single","multiple"))
 d
}
subgroup_ids <- c("age","severity","symptom","levels")
subgroup_labels <- list(age=c("Age <50 years","Age >=50 years"),
 severity=c("Baseline ODI <40","Baseline ODI >=40"),
 symptom=c("Leg symptoms <=12 months or none","Leg symptoms >12 months"),
 levels=c("Single level","Multiple levels"))
build_subgroup <- function(id) {
 stopifnot(id %in% subgroup_ids)
 d<-read_exploratory_data(); sg<-paste0("subgroup_",id)
 cv<-if(id=="levels")setdiff(COVS,"multilevel") else COVS
 rhs<-paste(c(paste0("treatment * ",sg),cv,"ct1","ct2","treatment:ct1","treatment:ct2"),collapse=" + ")
 form<-bf(as.formula(paste("odi_3m_zib ~",rhs)),zi~treatment+odi_baseline_z)
 p<-c(prior(normal(0,1),class=b),prior(normal(0,3),class=Intercept),
      prior(gamma(2,.1),class=phi),prior(normal(0,1.5),class=Intercept,dpar=zi),
      prior(normal(0,1),class=b,dpar=zi))
 list(id=id,sg=sg,target=d,data=d[!is.na(d$odi_3m),],formula=form,prior=p)
}
fit_subgroup <- function(id) {
 b<-build_subgroup(id); ncores<-as.integer(Sys.getenv("ENDO_CHAIN_CORES","4"))
 fn<-file.path(ROOT,"03_models",paste0("exploratory_subgroup_",id,"_calendar"))
 fit<-brm(b$formula,data=b$data,family=zero_inflated_beta(),prior=b$prior,
  chains=4,iter=2000,warmup=1000,cores=ncores,seed=SEED,backend="cmdstanr",
  control=list(adapt_delta=.95,max_treedepth=12),refresh=500,
  file=fn,file_refit="on_change")
 dg<-diagnostics(fit,paste0("exploratory_subgroup_",id,"_initial"),12)
 write_csv(dg,paste0("08_qa/exploratory_subgroup_",id,"_sampling_initial.csv"))
 if(!dg$pass) {
  saveRDS(fit,paste0(fn,"_initial.rds"))
  fit<-update(fit,iter=4000,warmup=2000,cores=ncores,seed=SEED,
   control=list(adapt_delta=.99,max_treedepth=14),file=NULL,refresh=1000)
  saveRDS(fit,paste0(fn,".rds"))
  dg<-diagnostics(fit,paste0("exploratory_subgroup_",id),14)
 }
 write_csv(dg,paste0("08_qa/exploratory_subgroup_",id,"_sampling.csv"))
 if(!dg$pass)stop("Unresolved subgroup sampling diagnostics: ",id)
 lev<-levels(b$target[[b$sg]]); draws<-list(); rows<-list()
 for(k in seq_along(lev)) {
  nd<-b$target[b$target[[b$sg]]==lev[k],]
  draws[[k]]<-gcomp(fit,nd,scale_factor=100)
  sm<-summarize_effect(draws[[k]]$delta)
  # No subgroup NI or superiority probabilities are used as formal tests.
  sm<-sm[c("mean","median","sd","lower","upper","draws")]
  counts<-unlist(lapply(c("MSD","ELD"),function(a)c(
   sum(nd$treatment==a),sum(nd$treatment==a & !is.na(nd$odi_3m)))))
  names(counts)<-c("MSD_total","MSD_observed","ELD_total","ELD_observed")
  rows[[k]]<-cbind(data.frame(subgroup=id,level=lev[k],label=subgroup_labels[[id]][k]),
                  as.data.frame(as.list(counts)),sm)
 }
 diff<-draws[[2]]$delta-draws[[1]]$delta
 ints<-summarize_effect(diff)[c("mean","median","sd","lower","upper","draws")]
 # Extract the named interaction, never depending on the order of records.
 dr<-as_draws_df(fit)
 nm<-paste0("b_treatmentELD:",b$sg,lev[2])
 stopifnot(nm %in% names(dr),length(dr[[nm]])==length(diff))
 ints<-cbind(data.frame(subgroup=id,contrast=paste(lev[2],"minus",lev[1])),ints)
 model_int<-cbind(data.frame(subgroup=id,coefficient=nm,scale="mean-model logit"),
  summarize_effect(dr[[nm]])[c("mean","median","sd","lower","upper","draws")])
 write_csv(do.call(rbind,rows),paste0("04_results/exploratory_subgroup_",id,".csv"))
 write_csv(ints,paste0("04_results/exploratory_subgroup_",id,"_contrast.csv"))
 write_csv(model_int,paste0("04_results/exploratory_subgroup_",id,"_model_interaction.csv"))
 saveRDS(list(levels=lev,level_draws=draws,odi_scale_difference=diff,
  model_scale_interaction=dr[[nm]],formula=b$formula,prior=b$prior),
  file.path(ROOT,"03_models",paste0("exploratory_subgroup_",id,"_draws.rds")))
 cat("Subgroup completed:",id,"\n")
}
run_subgroups <- function() {
 seen<-character()
 for(id in subgroup_ids) {
  b<-build_subgroup(id)
  fn<-cmdstanr::write_stan_file(make_stancode(b$formula,data=b$data,
       family=zero_inflated_beta(),prior=b$prior))
  if(!fn %in% seen) {cmdstanr::cmdstan_model(fn,quiet=TRUE);seen<-c(seen,fn)}
 }
 write_csv(data.frame(id=subgroup_ids),"08_qa/exploratory_subgroups_manifest.csv")
 Sys.setenv(ENDO_MANIFEST="exploratory_subgroups_manifest.csv",
  ENDO_SCRIPT="22_restore_exploratory.R",ENDO_BATCH="exploratory_subgroups")
 status<-system2("python3",shQuote(file.path(Sys.getenv("ENDO_CODE_DIR"),"04_model_scheduler.py")))
 if(status!=0L)stop("Subgroup worker failure. See private worker logs.")
 for(suffix in c("","_contrast","_model_interaction")) {
  vals<-lapply(subgroup_ids,function(id)read.csv(file.path(ROOT,"04_results",
                    paste0("exploratory_subgroup_",id,suffix,".csv"))))
  write_csv(do.call(rbind,vals),paste0("04_results/exploratory_subgroups",suffix,".csv"))
 }
}
run_forest <- function() {
 d<-read_exploratory_data(); d<-d[!is.na(d$odi_3m),]
 # Raw continuous covariates retain the original forest implementation.
 vars<-sub("_z$","",COVS)
 X<-model.matrix(reformulate(c(vars,"ct1","ct2")),data=d)[,-1,drop=FALSE]
 stopifnot(all(is.finite(X)),!any(grepl("treatment|odi_3m",colnames(X))),
  all(c("ct1","ct2") %in% colnames(X)))
 W<-as.numeric(d$treatment=="ELD");Y<-as.numeric(d$odi_3m)
 set.seed(SEED)
 cf<-causal_forest(X=X,Y=Y,W=W,num.trees=4000,min.node.size=5,
  honesty=TRUE,tune.parameters="all",num.threads=14,seed=SEED,
  compute.oob.predictions=TRUE)
 saveRDS(cf,file.path(ROOT,"03_models/exploratory_causal_forest_calendar.rds"))
 ate<-average_treatment_effect(cf,target.sample="all",method="AIPW")
 cate<- -as.numeric(predict(cf)$predictions)
 stopifnot(length(cate)==nrow(d),all(is.finite(cate)),
   all(is.finite(cf$Y.hat)),all(is.finite(cf$W.hat)),
   all(cf$W.hat>0 & cf$W.hat<1))
 summary<-data.frame(method="Causal forest",mean= -unname(ate[1]),
   se=unname(ate[2]),lower= -unname(ate[1])-qnorm(.975)*unname(ate[2]),
   upper= -unname(ate[1])+qnorm(.975)*unname(ate[2]),n=nrow(d))
 write_csv(summary,"04_results/exploratory_causal_forest_effect.csv")
 cal<-test_calibration(cf,vcov.type="HC3")
 mat_table<-function(x) data.frame(term=rownames(x),estimate=x[,1],se=x[,2],
                                  statistic=x[,3],p_value=x[,4],row.names=NULL)
 write_csv(mat_table(cal),"04_results/exploratory_causal_forest_calibration.csv")
 blp_vars<-c("age","odi_baseline","symptom_duration_leg")
 stopifnot(all(blp_vars %in% colnames(X)))
 blp<-best_linear_projection(cf,A=X[,blp_vars,drop=FALSE],vcov.type="HC3")
 bt<-mat_table(blp);bt$estimate<- -bt$estimate;bt$statistic<- -bt$statistic
 # Two-sided BLP p-values are unchanged by reversing the contrast.
 write_csv(bt,"04_results/exploratory_causal_forest_projection.csv")
 vi<-data.frame(variable=colnames(X),importance=as.numeric(variable_importance(cf)))
 vi<-vi[order(vi$importance,decreasing=TRUE),]
 write_csv(vi,"04_results/exploratory_causal_forest_importance.csv")
 quant<-data.frame(probability=c(0,.025,.25,.5,.75,.975,1),
   odi_contrast=as.numeric(quantile(cate,c(0,.025,.25,.5,.75,.975,1))))
 write_csv(quant,"04_results/exploratory_causal_forest_conditional_distribution.csv")
 overlap<-do.call(rbind,lapply(0:1,function(a)data.frame(arm=c("MSD","ELD")[a+1],
  n=sum(W==a),minimum=min(cf$W.hat[W==a]),median=median(cf$W.hat[W==a]),
  maximum=max(cf$W.hat[W==a]),below_005=sum(cf$W.hat[W==a]<.05),
  above_095=sum(cf$W.hat[W==a]>.95))))
 write_csv(overlap,"08_qa/exploratory_causal_forest_overlap.csv")
 saveRDS(list(row_index=which(!is.na(read_exploratory_data()$odi_3m)),X=X,
  cate=cate,Y_hat=cf$Y.hat,W_hat=cf$W.hat,tuning=cf$tuning.output,
  calibration=cal,projection=blp),file.path(ROOT,"03_models/exploratory_causal_forest_details.rds"))
 writeLines(c(capture.output(cf$tuning.output),capture.output(sessionInfo())),
   file.path(ROOT,"10_logs/exploratory_causal_forest_session.txt"))
 cat("Causal forest completed.\n")
}

# Weighting reconstruction. Propensities use all cohort procedures, without
# outcome or treatment-derived predictors. Outcome means use observed ODI only.
weight_estimate <- function(X,W,Y) {
 warnings<-character()
 fit<-withCallingHandlers(glm.fit(x=X,y=W,family=binomial(),
  control=glm.control(maxit=100)),warning=function(w) {
   warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
 e<-fit$fitted.values
 if(!fit$converged || any(!is.finite(e)) || length(unique(W))!=2L)
   stop("Logistic propensity fit is not estimable or did not converge.")
 assigned<-ifelse(W==1,e,1-e)
 if(any(assigned<=0))stop("Zero probability for an observed treatment.")
 sw<-ifelse(W==1,mean(W),1-mean(W))/assigned
 cap<-unname(quantile(sw,.99,type=7))
 weights<-list(capped_iptw=pmin(sw,cap),overlap=ifelse(W==1,1-e,e))
 obs<-!is.na(Y)
 means<-sapply(weights,function(w)sapply(0:1,function(a) {
  take<-obs & W==a
  if(sum(w[take])<=0)stop("No observed outcome weight in an arm.")
  weighted.mean(Y[take],w[take])
 }))
 list(delta=means[1,]-means[2,],means=means,weights=weights,e=e,
      cap=cap,uncapped=sw,fit=fit,warnings=unique(warnings))
}
weight_bootstrap_batch <- function(jobs,design,W,Y,seeds) {
 lapply(jobs,function(k) {
  set.seed(seeds[k]);ix<-sample.int(length(W),replace=TRUE)
  tryCatch({z<-weight_estimate(design[ix,,drop=FALSE],W[ix],Y[ix])
   data.frame(replicate=k,capped_iptw=z$delta[1],overlap=z$delta[2],
    warnings=paste(z$warnings,collapse=" | "),error="",row.names=NULL)
  },error=function(e)data.frame(replicate=k,capped_iptw=NA_real_,overlap=NA_real_,
                               warnings="",error=conditionMessage(e)))
 })
}
run_weighting <- function() {
 d<-read_exploratory_data()
 X<-model.matrix(reformulate(c(COVS,"ct1","ct2")),d)
 W<-as.numeric(d$treatment=="ELD");Y<-as.numeric(d$odi_3m)
 stopifnot(all(is.finite(X)),!any(grepl("treatment|odi_3m",colnames(X))))
 z<-weight_estimate(X,W,Y)
 saveRDS(list(X=X,W=W,Y=Y,estimate=z),file.path(ROOT,"03_models/exploratory_weighting_calendar.rds"))
 B<-2000L;set.seed(SEED);seeds<-sample.int(.Machine$integer.max,B)
 cl<-parallel::makePSOCKcluster(14L);on.exit(parallel::stopCluster(cl),add=TRUE)
 parallel::clusterExport(cl,c("weight_estimate","weight_bootstrap_batch"),envir=environment())
 batches<-split(seq_len(B),rep(seq_len(14L),length.out=B))
 raw<-parallel::parLapply(cl=cl,X=batches,fun=weight_bootstrap_batch,
                         design=X,W=W,Y=Y,seeds=seeds)
 boot<-do.call(rbind,unlist(raw,recursive=FALSE));boot<-boot[order(boot$replicate),]
 write_csv(boot,"08_qa/exploratory_weighting_bootstrap.csv")
 saveRDS(list(seeds=seeds,bootstrap=boot),file.path(ROOT,"03_models/exploratory_weighting_bootstrap.rds"))
 if(any(nzchar(boot$error)) || any(!is.finite(as.matrix(boot[c("capped_iptw","overlap")]))))
   stop("Weighting bootstrap failures saved. No final interval reported.")
 res<-do.call(rbind,lapply(names(z$weights),function(nm)data.frame(method=nm,
  mean=unname(z$delta[nm]),se=sd(boot[[nm]]),
  lower=unname(quantile(boot[[nm]],.025)),upper=unname(quantile(boot[[nm]],.975)),
  bootstrap_replicates=B,n_propensity=length(W),n_outcome=sum(!is.na(Y)))))
 write_csv(res,"04_results/exploratory_weighting_effects.csv")
 wm<-function(x,w)sum(w*x)/sum(w)
 # A fixed unweighted pooled SD makes pre/post weighting SMDs comparable.
 bal<-list();diags<-list();obs<-!is.na(Y)
 for(pop in c("full_cohort","observed_ODI")) {
  take<-if(pop=="full_cohort")rep(TRUE,length(W)) else obs
  for(nm in c("unweighted",names(z$weights))) {
   w<-if(nm=="unweighted")rep(1,length(W)) else z$weights[[nm]]
   for(a in 0:1) {
    ii<-take & W==a
    diags[[length(diags)+1]]<-data.frame(population=pop,method=nm,arm=c("MSD","ELD")[a+1],
     n=sum(ii),ESS=sum(w[ii])^2/sum(w[ii]^2),weight_min=min(w[ii]),
     weight_median=median(w[ii]),weight_max=max(w[ii]),
     propensity_min=min(z$e[ii]),propensity_median=median(z$e[ii]),propensity_max=max(z$e[ii]))
   }
   for(v in setdiff(colnames(X),"(Intercept)")) {
    a0<-take & W==0;a1<-take & W==1
    s<-sqrt((var(X[a0,v])+var(X[a1,v]))/2)
    difference<-wm(X[a1,v],w[a1])-wm(X[a0,v],w[a0])
    smd<-if(s>0)difference/s else if(difference==0)0 else NA_real_
    bal[[length(bal)+1]]<-data.frame(population=pop,method=nm,variable=v,smd=smd)
   }
  }
 }
 write_csv(do.call(rbind,diags),"08_qa/exploratory_weighting_diagnostics.csv")
 write_csv(do.call(rbind,bal),"08_qa/exploratory_weighting_balance.csv")
 warnings<-data.frame(source=c("cohort propensity","bootstrap propensity"),
   count=c(length(z$warnings),sum(nzchar(boot$warnings))),
   detail=c(paste(z$warnings,collapse=" | "),paste(unique(boot$warnings[nzchar(boot$warnings)]),collapse=" | ")))
 write_csv(warnings,"08_qa/exploratory_weighting_warnings.csv")
 write_csv(data.frame(cap=z$cap,percentile=.99,uncapped_max=max(z$uncapped),
  convergence=z$fit$converged,iterations=z$fit$iter,rank=z$fit$rank,columns=ncol(X)),
  "08_qa/exploratory_weighting_specification.csv")
 writeLines(c("Reconstruction: pooled 99th-percentile capped stabilised weights and overlap weights.",
  "Intervals: percentile bootstrap, propensity and cap refit per procedure resample.",
  capture.output(sessionInfo())),file.path(ROOT,"10_logs/exploratory_weighting_session.txt"))
 cat("Weighting completed.\n")
}

args<-commandArgs(trailingOnly=TRUE)
if(Sys.getenv("ENDO_MODEL_WORKER","")=="1") {
 for(id in args)fit_subgroup(id)
} else {
 modes<-if(length(args))args else c("subgroups","forest","weighting")
 stopifnot(all(modes %in% c("subgroups","forest","weighting")))
 for(mode in modes)switch(mode,subgroups=run_subgroups(),forest=run_forest(),weighting=run_weighting())
 writeLines(capture.output(sessionInfo()),file.path(ROOT,"10_logs/exploratory_restoration_session.txt"))
}
