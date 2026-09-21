# Purpose: Define primary ODI and perioperative models and target populations.
# Inputs: Shared configuration and the prepared three-month cohort.
# Outputs: primary_ids and build_primary() for fitting and posterior summaries.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
options(cmdstanr_write_stan_file_dir = file.path(ROOT, "03_models/stan_cache"))
primary_ids <- c("odi_3m", "day_surgery", "los_postop", "pt_comp_any_3m")
build_primary <- function(id) {
  stopifnot(id %in% primary_ids)
  d <- readRDS(file.path(ROOT, "02_data/derived/cohort_3m.rds"))
  depth <- if (id == "pt_comp_any_3m") 14L else 12L
  adapt <- if (id == "odi_3m") .95 else .99
  margin <- if (id == "odi_3m") 7 else NA_real_
  if (id == "odi_3m") {
    form <- bf(as.formula(paste("odi_3m_zib ~", RHS_PRIMARY)), zi ~ treatment + odi_baseline_z)
    fam <- zero_inflated_beta()
    p <- pri_zib()
    dat <- d[!is.na(d$odi_3m), ]
  } else if (id == "los_postop") {
    dat <- d[!is.na(d$los_postop), ]
    dat$los_ordinal <- ordered(pmin(dat$los_postop, 3L))
    form <- bf(as.formula(paste("los_ordinal ~", RHS)))
    fam <- cumulative("logit")
    p <- c(prior(normal(0, 2), class = b), prior(normal(0, 4), class = Intercept))
  } else {
    dat <- d
    form <- bf(as.formula(paste(id, "~", RHS)))
    fam <- bernoulli()
    p <- pri_binary
  }
  list(
    id = id, data = dat, target = d, formula = form, family = fam, prior = p,
    control = list(adapt_delta = adapt, max_treedepth = depth), margin = margin
  )
}
