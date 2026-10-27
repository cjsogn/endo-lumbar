# Purpose: Fit alternative primary likelihoods, priors and adjustment sets.
# Inputs: Prepared three-month cohort and coded source fields for level/stenosis.
# Outputs: Sensitivity models, posterior contrasts and sampling diagnostics.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
options(cmdstanr_write_stan_file_dir = file.path(ROOT, "03_models/stan_cache"))
d <- readRDS(file.path(ROOT, "02_data/derived/cohort_3m.rds"))
a <- readRDS(file.path(ROOT, "02_data/source/df_all.rds"))
ix <- match(d$ForlopsID, a$ForlopsID)
stopifnot(!anyNA(ix))
for (v in c("level_L34", "level_L45", "level_L5S1", "stenosis_lateral", "stenosis_foraminal")) d[[v]] <- a[[v]][ix]
ids <- c(
  "prior_skeptical", "prior_diffuse", "gaussian", "student_nu5", "horseshoe",
  "restricted", "operated_level", "beta_sv", "pure_disc"
)
make_job <- function(id) {
  full <- d
  rhs <- RHS_PRIMARY
  p <- pri_zib()
  f <- zero_inflated_beta()
  sf <- 100
  if (id == "pure_disc") full <- full[full$stenosis_central == 0 & full$stenosis_lateral == 0 & full$stenosis_foraminal == 0, ]
  full$y <- full$odi_3m_zib
  if (id == "restricted") rhs <- paste(
    "treatment + age_z + sex + odi_baseline_z + prior_surgery_any",
    "+ ct1 + ct2 + treatment:ct1 + treatment:ct2"
  )
  if (id == "operated_level") rhs <- paste(rhs, "+ level_L34 + level_L45 + level_L5S1")
  if (id == "prior_skeptical") p <- pri_zib(.5)
  if (id == "prior_diffuse") p <- pri_zib(2)
  if (id == "horseshoe") p <- c(
    prior(horseshoe(df = 1), class = b), prior(normal(0, 3), class = Intercept),
    prior(gamma(2, .1), class = phi), prior(normal(0, 1.5), class = Intercept, dpar = zi), prior(normal(0, 1), class = b, dpar = zi)
  )
  form <- bf(as.formula(paste("y ~", rhs)), zi ~ treatment + odi_baseline_z)
  dat <- full[!is.na(full$y), ]
  if (id %in% c("gaussian", "student_nu5")) {
    full$y <- full$odi_3m
    dat <- full
    sf <- 1
    p <- pri_odi_gaussian
    form <- bf(as.formula(paste("y | mi() ~", rhs)))
    f <- gaussian()
    if (id == "student_nu5") {
      f <- student()
      form <- bf(as.formula(paste("y | mi() ~", rhs)), nu = 5)
    }
  }
  if (id == "beta_sv") {
    n <- sum(!is.na(full$odi_3m))
    full$y <- ((full$odi_3m / 100) * (n - 1) + .5) / n
    dat <- full
    form <- bf(as.formula(paste("y | mi() ~", rhs)))
    f <- Beta()
    p <- c(
      prior(normal(0, 1), class = b, coef = treatmentELD), prior(normal(0, .5), class = b),
      prior(normal(0, 3), class = Intercept), prior(gamma(2, .1), class = phi)
    )
  }
  list(id = id, formula = form, prior = p, family = f, data = dat, target = full, scale = sf)
}
run_job <- function(id, cores) {
  b <- make_job(id)
  dest <- file.path(ROOT, "03_models", paste0("primary_", id, "_calendar"))
  con <- file(file.path(ROOT, "10_logs", paste0("sensitivity_", id, ".log")), "wt")
  sink(con)
  sink(con, type = "message")
  on.exit(
    {
      sink(type = "message")
      sink()
      close(con)
    },
    add = TRUE
  )
  ctl <- if (id %in% c("horseshoe", "student_nu5")) list(adapt_delta = .99, max_treedepth = 14) else list(adapt_delta = .95, max_treedepth = 12)
  fit <- brm(b$formula,
    data = b$data, family = b$family, prior = b$prior, chains = 4, iter = 2000, warmup = 1000,
    cores = cores, seed = SEED, control = ctl, file = dest, file_refit = "on_change", refresh = 1000
  )
  dg <- diagnostics(fit, paste0("sensitivity_", id), ctl$max_treedepth)
  if (!dg$pass) {
    saveRDS(fit, paste0(dest, "_initial.rds"))
    fit <- update(fit,
      iter = 6000, warmup = 2000, cores = cores, seed = SEED,
      control = list(adapt_delta = .999, max_treedepth = 15), file = NULL, refresh = 2000
    )
    saveRDS(fit, paste0(dest, ".rds"))
    dg <- diagnostics(fit, paste0("sensitivity_", id), 15)
  }
  write_csv(dg, paste0("08_qa/sensitivity_", id, "_diagnostics.csv"))
  if (!dg$pass) stop(paste("Unresolved sampler diagnostics:", id))
  x <- gcomp(fit, b$target, TRUE, b$scale)
  # SV back transformation is affine. Differences require n/(n-1), not an
  # intercept correction. Keep this explicit rather than calling y*100 exact.
  if (id == "beta_sv") {
    n <- sum(!is.na(d$odi_3m))
    x$delta <- x$delta * n / (n - 1)
    for (arm in c("MSD", "ELD")) x[[arm]] <- (x[[arm]] * n - 50) / (n - 1)
  }
  saveRDS(x, file.path(ROOT, "04_results", paste0("sensitivity_", id, "_draws.rds")))
  write_csv(cbind(id = id, summarize_effect(x$delta, 7)), paste0("04_results/sensitivity_", id, "_summary.csv"))
  capture.output(prior_summary(fit), file = file.path(ROOT, "08_qa", paste0("sensitivity_", id, "_resolved_priors.txt")))
}
if (Sys.getenv("ENDO_LOAD_ONLY", "") != "1") {
  if (Sys.getenv("ENDO_MODEL_WORKER", "") == "1") {
    for (id in commandArgs(trailingOnly = TRUE)) run_job(id, as.integer(Sys.getenv("ENDO_CHAIN_CORES")))
  } else {
    manifest <- list()
    seen <- character()
    for (id in ids) {
      b <- make_job(id)
      code <- make_stancode(b$formula, data = b$data, family = b$family, prior = b$prior)
      fn <- cmdstanr::write_stan_file(code)
      if (!fn %in% seen) {
        cat("Compiling", id, "\n")
        cmdstanr::cmdstan_model(fn, quiet = TRUE)
        seen <- c(seen, fn)
      }
      manifest[[id]] <- data.frame(
        id = id, stan_file = fn, fit_rows = nrow(b$data), target_rows = nrow(b$target),
        formula = paste(deparse(b$formula$formula), collapse = " ")
      )
    }
    write_csv(do.call(rbind, manifest), "08_qa/sensitivity_models_manifest.csv")
    Sys.setenv(
      ENDO_MANIFEST = "sensitivity_models_manifest.csv", ENDO_SCRIPT = "06_primary_model_sensitivities.R",
      ENDO_BATCH = "sensitivity"
    )
    s <- system2("python3", shQuote(file.path(CODE_DIR, "04_model_scheduler.py")))
    if (s != 0) stop("Sensitivity worker failure. Read endpoint logs.")
  }
}
