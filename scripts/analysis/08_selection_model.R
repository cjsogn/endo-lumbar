# Purpose: Fit a joint outcome/response model across a fixed sensitivity grid.
# Inputs: Prepared three-month cohort and shared analysis settings.
# Outputs: Stan source, selection-model fits, posterior contrasts and diagnostics.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
rho_grid <- c(-.10, -.05, 0, .05, .10)
ids <- paste0("rho_", seq_along(rho_grid))
d <- readRDS(file.path(ROOT, "02_data/derived/cohort_3m.rds"))
# Both outcome and response regressions use this reduced baseline set,
# the calendar spline and both treatment-by-calendar interactions.
f <- ~ treatment + age_z + sex + bmi_z + odi_baseline_z + eq5d_baseline_z + nrs_back_baseline_z +
  nrs_leg_baseline_z + motor_deficit + depression_anxiety + chronic_pain + prior_surgery_any +
  sick_leave + disability + ct1 + ct2 + treatment:ct1 + treatment:ct2
X <- model.matrix(f, d)[, -1, drop = FALSE]
e <- m <- d
e$treatment <- factor("ELD", levels = c("MSD", "ELD"))
m$treatment <- factor("MSD", levels = c("MSD", "ELD"))
dx <- colMeans(model.matrix(f, e)[, -1, drop = FALSE] - model.matrix(f, m)[, -1, drop = FALSE])
obs <- which(!is.na(d$odi_3m))
mis <- which(is.na(d$odi_3m))
dat <- list(
  N = nrow(d), N_obs = length(obs), N_mis = length(mis), K = ncol(X), obs_idx = obs,
  mis_idx = mis, y_obs = d$odi_3m[obs], X = X, R = as.integer(!is.na(d$odi_3m)), dx = dx, rho = 0
)
stan <- '
data {
 int<lower=1> N; int<lower=1> N_obs; int<lower=0> N_mis; int<lower=1> K;
 array[N_obs] int obs_idx; array[N_mis] int mis_idx;
 vector[N_obs] y_obs; matrix[N,K] X; array[N] int<lower=0,upper=1> R;
 vector[K] dx; real rho;
}
parameters {
 vector[K] beta; real alpha_y; real<lower=0> sigma;
 vector[K] gamma; real alpha_r; vector[N_mis] y_mis;
}
model {
 vector[N] y_full;
 beta ~ normal(0,10); alpha_y ~ normal(30,20); sigma ~ student_t(3,0,15);
 gamma ~ normal(0,2); alpha_r ~ normal(0,5);
 y_full[obs_idx]=y_obs; y_full[mis_idx]=y_mis;
 y_full ~ normal(alpha_y+X*beta,sigma);
 R ~ bernoulli_logit(alpha_r+X*gamma+rho*y_full);
}
generated quantities {real delta=-dot_product(dx,beta);}
'
fn <- file.path(ROOT, "03_models/selection_calendar.stan")
run_job <- function(id) {
  k <- match(id, ids)
  rho <- rho_grid[k]
  dat$rho <- rho
  con <- file(file.path(ROOT, "10_logs", paste0("selection_", id, ".log")), "wt")
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
  mod <- cmdstanr::cmdstan_model(fn, quiet = TRUE)
  outdir <- file.path(ROOT, "03_models", paste0("selection_", id))
  dir.create(outdir, showWarnings = FALSE)
  fit <- mod$sample(
    data = dat, chains = 4, parallel_chains = as.integer(Sys.getenv("ENDO_CHAIN_CORES", "4")),
    iter_warmup = 1000, iter_sampling = 1000, seed = SEED, adapt_delta = .95, max_treedepth = 12, refresh = 1000,
    output_dir = outdir
  )
  assess <- function(fit, depth) {
    s <- fit$summary()
    np <- fit$sampler_diagnostics(format = "draws_matrix")
    q <- data.frame(
      id = id, rho = rho, rhat_max = max(s$rhat, na.rm = TRUE), ess_bulk_min = min(s$ess_bulk, na.rm = TRUE),
      ess_tail_min = min(s$ess_tail, na.rm = TRUE), divergences = sum(np[, "divergent__"]),
      max_depth_hits = sum(np[, "treedepth__"] >= depth)
    )
    q$pass <- with(q, rhat_max < 1.01 & ess_bulk_min >= 400 & ess_tail_min >= 400 & divergences == 0 & max_depth_hits == 0)
    q
  }
  dg <- assess(fit, 12)
  if (!dg$pass) {
    fit$save_object(file.path(outdir, "initial_fit.rds"))
    fit <- mod$sample(
      data = dat, chains = 4, parallel_chains = as.integer(Sys.getenv("ENDO_CHAIN_CORES", "4")),
      iter_warmup = 2000, iter_sampling = 4000, seed = SEED, adapt_delta = .999, max_treedepth = 15, refresh = 2000,
      output_dir = outdir
    )
    dg <- assess(fit, 15)
  }
  write_csv(dg, paste0("08_qa/selection_", id, "_diagnostics.csv"))
  if (!dg$pass) stop(paste("Unresolved selection diagnostics", id))
  fit$save_object(file.path(outdir, "fit.rds"))
  draws <- as.numeric(fit$draws("delta", format = "draws_matrix"))
  saveRDS(draws, file.path(ROOT, "04_results", paste0("selection_", id, "_draws.rds")))
  write_csv(cbind(id = id, rho = rho, summarize_effect(draws, 7)), paste0("04_results/selection_", id, "_summary.csv"))
}
if (Sys.getenv("ENDO_LOAD_ONLY", "") != "1") {
  if (Sys.getenv("ENDO_MODEL_WORKER", "") == "1") for (id in commandArgs(trailingOnly = TRUE)) run_job(id) else {
    writeLines(stan, fn)
    cmdstanr::cmdstan_model(fn, quiet = TRUE)
    saveRDS(list(data = dat, formula = f, columns = colnames(X), rho = rho_grid), file.path(ROOT, "02_data/derived/selection_calendar_data.rds"))
    write_csv(data.frame(id = ids), "08_qa/selection_manifest.csv")
    Sys.setenv(ENDO_MANIFEST = "selection_manifest.csv", ENDO_SCRIPT = "08_selection_model.R", ENDO_BATCH = "selection")
    s <- system2("python3", shQuote(file.path(CODE_DIR, "04_model_scheduler.py")))
    if (s != 0) stop("Selection-model worker failed.")
  }
}
