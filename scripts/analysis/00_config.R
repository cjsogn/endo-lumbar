# Purpose: Shared analysis settings, priors and effect-summary functions.
# Inputs: ENDO_CODE_DIR, ENDO_WORK_DIR and ENDO_LUMBAR_RAW_DATA environment variables.
# Outputs: Shared objects and functions; private workspace directories.

# All data, model objects, numerical output and logs stay outside this checkout.
suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(haven)
})
CODE_DIR <- normalizePath(Sys.getenv("ENDO_CODE_DIR"), mustWork = TRUE)
REPO_DIR <- normalizePath(file.path(CODE_DIR, "../.."), mustWork = TRUE)
work <- Sys.getenv("ENDO_WORK_DIR")
if (!nzchar(work)) stop("Set ENDO_WORK_DIR to a private directory outside the repository.")
work <- path.expand(work)
if (!grepl("^/|^[A-Za-z]:[/\\\\]", work)) stop("ENDO_WORK_DIR must be an absolute path.")
# Resolve symlinked ancestors before creating anything.
ancestor <- work
while (!dir.exists(ancestor)) {
  next_parent <- dirname(ancestor)
  if (identical(next_parent, ancestor)) stop("Cannot resolve the work directory.")
  ancestor <- next_parent
}
resolved <- normalizePath(ancestor, mustWork = TRUE)
if (identical(resolved, REPO_DIR) || startsWith(resolved, paste0(REPO_DIR, "/")))
  stop("The work directory must be outside the code repository.")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
ROOT <- normalizePath(work, mustWork = TRUE)
if (identical(ROOT, REPO_DIR) || startsWith(ROOT, paste0(REPO_DIR, "/")))
  stop("The work directory must be outside the code repository.")
for (p in c("02_data/source", "02_data/derived", "03_models/stan_cache", "04_results", "08_qa", "10_logs"))
  dir.create(file.path(ROOT, p), recursive = TRUE, showWarnings = FALSE)
RAW_DATA <- Sys.getenv("ENDO_LUMBAR_RAW_DATA")
options(mc.cores = 4L, brms.backend = "cmdstanr")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
SEED <- 20260204L
COVS <- c(
  "age_z", "sex", "bmi_z", "smoking", "education", "employed_baseline",
  "sick_leave", "disability", "analgesic_baseline", "odi_baseline_z", "eq5d_baseline_z",
  "nrs_back_baseline_z", "nrs_leg_baseline_z", "symptom_duration_back", "symptom_duration_leg",
  "motor_deficit", "asa_cat", "depression_anxiety", "chronic_pain", "spondylolisthesis",
  "scoliosis", "prior_surgery_any", "n_prior_surgeries_z", "multilevel",
  "prolapse_intraforaminal", "prolapse_extralateral", "stenosis_central"
)
CONT <- c(
  "age", "bmi", "odi_baseline", "eq5d_baseline", "nrs_back_baseline",
  "nrs_leg_baseline", "n_prior_surgeries"
)
RHS_BASE <- paste(c("treatment", COVS, "ct1", "ct2"), collapse = " + ")
RHS <- paste(RHS_BASE, "+ treatment:ct1 + treatment:ct2")
RHS_PRIMARY <- RHS
write_csv <- function(x, path) write.csv(x, file.path(ROOT, path), row.names = FALSE, na = "")
standardize <- function(d) {
  for (v in CONT) d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
  d
}
pri_zib <- function(treatment_sd = 1) c(
  set_prior(sprintf("normal(0, %s)", treatment_sd), class = "b", coef = "treatmentELD"),
  prior(normal(0, .5), class = b), prior(normal(0, 3), class = Intercept),
  prior(gamma(2, .1), class = phi), prior(normal(0, 1.5), class = Intercept, dpar = zi),
  prior(normal(0, 1), class = b, dpar = zi)
)
pri_binary <- c(
  prior(normal(0, 1), class = b, coef = treatmentELD),
  prior(normal(0, .5), class = b), prior(normal(0, 5), class = Intercept)
)
pri_eq <- c(
  prior(normal(0, .5), class = b, coef = treatmentELD),
  prior(normal(0, .5), class = b, coef = odi_baseline_z), prior(normal(0, 2), class = b),
  prior(normal(.5, .5), class = Intercept), prior(student_t(3, 0, .3), class = sigma)
)
pri_odi_gaussian <- c(
  prior(normal(0, 10), class = b, coef = treatmentELD),
  prior(normal(1, .5), class = b, coef = odi_baseline_z), prior(normal(0, 2), class = b),
  prior(normal(30, 20), class = Intercept), prior(student_t(3, 0, 15), class = sigma)
)
summarize_effect <- function(x, margin = NA_real_) data.frame(
  mean = mean(x), median = median(x), sd = sd(x), lower = unname(quantile(x, .025)),
  upper = unname(quantile(x, .975)), p_superiority = mean(x > 0),
  margin = margin, p_ni = if (is.na(margin)) NA_real_ else mean(x > -margin), draws = length(x)
)
gcomp <- function(fit, d, lower_better = TRUE, scale_factor = 1) {
  # All retained posterior draws, and all procedures in the stated target population.
  out <- list()
  for (a in c("MSD", "ELD")) {
    nd <- d
    nd$treatment <- factor(a, levels = c("MSD", "ELD"))
    p <- posterior_epred(fit, newdata = nd, allow_new_levels = FALSE)
    stopifnot(length(dim(p)) == 2L, ncol(p) == nrow(d))
    out[[a]] <- rowMeans(p) * scale_factor
  }
  out$delta <- if (lower_better) out$MSD - out$ELD else out$ELD - out$MSD
  as.data.frame(out)
}
ordinal_calendar_contrast <- function(fit, d) {
  # For cumulative-logit models, smaller eta implies shorter stay. Average
  # individual conditional log odds ratios over the stated calendar distribution.
  # Exponentiation gives a geometric mean conditional OR, not a marginal OR.
  eta <- list()
  for (a in c("MSD", "ELD")) {
    nd <- d
    nd$treatment <- factor(a, levels = c("MSD", "ELD"))
    eta[[a]] <- posterior_linpred(fit, newdata = nd, incl_thres = FALSE)
    stopifnot(length(dim(eta[[a]])) == 2L, ncol(eta[[a]]) == nrow(d))
  }
  delta <- rowMeans(eta$MSD - eta$ELD)
  # Independent coefficient calculation checks direction and both interactions.
  dr <- as_draws_df(fit)
  expected <- -dr$b_treatmentELD
  for (v in c("ct1", "ct2")) {
    nm <- paste0("b_treatmentELD:", v)
    if (nm %in% names(dr)) expected <- expected - dr[[nm]] * mean(d[[v]])
  }
  stopifnot(max(abs(delta - expected)) < 1e-10)
  data.frame(delta = delta)
}
diagnostics <- function(fit, id, max_depth = 12L) {
  s <- posterior::summarise_draws(posterior::as_draws_array(fit))
  # Constant fixed parameters have undefined R-hat and ESS, not failed convergence.
  relevant <- s[is.finite(s$rhat) & is.finite(s$ess_bulk) & is.finite(s$ess_tail), ]
  np <- nuts_params(fit)
  r <- data.frame(
    id = id, rhat_max = max(relevant$rhat), ess_bulk_min = min(relevant$ess_bulk),
    ess_tail_min = min(relevant$ess_tail), divergences = sum(np$Value[np$Parameter == "divergent__"]),
    max_depth_hits = sum(np$Value[np$Parameter == "treedepth__"] >= max_depth)
  )
  r$pass <- with(r, rhat_max < 1.01 & ess_bulk_min >= 400 & ess_tail_min >= 400 &
    divergences == 0 & max_depth_hits == 0)
  write_csv(as.data.frame(s), paste0("08_qa/diagnostics_", id, "_parameters.csv"))
  r
}
