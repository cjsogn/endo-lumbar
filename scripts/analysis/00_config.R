# =============================================================================
# ENDO-LUMBAR: Configuration
# Endoscopic vs microsurgical lumbar discectomy - NORspine registry study
#
# This file is sourced by every analysis script. It loads packages, sets the
# project root via `here::here()` (so the pipeline works from any clone of
# the repository), and defines shared parameters, priors, and helpers.
#
# The raw NORspine data file is not part of the repository and must be
# requested from the registry. Set the environment variable
# ENDO_LUMBAR_RAW_DATA to the path of the SPSS export, or edit
# `paths$data_raw` directly below.
# =============================================================================

# --- Packages ----------------------------------------------------------------
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

suppressPackageStartupMessages({
  library(here)         # Project-relative paths
  library(haven)        # SPSS data import
  library(tidyverse)    # Data wrangling and visualization
  library(brms)         # Bayesian regression models (Stan backend)
  if (requireNamespace("cmdstanr", quietly = TRUE) &&
      !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))) {
    library(cmdstanr)
  }
  library(rms)          # Restricted cubic splines
  library(loo)          # LOO-CV, Pareto k, model comparison, stacking
  library(bayesplot)    # PPC, trace plots, posterior visualization
  library(posterior)    # Posterior summaries, convergence diagnostics
  library(grf)          # Causal forests (generalized random forests)
  library(EValue)       # E-values for unmeasured confounding
  library(projpred)     # Projection predictive variable selection
  library(tableone)     # Table 1 generation
  library(cobalt)       # Covariate balance (Love plots, SMD)
  library(patchwork)    # Figure composition
  library(scales)       # Axis formatting
  library(gt)           # Publication tables
  library(gtsummary)    # Summary tables
})

# --- Paths -------------------------------------------------------------------
# `project_root` is resolved from the repository root via the `here` package,
# so the pipeline works from any clone of the repo.
project_root <- here::here()

# Raw NORspine export: set via environment variable or edit directly.
data_raw_path <- Sys.getenv(
  "ENDO_LUMBAR_RAW_DATA",
  unset = file.path(project_root, "data", "raw_norspine_export.sav")
)

paths <- list(
  data_raw    = data_raw_path,
  data_clean  = file.path(project_root, "data"),
  output      = project_root,
  tables      = file.path(project_root, "tables"),
  figures     = file.path(project_root, "figures"),
  models      = file.path(project_root, "results"),
  results     = file.path(project_root, "results"),
  diagnostics = file.path(project_root, "diagnostics")
)

# --- Computational settings --------------------------------------------------
n_cores <- 14L
options(mc.cores = n_cores)
# Use cmdstanr if available, otherwise fall back to rstan
if (requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))) {
  options(brms.backend = "cmdstanr")
  cat("  Backend: cmdstanr\n")
} else {
  options(brms.backend = "rstan")
  cat("  Backend: rstan (cmdstanr not available)\n")
}
options(brms.file_refit = "on_change")

# MCMC settings
mcmc_settings <- list(
  chains   = 4L,
  iter     = 2000L,
  warmup   = 1000L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  seed     = 20260204L
)

# Convergence thresholds
convergence_thresholds <- list(
  rhat_mandatory    = 1.01,
  rhat_warning      = 1.05,
  ess_bulk_mandatory = 400,
  ess_tail_mandatory = 400,
  ess_warning        = 200,
  max_divergent      = 0
)

# --- Non-inferiority margins --------------------------------------------------
ni_margins <- list(
  # Tier 1 (primary)
  odi         = 7,      # ODI points (0-100 scale)
  # Tier 2 (secondary continuous)
  nrs_pain    = 1.0,    # NRS 0-10 scale
  eq5d        = 0.05,   # EQ-5D index points
  # Tier 2 (secondary binary, risk difference scale)
  responder   = 0.10,   # 10 percentage points
  rtw         = 0.10,   # 10 percentage points
  analgesic   = 0.10,   # 10 percentage points
  satisfaction = 0.10,  # 10 percentage points
  gpe         = 0.10    # 10 percentage points
)

# NI probability threshold
ni_threshold <- 0.95

# --- Effect direction convention -----------------------------------------------
# Delta > 0 = endoscopic superior for ALL outcomes
# For "lower is better" (ODI, NRS, LOS, complications): Delta = mu_MSD - mu_ELD
# For "higher is better" (EQ5D, responder, RTW, day surgery): Delta = mu_ELD - mu_MSD
# NI: P(Delta > -margin | data) > 0.95
# Superiority: P(Delta > 0 | data) > 0.95

# --- Prior distributions -----------------------------------------------------
priors_reference <- list(
  treatment_mean = 0,
  treatment_sd   = 10,   # Reference prior
  baseline_mean  = 1,
  baseline_sd    = 0.5,
  covariate_mean = 0,
  covariate_sd   = 2,
  sigma_df       = 3,
  sigma_mu       = 0,
  sigma_sigma    = 15
)

priors_skeptical <- list(treatment_sd = 4)
priors_diffuse   <- list(treatment_sd = 25)

# Binary outcome priors
priors_binary <- list(
  treatment_mean = 0,
  treatment_sd   = 1,    # on log-odds scale
  covariate_mean = 0,
  covariate_sd   = 2
)

# --- Complication rate gating thresholds --------------------------------------
rate_gating <- list(
  formal_test = 0.10,   # >= 10%: formal superiority test
  caution     = 0.05,   # 5-10%: test with caution
  descriptive = 0.05    # < 5%: descriptive only
)

# --- Calendar time RCS specification ------------------------------------------
rcs_knots <- 3  # 3 knots at 10th, 50th, 90th percentiles

# --- Pattern-mixture delta grid -----------------------------------------------
delta_grid <- seq(-8, 8, by = 2)

# --- Subgroup cut-points ------------------------------------------------------
subgroup_cuts <- list(
  age      = 50,
  odi_base = 40
  # symptom_duration = median (computed from data)
  # levels = single vs multilevel (computed from data)
)

# --- Causal forest settings ---------------------------------------------------
cf_settings <- list(
  num.trees     = 4000L,
  min.node.size = 5L,
  honesty       = TRUE,
  tune_forest   = TRUE
)

# --- Figure settings ---------------------------------------------------------
theme_pub <- theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70", linewidth = 0.5),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )
theme_set(theme_pub)

# Color palette for treatment groups
tx_colors <- c("ELD" = "#2166AC", "MSD" = "#B2182B")
tx_fills  <- c("ELD" = "#92C5DE", "MSD" = "#F4A582")

fig_width  <- 7
fig_height <- 5
fig_dpi    <- 300

# --- Helper functions --------------------------------------------------------

#' Compute G-computation ATE from brms posterior
#' @param fit brms model fit
#' @param newdata data frame with all patients
#' @param treatment_var name of treatment variable
#' @param outcome_type "continuous" or "binary"
#' @param lower_is_better logical, TRUE if lower outcome = better
#' @return data frame with ATE posterior draws
compute_gcomp_ate <- function(fit, newdata, treatment_var = "treatment",
                               outcome_type = "continuous",
                               lower_is_better = TRUE,
                               scale_factor = 1) {

  # Create counterfactual datasets
  nd_eld <- nd_msd <- newdata
  nd_eld[[treatment_var]] <- "ELD"
  nd_msd[[treatment_var]] <- "MSD"

  if (outcome_type == "continuous") {
    # Predicted values under each treatment (posterior draws)
    pred_eld <- posterior_epred(fit, newdata = nd_eld, allow_new_levels = TRUE)
    pred_msd <- posterior_epred(fit, newdata = nd_msd, allow_new_levels = TRUE)

    # Mean over patients for each draw
    # scale_factor converts (0,1) predictions back to original scale
    # (100 for ODI, 10 for NRS, 1 for Gaussian outcomes)
    mu_eld <- rowMeans(pred_eld) * scale_factor
    mu_msd <- rowMeans(pred_msd) * scale_factor

    # Delta convention: positive = ELD superior
    if (lower_is_better) {
      ate <- mu_msd - mu_eld  # MSD - ELD (positive if ELD lower = better)
    } else {
      ate <- mu_eld - mu_msd  # ELD - MSD (positive if ELD higher = better)
    }

  } else if (outcome_type == "binary") {
    pred_eld <- posterior_epred(fit, newdata = nd_eld, allow_new_levels = TRUE)
    pred_msd <- posterior_epred(fit, newdata = nd_msd, allow_new_levels = TRUE)

    p_eld <- rowMeans(pred_eld)
    p_msd <- rowMeans(pred_msd)

    if (lower_is_better) {
      ate <- p_msd - p_eld  # Risk difference: positive if ELD has lower risk
    } else {
      ate <- p_eld - p_msd  # Risk difference: positive if ELD has higher rate
    }
  }

  tibble(ate = ate)
}

#' Summarize ATE posterior with NI and superiority probabilities
#' @param ate_draws vector of ATE posterior draws
#' @param ni_margin non-inferiority margin (positive value)
#' @return named list with summary statistics
summarize_ate <- function(ate_draws, ni_margin = NULL) {
  result <- list(
    mean     = mean(ate_draws),
    median   = median(ate_draws),
    sd       = sd(ate_draws),
    cri_lo   = quantile(ate_draws, 0.025),
    cri_hi   = quantile(ate_draws, 0.975),
    p_superiority = mean(ate_draws > 0)
  )

  if (!is.null(ni_margin)) {
    result$p_ni <- mean(ate_draws > -ni_margin)
    result$ni_margin <- ni_margin
    result$ni_conclusion <- result$p_ni > ni_threshold
  }

  result
}

#' Check MCMC convergence against SAP thresholds
#' @param fit brms model fit
#' @return list with convergence status and details
check_convergence <- function(fit) {
  rhats <- brms::rhat(fit)
  np <- nuts_params(fit)

  s <- posterior::summarise_draws(posterior::as_draws(fit))

  list(
    rhat_max       = max(rhats, na.rm = TRUE),
    rhat_ok        = all(rhats < convergence_thresholds$rhat_mandatory, na.rm = TRUE),
    ess_bulk_min   = min(s$ess_bulk, na.rm = TRUE),
    ess_bulk_ok    = all(s$ess_bulk > convergence_thresholds$ess_bulk_mandatory, na.rm = TRUE),
    ess_tail_min   = min(s$ess_tail, na.rm = TRUE),
    ess_tail_ok    = all(s$ess_tail > convergence_thresholds$ess_tail_mandatory, na.rm = TRUE),
    n_divergent    = sum(subset(np, Parameter == "divergent__")$Value),
    divergent_ok   = sum(subset(np, Parameter == "divergent__")$Value) == 0,
    all_ok         = all(rhats < convergence_thresholds$rhat_mandatory, na.rm = TRUE) &&
                     all(s$ess_bulk > convergence_thresholds$ess_bulk_mandatory, na.rm = TRUE) &&
                     all(s$ess_tail > convergence_thresholds$ess_tail_mandatory, na.rm = TRUE) &&
                     sum(subset(np, Parameter == "divergent__")$Value) == 0
  )
}

#' Compute EQ-5D-5L index value using Norwegian value set
#' @param mobility integer 1-5
#' @param selfcare integer 1-5
#' @param usual integer 1-5
#' @param pain integer 1-5
#' @param anxiety integer 1-5
#' @return numeric EQ-5D index value
compute_eq5d_index <- function(mobility, selfcare, usual, pain, anxiety) {
  # Use the crosswalk value set (EQ-5D-5L -> 3L crosswalk)
  # If eq5d package is available, use it; otherwise use the pre-computed values
  # For now, the registry provides pre-computed EQ5DV3 values
  # This function is a placeholder; we use the registry-provided index
  NA_real_
}

#' Save figure with standard settings
save_fig <- function(plot, filename, width = fig_width, height = fig_height,
                     dpi = fig_dpi, path = paths$figures) {
  ggsave(
    filename = filename,
    plot = plot,
    path = path,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

# =============================================================================
# CENTRALIZED COVARIATE SPECIFICATION
# =============================================================================

# Full covariate set for outcome models (z-scored continuous, original binary/factor)
cov_string <- paste(
  "age_z + sex + bmi_z + smoking +",
  "education + employed_baseline + sick_leave + disability + analgesic_baseline +",
  "odi_baseline_z + eq5d_baseline_z +",
  "nrs_back_baseline_z + nrs_leg_baseline_z +",
  "symptom_duration_back + symptom_duration_leg + motor_deficit +",
  "asa_cat + depression_anxiety + chronic_pain +",
  "spondylolisthesis + scoliosis +",
  "prior_surgery_any + n_prior_surgeries_z + multilevel +",
  "prolapse_intraforaminal + prolapse_extralateral + stenosis_central"
)

# Covariate names in original (un-z-scored) form for falsification, Table 1, etc.
all_model_covs <- c(
  "age", "sex", "bmi", "smoking", "education", "employed_baseline",
  "sick_leave", "disability", "analgesic_baseline",
  "odi_baseline", "eq5d_baseline",
  "nrs_back_baseline", "nrs_leg_baseline",
  "symptom_duration_back", "symptom_duration_leg", "motor_deficit",
  "asa_cat", "depression_anxiety", "chronic_pain",
  "spondylolisthesis", "scoliosis",
  "prior_surgery_any", "n_prior_surgeries", "multilevel",
  "prolapse_intraforaminal", "prolapse_extralateral", "stenosis_central"
)

# Restricted covariate set for sensitivity analysis
cov_string_restricted <- "age_z + sex + odi_baseline_z + prior_surgery_any"

# Covariates for propensity score model (same set, un-z-scored)
ps_covariates <- c(
  "age", "sex", "bmi", "smoking", "education", "employed_baseline",
  "sick_leave", "disability", "analgesic_baseline",
  "odi_baseline", "eq5d_baseline",
  "nrs_back_baseline", "nrs_leg_baseline",
  "symptom_duration_back", "symptom_duration_leg", "motor_deficit",
  "asa_cat", "depression_anxiety", "chronic_pain",
  "spondylolisthesis", "scoliosis",
  "prior_surgery_any", "n_prior_surgeries", "multilevel",
  "prolapse_intraforaminal", "prolapse_extralateral", "stenosis_central"
)

# Continuous covariates that need z-scoring
covs_to_zscore <- c("age", "bmi", "odi_baseline", "eq5d_baseline",
                     "nrs_back_baseline", "nrs_leg_baseline", "n_prior_surgeries")

#' Standardize continuous covariates for MCMC sampling
#' @param d data frame
#' @return data frame with _z columns added
standardize_covs <- function(d) {
  d %>% mutate(
    age_z = scale(age)[,1],
    bmi_z = scale(bmi)[,1],
    odi_baseline_z = scale(odi_baseline)[,1],
    eq5d_baseline_z = scale(eq5d_baseline)[,1],
    nrs_back_baseline_z = scale(nrs_back_baseline)[,1],
    nrs_leg_baseline_z = scale(nrs_leg_baseline)[,1],
    n_prior_surgeries_z = scale(n_prior_surgeries)[,1]
  )
}

# Standard priors for continuous outcomes
priors_continuous <- c(
  prior(normal(0, 10), class = "b", coef = "treatmentELD"),
  prior(normal(1, 0.5), class = "b", coef = "odi_baseline_z"),
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 15), class = "sigma"),
  prior(normal(30, 20), class = "Intercept")
)

# Standard priors for binary outcomes (Normal(0, 0.5) for covariate regularization)
priors_binary_standard <- c(
  prior(normal(0, 1), class = "b", coef = "treatmentELD"),
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 5), class = "Intercept")
)

# Standard priors for beta regression (logit link)
priors_beta <- c(
  prior(normal(0, 1), class = "b", coef = "treatmentELD"),
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 3), class = "Intercept"),
  prior(gamma(2, 0.1), class = "phi")
)

# Prior sensitivity variants for beta (logit scale)
priors_beta_skeptical_sd  <- 0.5   # N(0, 0.5) on logit
priors_beta_reference_sd  <- 1.0   # N(0, 1) on logit
priors_beta_diffuse_sd    <- 2.0   # N(0, 2) on logit

# Standard priors for zero-inflated beta regression (ZIB, primary model)
priors_zib <- c(
  prior(normal(0, 1), class = "b", coef = "treatmentELD"),
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 3), class = "Intercept"),
  prior(gamma(2, 0.1), class = "phi"),
  prior(normal(0, 1.5), class = "Intercept", dpar = "zi"),
  prior(normal(0, 1), class = "b", dpar = "zi")
)

# Prior sensitivity variants for ZIB (logit scale, mu component)
priors_zib_skeptical_sd  <- 0.5   # N(0, 0.5) on logit
priors_zib_reference_sd  <- 1.0   # N(0, 1) on logit
priors_zib_diffuse_sd    <- 2.0   # N(0, 2) on logit

#' Transform a bounded outcome to (0,1) for beta regression
#' Uses Smithson & Verkuilen (2006): y' = (y * (n-1) + 0.5) / n
#' Preserves NAs for mi() compatibility
#' @param y numeric vector (may contain NAs)
#' @param upper upper bound of the scale (100 for ODI, 10 for NRS)
#' @return transformed vector on (0,1) open interval
transform_for_beta <- function(y, upper) {
  y_01 <- y / upper
  n <- sum(!is.na(y_01))
  (y_01 * (n - 1) + 0.5) / n
}

#' Transform a bounded outcome for zero-inflated beta regression (ZIB)
#' Values at 0 are preserved (handled by the zero-inflation component).
#' Values at the upper bound are clamped to 1 - 1e-6.
#' Preserves NAs for mi() compatibility.
#' @param y numeric vector (may contain NAs)
#' @param upper upper bound of the scale (100 for ODI, 10 for NRS)
#' @return transformed vector on [0, 1) interval
transform_for_zib <- function(y, upper) {
  y_01 <- y / upper
  pmin(y_01, 1 - 1e-6)
}

cat("Configuration loaded successfully.\n")
cat(sprintf("  Data: %s\n", paths$data_raw))
cat(sprintf("  Cores: %d\n", n_cores))
cat(sprintf("  MCMC: %d chains, %d iterations (%d warmup)\n",
            mcmc_settings$chains, mcmc_settings$iter, mcmc_settings$warmup))
cat(sprintf("  Primary NI margin: %d ODI points\n", ni_margins$odi))
cat(sprintf("  NI threshold: P > %.2f\n", ni_threshold))
cat(sprintf("  Model covariates: %d variables\n", length(all_model_covs)))
