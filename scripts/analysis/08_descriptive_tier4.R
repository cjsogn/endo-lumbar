# =============================================================================
# ENDO-LUMBAR: 08 Descriptive Outcomes (Tier 4)
# Operating time and Negative control
# SAP Sections 6.5, 20.1
# =============================================================================

source("/Users/cjsogn/endo_studies/lumbar/analysis/scripts/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
df_disc_12m <- readRDS(file.path(paths$data_clean, "df_disc_12m_eligible.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Tier 4: Descriptive Outcomes ===\n")

# Prepare data (standardize_covs and cov_string from 00_config.R)
df_disc <- standardize_covs(df_disc)

# =============================================================================
# 1. OPERATING TIME (SAP Section 6.5)
# =============================================================================

cat("\n--- Operating Time ---\n")
cat("Note: Descriptive only. No directional hypothesis (endo may be longer).\n")

# Descriptive statistics
op_time_desc <- df_disc %>%
  group_by(treatment) %>%
  dplyr::summarise(
    n = sum(!is.na(operating_time)),
    mean = mean(operating_time, na.rm = TRUE),
    sd = sd(operating_time, na.rm = TRUE),
    median = median(operating_time, na.rm = TRUE),
    q25 = quantile(operating_time, 0.25, na.rm = TRUE),
    q75 = quantile(operating_time, 0.75, na.rm = TRUE),
    min = min(operating_time, na.rm = TRUE),
    max = max(operating_time, na.rm = TRUE),
    .groups = "drop"
  )
print(op_time_desc)

# Bayesian model for adjusted difference
fit_optime <- brm(
  bf(as.formula(paste("operating_time | mi() ~ treatment +", cov_string))),
  data = df_disc,
  family = gaussian(),
  prior = c(
    prior(normal(0, 30), class = "b", coef = "treatmentELD"),
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 30), class = "sigma"),
    prior(normal(60, 30), class = "Intercept")
  ),
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.95),
  file = file.path(paths$models, "fit_tier4_optime"),
  file_refit = "on_change"
)

# G-computation
ate_optime <- compute_gcomp_ate(fit_optime, df_disc, "treatment", "continuous",
                                 lower_is_better = TRUE)
ate_optime_summary <- summarize_ate(ate_optime$ate)
cat(sprintf("\nAdjusted difference (MSD - ELD): %.1f min (95%% CrI: [%.1f, %.1f])\n",
            ate_optime_summary$mean, ate_optime_summary$cri_lo, ate_optime_summary$cri_hi))

# =============================================================================
# 2. NEGATIVE CONTROL: EQ-5D ANXIETY/DEPRESSION (SAP Section 20.1)
# =============================================================================

cat("\n--- Negative Control: EQ-5D Anxiety/Depression ---\n")
cat("Purpose: Diagnostic for residual confounding.\n")
cat("Expectation: No treatment effect (surgical technique should not affect anxiety/depression).\n")

# 3-month analysis
cat("\n  3-month:\n")
cat(sprintf("  Observed: %d/%d (%.1f%%)\n",
            sum(!is.na(df_disc$eq5d_anxiety_3m)), nrow(df_disc),
            100 * mean(!is.na(df_disc$eq5d_anxiety_3m))))

fit_negcontrol_3m <- brm(
  bf(as.formula(paste("eq5d_anxiety_3m | mi() ~ treatment +", cov_string))),
  data = df_disc,
  family = gaussian(),
  prior = c(
    prior(normal(0, 2), class = "b", coef = "treatmentELD"),
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 3), class = "sigma"),
    prior(normal(2, 2), class = "Intercept")
  ),
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.95),
  file = file.path(paths$models, "fit_tier4_negcontrol_3m"),
  file_refit = "on_change"
)

ate_neg3m <- compute_gcomp_ate(fit_negcontrol_3m, df_disc, "treatment", "continuous",
                                lower_is_better = TRUE)
ate_neg3m_summary <- summarize_ate(ate_neg3m$ate)
p_nonzero_3m <- 2 * min(mean(ate_neg3m$ate > 0), mean(ate_neg3m$ate < 0))

cat(sprintf("  ATE: %.3f (95%% CrI: [%.3f, %.3f])\n",
            ate_neg3m_summary$mean, ate_neg3m_summary$cri_lo, ate_neg3m_summary$cri_hi))
cat(sprintf("  P(beta != 0): %.3f\n", 1 - p_nonzero_3m))

# Check if 95% CrI excludes 0 (concerns about confounding)
if (ate_neg3m_summary$cri_lo > 0 | ate_neg3m_summary$cri_hi < 0) {
  cat("  WARNING: 95% CrI excludes 0. Raises concern about residual confounding.\n")
  cat("  Caveat: indirect pathway possible (faster recovery -> less anxiety).\n")
} else {
  cat("  95% CrI includes 0. No evidence of residual confounding from this test.\n")
}

# 12-month analysis (restricted to 12m-eligible patients)
cat("\n  12-month [12m-eligible subset]:\n")

# Standardize covariates in 12m-eligible dataset
df_disc_12m <- standardize_covs(df_disc_12m)

cat(sprintf("  N (12m-eligible): %d, Observed: %d (%.1f%%)\n",
            nrow(df_disc_12m),
            sum(!is.na(df_disc_12m$eq5d_anxiety_12m)),
            100 * mean(!is.na(df_disc_12m$eq5d_anxiety_12m))))

fit_negcontrol_12m <- brm(
  bf(as.formula(paste("eq5d_anxiety_12m | mi() ~ treatment +", cov_string))),
  data = df_disc_12m,
  family = gaussian(),
  prior = c(
    prior(normal(0, 2), class = "b", coef = "treatmentELD"),
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 3), class = "sigma"),
    prior(normal(2, 2), class = "Intercept")
  ),
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.95),
  file = file.path(paths$models, "fit_tier4_negcontrol_12m"),
  file_refit = "on_change"
)

ate_neg12m <- compute_gcomp_ate(fit_negcontrol_12m, df_disc_12m, "treatment", "continuous",
                                 lower_is_better = TRUE)
ate_neg12m_summary <- summarize_ate(ate_neg12m$ate)

cat(sprintf("  ATE: %.3f (95%% CrI: [%.3f, %.3f])\n",
            ate_neg12m_summary$mean, ate_neg12m_summary$cri_lo, ate_neg12m_summary$cri_hi))

# =============================================================================
# SAVE
# =============================================================================

tier4_results <- list(
  operating_time = list(
    descriptive = op_time_desc,
    ate_summary = ate_optime_summary,
    ate_draws = ate_optime
  ),
  negative_control_3m = list(
    ate_summary = ate_neg3m_summary,
    ate_draws = ate_neg3m,
    cri_excludes_zero = ate_neg3m_summary$cri_lo > 0 | ate_neg3m_summary$cri_hi < 0
  ),
  negative_control_12m = list(
    ate_summary = ate_neg12m_summary,
    ate_draws = ate_neg12m,
    cri_excludes_zero = ate_neg12m_summary$cri_lo > 0 | ate_neg12m_summary$cri_hi < 0
  )
)
saveRDS(tier4_results, file.path(paths$output, "tier4_results.rds"))

cat("\nTier 4 analysis complete.\n")
