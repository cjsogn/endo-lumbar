# =============================================================================
# ENDO-LUMBAR: 11 Model Specification Sensitivity
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Model Specification Sensitivity ===\n")

# Prepare data (standardize_covs, cov_string, cov_string_restricted from 00_config.R)
df_disc <- standardize_covs(df_disc)

# ZIB transformation for ODI (primary model)
df_disc$odi_3m_zib <- transform_for_zib(df_disc$odi_3m, upper = 100)
# Beta (SV) transformation for sensitivity analysis
df_disc$odi_3m_beta <- transform_for_beta(df_disc$odi_3m, upper = 100)

cov_string_full <- cov_string  # alias for clarity in this script

# Preserve full dataset for G-computation
df_disc_full <- df_disc

# Complete cases for ZIB models (mi() not supported for zero_inflated_beta)
df_cc <- df_disc %>% filter(!is.na(odi_3m))
df_cc$odi_3m_zib <- transform_for_zib(df_cc$odi_3m, upper = 100)
cat(sprintf("Complete cases for ZIB sensitivity: %d\n", nrow(df_cc)))

sensitivity_results <- list()

# =============================================================================
# Gaussian Likelihood
# =============================================================================

cat("\n--- Gaussian Likelihood ---\n")

fit_gaussian <- brm(
  bf(as.formula(paste("odi_3m | mi() ~ treatment +", cov_string_full))),
  data = df_disc,
  family = gaussian(),
  prior = priors_continuous,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta,
                 max_treedepth = mcmc_settings$max_treedepth),
  file = file.path(paths$models, "fit_sensitivity_gaussian"),
  file_refit = "on_change"
)

ate_gauss <- compute_gcomp_ate(fit_gaussian, df_disc, "treatment", "continuous",
                                lower_is_better = TRUE)
ate_gauss_summary <- summarize_ate(ate_gauss$ate, ni_margin = ni_margins$odi)
cat(sprintf("  ATE: %.2f (95%% CrI: [%.2f, %.2f]), P(NI): %.4f\n",
            ate_gauss_summary$mean, ate_gauss_summary$cri_lo, ate_gauss_summary$cri_hi,
            ate_gauss_summary$p_ni))

sensitivity_results$gaussian <- list(
  label = "Gaussian likelihood (original SAP)",
  ate_summary = ate_gauss_summary
)

# =============================================================================
# Beta Regression with SV Transform
# =============================================================================

cat("\n--- Beta Regression with SV Transform ---\n")

fit_beta_sv <- brm(
  bf(as.formula(paste("odi_3m_beta | mi() ~ treatment +", cov_string_full))),
  data = df_disc,
  family = Beta(),
  prior = priors_beta,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta,
                 max_treedepth = mcmc_settings$max_treedepth),
  file = file.path(paths$models, "fit_sensitivity_beta_sv"),
  file_refit = "on_change"
)

ate_beta <- compute_gcomp_ate(fit_beta_sv, df_disc, "treatment", "continuous",
                               lower_is_better = TRUE, scale_factor = 100)
ate_beta_summary <- summarize_ate(ate_beta$ate, ni_margin = ni_margins$odi)
cat(sprintf("  ATE: %.2f (95%% CrI: [%.2f, %.2f]), P(NI): %.4f\n",
            ate_beta_summary$mean, ate_beta_summary$cri_lo, ate_beta_summary$cri_hi,
            ate_beta_summary$p_ni))

sensitivity_results$beta_sv <- list(
  label = "Beta regression (SV transform)",
  ate_summary = ate_beta_summary
)

# =============================================================================
# Student-t Likelihood (estimated nu)
# =============================================================================

cat("\n--- Student-t Likelihood ---\n")

fit_studentt <- brm(
  bf(as.formula(paste("odi_3m | mi() ~ treatment +", cov_string_full))),
  data = df_disc,
  family = student(),
  prior = c(
    prior(normal(0, 10), class = "b", coef = "treatmentELD"),
    prior(normal(1, 0.5), class = "b", coef = "odi_baseline_z"),
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 15), class = "sigma"),
    prior(gamma(2, 0.1), class = "nu"),
    prior(normal(30, 20), class = "Intercept")
  ),
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = file.path(paths$models, "fit_sensitivity_studentt"),
  file_refit = "on_change"
)

ate_st <- compute_gcomp_ate(fit_studentt, df_disc, "treatment", "continuous",
                             lower_is_better = TRUE)
ate_st_summary <- summarize_ate(ate_st$ate, ni_margin = ni_margins$odi)
cat(sprintf("  ATE: %.2f (95%% CrI: [%.2f, %.2f]), P(NI): %.4f\n",
            ate_st_summary$mean, ate_st_summary$cri_lo, ate_st_summary$cri_hi,
            ate_st_summary$p_ni))

sensitivity_results$student_t <- list(
  label = "Student-t (nu = 5)",
  ate_summary = ate_st_summary
)

# =============================================================================
# Horseshoe prior sensitivity (ZIB model, all mu coefficients)
# =============================================================================

cat("\n--- Horseshoe Priors ---\n")

fit_hs <- brm(
  bf(as.formula(paste("odi_3m_zib ~ treatment +", cov_string_full)),
     zi ~ treatment + odi_baseline_z),
  data = df_cc,
  family = zero_inflated_beta(),
  prior = c(
    prior(horseshoe(df = 1), class = "b"),
    prior(normal(0, 3), class = "Intercept"),
    prior(gamma(2, 0.1), class = "phi"),
    prior(normal(0, 1.5), class = "Intercept", dpar = "zi"),
    prior(normal(0, 1), class = "b", dpar = "zi")
  ),
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  file = file.path(paths$models, "fit_sensitivity_horseshoe_zib"),
  file_refit = "on_change"
)

ate_hs <- compute_gcomp_ate(fit_hs, df_disc_full, "treatment", "continuous",
                             lower_is_better = TRUE, scale_factor = 100)
ate_hs_summary <- summarize_ate(ate_hs$ate, ni_margin = ni_margins$odi)
cat(sprintf("  ATE: %.2f (95%% CrI: [%.2f, %.2f]), P(NI): %.4f\n",
            ate_hs_summary$mean, ate_hs_summary$cri_lo, ate_hs_summary$cri_hi,
            ate_hs_summary$p_ni))

sensitivity_results$horseshoe <- list(
  label = "Horseshoe priors (all coefficients)",
  ate_summary = ate_hs_summary
)

# =============================================================================
# Restricted Covariate Set (ZIB model)
# =============================================================================

cat("\n--- Restricted Covariate Set ---\n")
cat("  Covariates: age, sex, baseline ODI, prior surgery\n")
# cov_string_restricted defined in 00_config.R

fit_restricted <- brm(
  bf(as.formula(paste("odi_3m_zib ~ treatment +", cov_string_restricted)),
     zi ~ treatment + odi_baseline_z),
  data = df_cc,
  family = zero_inflated_beta(),
  prior = priors_zib,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta),
  file = file.path(paths$models, "fit_sensitivity_restricted_zib"),
  file_refit = "on_change"
)


ate_restr <- compute_gcomp_ate(fit_restricted, df_disc_full, "treatment", "continuous",
                                lower_is_better = TRUE, scale_factor = 100)
ate_restr_summary <- summarize_ate(ate_restr$ate, ni_margin = ni_margins$odi)
cat(sprintf("  ATE: %.2f (95%% CrI: [%.2f, %.2f]), P(NI): %.4f\n",
            ate_restr_summary$mean, ate_restr_summary$cri_lo, ate_restr_summary$cri_hi,
            ate_restr_summary$p_ni))

sensitivity_results$restricted <- list(
  label = "Restricted covariate set",
  ate_summary = ate_restr_summary
)

# =============================================================================
# Pure Disc Herniation (no stenosis)
# =============================================================================

cat("\n--- Pure Disc Herniation (no stenosis) ---\n")

# Full pure disc subset for G-computation (all patients with covariate data)
df_pure_disc_full <- df_disc_full %>%
  filter(stenosis_central == 0 & stenosis_lateral == 0 & stenosis_foraminal == 0)

# Complete-case pure disc subset for model fitting
df_pure_disc <- df_cc %>%
  filter(stenosis_central == 0 & stenosis_lateral == 0 & stenosis_foraminal == 0)
cat(sprintf("  Pure disc herniation: %d complete / %d total (ELD=%d, MSD=%d)\n",
            nrow(df_pure_disc), nrow(df_pure_disc_full),
            sum(df_pure_disc$treatment == "ELD"),
            sum(df_pure_disc$treatment == "MSD")))

fit_pure <- brm(
  bf(as.formula(paste("odi_3m_zib ~ treatment +", cov_string_full)),
     zi ~ treatment + odi_baseline_z),
  data = df_pure_disc,
  family = zero_inflated_beta(),
  prior = priors_zib,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta),
  file = file.path(paths$models, "fit_sensitivity_pure_disc_zib"),
  file_refit = "on_change"
)


ate_pure <- compute_gcomp_ate(fit_pure, df_pure_disc_full, "treatment", "continuous",
                               lower_is_better = TRUE, scale_factor = 100)
ate_pure_summary <- summarize_ate(ate_pure$ate, ni_margin = ni_margins$odi)
cat(sprintf("  ATE: %.2f (95%% CrI: [%.2f, %.2f]), P(NI): %.4f\n",
            ate_pure_summary$mean, ate_pure_summary$cri_lo, ate_pure_summary$cri_hi,
            ate_pure_summary$p_ni))

sensitivity_results$pure_disc <- list(
  label = "Pure disc herniation (no stenosis)",
  ate_summary = ate_pure_summary
)

# =============================================================================
# Treatment x Calendar Time Interaction
# =============================================================================

cat("\n--- Treatment x Calendar Time Interaction ---\n")

# Add calendar_time_z to both df_cc (for fitting) and df_disc_full (for G-computation)
df_cc <- df_cc %>%
  mutate(calendar_time_z = scale(calendar_time)[,1])
df_disc_full <- df_disc_full %>%
  mutate(calendar_time_z = scale(calendar_time)[,1])

cov_interaction <- paste(cov_string_full, "+ calendar_time_z + treatment:calendar_time_z")

fit_caltime <- brm(
  bf(as.formula(paste("odi_3m_zib ~ treatment +", cov_interaction)),
     zi ~ treatment + odi_baseline_z),
  data = df_cc,
  family = zero_inflated_beta(),
  prior = c(
    prior(normal(0, 1), class = "b", coef = "treatmentELD"),
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 3), class = "Intercept"),
    prior(gamma(2, 0.1), class = "phi"),
    prior(normal(0, 1.5), class = "Intercept", dpar = "zi"),
    prior(normal(0, 1), class = "b", dpar = "zi")
  ),
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta),
  file = file.path(paths$models, "fit_sensitivity_caltime_zib"),
  file_refit = "on_change"
)

# Report interaction coefficient
interaction_draws <- as_draws_df(fit_caltime)$`b_treatmentELD:calendar_time_z`
interaction_mean <- mean(interaction_draws)
interaction_cri <- quantile(interaction_draws, c(0.025, 0.975))
p_positive <- mean(interaction_draws > 0)
cat(sprintf("  Interaction (treatment x calendar_time_z): %.2f [%.2f, %.2f]\n",
            interaction_mean, interaction_cri[1], interaction_cri[2]))
cat(sprintf("  P(interaction > 0) = %.3f, P(interaction < 0) = %.3f\n",
            p_positive, 1 - p_positive))
cat(sprintf("  %s\n",
            ifelse(interaction_cri[1] < 0 & interaction_cri[2] > 0,
                   "95% CrI includes zero: no evidence of time-varying treatment effect",
                   "95% CrI excludes zero: evidence of time-varying treatment effect")))


ate_caltime <- compute_gcomp_ate(fit_caltime, df_disc_full, "treatment", "continuous",
                                  lower_is_better = TRUE, scale_factor = 100)
ate_caltime_summary <- summarize_ate(ate_caltime$ate, ni_margin = ni_margins$odi)
cat(sprintf("  ATE (marginal): %.2f (95%% CrI: [%.2f, %.2f]), P(NI): %.4f\n",
            ate_caltime_summary$mean, ate_caltime_summary$cri_lo,
            ate_caltime_summary$cri_hi, ate_caltime_summary$p_ni))

sensitivity_results$caltime_interaction <- list(
  label = "Treatment x calendar time interaction",
  ate_summary = ate_caltime_summary,
  interaction_mean = interaction_mean,
  interaction_cri = interaction_cri
)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n=== Model Sensitivity Summary ===\n")
sens_table <- tibble(
  Analysis = c("Primary (ZIB regression)",
               map_chr(sensitivity_results, ~.x$label)),
  ATE = c(primary_results$ate_summary$mean,
          map_dbl(sensitivity_results, ~.x$ate_summary$mean)),
  CrI_lo = c(primary_results$ate_summary$cri_lo,
             map_dbl(sensitivity_results, ~.x$ate_summary$cri_lo)),
  CrI_hi = c(primary_results$ate_summary$cri_hi,
             map_dbl(sensitivity_results, ~.x$ate_summary$cri_hi)),
  P_NI = c(primary_results$ate_summary$p_ni,
           map_dbl(sensitivity_results, ~.x$ate_summary$p_ni)),
  NI_Conclusion = ifelse(P_NI > 0.95, "NI", "Not NI")
)

print(sens_table)
write.csv(sens_table, file.path(paths$tables, "model_sensitivity_summary.csv"),
          row.names = FALSE)

saveRDS(sensitivity_results, file.path(paths$output, "model_sensitivity_results.rds"))

cat("\nModel sensitivity analysis complete.\n")
