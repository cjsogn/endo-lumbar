# =============================================================================
# ENDO-LUMBAR: 06 Secondary Effectiveness Outcomes (Tier 2)
# SAP Section 15.1, 6.3
# =============================================================================

source("/Users/cjsogn/endo_studies/lumbar/analysis/scripts/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
df_disc_12m <- readRDS(file.path(paths$data_clean, "df_disc_12m_eligible.rds"))
scaling_params <- readRDS(file.path(paths$data_clean, "scaling_params.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Tier 2: Secondary Effectiveness Outcomes ===\n")
cat(sprintf("Full sample (3m outcomes): N=%d (ELD=%d, MSD=%d)\n",
            nrow(df_disc), sum(df_disc$treatment == "ELD"),
            sum(df_disc$treatment == "MSD")))
cat(sprintf("12m-eligible (surgery <= 2024-12-31): N=%d (ELD=%d, MSD=%d)\n",
            nrow(df_disc_12m), sum(df_disc_12m$treatment == "ELD"),
            sum(df_disc_12m$treatment == "MSD")))

# Standardize covariates in both datasets (function from 00_config.R)
df_disc <- standardize_covs(df_disc)
df_disc_12m <- standardize_covs(df_disc_12m)

# Covariate string and priors from 00_config.R (cov_string, priors_continuous)

# Common priors for binary outcomes (SAP Section 15.2, 24.1-24.2)
# Cauchy (student-t df=1) prior provides adaptive shrinkage for low-EPV settings:
# spike near zero shrinks noise covariates; heavy tails preserve strong confounders.
# This approximates horseshoe behavior while allowing coefficient-specific treatment prior.
# (brms does not allow mixing horseshoe special priors with coefficient-specific priors.)
priors_binary <- c(
  prior(normal(0, 1), class = "b", coef = "treatmentELD"),
  prior(student_t(1, 0, 1), class = "b"),
  prior(normal(0, 5), class = "Intercept")
)

# =============================================================================
# Helper: Fit and summarize a Tier 2 outcome
# =============================================================================

fit_tier2 <- function(outcome_var, outcome_label, family_type = "gaussian",
                      ni_margin, lower_is_better = TRUE, timepoint = "3m",
                      beta_upper = NULL, zi_baseline_var = NULL) {

  # Select appropriate dataset: 12m-eligible for 12m outcomes
  if (timepoint == "12m") {
    d <- df_disc_12m
    cat(sprintf("\n--- %s [12m-eligible subset] ---\n", outcome_label))
  } else {
    d <- df_disc
    cat(sprintf("\n--- %s ---\n", outcome_label))
  }

  cat(sprintf("  N: %d, Observed: %d (%.1f%%)\n",
              nrow(d), sum(!is.na(d[[outcome_var]])),
              100 * mean(!is.na(d[[outcome_var]]))))

  # Determine scale_factor for G-computation back-transformation
  sf <- 1

  if (family_type == "gaussian") {
    fml <- bf(as.formula(paste(outcome_var, "| mi() ~ treatment +", cov_string)))
    priors_use <- priors_continuous
    fam <- gaussian()
  } else if (family_type == "zoib") {
    # Transform outcome to [0,1) for ZOIB
    zoib_var <- paste0(outcome_var, "_zoib")
    d[[zoib_var]] <- transform_for_zoib(d[[outcome_var]], upper = beta_upper)
    sf <- beta_upper
    n_zeros <- sum(d[[zoib_var]] == 0, na.rm = TRUE)
    cat(sprintf("  ZOIB transform: [0, %d] -> [0, %.4f], zeros: %d (%.1f%%)\n",
                beta_upper, max(d[[zoib_var]], na.rm = TRUE),
                n_zeros, 100 * n_zeros / sum(!is.na(d[[zoib_var]]))))
    zi_formula <- if (!is.null(zi_baseline_var)) {
      paste("zi ~ treatment +", zi_baseline_var)
    } else {
      "zi ~ treatment"
    }
    # ZOIB does not support mi() in brms; use complete cases
    d <- d[!is.na(d[[outcome_var]]), ]
    cat(sprintf("  ZOIB complete cases: %d\n", nrow(d)))
    fml <- bf(as.formula(paste(zoib_var, "~ treatment +", cov_string)),
              as.formula(zi_formula))
    priors_use <- priors_zoib
    fam <- zero_inflated_beta()
  } else if (family_type == "beta") {
    # Transform outcome to (0,1) for beta regression
    beta_var <- paste0(outcome_var, "_beta")
    d[[beta_var]] <- transform_for_beta(d[[outcome_var]], upper = beta_upper)
    sf <- beta_upper
    cat(sprintf("  Beta transform: [0, %d] -> (%.4f, %.4f)\n",
                beta_upper, min(d[[beta_var]], na.rm = TRUE),
                max(d[[beta_var]], na.rm = TRUE)))
    fml <- bf(as.formula(paste(beta_var, "| mi() ~ treatment +", cov_string)))
    priors_use <- priors_beta
    fam <- Beta()
  } else {
    fml <- bf(as.formula(paste(outcome_var, "~ treatment +", cov_string)))
    priors_use <- priors_binary
    fam <- bernoulli()
  }

  model_name <- paste0("fit_tier2_", gsub("[^a-zA-Z0-9]", "_", outcome_var))
  if (family_type == "zoib") model_name <- paste0(model_name, "_zoib")
  if (family_type == "beta") model_name <- paste0(model_name, "_beta")

  fit <- brm(
    formula = fml,
    data = d,
    family = fam,
    prior = priors_use,
    chains = mcmc_settings$chains,
    iter = mcmc_settings$iter,
    warmup = mcmc_settings$warmup,
    cores = min(mcmc_settings$chains, n_cores),
    seed = mcmc_settings$seed,
    control = list(adapt_delta = mcmc_settings$adapt_delta,
                   max_treedepth = mcmc_settings$max_treedepth),
    file = file.path(paths$models, model_name),
    file_refit = "on_change"
  )

  # Check convergence
  conv <- check_convergence(fit)
  cat(sprintf("  Convergence: %s (Rhat max=%.4f, ESS_bulk min=%.0f)\n",
              ifelse(conv$all_ok, "PASS", "ISSUE"), conv$rhat_max, conv$ess_bulk_min))

  # G-computation
  otype <- ifelse(family_type %in% c("gaussian", "beta", "zoib"), "continuous", "binary")
  ate_draws <- compute_gcomp_ate(
    fit = fit,
    newdata = d,
    treatment_var = "treatment",
    outcome_type = otype,
    lower_is_better = lower_is_better,
    scale_factor = sf
  )

  ate_summary <- summarize_ate(ate_draws$ate, ni_margin = ni_margin)
  cat(sprintf("  ATE: %.3f (95%% CrI: [%.3f, %.3f])\n",
              ate_summary$mean, ate_summary$cri_lo, ate_summary$cri_hi))
  cat(sprintf("  P(NI, margin=%.2f): %.4f\n", ni_margin, ate_summary$p_ni))
  cat(sprintf("  P(Superiority): %.4f\n", ate_summary$p_superiority))

  list(
    outcome = outcome_var,
    label = outcome_label,
    timepoint = timepoint,
    n_sample = nrow(d),
    family = family_type,
    fit = fit,
    ate_draws = ate_draws,
    ate_summary = ate_summary,
    convergence = conv
  )
}

# =============================================================================
# FIT ALL TIER 2 OUTCOMES
# =============================================================================

tier2_results <- list()

# --- Continuous outcomes ---
# ZOIB regression for bounded outcomes (ODI 0-100, NRS 0-10) with zero-inflation
# Gaussian retained for EQ-5D (can have negative values with Norwegian value set)

# ODI 12 months (uses 12m-eligible subset)
tier2_results$odi_12m <- fit_tier2(
  "odi_12m", "ODI 12 months", "zoib",
  ni_margin = ni_margins$odi, lower_is_better = TRUE, timepoint = "12m",
  beta_upper = 100, zi_baseline_var = "odi_baseline_z"
)

# NRS back pain 3 months
tier2_results$nrs_back_3m <- fit_tier2(
  "nrs_back_3m", "NRS back pain 3 months", "zoib",
  ni_margin = ni_margins$nrs_pain, lower_is_better = TRUE, timepoint = "3m",
  beta_upper = 10, zi_baseline_var = "nrs_back_baseline_z"
)

# NRS back pain 12 months (uses 12m-eligible subset)
tier2_results$nrs_back_12m <- fit_tier2(
  "nrs_back_12m", "NRS back pain 12 months", "zoib",
  ni_margin = ni_margins$nrs_pain, lower_is_better = TRUE, timepoint = "12m",
  beta_upper = 10, zi_baseline_var = "nrs_back_baseline_z"
)

# NRS leg pain 3 months
tier2_results$nrs_leg_3m <- fit_tier2(
  "nrs_leg_3m", "NRS leg pain 3 months", "zoib",
  ni_margin = ni_margins$nrs_pain, lower_is_better = TRUE, timepoint = "3m",
  beta_upper = 10, zi_baseline_var = "nrs_leg_baseline_z"
)

# NRS leg pain 12 months (uses 12m-eligible subset)
tier2_results$nrs_leg_12m <- fit_tier2(
  "nrs_leg_12m", "NRS leg pain 12 months", "zoib",
  ni_margin = ni_margins$nrs_pain, lower_is_better = TRUE, timepoint = "12m",
  beta_upper = 10, zi_baseline_var = "nrs_leg_baseline_z"
)

# EQ-5D 3 months (Gaussian: EQ-5D can have negative values)
tier2_results$eq5d_3m <- fit_tier2(
  "eq5d_3m", "EQ-5D 3 months", "gaussian",
  ni_margin = ni_margins$eq5d, lower_is_better = FALSE, timepoint = "3m"
)

# EQ-5D 12 months (Gaussian: EQ-5D can have negative values)
tier2_results$eq5d_12m <- fit_tier2(
  "eq5d_12m", "EQ-5D 12 months", "gaussian",
  ni_margin = ni_margins$eq5d, lower_is_better = FALSE, timepoint = "12m"
)

# --- Binary outcomes ---

# Responder 3 months
tier2_results$responder_3m <- fit_tier2(
  "responder_3m", "Responder 3 months (>=30% or >=10pt)", "bernoulli",
  ni_margin = ni_margins$responder, lower_is_better = FALSE, timepoint = "3m"
)

# Return to work 3 months (excluding retired)
tier2_results$rtw_3m <- fit_tier2(
  "rtw_3m", "Return to work 3 months", "bernoulli",
  ni_margin = ni_margins$rtw, lower_is_better = FALSE, timepoint = "3m"
)

# Return to work 12 months (uses 12m-eligible subset)
tier2_results$rtw_12m <- fit_tier2(
  "rtw_12m", "Return to work 12 months", "bernoulli",
  ni_margin = ni_margins$rtw, lower_is_better = FALSE, timepoint = "12m"
)

# Analgesic use 3 months (lower = better: not using analgesics)
tier2_results$analgesic_3m <- fit_tier2(
  "analgesic_3m", "Analgesic use 3 months", "bernoulli",
  ni_margin = ni_margins$analgesic, lower_is_better = TRUE, timepoint = "3m"
)

# Analgesic use 12 months (uses 12m-eligible subset)
tier2_results$analgesic_12m <- fit_tier2(
  "analgesic_12m", "Analgesic use 12 months", "bernoulli",
  ni_margin = ni_margins$analgesic, lower_is_better = TRUE, timepoint = "12m"
)

# Patient satisfaction 3 months
tier2_results$satisfied_3m <- fit_tier2(
  "satisfied_3m", "Satisfaction 3 months", "bernoulli",
  ni_margin = ni_margins$satisfaction, lower_is_better = FALSE, timepoint = "3m"
)

# Patient satisfaction 12 months (uses 12m-eligible subset)
tier2_results$satisfied_12m <- fit_tier2(
  "satisfied_12m", "Satisfaction 12 months", "bernoulli",
  ni_margin = ni_margins$satisfaction, lower_is_better = FALSE, timepoint = "12m"
)

# GPE 3 months
tier2_results$gpe_3m <- fit_tier2(
  "gpe_success_3m", "GPE success 3 months", "bernoulli",
  ni_margin = ni_margins$gpe, lower_is_better = FALSE, timepoint = "3m"
)

# GPE 12 months (uses 12m-eligible subset)
tier2_results$gpe_12m <- fit_tier2(
  "gpe_success_12m", "GPE success 12 months", "bernoulli",
  ni_margin = ni_margins$gpe, lower_is_better = FALSE, timepoint = "12m"
)

# =============================================================================
# SUMMARY TABLE (Table 3)
# =============================================================================

cat("\n=== Tier 2 Summary Table ===\n")
tier2_table <- map_dfr(tier2_results, function(r) {
  tibble(
    Outcome = r$label,
    Timepoint = r$timepoint,
    N = r$n_sample,
    ATE = r$ate_summary$mean,
    CrI_lo = r$ate_summary$cri_lo,
    CrI_hi = r$ate_summary$cri_hi,
    P_NI = r$ate_summary$p_ni,
    P_Superiority = r$ate_summary$p_superiority,
    NI_Margin = r$ate_summary$ni_margin,
    NI_Conclusion = case_when(
      r$ate_summary$p_ni > 0.95 ~ "NI demonstrated",
      r$ate_summary$p_ni < 0.80 ~ "Concern",
      TRUE ~ "Inconclusive"
    ),
    Convergence = ifelse(r$convergence$all_ok, "OK", "Issue")
  )
})

print(tier2_table, n = 20)
write.csv(tier2_table, file.path(paths$tables, "table3_tier2_effectiveness.csv"),
          row.names = FALSE)

# Save all Tier 2 results
saveRDS(tier2_results, file.path(paths$output, "tier2_results.rds"))

cat("\nTier 2 analysis complete.\n")
