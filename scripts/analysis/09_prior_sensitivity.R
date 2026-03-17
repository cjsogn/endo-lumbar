# =============================================================================
# ENDO-LUMBAR: 09 Prior Sensitivity Analysis
# SAP Section 17
# =============================================================================

source("/Users/cjsogn/endo_studies/lumbar/analysis/scripts/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Prior Sensitivity Analysis (SAP Section 17) ===\n")

# Prepare data (standardize_covs and cov_string from 00_config.R)
df_disc <- standardize_covs(df_disc)

# ZOIB regression: transform ODI to [0,1)
df_disc$odi_3m_zoib <- transform_for_zoib(df_disc$odi_3m, upper = 100)

# Filter to complete cases (ZOIB does not support mi())
df_disc <- df_disc %>% filter(!is.na(odi_3m))
cat(sprintf("Complete cases for prior sensitivity: %d\n", nrow(df_disc)))

model_formula <- bf(
  as.formula(paste("odi_3m_zoib ~ treatment +", cov_string)),
  zi ~ treatment + odi_baseline_z
)

# =============================================================================
# 17.1 Prior Predictive Checks (before fitting data)
# =============================================================================

cat("\n--- 17.1 Prior Predictive Checks ---\n")

# Prior predictive check on logit scale (ZOIB mu component)
prior_pred_results <- list()
for (prior_label in c("Skeptical", "Reference", "Diffuse")) {
  prior_sd <- switch(prior_label,
    "Skeptical" = priors_zoib_skeptical_sd,
    "Reference" = priors_zoib_reference_sd,
    "Diffuse"   = priors_zoib_diffuse_sd
  )

  cat(sprintf("\n  %s prior: N(0, %.1f) on logit scale\n", prior_label, prior_sd))

  # Simulate from prior on logit scale, then approximate effect on ODI scale
  # At mean ODI ≈ 17/100 = 0.17, derivative of inverse-logit ≈ 0.14
  set.seed(mcmc_settings$seed)
  prior_draws_logit <- rnorm(1000, mean = 0, sd = prior_sd)
  # Approximate ODI-scale effect via delta method at the mean
  mean_mu <- mean(df_disc$odi_3m_zoib[df_disc$odi_3m_zoib > 0], na.rm = TRUE)
  deriv <- mean_mu * (1 - mean_mu)  # derivative of inv_logit at logit(mean_mu)
  prior_draws_odi <- prior_draws_logit * deriv * 100  # convert to ODI points

  cat(sprintf("    Approx 95%% prior mass on ODI scale: [%.1f, %.1f]\n",
              quantile(prior_draws_odi, 0.025), quantile(prior_draws_odi, 0.975)))

  prior_pred_results[[prior_label]] <- tibble(
    prior = prior_label,
    prior_sd = prior_sd,
    ate = prior_draws_odi
  )
}

prior_pred_df <- bind_rows(prior_pred_results)

p_prior_pred <- ggplot(prior_pred_df, aes(x = ate, fill = prior, color = prior)) +
  geom_density(alpha = 0.3, linewidth = 0.7) +
  geom_vline(xintercept = c(-ni_margins$odi, ni_margins$odi),
             linetype = "dashed", color = "red") +
  labs(
    x = "Approximate Treatment Effect (ODI points)",
    y = "Prior Density",
    title = "Prior Predictive Check: Treatment Effect Priors (ZOIB Model)",
    subtitle = "Dashed lines = NI margins; effects approximated via delta method"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  coord_cartesian(xlim = c(-20, 20))

save_fig(p_prior_pred, "prior_predictive_check.png")

# =============================================================================
# 17. FIT THREE PRIOR SPECIFICATIONS
# =============================================================================

prior_results <- list()

for (prior_label in c("Skeptical", "Reference", "Diffuse")) {
  prior_sd <- switch(prior_label,
    "Skeptical" = priors_zoib_skeptical_sd,
    "Reference" = priors_zoib_reference_sd,
    "Diffuse"   = priors_zoib_diffuse_sd
  )

  cat(sprintf("\n--- Fitting %s prior: N(0, %.1f) on logit ---\n", prior_label, prior_sd))

  model_priors <- c(
    set_prior(sprintf("normal(0, %.1f)", prior_sd), class = "b", coef = "treatmentELD"),
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 3), class = "Intercept"),
    prior(gamma(2, 0.1), class = "phi"),
    prior(normal(0, 1.5), class = "Intercept", dpar = "zi"),
    prior(normal(0, 1), class = "b", dpar = "zi")
  )

  model_name <- paste0("fit_prior_zoib_", tolower(prior_label))

  fit <- brm(
    formula = model_formula,
    data = df_disc,
    family = zero_inflated_beta(),
    prior = model_priors,
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

  # G-computation (scale_factor = 100 for ZOIB -> ODI)
  ate_draws <- compute_gcomp_ate(fit, df_disc, "treatment", "continuous",
                                  lower_is_better = TRUE, scale_factor = 100)
  ate_summary <- summarize_ate(ate_draws$ate, ni_margin = ni_margins$odi)

  # Prior-to-posterior update ratio (SAP Section 17.2)
  trt_posterior_sd <- sd(posterior::as_draws_df(fit)$b_treatmentELD)
  update_ratio <- trt_posterior_sd / prior_sd

  cat(sprintf("  ATE: %.2f (95%% CrI: [%.2f, %.2f])\n",
              ate_summary$mean, ate_summary$cri_lo, ate_summary$cri_hi))
  cat(sprintf("  P(NI): %.4f\n", ate_summary$p_ni))
  cat(sprintf("  Posterior SD / Prior SD: %.3f (target < 0.5)\n", update_ratio))
  if (update_ratio > 0.8) {
    cat("  WARNING: Prior dominance detected (ratio > 0.8)\n")
  }

  prior_results[[prior_label]] <- list(
    prior_label = prior_label,
    prior_sd = prior_sd,
    fit = fit,
    ate_draws = ate_draws,
    ate_summary = ate_summary,
    update_ratio = update_ratio,
    posterior_sd = trt_posterior_sd
  )
}

# =============================================================================
# 17.3 COMPARISON TABLE (Table 5)
# =============================================================================

cat("\n=== Prior Sensitivity Summary (Table 5) ===\n")

prior_table <- map_dfr(prior_results, function(r) {
  tibble(
    Prior = sprintf("%s N(0,%.1f)", r$prior_label, r$prior_sd),
    ATE = r$ate_summary$mean,
    CrI_lo = r$ate_summary$cri_lo,
    CrI_hi = r$ate_summary$cri_hi,
    P_NI = r$ate_summary$p_ni,
    P_Superiority = r$ate_summary$p_superiority,
    Post_SD = r$posterior_sd,
    Update_Ratio = r$update_ratio,
    NI_Conclusion = ifelse(r$ate_summary$p_ni > 0.95, "NI", "Not NI")
  )
})
print(prior_table)

# Check if NI conclusion changes between priors
conclusions <- prior_table$NI_Conclusion
if (length(unique(conclusions)) > 1) {
  cat("\nIMPORTANT: NI conclusion CHANGES between prior specifications.\n")
  cat("This indicates sensitivity to prior choice.\n")
} else {
  cat(sprintf("\nNI conclusion is CONSISTENT across all priors: %s\n", conclusions[1]))
}

write.csv(prior_table, file.path(paths$tables, "table5_prior_sensitivity.csv"),
          row.names = FALSE)

# Posterior comparison plot
posterior_df <- map_dfr(prior_results, function(r) {
  tibble(prior = r$prior_label, ate = r$ate_draws$ate)
})

p_prior_comp <- ggplot(posterior_df, aes(x = ate, fill = prior, color = prior)) +
  geom_density(alpha = 0.3, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey30") +
  geom_vline(xintercept = -ni_margins$odi, linetype = "dashed", color = "red") +
  labs(
    x = expression(paste(Delta, " (ODI points, positive = ELD superior)")),
    y = "Posterior Density",
    title = "Prior Sensitivity: Posterior ATE Distributions",
    subtitle = sprintf("NI margin = %d ODI points", ni_margins$odi)
  ) +
  scale_fill_brewer(palette = "Set2", name = "Prior") +
  scale_color_brewer(palette = "Set2", name = "Prior")

save_fig(p_prior_comp, "prior_sensitivity_posteriors.png")

saveRDS(prior_results, file.path(paths$output, "prior_sensitivity_results.rds"))

cat("\nPrior sensitivity analysis complete.\n")
