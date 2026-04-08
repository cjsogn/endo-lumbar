# =============================================================================
# ENDO-LUMBAR: 05 MCMC Diagnostics and Model Assessment
# Primary model: ZIB regression for ODI 3 months
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))
fit <- primary_results$fit

# Recreate standardized variables and ZIB-transformed outcome
df_disc <- standardize_covs(df_disc)
df_disc$odi_3m_zib <- transform_for_zib(df_disc$odi_3m, upper = 100)

cat("=== MCMC Diagnostics for Primary Model (ZIB Regression) ===\n")

# =============================================================================
# Convergence Diagnostics
# =============================================================================

cat("\n--- Convergence ---\n")

# Trace plots for key parameters
p_trace <- mcmc_trace(fit, pars = c("b_treatmentELD", "b_odi_baseline_z", "phi",
                                      "b_Intercept", "b_zi_Intercept",
                                      "b_zi_treatmentELD"),
                       facet_args = list(ncol = 1))
save_fig(p_trace, "trace_primary.png", width = 8, height = 8,
         path = paths$diagnostics)

# Rhat plot
rhat_vals <- brms::rhat(fit)
p_rhat <- mcmc_rhat(rhat_vals) +
  labs(title = "Rhat Diagnostics (ZIB Model)")
save_fig(p_rhat, "rhat_primary.png", width = 7, height = 5,
         path = paths$diagnostics)

# ESS plots
neff_vals <- neff_ratio(fit)
p_neff <- mcmc_neff(neff_vals) +
  labs(title = "Effective Sample Size Ratios")
save_fig(p_neff, "neff_primary.png", width = 7, height = 5,
         path = paths$diagnostics)

# Pairs plot for key parameters (check for correlations/funnels)
p_pairs <- mcmc_pairs(fit, pars = c("b_treatmentELD", "b_odi_baseline_z", "phi",
                                      "b_zi_treatmentELD"),
                       off_diag_args = list(size = 0.5, alpha = 0.3))
save_fig(p_pairs, "pairs_primary.png", width = 8, height = 8,
         path = paths$diagnostics)

# Detailed convergence summary
conv_summary <- posterior::summarise_draws(posterior::as_draws(fit))
cat("Parameters with Rhat > 1.01:\n")
bad_rhat <- conv_summary %>% filter(rhat > 1.01)
if (nrow(bad_rhat) > 0) {
  print(bad_rhat %>% dplyr::select(variable, rhat, ess_bulk, ess_tail))
} else {
  cat("  None (all PASS)\n")
}

cat("\nParameters with ESS_bulk < 400:\n")
bad_ess <- conv_summary %>% filter(ess_bulk < 400)
if (nrow(bad_ess) > 0) {
  print(bad_ess %>% dplyr::select(variable, rhat, ess_bulk, ess_tail))
} else {
  cat("  None (all PASS)\n")
}

# =============================================================================
# Posterior Predictive Checks (PPC)
# =============================================================================

cat("\n--- Posterior Predictive Checks ---\n")

# PPC over patients with observed outcomes
observed_idx <- which(!is.na(df_disc$odi_3m))
y_obs_zib <- df_disc$odi_3m_zib[observed_idx]
y_obs_odi <- df_disc$odi_3m[observed_idx]

# Posterior predictions on [0,1) scale, convert to ODI for visualization
yrep_zib <- posterior_predict(fit, newdata = df_disc[observed_idx, ], ndraws = 100)

# Convert predictions and observed to ODI scale for PPC display
yrep_odi <- yrep_zib * 100
# Clamp to valid ODI range [0, 100]
yrep_odi <- pmin(pmax(yrep_odi, 0), 100)

p_ppc_dens <- ppc_dens_overlay(y_obs_odi, yrep_odi) +
  coord_cartesian(xlim = c(0, max(y_obs_odi) + 5)) +
  labs(title = "PPC: Density Overlay (ODI 3 months, ZIB Model)",
       subtitle = "100 posterior predictive draws vs observed",
       x = "ODI Score")
save_fig(p_ppc_dens, "ppc_density_odi3m.png", width = 7, height = 5,
         path = paths$diagnostics)

# Test statistics on ODI scale
p_ppc_stat_mean <- ppc_stat(y_obs_odi, yrep_odi, stat = "mean") +
  labs(title = "PPC: Mean")
p_ppc_stat_sd <- ppc_stat(y_obs_odi, yrep_odi, stat = "sd") +
  labs(title = "PPC: SD")
p_ppc_stat_median <- ppc_stat(y_obs_odi, yrep_odi, stat = "median") +
  labs(title = "PPC: Median")

p_ppc_stats <- p_ppc_stat_mean + p_ppc_stat_sd + p_ppc_stat_median +
  plot_layout(ncol = 3) +
  plot_annotation(title = "PPC: Test Statistics (ZIB Model)")
save_fig(p_ppc_stats, "ppc_test_stats_odi3m.png", width = 12, height = 4,
         path = paths$diagnostics)

# Bayesian p-values (target: 0.05-0.95)
calc_bp <- function(stat_fn) {
  obs_stat <- stat_fn(y_obs_odi)
  rep_stats <- apply(yrep_odi, 1, stat_fn)
  mean(rep_stats >= obs_stat)
}

bp_mean <- calc_bp(mean)
bp_sd <- calc_bp(sd)
bp_median <- calc_bp(median)
bp_q10 <- calc_bp(function(x) quantile(x, 0.1))
bp_q90 <- calc_bp(function(x) quantile(x, 0.9))

# Skewness
skewness <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  sum(((x - m) / s)^3) / n
}
bp_skew <- calc_bp(skewness)

ppc_summary <- tibble(
  statistic = c("Mean", "SD", "Median", "Q10", "Q90", "Skewness"),
  bayesian_p = c(bp_mean, bp_sd, bp_median, bp_q10, bp_q90, bp_skew),
  pass = bayesian_p > 0.05 & bayesian_p < 0.95
)
cat("\nBayesian p-values (target: 0.05-0.95):\n")
print(ppc_summary)

# =============================================================================
# LOO-CV
# =============================================================================

cat("\n--- LOO-CV ---\n")

# LOO over observed outcomes
df_obs <- df_disc[!is.na(df_disc$odi_3m), ]
loo_primary <- tryCatch({
  loo(fit, newdata = df_obs, cores = n_cores)
}, error = function(e) {
  cat(sprintf("LOO-CV with newdata failed: %s\n", e$message))
  cat("Attempting LOO on full model...\n")
  tryCatch(loo(fit, cores = n_cores), error = function(e2) {
    cat(sprintf("LOO-CV fallback also failed: %s\n", e2$message))
    cat("Skipping LOO-CV due to mi() model structure.\n")
    NULL
  })
})

if (!is.null(loo_primary)) {
  cat("\nLOO-CV Summary:\n")
  print(loo_primary)
} else {
  cat("\nLOO-CV could not be computed for this mi() model.\n")
}

# Pareto k diagnostics
if (!is.null(loo_primary)) {
  k_vals <- loo_primary$diagnostics$pareto_k
  n_high_k <- sum(k_vals > 0.7, na.rm = TRUE)
  n_very_high_k <- sum(k_vals > 1.0, na.rm = TRUE)
  pct_high_k <- 100 * n_high_k / sum(!is.na(k_vals))

  cat(sprintf("\nPareto k > 0.7: %d (%.1f%%)\n", n_high_k, pct_high_k))
  cat(sprintf("Pareto k > 1.0: %d\n", n_very_high_k))

  if (pct_high_k > 5) {
    cat("WARNING: >5% observations with k > 0.7. Consider moment matching or full LOO.\n")
  }

  p_loo_k <- plot(loo_primary, label_points = TRUE)
  save_fig(p_loo_k, "loo_pareto_k_primary.png", width = 7, height = 5,
           path = paths$diagnostics)
} else {
  n_high_k <- NA
  n_very_high_k <- 0
  pct_high_k <- NA
}

# =============================================================================
# Residual Diagnostics
# =============================================================================

cat("\n--- Residual Diagnostics ---\n")

# Posterior mean fitted values (on response scale [0,1)) and convert to ODI
fitted_zib <- fitted(fit, newdata = df_disc[observed_idx, ])[, "Estimate"]
fitted_odi <- fitted_zib * 100  # Convert to ODI scale

# Raw residuals on ODI scale
residuals_raw <- y_obs_odi - fitted_odi
residuals_std <- residuals_raw / sd(residuals_raw)

resid_df <- tibble(
  fitted = fitted_odi,
  residual = residuals_std,
  odi_baseline = df_disc$odi_baseline[observed_idx],
  age = df_disc$age[observed_idx],
  treatment = df_disc$treatment[observed_idx]
)

# Residuals vs fitted
p_resid_fitted <- ggplot(resid_df, aes(x = fitted, y = residual)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_hline(yintercept = 0, color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "blue", linewidth = 0.5) +
  labs(x = "Fitted Values (ODI)", y = "Standardized Residuals",
       title = "Residuals vs Fitted Values")

# Q-Q plot
p_qq <- ggplot(resid_df, aes(sample = residual)) +
  stat_qq(alpha = 0.3, size = 1) +
  stat_qq_line(color = "red") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles",
       title = "Q-Q Plot of Residuals")

# Residuals vs key covariates
p_resid_odi <- ggplot(resid_df, aes(x = odi_baseline, y = residual)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "blue", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "Baseline ODI", y = "Standardized Residuals",
       title = "Residuals vs Baseline ODI")

p_resid_age <- ggplot(resid_df, aes(x = age, y = residual)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "blue", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "Age", y = "Standardized Residuals",
       title = "Residuals vs Age")

p_residuals <- (p_resid_fitted | p_qq) / (p_resid_odi | p_resid_age) +
  plot_annotation(title = "Residual Diagnostics: Primary Model (ZIB Regression)")
save_fig(p_residuals, "residual_diagnostics_primary.png", width = 10, height = 8,
         path = paths$diagnostics)

# =============================================================================
# Pseudo R²
# =============================================================================

cat("\n--- Pseudo R² ---\n")

# For ZIB, compute pseudo R² from posterior predictions
r2_vals <- tryCatch({
  bayes_R2(fit, newdata = df_obs)
}, error = function(e) {
  # Manual pseudo R²: var(predicted) / var(observed)
  cat(sprintf("Note: bayes_R2 fallback (%s)\n", e$message))
  pred <- posterior_epred(fit, newdata = df_obs)
  var_pred <- apply(pred, 1, var)
  var_obs <- var(y_obs_zib)
  var_pred / var_obs
})
r2_clean <- r2_vals[!is.na(r2_vals) & !is.nan(r2_vals)]

if (length(r2_clean) > 0) {
  cat(sprintf("Pseudo R²: %.3f (95%% CrI: [%.3f, %.3f])\n",
              mean(r2_clean), quantile(r2_clean, 0.025), quantile(r2_clean, 0.975)))
  if (mean(r2_clean) < 0.10) {
    cat("WARNING: R² < 0.10. Review covariate specification (SAP threshold: R² > 0.15).\n")
  } else if (mean(r2_clean) < 0.15) {
    cat("NOTE: R² between 0.10 and 0.15. Below SAP expected range but acceptable.\n")
  }
} else {
  cat("Pseudo R² could not be computed.\n")
  r2_clean <- NA
}

# =============================================================================
# Influential Observations
# =============================================================================

cat("\n--- Influential Observations ---\n")

if (n_very_high_k > 0) {
  influential_idx <- which(k_vals > 1.0)
  cat(sprintf("Observations with Pareto k > 1.0: %d\n", length(influential_idx)))
  cat("Characteristics:\n")
  print(df_disc[observed_idx[influential_idx],
                c("treatment", "age", "sex", "odi_baseline", "odi_3m")])
  cat("\nSensitivity analysis excluding these cases will be conducted.\n")
} else {
  cat("No observations with Pareto k > 1.0.\n")
}

# =============================================================================
# SAVE DIAGNOSTICS SUMMARY
# =============================================================================

diagnostics_summary <- list(
  convergence = primary_results$convergence,
  ppc = ppc_summary,
  loo = loo_primary,
  pareto_k = list(n_high = n_high_k, n_very_high = n_very_high_k, pct_high = pct_high_k),
  r2 = list(mean = mean(r2_clean), cri = quantile(r2_clean, c(0.025, 0.975)))
)
saveRDS(diagnostics_summary, file.path(paths$output, "diagnostics_primary.rds"))

cat("\nDiagnostics complete.\n")
