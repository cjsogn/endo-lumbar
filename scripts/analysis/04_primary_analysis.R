# =============================================================================
# ENDO-LUMBAR: 04 Primary Analysis - Bayesian G-Computation
# ODI at 3 months, Disc Herniation Population
# =============================================================================

source("/Users/cjsogn/ENDO_LUMBAR/scripts/analysis/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
var_meta <- readRDS(file.path(paths$data_clean, "var_meta.rds"))

cat("=== Primary Analysis: Bayesian G-Computation ===\n")
cat("Outcome: ODI at 3 months | Population: Disc herniation\n")
cat(sprintf("N = %d (ELD = %d, MSD = %d)\n",
            nrow(df_disc), sum(df_disc$treatment == "ELD"),
            sum(df_disc$treatment == "MSD")))
cat(sprintf("ODI 3m observed: %d (%.1f%%)\n",
            sum(!is.na(df_disc$odi_3m)),
            100 * mean(!is.na(df_disc$odi_3m))))

# =============================================================================
# 1. PREPARE DATA FOR BRMS
# =============================================================================

# Standardize continuous covariates for better MCMC sampling
df_disc <- standardize_covs(df_disc)

# Save scaling parameters for back-transformation
scaling_params <- list(
  age = c(mean = mean(df_disc$age), sd = sd(df_disc$age)),
  bmi = c(mean = mean(df_disc$bmi), sd = sd(df_disc$bmi)),
  odi_baseline = c(mean = mean(df_disc$odi_baseline, na.rm=TRUE),
                    sd = sd(df_disc$odi_baseline, na.rm=TRUE)),
  eq5d_baseline = c(mean = mean(df_disc$eq5d_baseline, na.rm=TRUE),
                     sd = sd(df_disc$eq5d_baseline, na.rm=TRUE)),
  nrs_back_baseline = c(mean = mean(df_disc$nrs_back_baseline, na.rm=TRUE),
                         sd = sd(df_disc$nrs_back_baseline, na.rm=TRUE)),
  nrs_leg_baseline = c(mean = mean(df_disc$nrs_leg_baseline, na.rm=TRUE),
                        sd = sd(df_disc$nrs_leg_baseline, na.rm=TRUE))
)
saveRDS(scaling_params, file.path(paths$data_clean, "scaling_params.rds"))

# =============================================================================
# 2. ZIB REGRESSION
# =============================================================================

# Transform ODI from [0, 100] to [0, 1) for zero-inflated beta regression
cat("\nTransforming ODI 3m for ZIB regression...\n")
df_disc$odi_3m_zib <- transform_for_zib(df_disc$odi_3m, upper = 100)
cat(sprintf("  ODI 3m ZIB range: [%.4f, %.4f] (observed only)\n",
            min(df_disc$odi_3m_zib, na.rm = TRUE),
            max(df_disc$odi_3m_zib, na.rm = TRUE)))
cat(sprintf("  Exact zeros: %d (%.1f%%)\n",
            sum(df_disc$odi_3m_zib == 0, na.rm = TRUE),
            100 * mean(df_disc$odi_3m_zib == 0, na.rm = TRUE)))

# =============================================================================
# 3. MODEL SPECIFICATION
# =============================================================================

# ZIB complete cases (mi() not supported); missing data handled in Script 10
# zi submodel includes treatment + baseline ODI

# Preserve full dataset for G-computation
df_disc_full <- df_disc

# Filter to complete cases for ZIB model fitting
df_disc <- df_disc %>% filter(!is.na(odi_3m))
cat(sprintf("Complete cases for ZIB model fitting: %d\n", nrow(df_disc)))
cat(sprintf("Full sample for G-computation: %d\n", nrow(df_disc_full)))

model_formula <- bf(
  as.formula(paste("odi_3m_zib ~ treatment +", cov_string)),
  zi ~ treatment + odi_baseline_z
)

# =============================================================================
# 4. PRIOR SPECIFICATION
# =============================================================================

# ZIB priors from 00_config.R (logit link scale)
model_priors <- priors_zib

# =============================================================================
# 5. FIT MODEL
# =============================================================================

cat("\nFitting primary Bayesian ZIB regression model...\n")
cat(sprintf("MCMC: %d chains, %d iterations (%d warmup)\n",
            mcmc_settings$chains, mcmc_settings$iter, mcmc_settings$warmup))
cat(sprintf("adapt_delta: %.2f, max_treedepth: %d\n",
            mcmc_settings$adapt_delta, mcmc_settings$max_treedepth))

fit_primary <- brm(
  formula = model_formula,
  data = df_disc,
  family = zero_inflated_beta(),
  prior = model_priors,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(
    adapt_delta = mcmc_settings$adapt_delta,
    max_treedepth = mcmc_settings$max_treedepth
  ),
  file = file.path(paths$models, "fit_primary_odi3m_zib"),
  file_refit = "on_change"
)

# =============================================================================
# 6. CONVERGENCE CHECK
# =============================================================================

cat("\n=== Convergence Diagnostics ===\n")
conv <- check_convergence(fit_primary)
cat(sprintf("Rhat max: %.4f (threshold: %.2f) -> %s\n",
            conv$rhat_max, convergence_thresholds$rhat_mandatory,
            ifelse(conv$rhat_ok, "PASS", "FAIL")))
cat(sprintf("Bulk ESS min: %.0f (threshold: %d) -> %s\n",
            conv$ess_bulk_min, convergence_thresholds$ess_bulk_mandatory,
            ifelse(conv$ess_bulk_ok, "PASS", "FAIL")))
cat(sprintf("Tail ESS min: %.0f (threshold: %d) -> %s\n",
            conv$ess_tail_min, convergence_thresholds$ess_tail_mandatory,
            ifelse(conv$ess_tail_ok, "PASS", "FAIL")))
cat(sprintf("Divergent transitions: %d (threshold: %d) -> %s\n",
            conv$n_divergent, convergence_thresholds$max_divergent,
            ifelse(conv$divergent_ok, "PASS", "FAIL")))

if (!conv$all_ok) {
  cat("\nWARNING: Convergence criteria not fully met. See diagnostics.\n")
}

# Print model summary
cat("\n=== Model Summary ===\n")
print(summary(fit_primary))

# =============================================================================
# 7. G-COMPUTATION
# =============================================================================

cat("\n=== G-Computation: Average Treatment Effect ===\n")

# Counterfactual predictions over full sample; scale_factor = 100 converts [0,1) to ODI

ate_draws <- compute_gcomp_ate(
  fit = fit_primary,
  newdata = df_disc_full,
  treatment_var = "treatment",
  outcome_type = "continuous",
  lower_is_better = TRUE,  # Lower ODI = better
  scale_factor = 100       # ZIB [0,1) -> ODI [0,100]
)

cat(sprintf("G-computation sample: N = %d (ELD = %d, MSD = %d)\n",
            nrow(df_disc_full), sum(df_disc_full$treatment == "ELD"),
            sum(df_disc_full$treatment == "MSD")))

# Summarize ATE
ate_summary <- summarize_ate(ate_draws$ate, ni_margin = ni_margins$odi)

cat(sprintf("\nATE (MSD - ELD): %.2f (95%% CrI: [%.2f, %.2f])\n",
            ate_summary$mean, ate_summary$cri_lo, ate_summary$cri_hi))
cat(sprintf("P(NI): P(Delta > -%.0f) = %.4f %s\n",
            ni_margins$odi, ate_summary$p_ni,
            ifelse(ate_summary$ni_conclusion, "[NON-INFERIOR]", "[NOT DEMONSTRATED]")))
cat(sprintf("P(Superiority): P(Delta > 0) = %.4f\n", ate_summary$p_superiority))

# =============================================================================
# 8. ROPE ANALYSIS
# =============================================================================

# ROPE = [-7, 7] ODI points
p_rope <- mean(ate_draws$ate > -ni_margins$odi & ate_draws$ate < ni_margins$odi)
p_inferior <- mean(ate_draws$ate <= -ni_margins$odi)
p_superior <- mean(ate_draws$ate >= ni_margins$odi)

cat(sprintf("\nROPE Analysis (±%d ODI):\n", ni_margins$odi))
cat(sprintf("  P(inferior, Delta <= -%d): %.4f\n", ni_margins$odi, p_inferior))
cat(sprintf("  P(equivalent, -%d < Delta < %d): %.4f\n",
            ni_margins$odi, ni_margins$odi, p_rope))
cat(sprintf("  P(superior, Delta >= %d): %.4f\n", ni_margins$odi, p_superior))

# =============================================================================
# 9. MARGIN SENSITIVITY CURVE
# =============================================================================

margins_grid <- seq(1, 15, by = 0.5)
margin_sensitivity <- tibble(
  margin = margins_grid,
  p_ni = map_dbl(margins_grid, ~mean(ate_draws$ate > -.x))
)

p_margin <- ggplot(margin_sensitivity, aes(x = margin, y = p_ni)) +
  geom_line(linewidth = 0.8, color = tx_colors["ELD"]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = ni_margins$odi, linetype = "dashed", color = "red") +
  annotate("text", x = ni_margins$odi + 0.3, y = 0.5, label = "Pre-specified\nmargin",
           hjust = 0, size = 3, color = "red") +
  annotate("text", x = 1, y = 0.96, label = "P = 0.95 threshold",
           hjust = 0, size = 3, color = "grey50") +
  labs(
    x = "Non-inferiority Margin (ODI points)",
    y = "P(Non-inferiority)",
    title = "Margin Sensitivity Curve: ODI at 3 Months",
    subtitle = "Disc herniation: ELD vs MSD"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))

save_fig(p_margin, "margin_sensitivity_odi3m.png")

# =============================================================================
# 10. POSTERIOR DENSITY PLOT
# =============================================================================

p_posterior <- ggplot(ate_draws, aes(x = ate)) +
  geom_density(fill = tx_fills["ELD"], alpha = 0.6, color = tx_colors["ELD"],
               linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey30") +
  geom_vline(xintercept = -ni_margins$odi, linetype = "dashed", color = "red",
             linewidth = 0.8) +
  annotate("text", x = -ni_margins$odi - 0.3, y = 0, label = paste0("NI margin = -", ni_margins$odi),
           hjust = 1, vjust = -0.5, size = 3, color = "red") +
  labs(
    x = expression(paste(Delta, " = ", mu[MSD] - mu[ELD], " (ODI points)")),
    y = "Posterior Density",
    title = "Posterior Distribution: Treatment Effect on ODI at 3 Months",
    subtitle = sprintf("ATE = %.1f (95%% CrI: [%.1f, %.1f]); P(NI) = %.3f",
                        ate_summary$mean, ate_summary$cri_lo, ate_summary$cri_hi,
                        ate_summary$p_ni)
  )

save_fig(p_posterior, "posterior_ate_odi3m.png")

# =============================================================================
# 11. SAVE PRIMARY RESULTS
# =============================================================================

primary_results <- list(
  fit = fit_primary,
  ate_draws = ate_draws,
  ate_summary = ate_summary,
  convergence = conv,
  rope = list(p_inferior = p_inferior, p_rope = p_rope, p_superior = p_superior),
  margin_sensitivity = margin_sensitivity,
  scaling_params = scaling_params,
  ni_conclusion = ate_summary$ni_conclusion,
  # Gating: if NI demonstrated, proceed to Tier 3

  gate_open = ate_summary$ni_conclusion
)

saveRDS(primary_results, file.path(paths$output, "primary_results.rds"))

cat("\n=== Primary Analysis Summary ===\n")
cat(sprintf("ATE: %.2f ODI points (95%% CrI: [%.2f, %.2f])\n",
            ate_summary$mean, ate_summary$cri_lo, ate_summary$cri_hi))
cat(sprintf("P(NI, margin=%d): %.4f -> %s\n",
            ni_margins$odi, ate_summary$p_ni,
            ifelse(ate_summary$ni_conclusion, "NON-INFERIORITY DEMONSTRATED",
                   "NON-INFERIORITY NOT DEMONSTRATED")))
cat(sprintf("P(Superiority): %.4f\n", ate_summary$p_superiority))
cat(sprintf("Gated testing: %s\n",
            ifelse(primary_results$gate_open, "GATE OPEN - proceed to Tier 3",
                   "GATE CLOSED - Tier 3 descriptive only")))

cat("\nPrimary analysis complete.\n")
