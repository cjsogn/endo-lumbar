# =============================================================================
# ENDO-LUMBAR: 14 Stenosis Population Analysis (Exploratory)
# =============================================================================

source("/Users/cjsogn/ENDO_LUMBAR/scripts/analysis/00_config.R")

df_sten <- readRDS(file.path(paths$data_clean, "df_sten_imp.rds"))

cat("=== Stenosis Population Analysis (Exploratory) ===\n")

n_eld_sten <- sum(df_sten$treatment == "ELD")
n_msd_sten <- sum(df_sten$treatment == "MSD")
cat(sprintf("Stenosis population: %d (ELD=%d, MSD=%d)\n",
            nrow(df_sten), n_eld_sten, n_msd_sten))

# =============================================================================
# FEASIBILITY CHECK
# =============================================================================

if (n_eld_sten < 20) {
  cat("\nELD stenosis group < 20. Only unadjusted comparisons per SAP.\n")

  # Unadjusted comparison
  cat("\n--- Unadjusted ODI 3 months ---\n")
  sten_summary <- df_sten %>%
    group_by(treatment) %>%
    dplyr::summarise(
      n = n(),
      odi_3m_n = sum(!is.na(odi_3m)),
      odi_3m_mean = mean(odi_3m, na.rm = TRUE),
      odi_3m_sd = sd(odi_3m, na.rm = TRUE),
      odi_baseline_mean = mean(odi_baseline, na.rm = TRUE),
      nrs_leg_3m_mean = mean(nrs_leg_3m, na.rm = TRUE),
      eq5d_3m_mean = mean(eq5d_3m, na.rm = TRUE),
      day_surgery_pct = 100 * mean(day_surgery == 1, na.rm = TRUE),
      los_mean = mean(los_total, na.rm = TRUE),
      op_time_mean = mean(operating_time, na.rm = TRUE),
      .groups = "drop"
    )
  print(sten_summary)

  write.csv(sten_summary,
            file.path(paths$tables, "table7_stenosis_unadjusted.csv"),
            row.names = FALSE)

  cat(sprintf("\nNote: Unadjusted comparisons only (%d ELD patients).\n", n_eld_sten))

  stenosis_results <- list(
    feasible = FALSE,
    n_eld = n_eld_sten,
    summary = sten_summary
  )

} else {
  cat("\nELD stenosis group >= 20. Full adjusted model.\n")

  # Prepare data (standardize_covs and cov_string from 00_config.R)
  df_sten <- standardize_covs(df_sten)

  fit_sten <- brm(
    bf(as.formula(paste("odi_3m | mi() ~ treatment +", cov_string))),
    data = df_sten,
    family = gaussian(),
    prior = c(
      prior(normal(0, 10), class = "b", coef = "treatmentELD"),
      prior(normal(1, 0.5), class = "b", coef = "odi_baseline_z"),
      prior(normal(0, 2), class = "b"),
      prior(student_t(3, 0, 15), class = "sigma"),
      prior(normal(30, 20), class = "Intercept")
    ),
    chains = mcmc_settings$chains,
    iter = mcmc_settings$iter,
    warmup = mcmc_settings$warmup,
    cores = min(mcmc_settings$chains, n_cores),
    seed = mcmc_settings$seed,
    control = list(adapt_delta = 0.95),
    file = file.path(paths$models, "fit_stenosis_primary"),
    file_refit = "on_change"
  )

  ate_sten <- compute_gcomp_ate(fit_sten, df_sten, "treatment", "continuous",
                                 lower_is_better = TRUE)
  ate_sten_summary <- summarize_ate(ate_sten$ate, ni_margin = ni_margins$odi)

  cat(sprintf("\nATE: %.2f (95%% CrI: [%.2f, %.2f])\n",
              ate_sten_summary$mean, ate_sten_summary$cri_lo, ate_sten_summary$cri_hi))
  cat(sprintf("P(NI): %.4f\n", ate_sten_summary$p_ni))

  # Prior dominance check
  trt_posterior_sd <- sd(posterior::as_draws_df(fit_sten)$b_treatmentELD)
  update_ratio <- trt_posterior_sd / priors_reference$treatment_sd
  cat(sprintf("Prior-to-posterior update ratio: %.3f\n", update_ratio))
  if (update_ratio > 0.8) {
    cat("WARNING: Prior dominance (ratio > 0.8). Interpret with caution.\n")
  }

  stenosis_results <- list(
    feasible = TRUE,
    n_eld = n_eld_sten,
    fit = fit_sten,
    ate_summary = ate_sten_summary,
    update_ratio = update_ratio
  )
}

saveRDS(stenosis_results, file.path(paths$output, "stenosis_results.rds"))

cat("\nStenosis analysis complete.\n")
