# =============================================================================
# ENDO-LUMBAR: 26 LOO-IC Model Comparison for Primary Outcome
# Formal model selection via LOO-CV (PSIS-LOO) for ODI at 3 months
# =============================================================================

source("/Users/cjsogn/ENDO_LUMBAR/scripts/analysis/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
df_disc <- standardize_covs(df_disc)

# Complete cases only (all models must use the same observations for LOO comparison)
df_cc <- df_disc %>% filter(!is.na(odi_3m))
cat(sprintf("Complete cases for LOO comparison: %d\n", nrow(df_cc)))

# Transform outcomes
df_cc$odi_3m_zoib <- transform_for_zoib(df_cc$odi_3m, upper = 100)
df_cc$odi_3m_beta <- transform_for_beta(df_cc$odi_3m, upper = 100)

cat(sprintf("  Exact zeros (ODI=0): %d (%.1f%%)\n",
    sum(df_cc$odi_3m_zoib == 0), 100 * mean(df_cc$odi_3m_zoib == 0)))

# =============================================================================
# FIT CANDIDATE MODELS (cached via file=)
# =============================================================================

model_dir <- file.path(paths$models, "loo_comparison")
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

common_mcmc <- list(
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

cat("\n=== Fitting candidate models ===\n")

# 1. ZIB with zi ~ treatment + baseline (primary specification)
cat("\n--- 1. ZIB (zi ~ treatment + baseline) ---\n")
fit_zib <- brm(
  bf(as.formula(paste("odi_3m_zoib ~ treatment +", cov_string)),
     zi ~ treatment + odi_baseline_z),
  data = df_cc,
  family = zero_inflated_beta(),
  prior = priors_zoib,
  chains = common_mcmc$chains, iter = common_mcmc$iter,
  warmup = common_mcmc$warmup, cores = common_mcmc$cores,
  seed = common_mcmc$seed, control = common_mcmc$control,
  file = file.path(model_dir, "loo_zib_full"),
  file_refit = "on_change"
)

# 2. ZIB with zi ~ 1 (intercept-only zero-inflation)
cat("\n--- 2. ZIB (zi ~ 1) ---\n")
fit_zib_int <- brm(
  bf(as.formula(paste("odi_3m_zoib ~ treatment +", cov_string)),
     zi ~ 1),
  data = df_cc,
  family = zero_inflated_beta(),
  prior = c(
    prior(normal(0, 1), class = "b", coef = "treatmentELD"),
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 3), class = "Intercept"),
    prior(gamma(2, 0.1), class = "phi"),
    prior(normal(0, 1.5), class = "Intercept", dpar = "zi")
  ),
  chains = common_mcmc$chains, iter = common_mcmc$iter,
  warmup = common_mcmc$warmup, cores = common_mcmc$cores,
  seed = common_mcmc$seed, control = common_mcmc$control,
  file = file.path(model_dir, "loo_zib_intercept"),
  file_refit = "on_change"
)

# 3. Gaussian
cat("\n--- 3. Gaussian ---\n")
fit_gaussian <- brm(
  bf(as.formula(paste("odi_3m ~ treatment +", cov_string))),
  data = df_cc,
  family = gaussian(),
  prior = priors_continuous,
  chains = common_mcmc$chains, iter = common_mcmc$iter,
  warmup = common_mcmc$warmup, cores = common_mcmc$cores,
  seed = common_mcmc$seed, control = common_mcmc$control,
  file = file.path(model_dir, "loo_gaussian"),
  file_refit = "on_change"
)

# 4. Beta with Smithson-Verkuilen transform
cat("\n--- 4. Beta (SV transform) ---\n")
fit_beta_sv <- brm(
  bf(as.formula(paste("odi_3m_beta ~ treatment +", cov_string))),
  data = df_cc,
  family = Beta(),
  prior = priors_beta,
  chains = common_mcmc$chains, iter = common_mcmc$iter,
  warmup = common_mcmc$warmup, cores = common_mcmc$cores,
  seed = common_mcmc$seed, control = common_mcmc$control,
  file = file.path(model_dir, "loo_beta_sv"),
  file_refit = "on_change"
)

# 5. Student-t (estimated nu)
cat("\n--- 5. Student-t ---\n")
fit_studentt <- brm(
  bf(as.formula(paste("odi_3m ~ treatment +", cov_string))),
  data = df_cc,
  family = student(),
  prior = c(
    prior(normal(0, 10), class = "b", coef = "treatmentELD"),
    prior(normal(1, 0.5), class = "b", coef = "odi_baseline_z"),
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 15), class = "sigma"),
    prior(gamma(2, 0.1), class = "nu"),
    prior(normal(30, 20), class = "Intercept")
  ),
  chains = common_mcmc$chains, iter = common_mcmc$iter,
  warmup = common_mcmc$warmup, cores = common_mcmc$cores,
  seed = common_mcmc$seed,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = file.path(model_dir, "loo_studentt"),
  file_refit = "on_change"
)

# 6. Hurdle lognormal
cat("\n--- 6. Hurdle lognormal ---\n")
# Transform: add small constant to avoid log(0) for non-zero values
df_cc$odi_3m_pos <- ifelse(df_cc$odi_3m == 0, 0, df_cc$odi_3m)
fit_hurdle_ln <- tryCatch({
  brm(
    bf(as.formula(paste("odi_3m_pos ~ treatment +", cov_string)),
       hu ~ treatment + odi_baseline_z),
    data = df_cc,
    family = hurdle_lognormal(),
    prior = c(
      prior(normal(0, 1), class = "b", coef = "treatmentELD"),
      prior(normal(0, 0.5), class = "b"),
      prior(normal(3, 1), class = "Intercept"),
      prior(student_t(3, 0, 1), class = "sigma"),
      prior(normal(0, 1.5), class = "Intercept", dpar = "hu"),
      prior(normal(0, 1), class = "b", dpar = "hu")
    ),
    chains = common_mcmc$chains, iter = common_mcmc$iter,
    warmup = common_mcmc$warmup, cores = common_mcmc$cores,
    seed = common_mcmc$seed,
    control = list(adapt_delta = 0.99, max_treedepth = 14),
    file = file.path(model_dir, "loo_hurdle_lognormal"),
    file_refit = "on_change"
  )
}, error = function(e) {
  cat(sprintf("  Hurdle lognormal failed: %s\n", e$message))
  NULL
})

# =============================================================================
# COMPUTE LOO-IC FOR EACH MODEL
# =============================================================================

cat("\n=== Computing LOO-IC ===\n")

# NOTE: Models on different response scales (Gaussian/Student-t on [0,100] vs
# Beta/ZIB on (0,1)) produce ELPD values on different log-likelihood scales.
# Direct loo_compare() across these families is only valid after Jacobian
# correction: Beta/ZIB log-likelihoods must be adjusted by adding
# log(1/100) = -log(100) per observation to account for the division-by-100
# transformation. Without this, Gaussian models appear artificially favored.
# Within-family comparisons (ZIB vs ZIB, Gaussian vs Student-t) are always valid.

models <- list(
  "ZIB (zi~tx+baseline)" = fit_zib,
  "ZIB (zi~1)"           = fit_zib_int,
  "Gaussian"             = fit_gaussian,
  "Beta (SV)"            = fit_beta_sv,
  "Student-t"            = fit_studentt
)
if (!is.null(fit_hurdle_ln)) {
  models[["Hurdle lognormal"]] <- fit_hurdle_ln
}

loo_list <- list()
for (nm in names(models)) {
  cat(sprintf("  Computing LOO for %s... ", nm))
  loo_list[[nm]] <- tryCatch({
    l <- loo(models[[nm]], cores = n_cores, moment_match = TRUE)
    k_vals <- l$diagnostics$pareto_k
    n_bad <- sum(k_vals > 0.7, na.rm = TRUE)
    cat(sprintf("ELPD = %.1f (SE = %.1f), Pareto k>0.7: %d\n",
                l$estimates["elpd_loo", "Estimate"],
                l$estimates["elpd_loo", "SE"],
                n_bad))
    l
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message))
    # Try without moment matching
    tryCatch({
      l <- loo(models[[nm]], cores = n_cores)
      cat(sprintf("  (without moment_match) ELPD = %.1f\n",
                  l$estimates["elpd_loo", "Estimate"]))
      l
    }, error = function(e2) {
      cat(sprintf("  Also failed without moment_match: %s\n", e2$message))
      NULL
    })
  })
}

# Remove NULLs
loo_list <- loo_list[!sapply(loo_list, is.null)]

# =============================================================================
# LOO COMPARISON TABLE
# =============================================================================

cat("\n=== LOO Model Comparison ===\n")
comp <- loo_compare(loo_list)
print(comp)

# Save comparison as CSV
comp_df <- as.data.frame(comp)
comp_df$model <- rownames(comp_df)
comp_df <- comp_df[, c("model", "elpd_diff", "se_diff", "elpd_loo", "se_elpd_loo",
                         "p_loo", "se_p_loo", "looic", "se_looic")]
write.csv(comp_df, file.path(paths$tables, "loo_model_comparison.csv"),
          row.names = FALSE)
cat("\nSaved: loo_model_comparison.csv\n")

# =============================================================================
# PARETO K DIAGNOSTICS FOR SELECTED MODEL (ZIB)
# =============================================================================

cat("\n=== Pareto k diagnostics for selected ZIB model ===\n")
if ("ZIB (zi~tx+baseline)" %in% names(loo_list)) {
  loo_zib <- loo_list[["ZIB (zi~tx+baseline)"]]
  k_vals <- loo_zib$diagnostics$pareto_k

  cat(sprintf("  N observations: %d\n", length(k_vals)))
  cat(sprintf("  Pareto k > 0.5: %d (%.1f%%)\n",
      sum(k_vals > 0.5), 100 * mean(k_vals > 0.5)))
  cat(sprintf("  Pareto k > 0.7: %d (%.1f%%)\n",
      sum(k_vals > 0.7), 100 * mean(k_vals > 0.7)))
  cat(sprintf("  Pareto k > 1.0: %d (%.1f%%)\n",
      sum(k_vals > 1.0), 100 * mean(k_vals > 1.0)))
  cat(sprintf("  Max Pareto k: %.3f\n", max(k_vals)))

  # Save Pareto k plot
  p_k <- plot(loo_zib, label_points = TRUE) +
    labs(title = "Pareto k Diagnostics for Primary ZIB Model",
         subtitle = sprintf("ODI at 3 months (N = %d). k > 0.7 in %d observations (%.1f%%)",
                            length(k_vals), sum(k_vals > 0.7),
                            100 * mean(k_vals > 0.7)))
  ggsave(file.path(paths$diagnostics, "loo_pareto_k_primary.png"), p_k,
         width = 8, height = 5, dpi = 300, bg = "white")
  cat("  Saved: loo_pareto_k_primary.png\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat("Best model: ", rownames(comp)[1], "\n")
if (nrow(comp) > 1) {
  cat(sprintf("ELPD difference (2nd best vs best): %.1f (SE: %.1f)\n",
      comp[2, "elpd_diff"],
      comp[2, "se_diff"]))
}

cat("\nCAVEAT: Cross-family comparison (Gaussian vs Beta/ZIB) requires Jacobian\n")
cat("correction. Within-family comparisons (e.g., ZIB vs ZIB) are directly valid.\n")

cat("\nLOO model comparison complete.\n")
