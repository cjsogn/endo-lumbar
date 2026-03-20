# =============================================================================
# ENDO-LUMBAR: 10 Missing Data Sensitivity Analyses
# =============================================================================

source("/Users/cjsogn/ENDO_LUMBAR/scripts/analysis/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Missing Data Sensitivity Analyses ===\n")

# Prepare data (standardize_covs and cov_string from 00_config.R)
df_disc <- standardize_covs(df_disc)

# ZIB transformation for ODI
df_disc$odi_3m_zib <- transform_for_zib(df_disc$odi_3m, upper = 100)

# ZIB complete-case model used as base for all sensitivity analyses

# =============================================================================
# PATTERN-MIXTURE MODEL (delta-adjustment)
# =============================================================================

cat("\n--- Pattern-Mixture Model ---\n")
cat(sprintf("Delta grid: %s ODI points\n", paste(delta_grid, collapse = ", ")))
cat("Delta adjusts predictions for patients with missing outcomes.\n")

# Get primary model and identify missing patients
missing_idx <- which(is.na(df_disc$odi_3m))
observed_idx <- which(!is.na(df_disc$odi_3m))
n_missing <- length(missing_idx)
cat(sprintf("Missing outcomes: %d/%d\n", n_missing, nrow(df_disc)))

fit_primary <- primary_results$fit

# Identify missing patients by treatment group
eld_missing <- missing_idx[df_disc$treatment[missing_idx] == "ELD"]
msd_missing <- missing_idx[df_disc$treatment[missing_idx] == "MSD"]
cat(sprintf("  ELD missing: %d, MSD missing: %d\n",
            length(eld_missing), length(msd_missing)))

# Compute counterfactual predictions ONCE (reused across all deltas)
nd_eld <- nd_msd <- df_disc
nd_eld$treatment <- factor("ELD", levels = levels(df_disc$treatment))
nd_msd$treatment <- factor("MSD", levels = levels(df_disc$treatment))

cat("  Computing base counterfactual predictions...\n")
pred_eld_base <- posterior_epred(fit_primary, newdata = nd_eld, allow_new_levels = TRUE)
pred_msd_base <- posterior_epred(fit_primary, newdata = nd_msd, allow_new_levels = TRUE)

# Scale predictions from [0,1) ZIB scale to ODI points for delta shifts
pred_eld_base <- pred_eld_base * 100
pred_msd_base <- pred_msd_base * 100

# For each delta: shift factual predictions for missing patients and recompute ATE
pm_results <- map_dfr(delta_grid, function(delta) {
  cat(sprintf("  delta = %+d: ", delta))

  pred_eld <- pred_eld_base
  pred_msd <- pred_msd_base

  # Apply delta shift to the FACTUAL treatment prediction for missing patients
  # ELD-missing patients: shift their ELD (factual) prediction
  # MSD-missing patients: shift their MSD (factual) prediction
  if (length(eld_missing) > 0) {
    pred_eld[, eld_missing] <- pred_eld[, eld_missing] + delta
  }
  if (length(msd_missing) > 0) {
    pred_msd[, msd_missing] <- pred_msd[, msd_missing] + delta
  }

  # Compute ATE (MSD - ELD, positive = ELD superior)
  mu_eld <- rowMeans(pred_eld)
  mu_msd <- rowMeans(pred_msd)
  ate <- mu_msd - mu_eld

  p_ni <- mean(ate > -ni_margins$odi)
  cat(sprintf("ATE = %.2f, P(NI) = %.4f\n", mean(ate), p_ni))

  tibble(
    delta = delta,
    ate_mean = mean(ate),
    ate_cri_lo = quantile(ate, 0.025),
    ate_cri_hi = quantile(ate, 0.975),
    p_ni = p_ni
  )
})

cat("\nPattern-mixture results:\n")
print(pm_results)

# Plot P(NI) vs delta
p_pm <- ggplot(pm_results, aes(x = delta, y = p_ni)) +
  geom_line(linewidth = 0.8, color = tx_colors["ELD"]) +
  geom_point(size = 2, color = tx_colors["ELD"]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text", x = 0, y = 0.96, label = "P = 0.95", hjust = -0.1,
           color = "red", size = 3) +
  labs(
    x = expression(delta ~ "(ODI points shift for missing outcomes)"),
    y = "P(Non-inferiority)",
    title = "Pattern-Mixture Sensitivity: P(NI) vs Delta",
    subtitle = "Delta shifts predicted outcomes for LTFU patients"
  ) +
  scale_x_continuous(breaks = delta_grid)

save_fig(p_pm, "pattern_mixture_sensitivity.png")

# =============================================================================
# TIPPING-POINT ANALYSIS
# =============================================================================

cat("\n--- Tipping-Point Analysis ---\n")
cat("Finding delta* that reverses the NI conclusion.\n")

# Fine-grained delta search (positive deltas only; adversarial direction)
delta_fine <- seq(0, 25, by = 0.5)

# Reuse base predictions computed above for pattern-mixture
tp_results <- map_dfr(delta_fine, function(delta) {
  pred_eld <- pred_eld_base
  pred_msd <- pred_msd_base

  # Adversarial shift: ELD-missing get worse, MSD-missing get better
  if (length(eld_missing) > 0) {
    pred_eld[, eld_missing] <- pred_eld[, eld_missing] + delta
  }
  if (length(msd_missing) > 0) {
    pred_msd[, msd_missing] <- pred_msd[, msd_missing] - delta
  }

  mu_eld <- rowMeans(pred_eld)
  mu_msd <- rowMeans(pred_msd)
  ate <- mu_msd - mu_eld
  p_ni <- mean(ate > -ni_margins$odi)

  tibble(delta = delta, ate_mean = mean(ate), p_ni = p_ni)
})

# Find tipping point: smallest positive delta where P(NI) drops below 0.95
tipping_point <- tp_results %>%
  filter(p_ni < 0.95) %>%
  slice_min(delta, n = 1) %>%
  pull(delta)

if (length(tipping_point) > 0) {
  cat(sprintf("  Tipping point (delta*): %.1f ODI points\n", tipping_point))
} else {
  cat("  No tipping point found in range [0, 25]. NI conclusion is very resilient.\n")
}

p_tp <- ggplot(tp_results, aes(x = delta, y = p_ni)) +
  geom_line(linewidth = 0.8, color = tx_colors["ELD"]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  {if (length(tipping_point) > 0)
    geom_vline(xintercept = tipping_point, linetype = "dotted", color = "grey50")} +
  labs(
    x = expression(delta ~ "(adversarial differential shift, ODI points)"),
    y = "P(Non-inferiority)",
    title = "Tipping-Point Analysis",
    subtitle = "How large an adversarial shift is needed to reverse NI?"
  )

save_fig(p_tp, "tipping_point_odi3m.png")

# =============================================================================
# GAUSSIAN MI() COMPARISON
# =============================================================================

cat("\n--- Gaussian mi() Comparison (MAR with full sample) ---\n")

cat(sprintf("Full sample: %d, Complete cases (primary): %d\n",
            nrow(df_disc), sum(!is.na(df_disc$odi_3m))))

fit_gaussian_mi <- brm(
  bf(as.formula(paste("odi_3m | mi() ~ treatment +", cov_string))),
  data = df_disc,
  family = gaussian(),
  prior = priors_continuous,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta),
  file = file.path(paths$models, "fit_gaussian_mi_comparison"),
  file_refit = "on_change"
)

ate_gauss_mi <- compute_gcomp_ate(fit_gaussian_mi, df_disc, "treatment", "continuous",
                                   lower_is_better = TRUE)
ate_gauss_mi_summary <- summarize_ate(ate_gauss_mi$ate, ni_margin = ni_margins$odi)

cat(sprintf("  Gaussian mi() ATE: %.2f (95%% CrI: [%.2f, %.2f])\n",
            ate_gauss_mi_summary$mean, ate_gauss_mi_summary$cri_lo, ate_gauss_mi_summary$cri_hi))
cat(sprintf("  Gaussian mi() P(NI): %.4f\n", ate_gauss_mi_summary$p_ni))
cat(sprintf("  Primary ZIB CC ATE: %.2f, P(NI): %.4f\n",
            primary_results$ate_summary$mean, primary_results$ate_summary$p_ni))

# Use primary results as the CC reference
ate_cc_summary <- primary_results$ate_summary

# =============================================================================
# COVARIATE SUB-MODEL SENSITIVITY
# =============================================================================

cat("\n--- Covariate Sub-Model Sensitivity ---\n")

# Build var_meta from the covariate missingness table
miss_table <- read.csv(file.path(paths$tables, "covariate_missingness_disc.csv"),
                        stringsAsFactors = FALSE)
var_meta <- list(
  covs_submodel = miss_table$variable[miss_table$pct_missing >= 5 &
                                       miss_table$pct_missing < 50],
  covs_latent   = miss_table$variable[miss_table$pct_missing >= 50]
)
cat(sprintf("  Covariates 5-50%% missing: %d (%s)\n",
            length(var_meta$covs_submodel),
            ifelse(length(var_meta$covs_submodel) > 0,
                   paste(var_meta$covs_submodel, collapse = ", "), "none")))
cat(sprintf("  Covariates >50%% missing: %d (%s)\n",
            length(var_meta$covs_latent),
            ifelse(length(var_meta$covs_latent) > 0,
                   paste(var_meta$covs_latent, collapse = ", "), "none")))

# Check which covariates need sub-models
if (length(var_meta$covs_submodel) == 0 && length(var_meta$covs_latent) == 0) {
  cat("\nAll covariates <5% missing; simple imputation sufficient.\n")

  submodel_sensitivity <- list(
    needed = FALSE,
    covs_submodel = character(0),
    covs_latent = character(0),
    note = "All covariates <5% missing; simple imputation used"
  )

} else {
  cat(sprintf("\nCovariates requiring sub-models (5-50%%): %s\n",
              paste(var_meta$covs_submodel, collapse = ", ")))
  cat(sprintf("Covariates requiring latent modeling (>50%%): %s\n",
              paste(var_meta$covs_latent, collapse = ", ")))

  # Load raw (unimputed) data for jointly modeling missing covariates
  df_raw <- readRDS(file.path(paths$data_clean, "df_disc_raw.rds"))
  df_raw <- standardize_covs(df_raw)

  # Build multi-formula brms model with mi() sub-models
  submodel_formulas <- list()
  submodel_formulas[[1]] <- bf(as.formula(paste("odi_3m | mi() ~ treatment +", cov_string)))

  for (cov_name in var_meta$covs_submodel) {
    cov_z <- paste0(cov_name, "_z")
    if (cov_z %in% names(df_raw)) {
      # Continuous covariate: Gaussian sub-model
      adj_covs <- setdiff(strsplit(cov_string, " \\+ ")[[1]], cov_z)
      adj_string <- paste(adj_covs, collapse = " + ")
      submodel_formulas[[length(submodel_formulas) + 1]] <-
        bf(as.formula(paste(cov_z, "| mi() ~ treatment +", adj_string)),
           family = gaussian())
    }
  }

  if (length(submodel_formulas) > 1) {
    combined_formula <- Reduce("+", submodel_formulas)

    cat(sprintf("Fitting joint model with %d sub-formulas...\n",
                length(submodel_formulas)))

    fit_submodel <- brm(
      formula = combined_formula,
      data = df_raw,
      prior = priors_continuous,
      chains = mcmc_settings$chains,
      iter = mcmc_settings$iter,
      warmup = mcmc_settings$warmup,
      cores = min(mcmc_settings$chains, n_cores),
      seed = mcmc_settings$seed,
      control = list(adapt_delta = 0.99, max_treedepth = 14),
      file = file.path(paths$models, "fit_sensitivity_submodel"),
      file_refit = "on_change"
    )

    ate_sub <- compute_gcomp_ate(fit_submodel, df_raw, "treatment",
                                  "continuous", lower_is_better = TRUE)
    ate_sub_summary <- summarize_ate(ate_sub$ate, ni_margin = ni_margins$odi)

    cat(sprintf("  ATE: %.2f (95%% CrI: [%.2f, %.2f]), P(NI): %.4f\n",
                ate_sub_summary$mean, ate_sub_summary$cri_lo,
                ate_sub_summary$cri_hi, ate_sub_summary$p_ni))

    submodel_sensitivity <- list(
      needed = TRUE,
      covs_submodel = var_meta$covs_submodel,
      covs_latent = var_meta$covs_latent,
      ate_summary = ate_sub_summary
    )
  } else {
    cat("No continuous covariates with 5-50% missing. Sub-model not fitted.\n")
    submodel_sensitivity <- list(
      needed = FALSE,
      covs_submodel = var_meta$covs_submodel,
      covs_latent = var_meta$covs_latent,
      note = "No continuous covariates in 5-50% range"
    )
  }
}

# =============================================================================
# SELECTION MODEL (custom Stan)
# =============================================================================

cat("\n--- Selection Model ---\n")
cat("Joint model for outcome Y and response indicator R via custom Stan.\n")

# Write Stan model for selection model
stan_selection_model <- '
data {
  int<lower=0> N;           // total patients
  int<lower=0> N_obs;       // observed outcomes
  int<lower=0> N_mis;       // missing outcomes
  int<lower=1> K;           // number of covariates
  array[N_obs] int obs_idx; // indices of observed
  array[N_mis] int mis_idx; // indices of missing
  vector[N_obs] y_obs;      // observed outcomes
  matrix[N, K] X;           // design matrix (with treatment)
  array[N] int<lower=0, upper=1> R;  // response indicator
  real rho;                 // sensitivity parameter
}
parameters {
  vector[K] beta;           // outcome model coefficients
  real alpha_y;             // outcome intercept
  real<lower=0> sigma;      // outcome SD
  vector[K] gamma;          // selection model coefficients
  real alpha_r;             // selection intercept
  vector[N_mis] y_mis;      // latent missing outcomes
}
model {
  // Priors
  beta ~ normal(0, 10);
  alpha_y ~ normal(30, 20);
  sigma ~ student_t(3, 0, 15);
  gamma ~ normal(0, 2);
  alpha_r ~ normal(0, 5);

  // Outcome model
  vector[N] y_full;
  y_full[obs_idx] = y_obs;
  y_full[mis_idx] = y_mis;

  y_full ~ normal(alpha_y + X * beta, sigma);

  // Selection model: P(R=1 | X, Y) = inv_logit(alpha_r + X*gamma + rho*Y)
  for (i in 1:N) {
    R[i] ~ bernoulli_logit(alpha_r + dot_product(X[i], gamma) + rho * y_full[i]);
  }
}
'

# Only run selection model if cmdstanr is available
if (requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))) {

  stan_file <- file.path(paths$models, "selection_model.stan")
  writeLines(stan_selection_model, stan_file)

  # Prepare Stan data (reduced covariate set for computational tractability)
  X_mat <- model.matrix(
    ~ treatment + age_z + sex + bmi_z + odi_baseline_z + eq5d_baseline_z +
      nrs_back_baseline_z + nrs_leg_baseline_z + motor_deficit +
      depression_anxiety + chronic_pain + prior_surgery_any +
      sick_leave + disability,
    data = df_disc
  )[, -1]  # remove intercept

  obs_idx <- which(!is.na(df_disc$odi_3m))
  mis_idx <- which(is.na(df_disc$odi_3m))

  stan_data <- list(
    N = nrow(df_disc),
    N_obs = length(obs_idx),
    N_mis = length(mis_idx),
    K = ncol(X_mat),
    obs_idx = obs_idx,
    mis_idx = mis_idx,
    y_obs = df_disc$odi_3m[obs_idx],
    X = X_mat,
    R = as.integer(!is.na(df_disc$odi_3m)),
    rho = 0  # will be varied
  )

  # Compile model
  mod <- cmdstanr::cmdstan_model(stan_file)

  # Fit for range of rho values
  rho_grid <- c(-0.10, -0.05, 0, 0.05, 0.10)
  sel_results <- list()

  for (rho_val in rho_grid) {
    cat(sprintf("  Fitting selection model with rho = %.2f...\n", rho_val))
    stan_data$rho <- rho_val

    fit_sel <- mod$sample(
      data = stan_data,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      seed = mcmc_settings$seed,
      adapt_delta = 0.95,
      refresh = 0
    )

    # Extract treatment effect (first coefficient in beta)
    beta_draws <- fit_sel$draws("beta[1]", format = "draws_matrix")
    ate_sel <- -as.vector(beta_draws)  # Negate for our sign convention

    sel_results[[as.character(rho_val)]] <- tibble(
      rho = rho_val,
      ate_mean = mean(ate_sel),
      ate_cri_lo = quantile(ate_sel, 0.025),
      ate_cri_hi = quantile(ate_sel, 0.975),
      p_ni = mean(ate_sel > -ni_margins$odi)
    )
  }

  sel_table <- bind_rows(sel_results)
  cat("\nSelection model results:\n")
  print(sel_table)

} else {
  cat("  cmdstanr not available. Skipping selection model.\n")
  cat("  Install CmdStan to run this analysis.\n")
  sel_table <- NULL
}

# =============================================================================
# SUMMARY TABLE (Table 6)
# =============================================================================

cat("\n=== Missing Data Sensitivity Summary (Table 6) ===\n")

missing_sensitivity_table <- bind_rows(
  tibble(Analysis = "Primary (ZIB, complete cases)",
         ATE = primary_results$ate_summary$mean,
         CrI_lo = primary_results$ate_summary$cri_lo,
         CrI_hi = primary_results$ate_summary$cri_hi,
         P_NI = primary_results$ate_summary$p_ni),
  tibble(Analysis = "Gaussian mi() (full sample, MAR)",
         ATE = ate_gauss_mi_summary$mean,
         CrI_lo = ate_gauss_mi_summary$cri_lo,
         CrI_hi = ate_gauss_mi_summary$cri_hi,
         P_NI = ate_gauss_mi_summary$p_ni),
  pm_results %>%
    mutate(Analysis = sprintf("Pattern-mixture (delta=%+d)", delta)) %>%
    dplyr::select(Analysis, ATE = ate_mean, CrI_lo = ate_cri_lo,
                  CrI_hi = ate_cri_hi, P_NI = p_ni)
)

print(missing_sensitivity_table)
write.csv(missing_sensitivity_table,
          file.path(paths$tables, "table6_missing_sensitivity.csv"),
          row.names = FALSE)

# Save all results
missing_results <- list(
  pattern_mixture = pm_results,
  tipping_point = tp_results,
  gaussian_mi = list(ate_summary = ate_gauss_mi_summary),
  complete_case = list(ate_summary = ate_cc_summary),
  selection_model = sel_table,
  submodel_sensitivity = submodel_sensitivity
)
saveRDS(missing_results, file.path(paths$output, "missing_sensitivity_results.rds"))

cat("\nMissing data sensitivity analysis complete.\n")
