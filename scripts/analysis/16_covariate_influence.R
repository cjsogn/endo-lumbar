# =============================================================================
# ENDO-LUMBAR: 16 Covariate Influence and Additional Sensitivity Analyses
# Post-hoc: raw vs adjusted comparison, leave-one-out, overlap restriction
# =============================================================================

source("/Users/cjsogn/ENDO_LUMBAR/scripts/analysis/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Covariate Influence Analysis ===\n")

# =============================================================================
# DATA PREPARATION
# =============================================================================

# Prepare data (standardize_covs, cov_string, priors_continuous from 00_config.R)
df_disc <- standardize_covs(df_disc)

cov_full <- cov_string       # alias for clarity in this script
shared_priors <- priors_continuous

# =============================================================================
# PART 1: RAW vs ADJUSTED COMPARISON
# =============================================================================

cat("--- Part 1: Raw vs Adjusted Comparison ---\n")

# Unadjusted (treatment only)
fit_unadj <- brm(
  bf(odi_3m | mi() ~ treatment),
  data = df_disc,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "b"),
    prior(student_t(3, 0, 15), class = "sigma"),
    prior(normal(30, 20), class = "Intercept")
  ),
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.95),
  file = file.path(paths$models, "fit_unadjusted"),
  file_refit = "on_change"
)

ate_unadj <- compute_gcomp_ate(fit_unadj, df_disc, "treatment", "continuous",
                                lower_is_better = TRUE)
ate_unadj_sum <- summarize_ate(ate_unadj$ate, ni_margin = ni_margins$odi)
cat(sprintf("  Unadjusted: ATE = %.2f [%.2f, %.2f], P(NI) = %.3f\n",
            ate_unadj_sum$mean, ate_unadj_sum$cri_lo, ate_unadj_sum$cri_hi,
            ate_unadj_sum$p_ni))

# Fully adjusted (primary model - use cached)
cat(sprintf("  Fully adjusted: ATE = %.2f [%.2f, %.2f], P(NI) = %.3f\n",
            primary_results$ate_summary$mean, primary_results$ate_summary$cri_lo,
            primary_results$ate_summary$cri_hi, primary_results$ate_summary$p_ni))

# =============================================================================
# PART 2: LEAVE-ONE-OUT COVARIATE INFLUENCE
# =============================================================================

cat("\n--- Part 2: Leave-One-Out Covariate Analysis ---\n")
cat("Dropping key covariates one at a time to quantify their influence on ATE.\n")

# Covariates to drop (high-influence suspects)
drop_specs <- list(
  list(name = "age", drop = "age_z", label = "Without age"),
  list(name = "smoking", drop = "smoking", label = "Without smoking"),
  list(name = "employment", drop = "employed_baseline", label = "Without employment status"),
  list(name = "baseline_odi", drop = "odi_baseline_z", label = "Without baseline ODI"),
  list(name = "depression", drop = "depression_anxiety", label = "Without depression/anxiety"),
  list(name = "asa", drop = "asa_cat", label = "Without ASA class")
)

loo_cov_results <- list()

for (spec in drop_specs) {
  cat(sprintf("\n  %s:\n", spec$label))

  # Build reduced covariate string
  cov_terms <- strsplit(cov_full, "\\s*\\+\\s*")[[1]]
  cov_terms <- trimws(cov_terms)
  cov_reduced <- cov_terms[!cov_terms %in% spec$drop]
  cov_string_reduced <- paste(cov_reduced, collapse = " + ")

  # Determine priors
  has_odi_base <- any(grepl("odi_baseline_z", cov_reduced))

  priors_reduced <- c(
    prior(normal(0, 10), class = "b", coef = "treatmentELD"),
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 15), class = "sigma"),
    prior(normal(30, 20), class = "Intercept")
  )
  if (has_odi_base) {
    priors_reduced <- c(priors_reduced,
      prior(normal(1, 0.5), class = "b", coef = "odi_baseline_z")
    )
  }

  fit_reduced <- brm(
    bf(as.formula(paste("odi_3m | mi() ~ treatment +", cov_string_reduced))),
    data = df_disc,
    family = gaussian(),
    prior = priors_reduced,
    chains = mcmc_settings$chains,
    iter = mcmc_settings$iter,
    warmup = mcmc_settings$warmup,
    cores = min(mcmc_settings$chains, n_cores),
    seed = mcmc_settings$seed,
    control = list(adapt_delta = 0.95),
    file = file.path(paths$models, paste0("fit_loo_cov_", spec$name)),
    file_refit = "on_change"
  )

  ate_reduced <- compute_gcomp_ate(fit_reduced, df_disc, "treatment", "continuous",
                                    lower_is_better = TRUE)
  ate_reduced_sum <- summarize_ate(ate_reduced$ate, ni_margin = ni_margins$odi)

  shift <- ate_reduced_sum$mean - primary_results$ate_summary$mean
  cat(sprintf("    ATE = %.2f [%.2f, %.2f], P(NI) = %.3f\n",
              ate_reduced_sum$mean, ate_reduced_sum$cri_lo, ate_reduced_sum$cri_hi,
              ate_reduced_sum$p_ni))
  cat(sprintf("    Shift from primary: %+.2f (positive = more favorable for ELD)\n", shift))

  loo_cov_results[[spec$name]] <- tibble(
    specification = spec$label,
    ate_mean = ate_reduced_sum$mean,
    cri_lo = ate_reduced_sum$cri_lo,
    cri_hi = ate_reduced_sum$cri_hi,
    p_ni = ate_reduced_sum$p_ni,
    shift_from_primary = shift
  )
}

loo_table <- bind_rows(
  tibble(specification = "Unadjusted (treatment only)",
         ate_mean = ate_unadj_sum$mean, cri_lo = ate_unadj_sum$cri_lo,
         cri_hi = ate_unadj_sum$cri_hi, p_ni = ate_unadj_sum$p_ni,
         shift_from_primary = ate_unadj_sum$mean - primary_results$ate_summary$mean),
  bind_rows(loo_cov_results),
  tibble(specification = "Fully adjusted (primary)",
         ate_mean = primary_results$ate_summary$mean,
         cri_lo = primary_results$ate_summary$cri_lo,
         cri_hi = primary_results$ate_summary$cri_hi,
         p_ni = primary_results$ate_summary$p_ni,
         shift_from_primary = 0)
)

cat("\n=== Covariate Influence Summary ===\n")
print(loo_table, n = 20)

write.csv(loo_table,
          file.path(paths$tables, "table8_covariate_influence.csv"),
          row.names = FALSE)

# =============================================================================
# PART 3: PS TRIMMING SENSITIVITY
# =============================================================================

cat("\n--- Part 3: PS Trimming Sensitivity ---\n")
cat("  (Data restricted to overlap period in script 01; see eMethods 7)\n")

# PS trimming
cat("\n  PS trimming (0.025 < PS < 0.975):\n")

ps_formula <- as.formula(paste("treatment_num ~",
  "age_z + sex + bmi_z + smoking + education + employed_baseline +",
  "analgesic_baseline + odi_baseline_z + eq5d_baseline_z +",
  "nrs_back_baseline_z + nrs_leg_baseline_z +",
  "symptom_duration_back + symptom_duration_leg + motor_deficit +",
  "asa_cat + depression_anxiety + chronic_pain +",
  "spondylolisthesis + scoliosis +",
  "prior_surgery_any + n_prior_surgeries_z + multilevel +",
  "prolapse_intraforaminal + prolapse_extralateral + stenosis_central"))

ps_fit <- glm(ps_formula, data = df_disc, family = binomial)
df_disc$ps <- predict(ps_fit, type = "response")

df_trimmed <- df_disc %>% filter(ps > 0.025 & ps < 0.975)
cat(sprintf("  After trimming: %d (ELD=%d, MSD=%d), dropped %d\n",
            nrow(df_trimmed),
            sum(df_trimmed$treatment == "ELD"),
            sum(df_trimmed$treatment == "MSD"),
            nrow(df_disc) - nrow(df_trimmed)))

df_trimmed <- standardize_covs(df_trimmed)

fit_trimmed <- brm(
  bf(as.formula(paste("odi_3m | mi() ~ treatment +", cov_full))),
  data = df_trimmed,
  family = gaussian(),
  prior = shared_priors,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.95),
  file = file.path(paths$models, "fit_ps_trimmed"),
  file_refit = "on_change"
)

ate_trimmed <- compute_gcomp_ate(fit_trimmed, df_trimmed, "treatment", "continuous",
                                  lower_is_better = TRUE)
ate_trimmed_sum <- summarize_ate(ate_trimmed$ate, ni_margin = ni_margins$odi)

cat(sprintf("  PS-trimmed: ATE = %.2f [%.2f, %.2f], P(NI) = %.3f\n",
            ate_trimmed_sum$mean, ate_trimmed_sum$cri_lo,
            ate_trimmed_sum$cri_hi, ate_trimmed_sum$p_ni))

# =============================================================================
# PART 4: ADD MISSING SAP COVARIATES
# =============================================================================

cat("\n--- Part 4: Model with Additional SAP Covariates ---\n")
cat("Adding prolapse location and specific operative levels per SAP.\n")

# Extended covariate set = full model + specific operative levels
cov_extended <- paste(cov_string, "+ level_L45 + level_L5S1 + level_L34")

fit_extended <- brm(
  bf(as.formula(paste("odi_3m | mi() ~ treatment +", cov_extended))),
  data = df_disc,
  family = gaussian(),
  prior = shared_priors,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.95),
  file = file.path(paths$models, "fit_extended_covariates_nocal"),
  file_refit = "on_change"
)

ate_extended <- compute_gcomp_ate(fit_extended, df_disc, "treatment", "continuous",
                                   lower_is_better = TRUE)
ate_extended_sum <- summarize_ate(ate_extended$ate, ni_margin = ni_margins$odi)

cat(sprintf("  Extended covariates: ATE = %.2f [%.2f, %.2f], P(NI) = %.3f\n",
            ate_extended_sum$mean, ate_extended_sum$cri_lo,
            ate_extended_sum$cri_hi, ate_extended_sum$p_ni))

# =============================================================================
# PART 5: SUMMARY FIGURE
# =============================================================================

cat("\n--- Part 5: Summary Figure ---\n")

all_specs <- bind_rows(
  tibble(spec = "Unadjusted", group = "Reference",
         ate = ate_unadj_sum$mean, lo = ate_unadj_sum$cri_lo,
         hi = ate_unadj_sum$cri_hi, p_ni = ate_unadj_sum$p_ni),
  tibble(spec = "Fully adjusted (primary)", group = "Reference",
         ate = primary_results$ate_summary$mean,
         lo = primary_results$ate_summary$cri_lo,
         hi = primary_results$ate_summary$cri_hi,
         p_ni = primary_results$ate_summary$p_ni),
  map_dfr(loo_cov_results, function(r) {
    tibble(spec = r$specification, group = "Drop covariate",
           ate = r$ate_mean, lo = r$cri_lo, hi = r$cri_hi, p_ni = r$p_ni)
  }),
  tibble(spec = "PS trimming", group = "Overlap",
         ate = ate_trimmed_sum$mean, lo = ate_trimmed_sum$cri_lo,
         hi = ate_trimmed_sum$cri_hi, p_ni = ate_trimmed_sum$p_ni),
  tibble(spec = "Extended SAP covariates", group = "Covariates",
         ate = ate_extended_sum$mean, lo = ate_extended_sum$cri_lo,
         hi = ate_extended_sum$cri_hi, p_ni = ate_extended_sum$p_ni)
)

all_specs$spec <- factor(all_specs$spec, levels = rev(all_specs$spec))

p_influence <- ggplot(all_specs, aes(x = ate, y = spec, color = group)) +
  geom_vline(xintercept = 0, color = "grey30") +
  geom_vline(xintercept = -ni_margins$odi, linetype = "dashed", color = "red",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(
    "Reference" = "#333333",
    "Drop covariate" = "#4E79A7",
    "Overlap" = "#59A14F",
    "Covariates" = "#F28E2B"
  )) +
  annotate("text", x = -ni_margins$odi - 0.3, y = 0.5, label = "NI margin",
           color = "red", size = 2.5, hjust = 1) +
  labs(
    x = expression(paste(Delta, " (positive = ELD superior, ODI points)")),
    y = NULL,
    color = "Analysis type",
    title = "Covariate Influence on Treatment Effect Estimate",
    subtitle = "How does each model specification change the ATE?"
  ) +
  theme(legend.position = "bottom")

save_fig(p_influence, "covariate_influence_forest.png", width = 10, height = 7)

# =============================================================================
# SAVE
# =============================================================================

write.csv(all_specs %>% mutate(spec = as.character(spec)),
          file.path(paths$tables, "table8_covariate_influence.csv"),
          row.names = FALSE)

covariate_influence <- list(
  unadjusted = ate_unadj_sum,
  loo_covariates = loo_cov_results,
  ps_trimmed = ate_trimmed_sum,
  extended_covariates = ate_extended_sum,
  summary_table = all_specs
)
saveRDS(covariate_influence, file.path(paths$output, "covariate_influence_results.rds"))

cat("\n=== Covariate influence analysis complete ===\n")
