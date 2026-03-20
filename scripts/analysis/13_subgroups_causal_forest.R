# =============================================================================
# ENDO-LUMBAR: 13 Subgroup Analyses and Causal Forest
# =============================================================================

source("/Users/cjsogn/ENDO_LUMBAR/scripts/analysis/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Subgroup Analyses and Causal Forest (Exploratory) ===\n")

# Prepare data (standardize_covs and cov_string from 00_config.R)
df_disc <- standardize_covs(df_disc)

df_disc_full <- df_disc

# ZIB transformation for ODI; filter to complete cases (mi() not supported for ZIB)
df_disc <- df_disc %>% filter(!is.na(odi_3m))
df_disc$odi_3m_zib <- transform_for_zib(df_disc$odi_3m, upper = 100)
cat(sprintf("Complete cases for subgroup analysis: %d / %d total\n", nrow(df_disc), nrow(df_disc_full)))

# =============================================================================
# PRE-SPECIFIED SUBGROUPS
# =============================================================================

cat("\n--- Pre-Specified Subgroup Analyses ---\n")

# Define subgroups
# Create subgroup variables on both filtered and full datasets
df_disc <- df_disc %>%
  mutate(
    subgroup_age = ifelse(age < 50, "<50", ">=50"),
    subgroup_severity = ifelse(odi_baseline < 40, "ODI<40", "ODI>=40"),
    subgroup_symptom = ifelse(symptom_duration_leg <= 3, "Short", "Long"),
    subgroup_levels = ifelse(multilevel == 1, "Multilevel", "Single")
  )
df_disc_full <- df_disc_full %>%
  mutate(
    subgroup_age = ifelse(age < 50, "<50", ">=50"),
    subgroup_severity = ifelse(odi_baseline < 40, "ODI<40", "ODI>=40"),
    subgroup_symptom = ifelse(symptom_duration_leg <= 3, "Short", "Long"),
    subgroup_levels = ifelse(multilevel == 1, "Multilevel", "Single")
  )

subgroups <- list(
  list(var = "subgroup_age", label = "Age (<50 vs >=50)", cut = 50),
  list(var = "subgroup_severity", label = "Baseline severity (ODI <40 vs >=40)", cut = 40),
  list(var = "subgroup_symptom", label = "Symptom duration (Short vs Long)", cut = 3),
  list(var = "subgroup_levels", label = "Levels (Single vs Multi)", cut = NA)
)

subgroup_results <- list()

for (sg in subgroups) {
  cat(sprintf("\n  %s:\n", sg$label))

  # Add interaction term (ZIB model for ODI)
  interaction_formula <- bf(
    as.formula(paste("odi_3m_zib ~ treatment *", sg$var, "+", cov_string)),
    zi ~ treatment + odi_baseline_z
  )

  fit_sg <- brm(
    formula = interaction_formula,
    data = df_disc,
    family = zero_inflated_beta(),
    prior = c(
      prior(normal(0, 1), class = "b"),
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
    control = list(adapt_delta = 0.95),
    file = file.path(paths$models, paste0("fit_subgroup_zib_", sg$var)),
    file_refit = "on_change"
  )

  # Compute ATE in each subgroup (scale ZIB -> ODI)
  for (level in unique(df_disc[[sg$var]])) {
    df_sub_full <- df_disc_full %>% filter(.data[[sg$var]] == level)
    ate_sub <- compute_gcomp_ate(fit_sg, df_sub_full, "treatment", "continuous",
                                  lower_is_better = TRUE, scale_factor = 100)
    ate_sub_summary <- summarize_ate(ate_sub$ate, ni_margin = ni_margins$odi)
    cat(sprintf("    %s (n=%d): ATE=%.2f [%.2f, %.2f], P(NI)=%.3f\n",
                level, nrow(df_sub_full), ate_sub_summary$mean,
                ate_sub_summary$cri_lo, ate_sub_summary$cri_hi,
                ate_sub_summary$p_ni))
  }

  # Extract interaction coefficient (heterogeneity test)
  interaction_coef <- paste0("treatmentELD:", sg$var,
                              unique(df_disc[[sg$var]])[2])
  draws_int <- tryCatch(
    posterior::as_draws_df(fit_sg)[[paste0("b_", interaction_coef)]],
    error = function(e) NULL
  )

  if (!is.null(draws_int)) {
    p_het <- 2 * min(mean(draws_int > 0), mean(draws_int < 0))
    cat(sprintf("    Interaction P(heterogeneity): %.3f\n", 1 - p_het))
  }

  subgroup_results[[sg$var]] <- list(label = sg$label, fit = fit_sg)
}

# =============================================================================
# CAUSAL FOREST
# =============================================================================

cat("\n--- Causal Forest ---\n")
cat("Package: grf, num.trees=4000, honesty=TRUE\n")
cat("Note: Restricted to patients with observed ODI 3m (complete-case for grf).\n")

# Prepare complete-case data for causal forest
df_cc <- df_disc %>% filter(!is.na(odi_3m))
cat(sprintf("Complete cases for causal forest: %d\n", nrow(df_cc)))

# Covariate matrix (all baseline covariates, median/mode imputed)
cf_covariates <- c(
  "age", "bmi", "odi_baseline", "eq5d_baseline",
  "nrs_back_baseline", "nrs_leg_baseline",
  "symptom_duration_back", "symptom_duration_leg",
  "n_prior_surgeries"
)
cf_factors <- c("sex", "smoking", "education", "employed_baseline",
                "sick_leave", "disability", "analgesic_baseline",
                "asa_cat", "motor_deficit", "depression_anxiety", "chronic_pain",
                "spondylolisthesis", "scoliosis", "prior_surgery_any", "multilevel",
                "prolapse_intraforaminal", "prolapse_extralateral", "stenosis_central")

# Create numeric matrix for grf
# Ensure no NAs remain in covariates (grf does not handle missing values)
cf_vars <- c(cf_covariates, cf_factors)
cf_complete <- complete.cases(df_cc[, cf_vars, drop = FALSE])
if (sum(!cf_complete) > 0) {
  cat(sprintf("  Dropping %d rows with NA covariates for causal forest.\n",
              sum(!cf_complete)))
  df_cc <- df_cc[cf_complete, ]
}

X_cf <- model.matrix(
  as.formula(paste("~", paste(cf_vars, collapse = " + "))),
  data = df_cc
)[, -1]  # remove intercept

W <- as.numeric(df_cc$treatment == "ELD")
Y <- as.numeric(df_cc$odi_3m)

# Fit causal forest
set.seed(mcmc_settings$seed)
cf <- causal_forest(
  X = X_cf,
  Y = Y,
  W = W,
  num.trees = cf_settings$num.trees,
  min.node.size = cf_settings$min.node.size,
  honesty = cf_settings$honesty,
  tune.parameters = "all"
)

# Average treatment effect
ate_cf <- average_treatment_effect(cf, target.sample = "all")
cat(sprintf("\nCausal forest ATE: %.2f (SE: %.2f)\n", ate_cf[1], ate_cf[2]))

# Negate for sign convention: positive = ELD superior (lower ODI)
cat(sprintf("Adjusted for sign convention: %.2f\n", -ate_cf[1]))

# =============================================================================
# Calibration test
# =============================================================================

cal_test <- test_calibration(cf)
cat("\nCalibration test:\n")
print(cal_test)

# Column name from test_calibration is "Pr(>t)"
p_col <- grep("Pr", colnames(cal_test), value = TRUE)[1]
if (!is.null(p_col) && !is.na(p_col) && cal_test[2, p_col] < 0.05) {
  cat(sprintf("NOTE: Significant heterogeneity detected (p = %.4f).\n", cal_test[2, p_col]))
}

# =============================================================================
# Best linear projection
# =============================================================================

cat("\nBest Linear Projection:\n")
blp_vars <- c("age", "odi_baseline", "symptom_duration_leg")
blp_idx <- which(colnames(X_cf) %in% blp_vars)

if (length(blp_idx) > 0) {
  blp <- best_linear_projection(cf, A = X_cf[, blp_idx])
  print(blp)
}

# =============================================================================
# Variable Importance
# =============================================================================

varimp <- variable_importance(cf)
varimp_df <- tibble(
  variable = colnames(X_cf),
  importance = as.vector(varimp)
) %>%
  arrange(desc(importance))

cat("\nVariable importance (top 10):\n")
print(varimp_df %>% head(10))

p_varimp <- ggplot(varimp_df %>% head(15),
                    aes(x = importance, y = reorder(variable, importance))) +
  geom_col(fill = tx_fills["ELD"], color = tx_colors["ELD"], width = 0.6) +
  labs(
    x = "Variable Importance",
    y = NULL,
    title = "Causal Forest: Variable Importance",
    subtitle = "Top 15 variables modifying treatment effect"
  )

save_fig(p_varimp, "causal_forest_varimp.png", width = 7, height = 6)

# =============================================================================
# CATE Distribution
# =============================================================================

cate_preds <- predict(cf)
# Negate for our sign convention (MSD - ELD)
cate_values <- -cate_preds[, 1]

df_cc$cate <- cate_values

p_cate <- ggplot(df_cc, aes(x = cate)) +
  geom_histogram(aes(fill = treatment), alpha = 0.5, bins = 30, position = "identity") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey30") +
  geom_vline(xintercept = -ni_margins$odi, linetype = "dashed", color = "red") +
  scale_fill_manual(values = tx_fills) +
  labs(
    x = expression(paste("CATE (", Delta, ", positive = ELD superior)")),
    y = "Count",
    title = "Distribution of Conditional Average Treatment Effects",
    subtitle = "Causal Forest Predictions"
  )

save_fig(p_cate, "causal_forest_cate.png", width = 7, height = 5)

# =============================================================================
# SAVE
# =============================================================================

cf_results <- list(
  causal_forest = cf,
  ate = ate_cf,
  calibration = cal_test,
  variable_importance = varimp_df,
  cate = cate_values,
  subgroup_results = subgroup_results
)
saveRDS(cf_results, file.path(paths$output, "causal_forest_results.rds"))

# =============================================================================
# PROJECTION PREDICTIVE VARIABLE SELECTION (projpred)
# =============================================================================

cat("\n--- Projection Predictive Variable Selection ---\n")

# Refit on complete cases (projpred does not support mi())
df_cc_proj <- df_disc %>% filter(!is.na(odi_3m))
df_cc_proj <- standardize_covs(df_cc_proj)

cat(sprintf("Refitting on complete cases for projpred: n=%d\n", nrow(df_cc_proj)))

fit_cc_proj <- brm(
  bf(as.formula(paste("odi_3m ~ treatment +", cov_string))),
  data = df_cc_proj,
  family = gaussian(),
  prior = priors_continuous,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.95),
  file = file.path(paths$models, "fit_projpred_cc"),
  file_refit = "on_change"
)

projpred_results <- tryCatch({
  ref_model <- get_refmodel(fit_cc_proj)

  cat("Running cv_varsel (this may take a while)...\n")
  vs <- cv_varsel(
    ref_model,
    method = "forward",
    cv_method = "LOO",
    nterms_max = min(20, length(all_model_covs)),
    seed = mcmc_settings$seed,
    verbose = TRUE
  )

  n_suggest <- suggest_size(vs, stat = "elpd")
  cat(sprintf("\nSuggested number of covariates: %d (out of %d)\n",
              n_suggest, length(all_model_covs)))

  sol_path <- solution_terms(vs)
  cat("\nVariable selection order:\n")
  for (i in seq_along(sol_path)) {
    cat(sprintf("  %d. %s\n", i, sol_path[i]))
  }
  cat(sprintf("\nMinimal sufficient set (%d vars): %s\n",
              n_suggest, paste(sol_path[1:n_suggest], collapse = ", ")))

  treatment_in_set <- "treatmentELD" %in% sol_path[1:n_suggest]
  cat(sprintf("Treatment (ELD) in minimal set: %s\n",
              ifelse(treatment_in_set, "YES", "NO")))

  p_projpred <- plot(vs, stats = "elpd") +
    geom_vline(xintercept = n_suggest, linetype = "dashed", color = "red") +
    labs(
      title = "Projection Predictive Variable Selection",
      subtitle = sprintf("Suggested subset: %d variables (red line)", n_suggest),
      x = "Number of Covariates",
      y = "ELPD (relative to full model)"
    )
  save_fig(p_projpred, "projpred_elpd.png", width = 8, height = 5)

  list(
    vs = vs,
    n_suggest = n_suggest,
    solution_path = sol_path,
    treatment_retained = treatment_in_set
  )
}, error = function(e) {
  cat(sprintf("  projpred error: %s\n", conditionMessage(e)))
  cat("  Saving partial results without projpred.\n")
  list(error = conditionMessage(e))
})

saveRDS(projpred_results, file.path(paths$output, "projpred_results.rds"))

cf_results$projpred <- projpred_results
saveRDS(cf_results, file.path(paths$output, "causal_forest_results.rds"))

cat("\nSubgroup, causal forest, and projpred analysis complete.\n")
