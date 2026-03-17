# =============================================================================
# ENDO-LUMBAR: 12 Falsification Tests and E-Value
# SAP Sections 20, 21
# =============================================================================

source("/Users/cjsogn/endo_studies/lumbar/analysis/scripts/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Falsification Tests and E-Value ===\n")

# =============================================================================
# 20.2 COVARIATE FALSIFICATION TEST
# =============================================================================

cat("\n--- 20.2 Covariate Falsification Test ---\n")
cat("After conditioning on all other covariates, treatment should not predict\n")
cat("any single baseline covariate.\n")

# All model covariates from 00_config.R (all_model_covs)

# Test covariates: these should not predict treatment after conditioning on all others
test_covs <- c("sex", "bmi", "education")

falsification_results <- map_dfr(test_covs, function(test_var) {
  # Adjustment set: all model covariates except the tested one
  adj_vars <- setdiff(all_model_covs, test_var)

  # Full model (with test variable) and reduced model (without)
  fml_full <- as.formula(paste("treatment_num ~", test_var, "+",
                               paste(adj_vars, collapse = " + ")))
  fml_reduced <- as.formula(paste("treatment_num ~",
                                  paste(adj_vars, collapse = " + ")))

  fit_full <- glm(fml_full, data = df_disc, family = binomial)
  fit_reduced <- glm(fml_reduced, data = df_disc, family = binomial)

  # Likelihood ratio test (handles both numeric and factor variables)
  lr_test <- anova(fit_reduced, fit_full, test = "LRT")
  lr_stat <- abs(lr_test$Deviance[2])
  lr_df <- lr_test$Df[2]
  p_value <- lr_test$`Pr(>Chi)`[2]

  # Report coefficient for numeric/binary; max absolute coefficient for factors
  if (is.numeric(df_disc[[test_var]]) || is.integer(df_disc[[test_var]])) {
    coef_test <- summary(fit_full)$coefficients[test_var, ]
    coefficient <- coef_test["Estimate"]
    std_coef <- coefficient * sd(as.numeric(df_disc[[test_var]]), na.rm = TRUE)
  } else {
    # Factor variable: report maximum absolute coefficient among levels
    test_coef_names <- grep(paste0("^", test_var), names(coef(fit_full)), value = TRUE)
    coefficients_vec <- coef(fit_full)[test_coef_names]
    coefficient <- coefficients_vec[which.max(abs(coefficients_vec))]
    std_coef <- coefficient  # already on log-odds scale for factors
  }

  concern <- p_value < 0.05 | abs(std_coef) > 0.1

  cat(sprintf("  %s: LRT chi-sq=%.2f (df=%d), p=%.4f, |max coef|=%.3f",
              test_var, lr_stat, lr_df, p_value, abs(coefficient)))
  if (concern) cat(" -> CONCERN") else cat(" -> OK")
  cat("\n")

  tibble(
    covariate = test_var,
    lr_chi_sq = lr_stat,
    lr_df = lr_df,
    p_value = p_value,
    max_coefficient = coefficient,
    std_coefficient = std_coef,
    concern = concern
  )
})

n_failures <- sum(falsification_results$concern)
cat(sprintf("\nCovariates with concern (p<0.05 or |coef|>0.1): %d/%d\n",
            n_failures, nrow(falsification_results)))
if (n_failures > 1) {
  cat("WARNING: Multiple failures suggest model misspecification or residual confounding.\n")
}

write.csv(falsification_results,
          file.path(paths$tables, "falsification_test.csv"),
          row.names = FALSE)

# =============================================================================
# 21. E-VALUE FOR UNMEASURED CONFOUNDING
# =============================================================================

cat("\n--- 21. E-Value Analysis ---\n")

ate_mean <- primary_results$ate_summary$mean
ate_cri_lo <- primary_results$ate_summary$cri_lo
ate_cri_hi <- primary_results$ate_summary$cri_hi

# Get pooled SD of ODI at 3 months
pooled_sd <- sd(df_disc$odi_3m, na.rm = TRUE)
cat(sprintf("  Pooled SD of ODI 3m: %.1f\n", pooled_sd))

# Convert ATE to standardized mean difference
smd <- ate_mean / pooled_sd
cat(sprintf("  ATE = %.2f -> SMD = %.3f\n", ate_mean, smd))

# E-value computation using the EValue package
# For continuous outcomes: convert SMD to approximate RR
# RR = exp(0.91 * SMD) is the square-root conversion
rr_est <- exp(0.91 * abs(smd))

# Compute E-value directly: E = RR + sqrt(RR * (RR - 1))
compute_evalue <- function(rr) {
  rr <- max(rr, 1/rr)  # ensure RR >= 1
  rr + sqrt(rr * (rr - 1))
}

evalue_point_val <- compute_evalue(rr_est)

cat(sprintf("  Approximate RR: %.3f\n", rr_est))
cat(sprintf("  E-value (point estimate): %.2f\n", evalue_point_val))

# E-value for CrI bound closest to null
cri_bound <- ifelse(abs(ate_cri_lo) < abs(ate_cri_hi), ate_cri_lo, ate_cri_hi)
smd_bound <- cri_bound / pooled_sd
rr_bound <- exp(0.91 * abs(smd_bound))
evalue_bound_val <- compute_evalue(rr_bound)

cat(sprintf("  E-value (CrI bound closest to null): %.2f\n", evalue_bound_val))
cat(sprintf("  CrI bound used: %.2f (SMD = %.3f)\n", cri_bound, smd_bound))

# Interpretation
cat("\n  Interpretation:\n")
if (evalue_point_val > 2.0) {
  cat(sprintf("  E-value = %.2f > 2.0: Reassuring against unmeasured confounding.\n",
              evalue_point_val))
  cat("  An unmeasured confounder would need to more than double both\n")
  cat("  P(ELD) and P(poor outcome) to explain away the observed effect.\n")
} else {
  cat(sprintf("  E-value = %.2f <= 2.0: Moderate vulnerability to confounding.\n",
              evalue_point_val))
}

# E-values for binary Tier 3 outcomes
cat("\n  E-values for Tier 3 binary outcomes:\n")
tier3_results <- tryCatch(
  readRDS(file.path(paths$output, "tier3_results.rds")),
  error = function(e) NULL
)

if (!is.null(tier3_results)) {
  for (nm in names(tier3_results)) {
    r <- tier3_results[[nm]]
    if (!isTRUE(r$descriptive_only) && !is.null(r$ate_summary)) {
      rd <- r$ate_summary$mean
      # Get baseline risk (MSD event rate) from data for RD-to-RR conversion
      if (nm %in% names(df_disc)) {
        p0 <- mean(df_disc[[nm]][df_disc$treatment == "MSD"], na.rm = TRUE)
      } else {
        p0 <- 0.5  # conservative fallback if column not found
      }
      # Convert RD to RR using baseline risk: RR = (p0 + RD) / p0
      if (p0 > 0.01 && p0 < 0.99) {
        rr_tier3 <- (p0 + rd) / p0
      } else {
        rr_tier3 <- 1 + 2 * rd  # fallback for extreme rates (p0~0.5)
      }
      rr_tier3 <- max(rr_tier3, 1 / rr_tier3)  # ensure >= 1
      ev_tier3 <- compute_evalue(rr_tier3)
      cat(sprintf("  %s: RD = %.3f, p0(MSD) = %.3f, RR = %.2f, E-value = %.2f\n",
                  r$label, rd, p0, rr_tier3, ev_tier3))
    }
  }
}

# =============================================================================
# SAVE
# =============================================================================

evalue_results <- list(
  primary = list(
    ate_mean = ate_mean,
    pooled_sd = pooled_sd,
    smd = smd,
    rr_est = rr_est,
    evalue_point = evalue_point_val,
    evalue_bound = evalue_bound_val
  ),
  falsification = falsification_results
)
saveRDS(evalue_results, file.path(paths$output, "evalue_results.rds"))

cat("\nFalsification and E-value analysis complete.\n")
