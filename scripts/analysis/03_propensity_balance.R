# =============================================================================
# ENDO-LUMBAR: 03 Propensity Score and Covariate Balance (Descriptive)
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
var_meta <- readRDS(file.path(paths$data_clean, "var_meta.rds"))

cat("=== Propensity Score and Covariate Balance (Descriptive Only) ===\n")
cat("Note: PS analysis is descriptive only.\n\n")

# =============================================================================
# 14.1 Propensity Score Estimation
# =============================================================================

# Covariates for PS model from 00_config.R (ps_covariates)

# Build PS formula
ps_formula <- as.formula(
  paste("treatment_num ~", paste(ps_covariates, collapse = " + "))
)

# Fit logistic regression
ps_model <- glm(ps_formula, data = df_disc, family = binomial)
df_disc$ps <- predict(ps_model, type = "response")

cat(sprintf("PS model fitted: %d covariates\n", length(ps_covariates)))
cat(sprintf("PS range: [%.4f, %.4f]\n", min(df_disc$ps), max(df_disc$ps)))

# =============================================================================
# 14.2 Overlap Assessment
# =============================================================================

# PS distribution by treatment group
p_overlap <- ggplot(df_disc, aes(x = ps, fill = treatment, color = treatment)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  scale_fill_manual(values = tx_fills) +
  scale_color_manual(values = tx_colors) +
  labs(
    x = "Propensity Score P(ELD | W)",
    y = "Density",
    title = "Propensity Score Distributions",
    subtitle = "Disc herniation population (descriptive only)"
  ) +
  theme(legend.position = c(0.8, 0.8))

save_fig(p_overlap, "ps_overlap_density.png", width = 7, height = 5)

# Histogram version
p_overlap_hist <- ggplot(df_disc, aes(x = ps, fill = treatment)) +
  geom_histogram(
    data = filter(df_disc, treatment == "MSD"),
    aes(y = after_stat(density)),
    alpha = 0.5, bins = 30
  ) +
  geom_histogram(
    data = filter(df_disc, treatment == "ELD"),
    aes(y = -after_stat(density)),
    alpha = 0.5, bins = 30
  ) +
  scale_fill_manual(values = tx_fills) +
  labs(
    x = "Propensity Score P(ELD | W)",
    y = "Density (MSD above, ELD below)",
    title = "Propensity Score Mirror Histogram"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5)

save_fig(p_overlap_hist, "ps_overlap_mirror.png", width = 7, height = 5)

# Effective overlap (proportion with PS in [0.025, 0.975])
n_overlap <- sum(df_disc$ps >= 0.025 & df_disc$ps <= 0.975)
pct_overlap <- 100 * n_overlap / nrow(df_disc)
cat(sprintf("\nEffective overlap (PS in [0.025, 0.975]): %d/%d (%.1f%%)\n",
            n_overlap, nrow(df_disc), pct_overlap))

# Patients with extreme PS
n_extreme_low <- sum(df_disc$ps < 0.1)
n_extreme_high <- sum(df_disc$ps > 0.9)
cat(sprintf("Extreme PS: <0.1: %d, >0.9: %d\n", n_extreme_low, n_extreme_high))

if (n_extreme_high > 0) {
  cat("\nCharacteristics of patients with PS > 0.9:\n")
  extreme_high <- df_disc %>% filter(ps > 0.9)
  print(extreme_high %>% dplyr::select(treatment, age, sex, odi_baseline, calendar_time, ps))
}

# Flag poor overlap
if (pct_overlap < 90) {
  cat("\nWARNING: >10% of sample in PS violation zone.\n")
}

# =============================================================================
# 14.3 Covariate Balance (SMD)
# =============================================================================

# Compute unadjusted SMDs using cobalt
bal <- bal.tab(
  ps_formula,
  data = df_disc,
  treat = df_disc$treatment_num,
  binary = "std",
  continuous = "std",
  s.d.denom = "pooled"
)

cat("\n=== Covariate Balance (Standardized Mean Differences) ===\n")
print(bal)

# Extract SMD values for Love plot
smd_df <- bal$Balance %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  dplyr::rename(smd = Diff.Un) %>%
  filter(variable != "distance") %>%
  mutate(
    abs_smd = abs(smd),
    balance_status = case_when(
      abs_smd < 0.1 ~ "Good (|SMD| < 0.1)",
      abs_smd < 0.25 ~ "Moderate (0.1-0.25)",
      TRUE ~ "Substantial (|SMD| > 0.25)"
    )
  ) %>%
  arrange(desc(abs_smd))

cat("\nVariables with |SMD| > 0.1:\n")
print(smd_df %>% filter(abs_smd > 0.1) %>% dplyr::select(variable, smd, abs_smd))

# Love plot
p_love <- ggplot(smd_df, aes(x = abs_smd, y = reorder(variable, abs_smd))) +
  geom_point(aes(color = balance_status), size = 2.5) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.25, linetype = "dotted", color = "red") +
  scale_color_manual(
    values = c(
      "Good (|SMD| < 0.1)" = "#2166AC",
      "Moderate (0.1-0.25)" = "#F28E2B",
      "Substantial (|SMD| > 0.25)" = "#E15759"
    ),
    name = "Balance"
  ) +
  labs(
    x = "Absolute Standardized Mean Difference",
    y = NULL,
    title = "Covariate Balance: Love Plot (Unadjusted)",
    subtitle = "Disc herniation: ELD vs MSD"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))

save_fig(p_love, "love_plot_unadjusted.png", width = 8, height = 7)

# =============================================================================
# 14.3 Variance Ratios (continuous covariates)
# =============================================================================

cont_covs <- c("age", "bmi", "odi_baseline", "eq5d_baseline",
               "nrs_back_baseline", "nrs_leg_baseline",
               "n_prior_surgeries")

vr_df <- map_dfr(cont_covs, function(v) {
  var_eld <- var(df_disc[[v]][df_disc$treatment == "ELD"], na.rm = TRUE)
  var_msd <- var(df_disc[[v]][df_disc$treatment == "MSD"], na.rm = TRUE)
  tibble(
    variable = v,
    var_ELD = var_eld,
    var_MSD = var_msd,
    variance_ratio = var_eld / var_msd,
    acceptable = variance_ratio >= 0.5 & variance_ratio <= 2.0
  )
})

cat("\n=== Variance Ratios (continuous covariates) ===\n")
cat("Target: 0.5-2.0\n")
print(vr_df)

write.csv(vr_df, file.path(paths$tables, "variance_ratios_disc.csv"), row.names = FALSE)

# =============================================================================
# Save balance summary
# =============================================================================

balance_summary <- list(
  smd_df = smd_df,
  variance_ratios = vr_df,
  ps_summary = list(
    range = range(df_disc$ps),
    pct_overlap = pct_overlap,
    n_extreme_low = n_extreme_low,
    n_extreme_high = n_extreme_high
  )
)
saveRDS(balance_summary, file.path(paths$output, "balance_summary.rds"))

# Save PS-augmented data
saveRDS(df_disc, file.path(paths$data_clean, "df_disc_ps.rds"))

cat("\nPropensity score and balance analysis complete.\n")
