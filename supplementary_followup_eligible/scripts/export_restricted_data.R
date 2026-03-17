#!/usr/bin/env Rscript
# =============================================================================
# Export follow-up-eligible restricted samples for supplementary analysis
#
# 3m-eligible: surgery_date <= 2025-09-30 (at least 3 months before data cutoff)
# 12m-eligible: surgery_date <= 2024-12-31 (at least 12 months before data cutoff)
# Data cutoff: ~2025-12-31 (end of study period)
# =============================================================================

cat("=== Exporting Follow-Up Eligible Restricted Samples ===\n")

# Load the imputed analysis dataset (same source as main analysis)
paths <- list(data_clean = "/Users/cjsogn/endo_studies/lumbar/analysis/data")
df_full <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))

cat(sprintf("Full sample: N=%d (ELD=%d, MSD=%d)\n",
            nrow(df_full), sum(df_full$treatment == "ELD"),
            sum(df_full$treatment == "MSD")))
cat(sprintf("Surgery date range: %s to %s\n",
            as.character(min(df_full$surgery_date)),
            as.character(max(df_full$surgery_date))))

# Apply date restrictions
cutoff_3m  <- as.Date("2025-09-30")
cutoff_12m <- as.Date("2024-12-31")

df_3m  <- df_full[df_full$surgery_date <= cutoff_3m, ]
df_12m <- df_full[df_full$surgery_date <= cutoff_12m, ]

cat(sprintf("\n3m-eligible (surgery <= %s): N=%d (ELD=%d, MSD=%d)\n",
            as.character(cutoff_3m), nrow(df_3m),
            sum(df_3m$treatment == "ELD"), sum(df_3m$treatment == "MSD")))
cat(sprintf("12m-eligible (surgery <= %s): N=%d (ELD=%d, MSD=%d)\n",
            as.character(cutoff_12m), nrow(df_12m),
            sum(df_12m$treatment == "ELD"), sum(df_12m$treatment == "MSD")))

# Same variables as main analysis export (22a_export_data_for_pymc.R)
vars_to_export <- c(
  "treatment_num",
  # 27 covariates
  "age", "sex", "bmi", "smoking", "education", "employed_baseline",
  "sick_leave", "disability", "analgesic_baseline",
  "odi_baseline", "eq5d_baseline",
  "nrs_back_baseline", "nrs_leg_baseline",
  "symptom_duration_back", "symptom_duration_leg", "motor_deficit",
  "asa_cat", "depression_anxiety", "chronic_pain",
  "spondylolisthesis", "scoliosis",
  "prior_surgery_any", "n_prior_surgeries", "multilevel",
  "prolapse_intraforaminal", "prolapse_extralateral", "stenosis_central",
  # 17 outcomes
  "odi_3m", "odi_12m",
  "nrs_back_3m", "nrs_back_12m",
  "nrs_leg_3m", "nrs_leg_12m",
  "eq5d_3m", "eq5d_12m",
  "responder_3m",
  "rtw_3m", "rtw_12m",
  "analgesic_3m", "analgesic_12m",
  "satisfied_3m", "satisfied_12m",
  "gpe_success_3m", "gpe_success_12m"
)

# Check all variables exist
missing_vars <- setdiff(vars_to_export, colnames(df_full))
if (length(missing_vars) > 0) {
  stop("Missing variables: ", paste(missing_vars, collapse = ", "))
}

# Export
out_dir <- "/Users/cjsogn/ENDO_LUMBAR_FINAL/supplementary_followup_eligible/data"
write.csv(df_3m[, vars_to_export], file.path(out_dir, "df_3m_eligible.csv"),
          row.names = FALSE, na = "")
write.csv(df_12m[, vars_to_export], file.path(out_dir, "df_12m_eligible.csv"),
          row.names = FALSE, na = "")

# Report outcome missingness
cat("\n--- Outcome missingness: 3m-eligible sample ---\n")
outcomes_3m <- c("odi_3m", "nrs_back_3m", "nrs_leg_3m", "eq5d_3m",
                 "responder_3m", "rtw_3m", "analgesic_3m", "satisfied_3m",
                 "gpe_success_3m")
for (v in outcomes_3m) {
  n_miss <- sum(is.na(df_3m[[v]]))
  cat(sprintf("  %-25s: %d missing (%.1f%%)\n", v, n_miss,
              100 * n_miss / nrow(df_3m)))
}

cat("\n--- Outcome missingness: 12m-eligible sample ---\n")
outcomes_12m <- c("odi_12m", "nrs_back_12m", "nrs_leg_12m", "eq5d_12m",
                  "rtw_12m", "analgesic_12m", "satisfied_12m",
                  "gpe_success_12m")
for (v in outcomes_12m) {
  n_miss <- sum(is.na(df_12m[[v]]))
  cat(sprintf("  %-25s: %d missing (%.1f%%)\n", v, n_miss,
              100 * n_miss / nrow(df_12m)))
}

# Also report follow-up rates by treatment for comparison with main analysis
cat("\n--- Follow-up rates by treatment ---\n")
for (tx in c("ELD", "MSD")) {
  d3 <- df_3m[df_3m$treatment == tx, ]
  d12 <- df_12m[df_12m$treatment == tx, ]
  cat(sprintf("%s 3m-eligible (N=%d): ODI 3m observed = %d (%.1f%%)\n",
              tx, nrow(d3), sum(!is.na(d3$odi_3m)),
              100 * mean(!is.na(d3$odi_3m))))
  cat(sprintf("%s 12m-eligible (N=%d): ODI 12m observed = %d (%.1f%%)\n",
              tx, nrow(d12), sum(!is.na(d12$odi_12m)),
              100 * mean(!is.na(d12$odi_12m))))
}

cat("\nExport complete.\n")
