# =============================================================================
# ENDO-LUMBAR: 22a Export Data for PyMC Analysis
# Exports RDS datasets to CSV for Python/PyMC consumption
# =============================================================================

source("/Users/cjsogn/endo_studies/lumbar/analysis/scripts/00_config.R")

cat("=== Exporting Data for PyMC Analysis ===\n")

# Load both datasets
df_full <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
df_12m  <- readRDS(file.path(paths$data_clean, "df_disc_12m_eligible.rds"))

cat(sprintf("Full sample: N=%d (ELD=%d, MSD=%d)\n",
            nrow(df_full), sum(df_full$treatment == "ELD"),
            sum(df_full$treatment == "MSD")))
cat(sprintf("12m-eligible: N=%d (ELD=%d, MSD=%d)\n",
            nrow(df_12m), sum(df_12m$treatment == "ELD"),
            sum(df_12m$treatment == "MSD")))

# Variables to export
vars_to_export <- c(
  "treatment_num",
  # 27 covariates (matching all_model_covs from 00_config.R)
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
  stop("Missing variables in df_full: ", paste(missing_vars, collapse = ", "))
}
missing_vars_12m <- setdiff(vars_to_export, colnames(df_12m))
if (length(missing_vars_12m) > 0) {
  stop("Missing variables in df_12m: ", paste(missing_vars_12m, collapse = ", "))
}

# Ensure output directory exists
out_dir <- "/Users/cjsogn/ENDO_LUMBAR/data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Export CSVs
write.csv(df_full[, vars_to_export], file.path(out_dir, "df_disc_full.csv"),
          row.names = FALSE, na = "")
write.csv(df_12m[, vars_to_export], file.path(out_dir, "df_disc_12m.csv"),
          row.names = FALSE, na = "")

cat(sprintf("\nExported df_disc_full.csv: %d rows x %d cols\n",
            nrow(df_full), length(vars_to_export)))
cat(sprintf("Exported df_disc_12m.csv: %d rows x %d cols\n",
            nrow(df_12m), length(vars_to_export)))

# Report missingness for outcomes
cat("\nOutcome missingness (full sample):\n")
outcomes <- vars_to_export[grepl("_3m$|_12m$", vars_to_export)]
for (v in outcomes) {
  n_miss <- sum(is.na(df_full[[v]]))
  cat(sprintf("  %-25s: %d missing (%.1f%%)\n", v, n_miss,
              100 * n_miss / nrow(df_full)))
}

cat("\nExport complete.\n")
