# =============================================================================
# ENDO-LUMBAR: 02 Descriptive Statistics (Table 1)
# Baseline characteristics by treatment group
# SAP Section 25.1 (Table 1)
# =============================================================================

source("/Users/cjsogn/endo_studies/lumbar/analysis/scripts/00_config.R")

# Load data
df_disc <- readRDS(file.path(paths$data_clean, "df_disc.rds"))
df_sten <- readRDS(file.path(paths$data_clean, "df_sten.rds"))

# =============================================================================
# TABLE 1: Disc herniation population
# =============================================================================

cat("=== Table 1: Disc Herniation Population ===\n")

# Variables for Table 1
table1_vars <- c(
  # Demographics
  "age", "sex", "bmi", "smoking",
  # Socioeconomic
  "education", "employed_baseline", "sick_leave", "disability",
  "analgesic_baseline",
  # Symptoms
  "odi_baseline", "eq5d_baseline", "nrs_back_baseline", "nrs_leg_baseline",
  "symptom_duration_back", "symptom_duration_leg", "motor_deficit",
  # Comorbidities
  "asa_cat", "depression_anxiety", "chronic_pain",
  # Pathology
  "prolapse_intraforaminal", "prolapse_extralateral",
  "spondylolisthesis", "scoliosis",
  # Surgical
  "prior_surgery", "prior_surgery_any", "n_prior_surgeries",
  "n_levels", "multilevel",
  "level_L45", "level_L5S1", "level_L34",
  "approach",
  "surgery_year"
)

# Identify factor vs numeric variables
cat_vars <- c("sex", "smoking", "education", "asa_cat",
              "prior_surgery", "approach", "surgery_year")
nonnormal_vars <- c("n_prior_surgeries", "symptom_duration_back", "symptom_duration_leg",
                     "n_levels", "los_total", "operating_time")

# Ensure surgery_year is factor for table
df_disc <- df_disc %>% mutate(surgery_year = factor(surgery_year))
df_sten <- df_sten %>% mutate(surgery_year = factor(surgery_year))

# Create Table 1 using tableone
tab1_disc <- CreateTableOne(
  vars = table1_vars,
  strata = "treatment",
  data = df_disc,
  factorVars = cat_vars,
  test = FALSE,
  smd = TRUE
)

# Print with SMD
tab1_disc_print <- print(tab1_disc,
                          smd = TRUE,
                          showAllLevels = TRUE,
                          noSpaces = TRUE,
                          printToggle = FALSE,
                          nonnormal = nonnormal_vars)

cat("\n--- Disc Herniation: Baseline Characteristics ---\n")
print(tab1_disc_print, quote = FALSE)

# Save Table 1 as CSV
write.csv(tab1_disc_print,
          file.path(paths$tables, "table1_disc_herniation.csv"))

# =============================================================================
# Follow-up rates by treatment group
# =============================================================================

cat("\n=== Follow-up Rates ===\n")
fu_table <- df_disc %>%
  group_by(treatment) %>%
  dplyr::summarise(
    n = n(),
    fu_3m_n = sum(fu_3m_completed == 1, na.rm = TRUE),
    fu_3m_pct = 100 * mean(fu_3m_completed == 1, na.rm = TRUE),
    fu_12m_n = sum(fu_12m_completed == 1, na.rm = TRUE),
    fu_12m_pct = 100 * mean(fu_12m_completed == 1, na.rm = TRUE),
    odi_3m_n = sum(!is.na(odi_3m)),
    odi_3m_pct = 100 * mean(!is.na(odi_3m)),
    odi_12m_n = sum(!is.na(odi_12m)),
    odi_12m_pct = 100 * mean(!is.na(odi_12m)),
    .groups = "drop"
  )
print(fu_table)

write.csv(fu_table, file.path(paths$tables, "followup_rates_disc.csv"), row.names = FALSE)

# =============================================================================
# Outcome summary by treatment group (observed only)
# =============================================================================

cat("\n=== Outcome Summary (Observed) ===\n")
outcome_summary <- df_disc %>%
  group_by(treatment) %>%
  dplyr::summarise(
    # Primary
    odi_3m_mean = mean(odi_3m, na.rm = TRUE),
    odi_3m_sd = sd(odi_3m, na.rm = TRUE),
    odi_3m_n = sum(!is.na(odi_3m)),
    # Secondary continuous
    odi_12m_mean = mean(odi_12m, na.rm = TRUE),
    odi_12m_sd = sd(odi_12m, na.rm = TRUE),
    nrs_back_3m_mean = mean(nrs_back_3m, na.rm = TRUE),
    nrs_leg_3m_mean = mean(nrs_leg_3m, na.rm = TRUE),
    eq5d_3m_mean = mean(eq5d_3m, na.rm = TRUE),
    # Perioperative
    day_surgery_pct = 100 * mean(day_surgery == 1, na.rm = TRUE),
    los_mean = mean(los_total, na.rm = TRUE),
    los_median = median(los_total, na.rm = TRUE),
    op_time_mean = mean(operating_time, na.rm = TRUE),
    perop_comp_pct = 100 * mean(perop_comp_any == 1, na.rm = TRUE),
    .groups = "drop"
  )
print(outcome_summary)

write.csv(outcome_summary, file.path(paths$tables, "outcome_summary_disc.csv"), row.names = FALSE)

# =============================================================================
# LTFU analysis: compare responders vs non-responders at baseline
# =============================================================================

cat("\n=== LTFU Pattern Analysis ===\n")
df_disc <- df_disc %>%
  mutate(ltfu_3m = as.integer(is.na(odi_3m)))

ltfu_table <- CreateTableOne(
  vars = c("age", "sex", "bmi", "odi_baseline", "eq5d_baseline",
           "asa_cat", "depression_anxiety", "treatment"),
  strata = "ltfu_3m",
  data = df_disc,
  factorVars = c("sex", "asa_cat", "treatment"),
  test = TRUE
)
ltfu_print <- print(ltfu_table, smd = TRUE, printToggle = FALSE, test = TRUE)
cat("\nBaseline comparison: Responders (0) vs LTFU (1) at 3 months\n")
print(ltfu_print, quote = FALSE)

write.csv(ltfu_print, file.path(paths$tables, "ltfu_analysis_disc.csv"))

# =============================================================================
# TABLE 1: Stenosis population (exploratory)
# =============================================================================

cat("\n=== Table 1: Stenosis Population (Exploratory) ===\n")

# Only if there are enough patients
if (sum(df_sten$treatment == "ELD") >= 5) {
  tab1_sten <- CreateTableOne(
    vars = table1_vars,
    strata = "treatment",
    data = df_sten,
    factorVars = cat_vars,
    test = FALSE,
    smd = TRUE
  )
  tab1_sten_print <- print(tab1_sten, smd = TRUE, showAllLevels = TRUE,
                            printToggle = FALSE)
  print(tab1_sten_print, quote = FALSE)
  write.csv(tab1_sten_print,
            file.path(paths$tables, "table1_stenosis.csv"))
} else {
  cat("  ELD stenosis group too small for Table 1\n")
}

cat("\nTable 1 generation complete.\n")
