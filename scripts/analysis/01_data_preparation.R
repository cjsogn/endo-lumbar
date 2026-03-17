# =============================================================================
# ENDO-LUMBAR: 01 Data Preparation
# Load, clean, define populations, derive variables
# SAP Sections: 3, 10, 11
# =============================================================================

source("/Users/cjsogn/endo_studies/lumbar/analysis/scripts/00_config.R")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading data...\n")
raw <- read_sav(paths$data_raw, encoding = "latin1")
cat(sprintf("  Raw data: %d patients, %d variables\n", nrow(raw), ncol(raw)))

# =============================================================================
# 2. DEFINE TREATMENT GROUPS (SAP Section 3)
# =============================================================================
# OpMikroV3: 1=Mikroskopi, 2=Lupebriller, 3=Endoskopi, 0=Nei, 9=Ikke utfylt
# ELD = Endoscopic (OpMikroV3 == 3)
# MSD = Microsurgical (OpMikroV3 == 1, i.e., microscope-assisted)
# Exclude: Loupes (2), None (0), Not filled (9)

df <- raw %>%
  mutate(
    treatment = case_when(
      OpMikroV3 == 3 ~ "ELD",
      OpMikroV3 == 1 ~ "MSD",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(treatment))

cat(sprintf("  After treatment group filter: %d patients (ELD=%d, MSD=%d)\n",
            nrow(df), sum(df$treatment == "ELD"), sum(df$treatment == "MSD")))

# =============================================================================
# 3. DEFINE ANALYSIS POPULATIONS (SAP Section 3)
# =============================================================================
# HovedInngrepV2V3: 1=Prolaps kirurgi, 2=Midtlinje dekompresjon, 3=Laminektomi,
#                   5=Fusjonskirurgi, etc.
# ProlapsoprAlle: 1=prolapse surgery, 0=no
# LSSopr: 1=stenosis surgery, 0=no
# OpIndCauda: cauda equina indication

# Exclude cauda equina (SAP Section 3.2)
df <- df %>% filter(OpIndCauda != 1 | is.na(OpIndCauda))
cat(sprintf("  After cauda equina exclusion: %d\n", nrow(df)))

# Define populations
df <- df %>%
  mutate(
    # Disc herniation population: prolapse surgery
    pop_disc = as.integer(HovedInngrepV2V3 == 1),
    # Stenosis population: decompression (midline-preserving or laminectomy)
    pop_stenosis = as.integer(HovedInngrepV2V3 %in% c(2, 3)),
    # Exclude fusion, revision, osteotomy, disc replacement for primary comparisons
    eligible = as.integer(HovedInngrepV2V3 %in% c(1, 2, 3))
  )

# Restrict to eligible procedures (prolapse or decompression)
df <- df %>% filter(eligible == 1)
cat(sprintf("  After procedure eligibility filter: %d\n", nrow(df)))

# Population counts
cat(sprintf("  Disc herniation population: %d (ELD=%d, MSD=%d)\n",
            sum(df$pop_disc == 1),
            sum(df$pop_disc == 1 & df$treatment == "ELD"),
            sum(df$pop_disc == 1 & df$treatment == "MSD")))
cat(sprintf("  Stenosis population: %d (ELD=%d, MSD=%d)\n",
            sum(df$pop_stenosis == 1),
            sum(df$pop_stenosis == 1 & df$treatment == "ELD"),
            sum(df$pop_stenosis == 1 & df$treatment == "MSD")))

# =============================================================================
# 4. BASELINE COVARIATES (SAP Section 10)
# =============================================================================

df <- df %>%
  mutate(
    # --- Demographics ---
    age = Alder,
    sex = factor(Kjonn, levels = c(1, 2), labels = c("Male", "Female")),
    bmi = BMI,

    # Smoking: 0=No, 1=Current, 2=Former, 9=Not filled
    smoking = case_when(
      RokerV3 == 0 ~ "Never",
      RokerV3 == 1 ~ "Current",
      RokerV3 == 2 ~ "Former",
      RokerV3 == 9 ~ NA_character_
    ),
    smoking = factor(smoking, levels = c("Never", "Former", "Current")),

    # --- Socioeconomic ---
    education = case_when(
      Utd %in% c(1, 2, 3) ~ "Secondary_or_less",
      Utd %in% c(4, 5) ~ "Higher_education",
      Utd == 9 ~ NA_character_
    ),
    education = factor(education, levels = c("Secondary_or_less", "Higher_education")),

    # Employment status at baseline
    employed_baseline = case_when(
      ArbstatusPreV2V3 == 1 ~ 1L,           # Working (full/part time)
      ArbstatusPreV2V3 %in% c(3:9) ~ 0L,    # Not working
      TRUE ~ NA_integer_
    ),

    # Sick leave (preserve NA when employment status is unknown)
    sick_leave = case_when(
      is.na(ArbstatusPreV2V3) ~ NA_integer_,
      ArbstatusPreV2V3 %in% c(6, 7) ~ 1L,   # Fully or partially sick-listed
      TRUE ~ 0L
    ),

    # Disability pension applied/receiving (preserve NA when employment status is unknown)
    disability = case_when(
      is.na(ArbstatusPreV2V3) ~ NA_integer_,
      ArbstatusPreV2V3 == 9 ~ 1L,            # Disability pension
      TRUE ~ 0L
    ),

    # Analgesic use at baseline
    analgesic_baseline = case_when(
      SmStiPre == 1 ~ 1L,
      SmStiPre == 0 ~ 0L,
      SmStiPre == 9 ~ NA_integer_
    ),

    # --- Symptoms ---
    odi_baseline = OswTotPre,
    eq5d_baseline = EQ5DV3Pre,
    nrs_back_baseline = SmRyPre,
    nrs_leg_baseline = SmBePre,

    # Symptom duration (back/hip): ordinal 1-5
    # 1=None, 2=<3m, 3=3-12m, 4=12-24m, 5=>24m
    symptom_duration_back = case_when(
      SymptVarighRyggHof %in% 1:5 ~ as.integer(SymptVarighRyggHof),
      TRUE ~ NA_integer_
    ),

    # Symptom duration (radiating): ordinal 1-5
    symptom_duration_leg = case_when(
      SympVarighUtstr %in% 1:5 ~ as.integer(SympVarighUtstr),
      TRUE ~ NA_integer_
    ),

    # Motor deficit at baseline
    motor_deficit = case_when(
      OpIndParese == 1 ~ 1L,
      OpIndParese == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),

    # --- Comorbidities ---
    asa = case_when(
      ASA %in% 1:5 ~ as.integer(ASA),
      TRUE ~ NA_integer_
    ),
    asa_cat = factor(
      case_when(
        ASA == 1 ~ "ASA_1",
        ASA == 2 ~ "ASA_2",
        ASA >= 3 & ASA <= 5 ~ "ASA_3plus"
      ),
      levels = c("ASA_1", "ASA_2", "ASA_3plus")
    ),

    depression_anxiety = as.integer(SykdDepresjonAngst == 1),
    chronic_pain = as.integer(SykdGeneralisertSmSyndr == 1),

    # --- Pathology ---
    # Prolapse location categories
    prolapse_intraforaminal = as.integer(RfIntrforaminaltProl == 1),
    prolapse_extralateral = as.integer(RfEkstrLatProl == 1),

    # Stenosis types
    stenosis_central = as.integer(RfSentr == 1),
    stenosis_lateral = as.integer(RfLateral == 1),
    stenosis_foraminal = as.integer(RfForaminalSS == 1),

    # Spondylolisthesis (isthmic or degenerative)
    spondylolisthesis = as.integer(RfSpondtypeIsmisk == 1 | RfSpondtypeDegen == 1),

    # Scoliosis
    scoliosis = as.integer(RfDegskol == 1),

    # --- Surgical factors ---
    # Prior lumbar surgery
    prior_surgery = case_when(
      TidlOpr == 4 ~ "None",
      TidlOpr == 1 ~ "Same_level",
      TidlOpr == 2 ~ "Other_level",
      TidlOpr == 3 ~ "Both",
      TRUE ~ NA_character_
    ),
    prior_surgery = factor(prior_surgery, levels = c("None", "Same_level", "Other_level", "Both")),
    prior_surgery_any = as.integer(TidlOpr != 4),
    n_prior_surgeries = TidlOprAntall,

    # Surgical level(s)
    level_L45 = as.integer(OpL45 == 1),
    level_L5S1 = as.integer(OpL5S1 == 1),
    level_L34 = as.integer(OpL34 == 1),
    level_L23 = as.integer(OpL23 == 1),
    level_upper = as.integer(OpTh12L1 == 1 | OpL1L2 == 1 | OpL23 == 1),

    # Number of operated levels
    n_levels = as.integer(OpTh12L1 == 1) + as.integer(OpL1L2 == 1) +
               as.integer(OpL23 == 1) + as.integer(OpL34 == 1) +
               as.integer(OpL45 == 1) + as.integer(OpL5S1 == 1),
    multilevel = as.integer(n_levels > 1),

    # Calendar time (days since start of study period)
    surgery_date = as.Date(OpDato),
    calendar_time = as.numeric(surgery_date - as.Date("2022-01-01")),
    surgery_year = as.integer(format(surgery_date, "%Y")),

    # Surgical approach
    approach = case_when(
      OpTilgangV3 == 1 ~ "Midline",
      OpTilgangV3 == 3 ~ "Wiltse",
      OpTilgangV3 %in% c(2, 4) ~ "Other",
      TRUE ~ NA_character_
    ),
    approach = factor(approach, levels = c("Midline", "Wiltse", "Other")),

    # Admission category
    admission_elective = as.integer(OpKat == 1)
  )

# =============================================================================
# 5. OUTCOME VARIABLES
# =============================================================================

df <- df %>%
  mutate(
    # --- Tier 1: Primary outcome ---
    odi_3m = OswTot3mnd,

    # --- Tier 2: Secondary effectiveness ---
    odi_12m = OswTot12mnd,
    nrs_back_3m = SmRy3mnd,
    nrs_back_12m = SmRy12mnd,
    nrs_leg_3m = SmBe3mnd,
    nrs_leg_12m = SmBe12mnd,
    eq5d_3m = EQ5DV33mnd,
    eq5d_12m = EQ5DV312mnd,

    # Responder: >=30% improvement OR >=10 point improvement in ODI
    odi_change = OswTotPre - OswTot3mnd,
    odi_pct_change = ifelse(OswTotPre > 0, odi_change / OswTotPre * 100, NA_real_),
    responder_3m = case_when(
      is.na(OswTot3mnd) | is.na(OswTotPre) ~ NA_integer_,
      odi_pct_change >= 30 | odi_change >= 10 ~ 1L,
      TRUE ~ 0L
    ),

    # Return to work at 3 months
    # Working at 3m among those who were of working age (not retired)
    rtw_3m = case_when(
      is.na(Arbstatus3mndV2V3) ~ NA_integer_,
      Arbstatus3mndV2V3 == 1 ~ 1L,      # Full/part-time work
      Arbstatus3mndV2V3 == 4 ~ NA_integer_,  # Retired: exclude from RTW
      TRUE ~ 0L
    ),
    rtw_12m = case_when(
      is.na(Arbstatus12mndV2V3) ~ NA_integer_,
      Arbstatus12mndV2V3 == 1 ~ 1L,
      Arbstatus12mndV2V3 == 4 ~ NA_integer_,
      TRUE ~ 0L
    ),

    # Analgesic use at follow-up
    analgesic_3m = case_when(
      is.na(SmSti3mnd) ~ NA_integer_,
      SmSti3mnd == 1 ~ 1L,
      SmSti3mnd == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    analgesic_12m = case_when(
      is.na(SmSti12mnd) ~ NA_integer_,
      SmSti12mnd == 1 ~ 1L,
      SmSti12mnd == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),

    # Patient satisfaction (dichotomized: 1-2 = satisfied, 3-5 = not)
    # SAP: "Satisfied (very satisfied/satisfied) vs Not satisfied"
    satisfied_3m = case_when(
      is.na(Fornoyd3mnd) ~ NA_integer_,
      Fornoyd3mnd %in% c(1, 2) ~ 1L,
      Fornoyd3mnd %in% c(3, 4, 5) ~ 0L
    ),
    satisfied_12m = case_when(
      is.na(Fornoyd12mnd) ~ NA_integer_,
      Fornoyd12mnd %in% c(1, 2) ~ 1L,
      Fornoyd12mnd %in% c(3, 4, 5) ~ 0L
    ),

    # GPE (dichotomized: 1-2 = success, 3-7 = not)
    # SAP: "Success (much better/better) vs Not success"
    gpe_success_3m = case_when(
      is.na(Nytte3mnd) ~ NA_integer_,
      Nytte3mnd %in% c(1, 2) ~ 1L,
      Nytte3mnd %in% c(3, 4, 5, 6, 7) ~ 0L
    ),
    gpe_success_12m = case_when(
      is.na(Nytte12mnd) ~ NA_integer_,
      Nytte12mnd %in% c(1, 2) ~ 1L,
      Nytte12mnd %in% c(3, 4, 5, 6, 7) ~ 0L
    ),

    # --- Tier 3: Perioperative/Safety ---
    # Day surgery
    day_surgery = case_when(
      Dagkirurgi == 1 ~ 1L,
      Dagkirurgi == 0 ~ 0L,
      Dagkirurgi == 9 ~ NA_integer_
    ),

    # Length of stay
    los_total = Liggedogn,
    los_postop = LiggetidPostOp,

    # Peroperative complications (surgeon-reported composite)
    perop_comp_any = case_when(
      PeropKomp == 1 ~ 1L,
      PeropKomp == 0 ~ 0L,
      PeropKomp == 9 ~ NA_integer_
    ),
    perop_dural_tear = as.integer(PeropKompDura == 1),
    perop_nerve_injury = as.integer(PeropKompNerve == 1),
    perop_wrong_level = as.integer(PeropKompFeilnivSide == 1),
    perop_bleeding = as.integer(PeropKompTransfuBlodning == 1),
    perop_respiratory = as.integer(PeropKompResp == 1),
    perop_cardiovascular = as.integer(PeropKompKardio == 1),
    perop_anaphylaxis = as.integer(PeropKompAnafy == 1),
    perop_other = as.integer(PeropKompAnnet == 1),

    # Patient-reported complications at 3 months (composite)
    # Individual components are 0/1/NA (NA = no follow-up)
    pt_comp_any_3m = case_when(
      is.na(KpBlod3mnd) & is.na(KpDVT3mnd) & is.na(KpInfOverfla3mnd) &
        is.na(KpInfDyp3mnd) & is.na(KpLE3mnd) & is.na(KpLungebet3mnd) &
        is.na(KpUVI3mnd) ~ NA_integer_,
      KpBlod3mnd == 1 | KpDVT3mnd == 1 | KpInfOverfla3mnd == 1 |
        KpInfDyp3mnd == 1 | KpLE3mnd == 1 | KpLungebet3mnd == 1 |
        KpUVI3mnd == 1 ~ 1L,
      TRUE ~ 0L
    ),
    pt_comp_wound_superficial = case_when(
      is.na(KpInfOverfla3mnd) ~ NA_integer_, TRUE ~ as.integer(KpInfOverfla3mnd == 1)
    ),
    pt_comp_wound_deep = case_when(
      is.na(KpInfDyp3mnd) ~ NA_integer_, TRUE ~ as.integer(KpInfDyp3mnd == 1)
    ),
    pt_comp_uti = case_when(
      is.na(KpUVI3mnd) ~ NA_integer_, TRUE ~ as.integer(KpUVI3mnd == 1)
    ),
    pt_comp_pneumonia = case_when(
      is.na(KpLungebet3mnd) ~ NA_integer_, TRUE ~ as.integer(KpLungebet3mnd == 1)
    ),
    pt_comp_dvt = case_when(
      is.na(KpDVT3mnd) ~ NA_integer_, TRUE ~ as.integer(KpDVT3mnd == 1)
    ),
    pt_comp_pe = case_when(
      is.na(KpLE3mnd) ~ NA_integer_, TRUE ~ as.integer(KpLE3mnd == 1)
    ),
    pt_comp_bleeding = case_when(
      is.na(KpBlod3mnd) ~ NA_integer_, TRUE ~ as.integer(KpBlod3mnd == 1)
    ),

    # Reoperation during hospital stay
    reop_during_stay = case_when(
      ReopUnderOpph == 1 ~ 1L,
      ReopUnderOpph == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),

    # --- Tier 4: Descriptive ---
    operating_time = KnivtidTot,

    # Negative control: EQ-5D anxiety/depression dimension (1-5)
    eq5d_anxiety_baseline = case_when(
      EqangstV3Pre %in% 1:5 ~ as.integer(EqangstV3Pre),
      TRUE ~ NA_integer_
    ),
    eq5d_anxiety_3m = case_when(
      EqangstV33mnd %in% 1:5 ~ as.integer(EqangstV33mnd),
      TRUE ~ NA_integer_
    ),
    eq5d_anxiety_12m = case_when(
      EqangstV312mnd %in% 1:5 ~ as.integer(EqangstV312mnd),
      TRUE ~ NA_integer_
    ),

    # Follow-up status
    fu_3m_completed = as.integer(Ferdigstilt1b3mnd == 1),
    fu_12m_completed = as.integer(Ferdigstilt1b12mnd == 1),

    # Patient vital status
    deceased = as.integer(PasientDod == 1)
  )

# =============================================================================
# 6. TREATMENT VARIABLE CODING
# =============================================================================
df <- df %>%
  mutate(
    treatment = factor(treatment, levels = c("MSD", "ELD")),
    treatment_num = as.integer(treatment == "ELD")
  )

# =============================================================================
# 7. MISSING DATA ASSESSMENT (SAP Section 11)
# =============================================================================

# Define covariate list for missing assessment
covariates <- c(
  "age", "sex", "bmi", "smoking", "education", "employed_baseline",
  "sick_leave", "disability", "analgesic_baseline",
  "odi_baseline", "eq5d_baseline", "nrs_back_baseline", "nrs_leg_baseline",
  "symptom_duration_back", "symptom_duration_leg", "motor_deficit",
  "asa_cat", "depression_anxiety", "chronic_pain",
  "spondylolisthesis", "scoliosis",
  "prior_surgery", "prior_surgery_any", "n_prior_surgeries",
  "multilevel"
)

outcomes <- c(
  "odi_3m", "odi_12m",
  "nrs_back_3m", "nrs_back_12m", "nrs_leg_3m", "nrs_leg_12m",
  "eq5d_3m", "eq5d_12m",
  "responder_3m", "rtw_3m", "rtw_12m",
  "analgesic_3m", "analgesic_12m",
  "satisfied_3m", "satisfied_12m",
  "gpe_success_3m", "gpe_success_12m",
  "day_surgery", "los_total", "los_postop",
  "perop_comp_any", "pt_comp_any_3m",
  "reop_during_stay", "operating_time",
  "eq5d_anxiety_3m", "eq5d_anxiety_12m"
)

# Missing data summary
cat("\n=== Missing Data Assessment ===\n")
cat("\nCovariates:\n")
miss_cov <- df %>%
  dplyr::summarise(across(all_of(covariates), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    pct_missing = n_missing / nrow(df) * 100,
    category = case_when(
      pct_missing < 5 ~ "<5% (median/mode impute)",
      pct_missing < 50 ~ "5-50% (brms sub-model)",
      TRUE ~ ">50% (latent + informative prior)"
    )
  )
print(miss_cov, n = 30)

cat("\nOutcomes:\n")
miss_out <- df %>%
  dplyr::summarise(across(all_of(outcomes), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = n_missing / nrow(df) * 100)
print(miss_out, n = 30)

# =============================================================================
# 8. HANDLE MISSING COVARIATES (SAP Section 11.2-11.4)
# =============================================================================

# Identify missing rates for disc herniation population specifically
df_disc <- df %>% filter(pop_disc == 1)

cat("\n=== Missing Data in Disc Herniation Population ===\n")
miss_disc <- df_disc %>%
  dplyr::summarise(across(all_of(covariates), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    pct_missing = n_missing / nrow(df_disc) * 100,
    handling = case_when(
      pct_missing == 0 ~ "Complete",
      pct_missing < 5 ~ "Median/mode impute",
      pct_missing < 50 ~ "brms sub-model",
      TRUE ~ "Latent + prior"
    )
  ) %>%
  arrange(desc(pct_missing))
print(miss_disc, n = 30)

# Apply deterministic imputation for <5% missing covariates
# (Sensitivity analysis will jointly model these instead)
impute_median_mode <- function(x) {
  if (is.numeric(x)) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
  } else if (is.factor(x) || is.character(x)) {
    mode_val <- names(sort(table(x), decreasing = TRUE))[1]
    x[is.na(x)] <- mode_val
  }
  x
}

# Identify which covariates need simple imputation vs sub-models
# (Computed within each analysis population)

# =============================================================================
# 9. CREATE ANALYSIS DATASETS
# =============================================================================

# Disc herniation population (primary) — restricted to overlap period
df_disc <- df %>%
  filter(pop_disc == 1) %>%
  filter(surgery_date >= as.Date("2023-10-01")) %>%
  mutate(population = "disc_herniation")

# Stenosis population (exploratory) — restricted to overlap period
df_sten <- df %>%
  filter(pop_stenosis == 1) %>%
  filter(surgery_date >= as.Date("2023-10-01")) %>%
  mutate(population = "stenosis")

cat(sprintf("\n=== Analysis Populations ===\n"))
cat(sprintf("Disc herniation: %d (ELD=%d, MSD=%d)\n",
            nrow(df_disc),
            sum(df_disc$treatment == "ELD"),
            sum(df_disc$treatment == "MSD")))
cat(sprintf("  3m follow-up: %d (%.1f%%)\n",
            sum(df_disc$fu_3m_completed == 1, na.rm = TRUE),
            100 * mean(df_disc$fu_3m_completed == 1, na.rm = TRUE)))
cat(sprintf("  12m follow-up: %d (%.1f%%)\n",
            sum(df_disc$fu_12m_completed == 1, na.rm = TRUE),
            100 * mean(df_disc$fu_12m_completed == 1, na.rm = TRUE)))

cat(sprintf("\nStenosis: %d (ELD=%d, MSD=%d)\n",
            nrow(df_sten),
            sum(df_sten$treatment == "ELD"),
            sum(df_sten$treatment == "MSD")))
cat(sprintf("  3m follow-up: %d (%.1f%%)\n",
            sum(df_sten$fu_3m_completed == 1, na.rm = TRUE),
            100 * mean(df_sten$fu_3m_completed == 1, na.rm = TRUE)))

# =============================================================================
# 10. APPLY SIMPLE IMPUTATION FOR <5% MISSING COVARIATES
# =============================================================================

# For the disc herniation analysis dataset
# First identify which covariates have <5% missing
disc_miss_pct <- df_disc %>%
  dplyr::summarise(across(all_of(covariates), ~mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "var", values_to = "pct") %>%
  deframe()

covs_impute_simple <- names(disc_miss_pct[disc_miss_pct > 0 & disc_miss_pct < 0.05])
covs_submodel      <- names(disc_miss_pct[disc_miss_pct >= 0.05 & disc_miss_pct < 0.50])
covs_latent        <- names(disc_miss_pct[disc_miss_pct >= 0.50])
covs_complete      <- names(disc_miss_pct[disc_miss_pct == 0])

cat("\n=== Covariate Missing Data Handling (Disc Herniation) ===\n")
cat("Complete (0%):", paste(covs_complete, collapse = ", "), "\n")
cat("Simple impute (<5%):", paste(covs_impute_simple, collapse = ", "), "\n")
cat("Sub-model (5-50%):", paste(covs_submodel, collapse = ", "), "\n")
cat("Latent (>50%):", paste(covs_latent, collapse = ", "), "\n")

# Save missingness table for manuscript
miss_disc_table <- df_disc %>%
  dplyr::summarise(across(all_of(covariates), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    n_total = nrow(df_disc),
    pct_missing = round(n_missing / n_total * 100, 1),
    handling = case_when(
      pct_missing == 0 ~ "Complete",
      pct_missing < 5 ~ "Median/mode impute",
      pct_missing < 50 ~ "brms sub-model",
      TRUE ~ "Latent + informative prior"
    )
  ) %>%
  arrange(desc(pct_missing))
write.csv(miss_disc_table, file.path(paths$tables, "covariate_missingness_disc.csv"),
          row.names = FALSE)

# Apply simple imputation
df_disc_imp <- df_disc
for (v in covs_impute_simple) {
  df_disc_imp[[v]] <- impute_median_mode(df_disc_imp[[v]])
}

# Keep unimputed version for sensitivity analysis
df_disc_raw <- df_disc

# Same for stenosis
df_sten_imp <- df_sten
sten_miss_pct <- df_sten %>%
  dplyr::summarise(across(all_of(covariates), ~mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "var", values_to = "pct") %>%
  deframe()
covs_impute_simple_sten <- names(sten_miss_pct[sten_miss_pct > 0 & sten_miss_pct < 0.05])
for (v in covs_impute_simple_sten) {
  df_sten_imp[[v]] <- impute_median_mode(df_sten_imp[[v]])
}

# =============================================================================
# 11. FEASIBILITY CHECK (SAP Section 9)
# =============================================================================

cat("\n=== Feasibility Assessment (SAP Section 9) ===\n")
eld_disc_n <- sum(df_disc$treatment == "ELD")
cat(sprintf("ELD disc herniation: %d (threshold: 50, 80%% assurance: 55)\n", eld_disc_n))
if (eld_disc_n >= 55) {
  cat("  PASS: Exceeds 80% assurance threshold\n")
} else if (eld_disc_n >= 50) {
  cat("  MARGINAL: Between proceed threshold (50) and assurance threshold (55)\n")
  cat("  Analysis proceeds with ~75% assurance\n")
} else {
  cat("  BELOW THRESHOLD: Consider extending extraction period\n")
}

# Stenosis
eld_sten_n <- sum(df_sten$treatment == "ELD")
cat(sprintf("\nELD stenosis: %d (threshold for full model: 20)\n", eld_sten_n))
if (eld_sten_n >= 20) {
  cat("  Full adjusted model feasible\n")
} else {
  cat("  Only unadjusted comparison feasible\n")
}

# =============================================================================
# 11b. CREATE 12-MONTH ELIGIBLE DATASETS
# =============================================================================
# Patients eligible for 12-month follow-up are those who either:
# (a) had surgery on or before Dec 31, 2024 (sufficient time to reach 12 months), OR
# (b) have confirmed 12-month questionnaire completion (Ferdigstilt1b12mnd == 1),
#     indicating the registry received their 12m response regardless of surgery date.
# Patients operated later without confirmed completion have structurally absent
# 12m data (not missing at random) and are excluded from 12m analyses.

df_disc_12m_eligible <- df_disc_imp %>%
  filter(surgery_date <= as.Date("2024-12-31") |
         haven::zap_labels(Ferdigstilt1b12mnd) == 1)

cat(sprintf("\n=== 12-Month Eligible Population (Disc Herniation) ===\n"))
cat(sprintf("  Eligible: %d (ELD=%d, MSD=%d)\n",
            nrow(df_disc_12m_eligible),
            sum(df_disc_12m_eligible$treatment == "ELD"),
            sum(df_disc_12m_eligible$treatment == "MSD")))
cat(sprintf("  12m questionnaire completed: %d (%.1f%%)\n",
            sum(df_disc_12m_eligible$fu_12m_completed == 1, na.rm = TRUE),
            100 * mean(df_disc_12m_eligible$fu_12m_completed == 1, na.rm = TRUE)))
cat(sprintf("  ODI 12m observed: %d (%.1f%%)\n",
            sum(!is.na(df_disc_12m_eligible$odi_12m)),
            100 * mean(!is.na(df_disc_12m_eligible$odi_12m))))

# =============================================================================
# 12. SAVE ANALYSIS DATASETS
# =============================================================================

# Save as RDS for R analysis
saveRDS(df, file.path(paths$data_clean, "df_all.rds"))
saveRDS(df_disc, file.path(paths$data_clean, "df_disc.rds"))
saveRDS(df_disc_imp, file.path(paths$data_clean, "df_disc_imp.rds"))
saveRDS(df_disc_raw, file.path(paths$data_clean, "df_disc_raw.rds"))
saveRDS(df_disc_12m_eligible, file.path(paths$data_clean, "df_disc_12m_eligible.rds"))
saveRDS(df_sten, file.path(paths$data_clean, "df_sten.rds"))
saveRDS(df_sten_imp, file.path(paths$data_clean, "df_sten_imp.rds"))

# Save variable metadata
var_meta <- list(
  covariates = covariates,
  outcomes = outcomes,
  covs_complete = covs_complete,
  covs_impute_simple = covs_impute_simple,
  covs_submodel = covs_submodel,
  covs_latent = covs_latent,
  ni_margins = ni_margins,
  treatment_levels = levels(df$treatment)
)
saveRDS(var_meta, file.path(paths$data_clean, "var_meta.rds"))

cat("\n=== Data preparation complete ===\n")
cat(sprintf("Files saved to: %s\n", paths$data_clean))
