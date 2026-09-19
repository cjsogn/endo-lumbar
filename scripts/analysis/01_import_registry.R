source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(tibble)})
if(!nzchar(RAW_DATA) || !file.exists(RAW_DATA)) stop("Set ENDO_LUMBAR_RAW_DATA to the private SPSS registry export.")
raw <- read_sav(RAW_DATA,encoding="latin1")
stopifnot(!anyNA(raw$ForlopsID),!anyDuplicated(raw$ForlopsID))
df <- raw %>%
  mutate(
    treatment = case_when(
      OpMikroV3 == 3 ~ "ELD",
      OpMikroV3 == 1 ~ "MSD",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(treatment))

cat(sprintf("  After treatment group filter: %d procedures (ELD=%d, MSD=%d)\n",
            nrow(df), sum(df$treatment == "ELD"), sum(df$treatment == "MSD")))

# =============================================================================
# 3. DEFINE ANALYSIS POPULATIONS
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
    # Exclude fusion, osteotomy and disc replacement procedure categories
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
# 4. BASELINE COVARIATES
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
    # 1 = working at follow-up, NA = retired at follow-up (excluded),
    # 0 = all other employment statuses at follow-up (sick leave, disability,
    # unemployed, student, homemaker). Baseline employment status is included
    # as a covariate in all models.
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
impute_median_mode <- function(x) {
  if (is.numeric(x)) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
  } else if (is.factor(x) || is.character(x)) {
    mode_val <- names(sort(table(x), decreasing = TRUE))[1]
    x[is.na(x)] <- mode_val
  }
  x
}


# Define the overlap cohort before imputing baseline covariates.
df_disc <- df %>% filter(pop_disc==1,surgery_date>=as.Date("2023-10-01")) %>%
 mutate(population="disc_herniation")
rates <- vapply(df_disc[covariates],function(x)mean(is.na(x)),numeric(1))
# This release implements the low-missingness baseline rule used in the study.
# Stop on data needing a different missing-covariate model.
if(any(rates>=.05)) stop("Baseline missingness exceeds the implemented median/mode rule. Review the data and analysis specification.")
df_disc_imp <- df_disc
for(v in names(rates)[rates>0 & rates<.05]) df_disc_imp[[v]] <- impute_median_mode(df_disc_imp[[v]])
saveRDS(df,file.path(ROOT,"02_data/source/df_all.rds"))
saveRDS(df_disc,file.path(ROOT,"02_data/source/df_disc_raw.rds"))
saveRDS(df_disc_imp,file.path(ROOT,"02_data/source/df_disc_imp.rds"))
write_csv(data.frame(variable=names(rates),missing_fraction=rates),"04_results/baseline_missingness.csv")
