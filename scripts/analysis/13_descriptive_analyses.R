# Purpose: Summarise baseline characteristics, case sequence and recorded reoperations.
# Inputs: Coded source data, prepared cohorts and the authorised SPSS export.
# Outputs: Descriptive summaries and a private procedure-level reoperation linkage.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(tableone)
})
df_disc <- readRDS(file.path(ROOT, "02_data/source/df_disc_raw.rds"))
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
  "prolapse_intraforaminal", "prolapse_extralateral", "stenosis_central",
  "spondylolisthesis", "scoliosis",
  # Surgical
  "prior_surgery", "prior_surgery_any", "n_prior_surgeries",
  "n_levels", "multilevel",
  "level_L45", "level_L5S1", "level_L34",
  "approach",
  "surgery_year"
)

# Identify factor vs numeric variables
cat_vars <- c(
  "sex", "smoking", "education", "asa_cat",
  "prior_surgery", "approach", "surgery_year"
)
nonnormal_vars <- c(
  "n_prior_surgeries", "symptom_duration_back", "symptom_duration_leg",
  "n_levels", "los_total", "operating_time"
)

# Ensure surgery_year is factor for table
df_disc <- df_disc %>% mutate(surgery_year = factor(surgery_year))


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
  nonnormal = nonnormal_vars
)

cat("\n--- Disc Herniation: Baseline Characteristics ---\n")
print(tab1_disc_print, quote = FALSE)

# Save Table 1 as CSV
write.csv(
  tab1_disc_print,
  file.path(ROOT, "04_results/table1_disc_herniation.csv")
)


wr <- function(x, f) write_csv(x, file.path("04_results", f))
d <- readRDS(file.path(ROOT, '02_data/derived/cohort_3m.rds'))
raw <- readRDS(file.path(ROOT, '02_data/source/df_all.rds'))
lc <- raw %>%
  filter(ForlopsID %in% d$ForlopsID, treatment == 'ELD') %>%
  arrange(surgery_date) %>%
  mutate(
    case_no = row_number(), operating_time = ifelse(operating_time >= 500, NA_real_, operating_time),
    pure_disc = !(stenosis_central == 1 | stenosis_lateral == 1 | stenosis_foraminal == 1)
  )
stopifnot(nrow(lc) == sum(d$treatment == "ELD"))
pure <- lc %>%
  filter(pure_disc) %>%
  mutate(quartile = ntile(case_no, 4))
q <- pure %>%
  group_by(quartile) %>%
  summarise(
    n = sum(!is.na(operating_time)), case_no = mean(case_no),
    median = median(operating_time, na.rm = TRUE), mean = mean(operating_time, na.rm = TRUE), .groups = 'drop'
  )
wr(q, 'learning_quartiles_verified.csv')
rho <- function(x, y, name) {
  v <- suppressWarnings(cor.test(x$case_no, x[[y]], method = 'spearman', exact = FALSE))

  data.frame(analysis = name, n = sum(!is.na(x[[y]])), rho = unname(v$estimate), p = v$p.value)
}
wr(rbind(
  rho(pure, 'operating_time', 'Pure-disc operating time'), rho(lc, 'operating_time', 'All ELD operating time'),
  rho(pure, 'odi_3m', 'Pure-disc ODI')
), 'learning_descriptive_verified.csv')
eld <- d[d$treatment == 'ELD', ]
eld$corridor <- ifelse(as.numeric(eld$OpTilgangV3) == 1, 'Inferred IL', 'Inferred TF')
stopifnot(all(as.numeric(eld$OpTilgangV3) %in% c(1, 3)))
eld$operating_time[eld$operating_time >= 500] <- NA_real_
cs <- lapply(c('odi_3m', 'nrs_leg_3m', 'operating_time'), function(y) {
  sub <- eld[!is.na(eld[[y]]), ]
  sub$corridor <- factor(sub$corridor, levels = c('Inferred IL', 'Inferred TF'))

  sub %>%
    group_by(corridor) %>%
    summarise(outcome = y, n = n(), mean = mean(.data[[y]]), median = median(.data[[y]]), .groups = 'drop')
})
wr(bind_rows(cs), 'corridor_descriptive_verified.csv')

d12 <- readRDS(file.path(ROOT, "02_data/derived/cohort_12m_date_eligible.rds"))
raw <- haven::read_sav(RAW_DATA, encoding = "latin1")
# Same-patient later surgery is linked across the full export. The same-level
# discectomy flag is a recurrence proxy, not adjudicated recurrent/residual disc.
link <- lapply(seq_len(nrow(d)), function(i) {
  subs <- raw[!is.na(raw$PasientID) & raw$PasientID == d$PasientID[i] &
    !is.na(raw$OpDato) & raw$OpDato > d$surgery_date[i] &
    raw$ForlopsID != d$ForlopsID[i], ]
  days <- as.numeric(as.Date(subs$OpDato) - d$surgery_date[i])
  proxy <- subs$TidlOpsammeNiv == 1 & subs$HovedInngrepV2V3 == 1
  earliest <- function(z) if (length(z)) min(z) else NA_real_
  data.frame(
    ForlopsID = d$ForlopsID[i], first_reop_any_days = earliest(days),
    first_rehern_days = earliest(days[which(proxy)]),
    first_other_reop_days = earliest(days[which(!proxy)]), n_reops = nrow(subs)
  )
})
link <- do.call(rbind, link)

link$treatment <- d$treatment
link$PasientID <- d$PasientID
saveRDS(link, file.path(ROOT, "02_data/derived/reoperation_linkage.rds"))
ev <- link[link$ForlopsID %in% d12$ForlopsID & !is.na(link$first_reop_any_days) &
  link$first_reop_any_days > 0 & link$first_reop_any_days <= 365, ]
re <- bind_rows(lapply(c('Recorded reoperation', 'Reherniation proxy', 'Without reherniation flag'), function(y) {
  bind_rows(lapply(c('ELD', 'MSD'), function(a) {
    s <- ev[ev$treatment == a, ]
    n <- sum(d12$treatment == a)
    k <- if (y == 'Recorded reoperation') nrow(s) else if (y == 'Reherniation proxy') sum(!is.na(s$first_rehern_days) & s$first_rehern_days > 0 & s$first_rehern_days <= 365) else sum(is.na(s$first_rehern_days))
    b <- binom.test(k, n)
    data.frame(outcome = y, treatment = a, k = k, n = n, mean = 100 * k / n, lower = 100 * b$conf.int[1], upper = 100 * b$conf.int[2])
  }))
}))
wr(re, 'reoperation_date_eligible_verified.csv')
