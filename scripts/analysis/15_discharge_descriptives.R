# Purpose: Reconcile planned admission type with actual postoperative stay.
# Inputs: Prepared cohort and the authorised SPSS export.
# Outputs: Aggregate cross-tabulations and historical MSD discharge summaries.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
suppressPackageStartupMessages(library(dplyr))
d <- readRDS(file.path(ROOT, "02_data/derived/cohort_3m.rds"))
stopifnot(!anyNA(d$day_surgery), !anyNA(d$los_postop))
cross <- d %>% mutate(
  planned_day_case_code = factor(day_surgery_registered, levels = 0:1, labels = c("No", "Yes")),
  actual_postoperative_nights = factor(pmin(los_postop, 3), levels = 0:3,
                                      labels = c("0", "1", "2", "3 or more"))
)
write_csv(as.data.frame(table(
  arm = cross$treatment, planned_day_case = cross$planned_day_case_code,
  postoperative_nights = cross$actual_postoperative_nights, useNA = "ifany"
)), "04_results/planned_day_case_by_actual_stay.csv")

summary <- bind_rows(lapply(levels(d$treatment), function(a) {
  x <- d[d$treatment == a, ]
  data.frame(
    arm = a, n = nrow(x), same_day_discharge = sum(x$day_surgery),
    same_day_percent = 100 * mean(x$day_surgery),
    date_stay_median = median(x$los_postop),
    date_stay_q1 = unname(quantile(x$los_postop, .25)),
    date_stay_q3 = unname(quantile(x$los_postop, .75)),
    reported_stay_missing = sum(is.na(x$los_postop_registered)),
    planned_yes_overnight = sum(x$day_surgery_registered == 1 & x$los_postop > 0, na.rm = TRUE),
    planned_no_same_day = sum(x$day_surgery_registered == 0 & x$los_postop == 0, na.rm = TRUE)
  )
}))
write_csv(summary, "04_results/discharge_summary.csv")

# Earlier procedures supply descriptive context only and do not enter adjusted models.
raw <- haven::read_sav(RAW_DATA, encoding = "latin1")
h <- raw[which(raw$OpMikroV3 %in% c(1, 3) & raw$HovedInngrepV2V3 == 1), ]
pre <- h[which(as.Date(h$OpDato) < as.Date("2023-10-01")), ]
write_csv(data.frame(
  preintroduction_before_cauda_exclusion = nrow(pre),
  preintroduction_cauda_exclusions = sum(pre$OpIndCauda == 1, na.rm = TRUE)
), "04_results/historical_context_eligibility.csv")
h <- h[which(h$OpIndCauda != 1 | is.na(h$OpIndCauda)), ]
operation <- as.Date(h$OpDato)
discharge <- as.Date(h$UtskrivelseDato, format = "%Y-%m-%d")
stopifnot(all(is.na(discharge) == is.na(h$UtskrivelseDato)))
nights <- as.numeric(discharge - operation)
stopifnot(all(nights[!is.na(nights)] >= 0))
h$same_day_discharge <- ifelse(is.na(nights), NA_integer_, as.integer(nights == 0))
h$arm <- ifelse(h$OpMikroV3 == 1, "MSD", "ELD")
h$period <- ifelse(operation < as.Date("2023-10-01"), "Before October 2023", format(operation, "%Y"))
h$day_flag <- ifelse(h$Dagkirurgi %in% 0:1, as.numeric(h$Dagkirurgi), NA_real_)
trend <- h %>% filter(arm == "MSD") %>% group_by(period) %>% summarise(
  n = n(), flag_observed = sum(!is.na(day_flag)), planned_yes = sum(day_flag == 1, na.rm = TRUE),
  dates_observed = sum(!is.na(same_day_discharge)), actual_same_day = sum(same_day_discharge, na.rm = TRUE),
  actual_same_day_percent = 100 * mean(same_day_discharge, na.rm = TRUE), .groups = "drop"
)
write_csv(trend, "04_results/msd_discharge_trend.csv")
