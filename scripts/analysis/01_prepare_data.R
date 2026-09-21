# Purpose: Prepare standardised cohorts, calendar terms and endpoint denominators.
# Inputs: Raw and baseline-imputed disc cohorts from 01_import_registry.R.
# Outputs: Derived cohorts, calendar specification, data dictionary and summaries.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
d <- readRDS(file.path(ROOT, "02_data/source/df_disc_imp.rds"))
raw <- readRDS(file.path(ROOT, "02_data/source/df_disc_raw.rds"))
stopifnot(identical(d$ForlopsID, raw$ForlopsID), !anyDuplicated(d$ForlopsID))
d <- standardize(d)
for (tp in c("3m", "12m")) {
  field <- if (tp == "3m") "EqangstV33mnd" else "EqangstV312mnd"
  z <- haven::zap_labels(d[[field]])
  d[[paste0("eq5d_anxiety_", tp)]] <- ifelse(z %in% 1:5, as.integer(z), NA_integer_)
}
stopifnot(all(complete.cases(d[COVS])))
d$calendar_time <- as.numeric(d$surgery_date - as.Date("2022-01-01"))
ct <- scale(d$calendar_time)
d$calendar_time_z <- as.numeric(ct)
b <- splines::ns(d$calendar_time_z, df = 2)
d$ct1 <- b[, 1]
d$ct2 <- b[, 2]
cal <- list(
  origin = as.Date("2022-01-01"), center = attr(ct, "scaled:center"),
  scale = attr(ct, "scaled:scale"), basis = b, attributes = attributes(b),
  convention = "ns(df=2): median interior knot, boundary knots at observed extremes"
)
cal$knot_dates <- cal$origin + cal$center +
  c(attr(b, "Boundary.knots")[1], attr(b, "knots"), attr(b, "Boundary.knots")[2]) * cal$scale
saveRDS(cal, file.path(ROOT, "02_data/derived/calendar_basis_specification.rds"))
write_csv(data.frame(
  type = c("lower boundary", "interior", "upper boundary"),
  date = as.character(cal$knot_dates)
), "04_results/calendar_knots.csv")
d$odi_3m_zib <- pmin(d$odi_3m / 100, 1 - 1e-6)
d12 <- standardize(d[d$surgery_date <= as.Date("2025-02-28"), ])
stopifnot(nrow(d12) > 0L, all(table(d12$treatment) > 0L))
saveRDS(d, file.path(ROOT, "02_data/derived/cohort_3m.rds"))
saveRDS(d12, file.path(ROOT, "02_data/derived/cohort_12m_date_eligible.rds"))
scaling <- do.call(rbind, lapply(list(full = d, eligible_12m = d12), function(x)
  data.frame(variable = CONT, center = sapply(x[CONT], mean), scale = sapply(x[CONT], sd))))
write_csv(
  data.frame(cohort = rep(c("full", "eligible_12m"), each = length(CONT)), scaling),
  "04_results/baseline_scaling.csv"
)
outcomes <- c(
  "odi_3m", "odi_12m", "nrs_back_3m", "nrs_back_12m", "nrs_leg_3m", "nrs_leg_12m",
  "eq5d_3m", "eq5d_12m", "responder_3m", "rtw_3m", "rtw_12m", "analgesic_3m", "analgesic_12m",
  "satisfied_3m", "satisfied_12m", "gpe_success_3m", "gpe_success_12m", "day_surgery",
  "los_postop", "perop_comp_any", "pt_comp_any_3m", "operating_time"
)
counts <- do.call(rbind, lapply(outcomes, function(y) {
  x <- if (grepl("12m$", y)) d12 else d
  do.call(rbind, lapply(levels(x$treatment), function(a) {
    z <- x[[y]][x$treatment == a]
    data.frame(outcome = y, arm = a, eligible = length(z), observed = sum(!is.na(z)), missing = sum(is.na(z)))
  }))
}))
write_csv(counts, "04_results/endpoint_denominators.csv")
codes <- haven::zap_labels(d$OpTilgangV3[d$treatment == "ELD"])
write_csv(as.data.frame(table(
  code = factor(codes, levels = sort(unique(c(codes, 0:4, 9)))),
  useNA = "always"
)), "04_results/eld_approach_codes.csv")
events <- do.call(rbind, lapply(c("perop_comp_any", "pt_comp_any_3m"), function(y)
  do.call(rbind, lapply(levels(d$treatment), function(a) {
    z <- d[[y]][d$treatment == a]
    n <- sum(!is.na(z))
    k <- sum(z == 1, na.rm = TRUE)
    ci <- binom.test(k, n)$conf.int
    data.frame(outcome = y, arm = a, events = k, n = n, proportion = k / n, exact95_lower = ci[1], exact95_upper = ci[2])
  }))))
write_csv(events, "04_results/complication_counts_exact_intervals.csv")
n <- sum(!is.na(d$pt_comp_any_3m))
k <- sum(d$pt_comp_any_3m == 1, na.rm = TRUE)
write_csv(
  data.frame(
    events = k, n = n, overall_rate = k / n, rule = "Overall observed cohort event rate",
    classification = if (k / n >= .10) "formal" else if (k / n >= .05) "test with caution" else "descriptive"
  ),
  "04_results/patient_complication_rate_rule.csv"
)
write_csv(
  data.frame(
    arm = levels(raw$treatment),
    baseline_odi_missing = sapply(levels(raw$treatment), function(a) sum(is.na(raw$odi_baseline[raw$treatment == a])))
  ),
  "04_results/eligibility_baseline_odi_discrepancy.csv"
)
dict <- data.frame(
  variable = names(raw),
  class = vapply(raw, function(x) paste(class(x), collapse = "/"), character(1)),
  # Exact matching prevents an absent label from matching value-level labels.
  label = vapply(raw, function(x) {
    z <- attr(x, "label", exact = TRUE)
    if (is.null(z)) return("")
    stopifnot(length(z) == 1L)
    as.character(z)
  }, character(1)),
  missing = vapply(raw, function(x) sum(is.na(x)), integer(1)), row.names = NULL
)
stopifnot(
  identical(names(dict), c("variable", "class", "label", "missing")),
  nrow(dict) == ncol(raw), !anyDuplicated(dict$variable)
)
write_csv(dict, "02_data/derived/data_dictionary.csv")
writeLines(
  c(capture.output(sessionInfo()), "No patient-level data are printed in logs."),
  file.path(ROOT, "10_logs/data_preparation_session.txt")
)
print(counts[counts$outcome %in% c("odi_3m", "odi_12m"), ])
print(events)
