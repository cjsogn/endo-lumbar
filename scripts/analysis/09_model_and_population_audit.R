# Purpose: Check standardisation populations, TMLE targeting and registry fields.
# Inputs: Prepared/source cohorts, Bayesian models and observed-case TMLE fits.
# Outputs: Combined estimates, population checks, overlap and field inventories.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
d <- readRDS(file.path(ROOT, "02_data/derived/cohort_3m.rds"))
d12 <- readRDS(file.path(ROOT, "02_data/derived/cohort_12m_date_eligible.rds"))
new <- read.csv(file.path(ROOT, "04_results/verified_calendar_primary_perioperative.csv"))
for (fn in list.files(file.path(ROOT, "04_results"), "_calendar_summary.csv$", full.names = TRUE)) {
  if (grepl("tipping_point|tmle_", basename(fn))) next
  z <- read.csv(fn)
  new <- rbind(new, z[names(new)])
}
stopifnot(!anyDuplicated(new$id))
write_csv(new, "04_results/calendar_all_bayesian_results.csv")
# Match the Bayesian standardisation population to each observed-case TMLE.
matched <- list()
targeting <- list()
for (y in c("odi_3m", "odi_12m", "day_surgery")) {
  z <- if (grepl("12m", y)) d12 else d
  fit <- readRDS(file.path(ROOT, "03_models", paste0(y, "_calendar.rds")))
  if (y == "odi_12m") {
    z$y <- pmin(z$odi_12m / 100, 1 - 1e-6)
    z$baseline_zi <- z$odi_baseline_z
  }
  obs <- !is.na(z[[y]])
  dr <- gcomp(fit, z[obs, ], y != "day_surgery", if (y == "day_surgery") 1 else 100)
  matched[[y]] <- cbind(id = y, target_n = sum(obs), summarize_effect(dr$delta, if (y == "day_surgery") NA else 7))
  f <- readRDS(file.path(ROOT, "03_models", paste0("tmle_", y, "_calendar.rds")))$fit
  sign <- if (y == "day_surgery") 1 else -1
  targeting[[y]] <- data.frame(
    id = y, initial_prediction_difference = sign * mean(f$Qinit$Q[, 2] - f$Qinit$Q[, 1]),
    targeted_difference = sign * mean(f$Qstar[, 2] - f$Qstar[, 1]),
    min_ps = min(f$g$g1W), max_ps = max(f$g$g1W), ps_lower_bound = f$gbound[1],
    n_below_lower_bound = sum(f$g$g1W < f$gbound[1])
  )
}
write_csv(do.call(rbind, matched), "04_results/bayesian_matched_observed_populations.csv")
write_csv(do.call(rbind, targeting), "04_results/tmle_targeting_step_diagnostics.csv")

# Baseline composition by actual 12-month response, without inferring attrition
# mechanisms from response proportions alone.
comp12 <- do.call(rbind, lapply(c("MSD", "ELD"), function(a)
  do.call(rbind, lapply(c(TRUE, FALSE), function(observed) {
    x <- d12[d12$treatment == a & (!is.na(d12$odi_12m)) == observed, ]
    data.frame(
      arm = a, odi_observed = observed, n = nrow(x), variable = CONT,
      mean = sapply(x[CONT], mean), sd = sapply(x[CONT], sd)
    )
  }))))
write_csv(comp12, "04_results/12m_baseline_by_response_and_arm.csv")
psfit <- glm(as.formula(paste("I(treatment == 'ELD') ~", paste(c(COVS, "ct1", "ct2"), collapse = " + "))),
  data = d, family = binomial()
)
saveRDS(psfit, file.path(ROOT, "03_models/calendar_propensity_diagnostic.rds"))
ps <- predict(psfit, type = "response")
write_csv(
  data.frame(treatment = d$treatment, propensity = ps, year = format(d$surgery_date, "%Y")),
  "04_results/calendar_propensity_values.csv"
)
support <- do.call(rbind, lapply(split(seq_len(nrow(d)), d$treatment), function(i)
  data.frame(
    n = length(i), min_ps = min(ps[i]), q05 = quantile(ps[i], .05), median = median(ps[i]),
    max_ps = max(ps[i]), below_005 = sum(ps[i] < .05), above_095 = sum(ps[i] > .95)
  )))
write_csv(data.frame(arm = rownames(support), support), "04_results/calendar_overlap_summary.csv")
write_csv(
  as.data.frame(table(arm = d$treatment, year = format(d$surgery_date, "%Y"))),
  "04_results/calendar_support_by_year.csv"
)

# Field availability is extracted from this study dataset, not inferred from a
# later public version of the registry form.
raw <- readRDS(file.path(ROOT, "02_data/source/df_disc_raw.rds"))
fields <- grep("^PeropKomp|^Kp.*3mnd$|ReopUnderOpph|^OpTilgangV3$|^OpMikroV3$", names(raw), value = TRUE)
inventory <- do.call(rbind, lapply(fields, function(v) {
  x <- raw[[v]]
  lab <- attr(x, "label", exact = TRUE)
  lev <- attr(x, "labels")
  data.frame(
    field = v, label = if (is.null(lab)) "" else lab,
    codes = if (is.null(lev)) "" else paste(names(lev), lev, sep = "=", collapse = " | "), missing = sum(is.na(x))
  )
}))
write_csv(inventory, "04_results/registry_complication_field_inventory.csv")
components <- do.call(rbind, lapply(fields, function(v) do.call(rbind, lapply(c("MSD", "ELD"), function(a) {
  x <- haven::zap_labels(raw[[v]][raw$treatment == a])
  data.frame(
    field = v, arm = a,
    n_nonmissing = sum(!is.na(x)), code1 = sum(x == 1, na.rm = TRUE), code0 = sum(x == 0, na.rm = TRUE)
  )
}))))
write_csv(components, "04_results/registry_complication_components_raw_codes.csv")

# Read every completed main-model diagnostic. Summaries do not replace PPC
# inspection or justify family selection by incomparable likelihood scales.
files <- list.files(file.path(ROOT, "08_qa"), "_calendar_diagnostics.csv$", full.names = TRUE)
dg <- do.call(rbind, lapply(files, read.csv))
stopifnot(all(dg$pass))
write_csv(dg, "08_qa/all_calendar_model_diagnostics.csv")
print(do.call(rbind, matched))
print(do.call(rbind, targeting))
