# Purpose: Assess primary ODI standardisation after a follow-up opportunity restriction.
# Inputs: Existing primary fit, prepared cohort, explicit reference date and date basis.
# Outputs: Restricted and full-population posterior summaries, counts and verification.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
reference_text <- Sys.getenv("ENDO_FOLLOWUP_REFERENCE")
reference_basis <- Sys.getenv("ENDO_FOLLOWUP_REFERENCE_BASIS")
if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", reference_text) ||
    is.na(as.Date(reference_text))) stop("Supply ENDO_FOLLOWUP_REFERENCE as YYYY-MM-DD.")
if (!reference_basis %in% c("SPSS file creation date", "confirmed registry cutoff"))
  stop("Specify whether the reference is an SPSS file creation date or confirmed registry cutoff.")
reference <- as.Date(reference_text)
# Subtract three calendar months, clamping the day at the target month's end.
month_start <- as.Date(format(reference, "%Y-%m-01"))
cutoff_month <- seq(month_start, by = "-3 months", length.out = 2L)[2]
next_month <- seq(cutoff_month, by = "month", length.out = 2L)[2]
day <- min(as.integer(format(reference, "%d")), as.numeric(next_month - cutoff_month))
operation_cutoff <- cutoff_month + day - 1L

d <- readRDS(file.path(ROOT, "02_data/derived/cohort_3m.rds"))
fit_path <- file.path(ROOT, "03_models/odi_3m_calendar.rds")
fit <- readRDS(fit_path)
eligible <- d$surgery_date <= operation_cutoff
stopifnot(!anyNA(eligible), any(eligible), all(table(d$treatment[eligible]) > 0L))
if (any(!is.na(d$odi_3m[!eligible])))
  stop("The restriction removes observed ODI values. This procedure cannot reuse the original fit; an appropriate refit is required.")

# The likelihood, baseline scaling, imputations, spline basis and priors stay fixed.
# Reusing draws is justified only when every observed fitting row is retained.
observed <- d[!is.na(d$odi_3m), ]
stopifnot(nrow(fit$data) == nrow(observed),
          identical(observed$ForlopsID, d$ForlopsID[eligible & !is.na(d$odi_3m)]),
          identical(fit$family$family, "zero_inflated_beta"))
# brms retains model variables, not the registry identifier. Verify every model
# input row by row, including the response, instead of relying on a missing ID.
if ("ForlopsID" %in% names(fit$data))
  stopifnot(identical(fit$data$ForlopsID, observed$ForlopsID))
variables <- unique(c("treatment", COVS, "ct1", "ct2", "odi_3m_zib"))
checks <- data.frame(variable = variables, matches_retained_fit = vapply(variables, function(v)
  isTRUE(all.equal(observed[[v]], fit$data[[v]], check.attributes = FALSE)), logical(1)))
stopifnot(all(checks$matches_retained_fit),
          all(COVS %in% all.vars(fit$formula$formula)),
          all(c("treatmentELD:ct1", "treatmentELD:ct2") %in% rownames(fixef(fit))))
write_csv(checks, "08_qa/followup_reference_model_inputs.csv")

predictions <- lapply(c("MSD", "ELD"), function(arm) {
  nd <- d
  nd$treatment <- factor(arm, levels = c("MSD", "ELD"))
  out <- posterior_epred(fit, newdata = nd, allow_new_levels = FALSE)
  stopifnot(length(dim(out)) == 2L, ncol(out) == nrow(d), all(is.finite(out)))
  out
})
individual_contrasts <- 100 * (predictions[[1]] - predictions[[2]])
draws <- list(
  full = rowMeans(individual_contrasts),
  restricted = rowMeans(individual_contrasts[, eligible, drop = FALSE])
)
# A separate restricted newdata call verifies the column-subset calculation.
independent <- gcomp(fit, d[eligible, ], lower_better = TRUE, scale_factor = 100)$delta
max_difference <- max(abs(draws$restricted - independent))
stopifnot(max_difference < 1e-10)

summaries <- do.call(rbind, lapply(names(draws), function(population)
  do.call(rbind, lapply(c(7, 5, 3), function(margin)
    cbind(population = population, n = if (population == "full") nrow(d) else sum(eligible),
          summarize_effect(draws[[population]], margin))))))
counts <- do.call(rbind, lapply(levels(d$treatment), function(arm) {
  i <- d$treatment == arm
  late <- i & !eligible
  data.frame(arm = arm, full = sum(i), retained = sum(i & eligible), excluded = sum(late),
    observed_full = sum(i & !is.na(d$odi_3m)),
    observed_retained = sum(i & eligible & !is.na(d$odi_3m)),
    missing_full = sum(i & is.na(d$odi_3m)),
    missing_retained = sum(i & eligible & is.na(d$odi_3m)),
    excluded_min_days = if (any(late)) min(as.numeric(reference - d$surgery_date[late])) else NA_real_,
    excluded_max_days = if (any(late)) max(as.numeric(reference - d$surgery_date[late])) else NA_real_)
}))
provenance <- data.frame(reference_date = reference_text, reference_basis = reference_basis,
  operation_cutoff = as.character(operation_cutoff), fitted_rows = nrow(fit$data),
  posterior_draws = length(draws$full), full_target_n = nrow(d), restricted_target_n = sum(eligible),
  model_md5 = unname(tools::md5sum(fit_path)), observed_fitting_rows_unchanged = TRUE,
  scaling_and_spline_unchanged = TRUE, separate_calculation_max_difference = max_difference)
write_csv(summaries, "04_results/odi_followup_reference_sensitivity.csv")
write_csv(counts, "04_results/followup_reference_counts.csv")
write_csv(provenance, "08_qa/followup_reference_provenance.csv")
saveRDS(draws, file.path(ROOT, "04_results/odi_followup_reference_draws.rds"))
depth <- if (inherits(fit$fit, "stanfit")) {
  fit$fit@stan_args[[1]]$control$max_treedepth
} else if (inherits(fit$fit, "CmdStanMCMC")) {
  fit$fit$metadata()$max_treedepth
} else stop("Cannot verify the retained fit's maximum tree-depth setting.")
stopifnot(length(depth) == 1L, is.finite(depth), depth > 0L)
dg <- diagnostics(fit, "odi_followup_reference_retained_fit", depth)
stopifnot(dg$pass)
write_csv(dg, "08_qa/followup_reference_retained_diagnostics.csv")
capture.output(list(formula = fit$formula, priors = prior_summary(fit),
  reference_date = reference_text, reference_basis = reference_basis,
  sampling = "All retained draws reused. No new sampling, imputation or refit."),
  file = file.path(ROOT, "08_qa/followup_reference_specification.txt"))
print(counts)
print(summaries)
cat("PASS: observed fitting data unchanged; restricted standardisation independently reproduced.\n")
