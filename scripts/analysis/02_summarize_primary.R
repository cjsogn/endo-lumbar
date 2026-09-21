# Purpose: Standardise primary/perioperative effects and evaluate ODI margins.
# Inputs: Prepared cohort, model definitions and fitted primary/perioperative models.
# Outputs: Posterior draws, effect summaries, margin checks and provenance records.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "02_primary_jobs.R"))
res <- diag <- audit <- list()
for (id in primary_ids) {
  b <- build_primary(id)
  fit <- readRDS(file.path(ROOT, "03_models", paste0(id, "_calendar.rds")))
  actual_depth <- fit$fit@stan_args[[1]]$control$max_treedepth
  diag[[id]] <- diagnostics(fit, id, actual_depth)
  stopifnot(diag[[id]]$pass)
  capture.output(list(
    formula = fit$formula, prior = prior_summary(fit),
    sampler = fit$fit@stan_args
  ), file = file.path(ROOT, "08_qa", paste0("verified_", id, "_specification.txt")))
  interaction_present <- all(c("treatmentELD:ct1", "treatmentELD:ct2") %in% rownames(fixef(fit)))
  stopifnot(interaction_present, all(COVS %in% all.vars(fit$formula$formula)))
  path <- file.path(ROOT, "03_models", paste0(id, "_calendar.rds"))
  response <- all.vars(fit$formula$formula)[1]
  audit[[id]] <- data.frame(
    id = id, source = path, fit_rows = nrow(fit$data),
    observed = sum(!is.na(fit$data[[response]])),
    target_n = if (id == "los_postop") nrow(b$data) else nrow(b$target),
    full_covariates = TRUE, calendar_interaction = interaction_present,
    source_md5 = unname(tools::md5sum(path))
  )
  if (id == "los_postop") draws <- ordinal_calendar_contrast(fit, b$data) else
    draws <- gcomp(fit, b$target, id != "day_surgery", if (id == "odi_3m") 100 else 1)
  res[[id]] <- cbind(id = id, summarize_effect(draws$delta, b$margin))
  saveRDS(draws, file.path(ROOT, "04_results", paste0(id, "_calendar_draws.rds")))
  if (id == "odi_3m") {
    write_csv(
      do.call(rbind, lapply(c(7, 5, 3), function(m) summarize_effect(draws$delta, m))),
      "04_results/primary_stricter_margins.csv"
    )
    saveRDS(
      gcomp(fit, b$target[!is.na(b$target$odi_3m), ], TRUE, 100),
      file.path(ROOT, "04_results/odi_3m_calendar_observed_target_draws.rds")
    )
  }
}
write_csv(do.call(rbind, res), "04_results/verified_calendar_primary_perioperative.csv")
write_csv(do.call(rbind, diag), "08_qa/verified_calendar_diagnostics.csv")
write_csv(do.call(rbind, audit), "08_qa/verified_calendar_provenance.csv")
