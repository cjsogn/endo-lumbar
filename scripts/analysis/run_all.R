# =============================================================================
# ENDO-LUMBAR: Master Run Script
# Executes all analyses in sequence with timing and error handling
# =============================================================================

cat("====================================================================\n")
cat("  ENDO-LUMBAR: Full Analysis Pipeline (Re-run)\n")
cat("  Endoscopic vs Microsurgical Lumbar Discectomy\n")
cat("  SAP Version 2.0 | March 2026\n")
cat("====================================================================\n\n")

start_time <- Sys.time()
script_dir <- "/Users/cjsogn/ENDO_LUMBAR/scripts"

scripts <- c(
  "01_data_preparation.R",
  "02_descriptive_table1.R",
  "03_propensity_balance.R",
  "04_primary_analysis.R",
  "05_mcmc_diagnostics.R",
  "06_secondary_effectiveness.R",
  "07_perioperative_superiority.R",
  "08_descriptive_tier4.R",
  "09_prior_sensitivity.R",
  "10_missing_data_sensitivity.R",
  "11_model_sensitivity.R",
  "12_falsification_evalue.R",
  "13_subgroups_causal_forest.R",
  "14_stenosis_exploratory.R",
  "15_tables_figures.R",
  "16_covariate_influence.R",
  "17_eld_approach_comparison.R",
  "18_learning_curve.R",
  "19_odi_spider_plot.R",
  "20_publication_figures.R"
)

results <- list()

for (i in seq_along(scripts)) {
  s <- scripts[i]
  path <- file.path(script_dir, s)

  if (!file.exists(path)) {
    cat(sprintf("\n[%d/%d] SKIP: %s (not found)\n", i, length(scripts), s))
    results[[s]] <- list(status = "SKIPPED", time = 0)
    next
  }

  cat(sprintf("\n[%d/%d] %s ... (%s)\n", i, length(scripts), s, Sys.time()))
  t0 <- Sys.time()

  status <- tryCatch({
    source(path, local = new.env(parent = globalenv()))
    "OK"
  }, error = function(e) {
    cat(sprintf("\n  [ERROR] %s: %s\n", s, conditionMessage(e)))
    paste("ERROR:", conditionMessage(e))
  })

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  [%s] %s (%.1f min)\n",
              ifelse(status == "OK", "DONE", "FAIL"), s, elapsed))

  results[[s]] <- list(status = status, time = elapsed)
}

# --- Summary ---
total <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

cat("\n====================================================================\n")
cat("  PIPELINE SUMMARY\n")
cat(sprintf("  Total time: %.1f min (%.1f hours)\n", total, total/60))
cat("====================================================================\n")

for (s in names(results)) {
  r <- results[[s]]
  cat(sprintf("  %-45s %s (%.1f min)\n", s, r$status, r$time))
}

n_ok <- sum(sapply(results, function(r) r$status == "OK"))
n_fail <- sum(sapply(results, function(r) grepl("^ERROR", r$status)))
cat(sprintf("\n  OK: %d | FAILED: %d | SKIPPED: %d\n",
            n_ok, n_fail, length(scripts) - n_ok - n_fail))

# Print key results if available
if (file.exists(file.path("/Users/cjsogn/endo_studies/lumbar/analysis/output/", "primary_results.rds"))) {
  primary <- readRDS(file.path("/Users/cjsogn/endo_studies/lumbar/analysis/output/", "primary_results.rds"))
  cat(sprintf("\n  PRIMARY RESULT:\n"))
  cat(sprintf("    ATE (ODI 3m): %.2f (95%% CrI: [%.2f, %.2f])\n",
              primary$ate_summary$mean, primary$ate_summary$cri_lo,
              primary$ate_summary$cri_hi))
  cat(sprintf("    P(NI): %.4f -> %s\n",
              primary$ate_summary$p_ni,
              ifelse(primary$ate_summary$ni_conclusion, "NON-INFERIOR", "NOT DEMONSTRATED")))
}
