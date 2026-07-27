# =============================================================================
# ENDO-LUMBAR: 17 Refit of 12-month outcomes under date-based eligibility
#
# The original 12-month eligible cohort was defined as
#     surgery_date <= 2024-12-31  OR  Ferdigstilt1b12mnd == 1
# The second clause admits patients operated after the cutoff only if they
# responded, which conditions the analysis set on the outcome being observed.
#
# This script replaces it with a purely calendar-based rule
#     surgery_date <= 2025-02-28
# (the latest surgery with an observed 12-month ODI is 2025-03-07) and refits
# every 12-month outcome model.
#
# It also refits the EQ-5D models with priors appropriate to the EQ-5D scale.
# The original run applied `priors_continuous`, which was written for ODI on a
# 0-100 scale (Intercept normal(30, 20), sigma student_t(3, 0, 15),
# odi_baseline_z normal(1, 0.5)) to an outcome bounded near [-0.6, 1].
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

ELIGIBILITY_CUTOFF <- as.Date("2025-02-28")

out_dir <- file.path(paths$results, "refit12m")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))

# --- New 12-month eligible cohort --------------------------------------------
df_12m <- df_disc %>% filter(surgery_date <= ELIGIBILITY_CUTOFF)

cat("=== 12-month eligible cohort (date-based) ===\n")
cat(sprintf("  Cutoff: %s\n", format(ELIGIBILITY_CUTOFF)))
cat(sprintf("  N = %d (ELD = %d, MSD = %d)\n", nrow(df_12m),
            sum(df_12m$treatment == "ELD"), sum(df_12m$treatment == "MSD")))
for (tx in c("ELD", "MSD")) {
  s <- df_12m[df_12m$treatment == tx, ]
  cat(sprintf("  %s: ODI 12m observed %d/%d (%.1f%%)\n", tx,
              sum(!is.na(s$odi_12m)), nrow(s), 100 * mean(!is.na(s$odi_12m))))
}
saveRDS(df_12m, file.path(paths$data_clean, "df_disc_12m_eligible_datebased.rds"))

df_12m   <- standardize_covs(df_12m)
df_disc  <- standardize_covs(df_disc)

# --- EQ-5D priors on the EQ-5D scale (as described in eMethods 2.4) ----------
priors_eq5d <- c(
  prior(normal(0, 0.5),               class = "b", coef = "treatmentELD"),
  prior(normal(0, 0.5),               class = "b", coef = "odi_baseline_z"),
  prior(normal(0, 2),                 class = "b"),
  prior(student_t(3, 0, 0.3),         class = "sigma"),
  prior(normal(0.5, 0.5),             class = "Intercept")
)

# --- Generic fitter ----------------------------------------------------------
fit_outcome <- function(d, outcome_var, label, family_type, ni_margin,
                        lower_is_better, beta_upper = NULL,
                        zi_baseline_var = NULL, priors_use = NULL,
                        tag = "") {

  cat(sprintf("\n--- %s%s ---\n", label, tag))
  d_full <- d
  sf <- 1

  if (family_type == "gaussian") {
    fml <- bf(as.formula(paste(outcome_var, "| mi() ~ treatment +", cov_string)))
    if (is.null(priors_use)) priors_use <- priors_continuous
    fam <- gaussian()
  } else if (family_type == "zib") {
    zib_var <- paste0(outcome_var, "_zib")
    d[[zib_var]] <- transform_for_zib(d[[outcome_var]], upper = beta_upper)
    sf <- beta_upper
    zi_formula <- if (!is.null(zi_baseline_var)) {
      paste("zi ~ treatment +", zi_baseline_var)
    } else "zi ~ treatment"
    d_full <- d
    d <- d[!is.na(d[[outcome_var]]), ]
    fml <- bf(as.formula(paste(zib_var, "~ treatment +", cov_string)),
              as.formula(zi_formula))
    priors_use <- priors_zib
    fam <- zero_inflated_beta()
  } else {
    d_full <- d
    d <- d[!is.na(d[[outcome_var]]), ]
    fml <- bf(as.formula(paste(outcome_var, "~ treatment +", cov_string)))
    priors_use <- priors_binary_standard
    fam <- bernoulli()
  }

  cat(sprintf("  G-comp sample %d, model rows %d\n", nrow(d_full), nrow(d)))

  fit <- brm(
    formula = fml, data = d, family = fam, prior = priors_use,
    chains = mcmc_settings$chains, iter = mcmc_settings$iter,
    warmup = mcmc_settings$warmup, cores = min(mcmc_settings$chains, n_cores),
    seed = mcmc_settings$seed,
    control = list(adapt_delta = mcmc_settings$adapt_delta,
                   max_treedepth = mcmc_settings$max_treedepth),
    file = file.path(out_dir, paste0("fit_", outcome_var, tag)),
    file_refit = "on_change"
  )

  conv <- check_convergence(fit)
  otype <- if (family_type %in% c("gaussian", "zib")) "continuous" else "binary"
  ate <- compute_gcomp_ate(fit, newdata = d_full, treatment_var = "treatment",
                           outcome_type = otype,
                           lower_is_better = lower_is_better, scale_factor = sf)
  s <- summarize_ate(ate$ate, ni_margin = ni_margin)

  cat(sprintf("  ATE %.4f [%.4f, %.4f]  P(NI)=%.4f  P(Sup)=%.4f  Rhat=%.4f ESS=%.0f div=%d\n",
              s$mean, s$cri_lo, s$cri_hi, s$p_ni, s$p_superiority,
              conv$rhat_max, conv$ess_bulk_min, conv$n_divergent))

  tibble(
    Outcome = label, N = nrow(d_full), N_observed = nrow(d),
    ATE = s$mean, CrI_lo = s$cri_lo, CrI_hi = s$cri_hi,
    P_NI = s$p_ni, P_Superiority = s$p_superiority, NI_Margin = ni_margin,
    NI_Conclusion = ifelse(s$ni_conclusion, "NI demonstrated", "NI not demonstrated"),
    Rhat_max = conv$rhat_max, ESS_bulk_min = conv$ess_bulk_min,
    Divergent = conv$n_divergent
  )
}

# =============================================================================
# 12-month outcomes under the date-based eligible cohort
# =============================================================================

res <- bind_rows(
  fit_outcome(df_12m, "odi_12m", "ODI 12 months", "zib",
              ni_margins$odi, TRUE, beta_upper = 100,
              zi_baseline_var = "odi_baseline_z"),
  fit_outcome(df_12m, "nrs_back_12m", "NRS back pain 12 months", "zib",
              ni_margins$nrs_pain, TRUE, beta_upper = 10,
              zi_baseline_var = "nrs_back_baseline_z"),
  fit_outcome(df_12m, "nrs_leg_12m", "NRS leg pain 12 months", "zib",
              ni_margins$nrs_pain, TRUE, beta_upper = 10,
              zi_baseline_var = "nrs_leg_baseline_z"),
  fit_outcome(df_12m, "eq5d_12m", "EQ-5D 12 months", "gaussian",
              ni_margins$eq5d, FALSE, priors_use = priors_eq5d),
  fit_outcome(df_12m, "rtw_12m", "Return to work 12 months", "bernoulli",
              ni_margins$rtw, FALSE),
  fit_outcome(df_12m, "analgesic_12m", "Analgesic use 12 months", "bernoulli",
              ni_margins$analgesic, TRUE),
  fit_outcome(df_12m, "satisfied_12m", "Satisfaction 12 months", "bernoulli",
              ni_margins$satisfaction, FALSE),
  fit_outcome(df_12m, "gpe_success_12m", "GPE success 12 months", "bernoulli",
              ni_margins$gpe, FALSE)
)

write.csv(res, file.path(paths$tables, "table3b_tier2_12m_datebased.csv"),
          row.names = FALSE)

# =============================================================================
# EQ-5D 3 months refitted with scale-appropriate priors (sensitivity)
# =============================================================================

res_eq3 <- fit_outcome(df_disc, "eq5d_3m", "EQ-5D 3 months", "gaussian",
                       ni_margins$eq5d, FALSE, priors_use = priors_eq5d,
                       tag = "_scaledpriors")
write.csv(res_eq3, file.path(paths$tables, "table_eq5d_3m_scaled_priors.csv"),
          row.names = FALSE)

# =============================================================================
# Reoperation rates under the same date-based cohort
# =============================================================================

cat("\n=== Reoperation denominators under date-based eligibility ===\n")
cat(sprintf("  12m-eligible (date-based): ELD %d, MSD %d\n",
            sum(df_12m$treatment == "ELD"), sum(df_12m$treatment == "MSD")))

cat("\n=== Refit complete ===\n")
print(as.data.frame(res), row.names = FALSE)
