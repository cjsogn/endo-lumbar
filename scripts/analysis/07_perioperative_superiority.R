# =============================================================================
# ENDO-LUMBAR: 07 Perioperative/Safety Outcomes (Tier 3 - Superiority)
# SAP Sections 6.4, 7, 15.2
# =============================================================================

source("/Users/cjsogn/endo_studies/lumbar/analysis/scripts/00_config.R")

df_disc <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))
primary_results <- readRDS(file.path(paths$output, "primary_results.rds"))

cat("=== Tier 3: Perioperative and Safety Outcomes (Superiority) ===\n")

# Check gating status (SAP Section 7)
gate_open <- primary_results$gate_open
cat(sprintf("Primary NI: %s -> Gate: %s\n",
            ifelse(gate_open, "DEMONSTRATED", "NOT DEMONSTRATED"),
            ifelse(gate_open, "OPEN (confirmatory)", "CLOSED (descriptive only)")))

# Prepare data (standardize_covs and cov_string from 00_config.R)
df_disc <- standardize_covs(df_disc)

# =============================================================================
# RATE-BASED GATING (SAP Section 6.4)
# =============================================================================

cat("\n=== Rate-Based Gating Assessment ===\n")

rate_check <- function(outcome_var, label) {
  observed <- df_disc[[outcome_var]]
  n_obs <- sum(!is.na(observed))
  n_events <- sum(observed == 1, na.rm = TRUE)
  rate <- n_events / n_obs
  cat(sprintf("  %s: %d/%d (%.1f%%)", label, n_events, n_obs, 100 * rate))

  gating <- case_when(
    rate >= 0.10 ~ "Formal test",
    rate >= 0.05 ~ "Test with caution",
    TRUE ~ "Descriptive only"
  )
  cat(sprintf(" -> %s\n", gating))

  # EPV check (SAP Section 24.1)
  n_params <- 25  # approximate number of model parameters
  epv <- min(n_events, n_obs - n_events) / n_params
  cat(sprintf("    EPV: %.1f", epv))
  if (epv < 10) cat(" (LOW: regularized priors needed)")
  cat("\n")

  list(rate = rate, gating = gating, n_events = n_events, n_obs = n_obs, epv = epv)
}

perop_gate <- rate_check("perop_comp_any", "Peroperative complications (surgeon)")
ptcomp_gate <- rate_check("pt_comp_any_3m", "Patient-reported complications (3m)")
reop_gate <- rate_check("reop_during_stay", "Reoperation during stay")
daysurg_gate <- rate_check("day_surgery", "Day surgery")

# =============================================================================
# FIT TIER 3 MODELS
# =============================================================================

tier3_results <- list()

# --- Day Surgery Rate ---
cat("\n--- Day Surgery Rate ---\n")

# Priors for binary outcomes (SAP Sections 15.2, 24.1-24.2)
# Cauchy (student-t df=1) prior provides adaptive shrinkage for low-EPV settings:
# spike near zero shrinks noise covariates; heavy tails preserve strong confounders.
# This approximates horseshoe behavior while allowing coefficient-specific treatment prior.
# (brms does not allow mixing horseshoe special priors with coefficient-specific priors.)
priors_binary <- c(
  prior(normal(0, 1), class = "b", coef = "treatmentELD"),
  prior(student_t(1, 0, 1), class = "b"),
  prior(normal(0, 5), class = "Intercept")
)

cat(sprintf("  EPV: %.1f — using horseshoe priors for adaptive regularization\n",
            daysurg_gate$epv))

fit_daysurg <- brm(
  bf(as.formula(paste("day_surgery ~ treatment +", cov_string))),
  data = df_disc,
  family = bernoulli(),
  prior = priors_binary,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = file.path(paths$models, "fit_tier3_day_surgery"),
  file_refit = "on_change"
)

ate_daysurg <- compute_gcomp_ate(fit_daysurg, df_disc, "treatment", "binary",
                                  lower_is_better = FALSE)
tier3_results$day_surgery <- list(
  label = "Day surgery rate",
  ate_draws = ate_daysurg,
  ate_summary = summarize_ate(ate_daysurg$ate),
  gating = daysurg_gate
)

# --- Length of Stay (Postoperative, Ordinal) ---
cat("\n--- Length of Stay (Postoperative, Ordinal Cumulative Model) ---\n")

# Using los_postop (surgery to discharge) rather than los_total (admission to discharge),
# because los_total includes preoperative admission days unrelated to the surgical technique.
# LOS is modeled as an ordinal outcome with 4 categories: 0 (day surgery), 1 (one night),
# 2 (two nights), 3+ (extended stay). This is appropriate because the data consists of
# discrete point masses with rare extreme outliers (up to 62 days). Continuous
# and count models are misspecified for this distribution. The ordinal cumulative (proportional
# odds) model estimates the adjusted odds of being in a lower LOS category for ELD vs MSD.
# cumulative() does not support mi(), so we use complete cases (< 5% missing).
df_los <- df_disc %>% filter(!is.na(los_postop))

# Create ordinal LOS variable: 0 = day surgery, 1 = one night, 2 = two nights, 3+ = extended
df_los <- df_los %>%
  mutate(
    los_ordinal = case_when(
      los_postop == 0 ~ 0L,
      los_postop == 1 ~ 1L,
      los_postop == 2 ~ 2L,
      los_postop >= 3 ~ 3L
    ),
    los_ordinal = ordered(los_ordinal)
  )

cat(sprintf("  Complete cases: %d/%d (%.1f%% missing)\n",
            nrow(df_los), nrow(df_disc), 100 * (1 - nrow(df_los) / nrow(df_disc))))
cat(sprintf("  Ordinal categories:\n"))
cat(sprintf("    0 (day surgery):  %d (%.1f%%)\n",
            sum(df_los$los_ordinal == 0), 100 * mean(df_los$los_ordinal == 0)))
cat(sprintf("    1 (one night):    %d (%.1f%%)\n",
            sum(df_los$los_ordinal == 1), 100 * mean(df_los$los_ordinal == 1)))
cat(sprintf("    2 (two nights):   %d (%.1f%%)\n",
            sum(df_los$los_ordinal == 2), 100 * mean(df_los$los_ordinal == 2)))
cat(sprintf("    3+ (extended):    %d (%.1f%%)\n",
            sum(df_los$los_ordinal == 3), 100 * mean(df_los$los_ordinal == 3)))

fit_los <- brm(
  bf(as.formula(paste("los_ordinal ~ treatment +", cov_string))),
  data = df_los,
  family = cumulative("logit"),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(normal(0, 4), class = "Intercept")
  ),
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = file.path(paths$models, "fit_tier3_los_ordinal"),
  file_refit = "on_change"
)

# Extract treatment coefficient (log-OR for proportional odds)
# In cumulative models, negative coefficient = lower category more likely
# Convention: we want positive = ELD superior (lower LOS), so negate
draws_los <- as_draws_df(fit_los)
log_or_draws <- -draws_los$b_treatmentELD  # negate: positive = ELD in lower category

log_or_summary <- list(
  mean = mean(log_or_draws),
  median = median(log_or_draws),
  sd = sd(log_or_draws),
  cri_lo = unname(quantile(log_or_draws, 0.025)),
  cri_hi = unname(quantile(log_or_draws, 0.975)),
  p_superiority = mean(log_or_draws > 0),
  or_mean = exp(mean(log_or_draws)),
  or_lo = exp(unname(quantile(log_or_draws, 0.025))),
  or_hi = exp(unname(quantile(log_or_draws, 0.975)))
)

cat(sprintf("  Log-OR (lower LOS): %.2f [%.2f, %.2f]\n",
            log_or_summary$mean, log_or_summary$cri_lo, log_or_summary$cri_hi))
cat(sprintf("  OR: %.2f [%.2f, %.2f]\n",
            log_or_summary$or_mean, log_or_summary$or_lo, log_or_summary$or_hi))
cat(sprintf("  P(Superiority): %.4f\n", log_or_summary$p_superiority))

# G-computation: predicted category probabilities under each treatment
nd_eld <- df_los %>% mutate(treatment = factor("ELD", levels = levels(treatment)))
nd_msd <- df_los %>% mutate(treatment = factor("MSD", levels = levels(treatment)))
pp_eld <- posterior_epred(fit_los, newdata = nd_eld)  # dims: draws x patients x categories
pp_msd <- posterior_epred(fit_los, newdata = nd_msd)

# Mean predicted probability per category
cat("  Predicted category probabilities:\n")
for (k in 1:4) {
  p_eld_k <- mean(rowMeans(pp_eld[, , k]))
  p_msd_k <- mean(rowMeans(pp_msd[, , k]))
  cat(sprintf("    Category %d: ELD=%.1f%%, MSD=%.1f%%\n", k - 1, 100 * p_eld_k, 100 * p_msd_k))
}

tier3_results$los <- list(
  label = "Length of stay (postop)",
  ate_draws = list(ate = log_or_draws),
  ate_summary = log_or_summary,
  gating = list(gating = "Ordinal"),
  model_note = "Ordinal cumulative model on los_postop (0/1/2/3+ days, complete cases)"
)

# --- Peroperative Complications (surgeon-reported composite) ---
cat("\n--- Peroperative Complications (Surgeon-Reported) ---\n")

if (perop_gate$gating != "Descriptive only") {
  cat(sprintf("  EPV: %.1f — using horseshoe priors for adaptive regularization\n",
              perop_gate$epv))

  fit_perop <- brm(
    bf(as.formula(paste("perop_comp_any ~ treatment +", cov_string))),
    data = df_disc,
    family = bernoulli(),
    prior = priors_binary,
    chains = mcmc_settings$chains,
    iter = mcmc_settings$iter,
    warmup = mcmc_settings$warmup,
    cores = min(mcmc_settings$chains, n_cores),
    seed = mcmc_settings$seed,
    control = list(adapt_delta = 0.99, max_treedepth = 14),
    file = file.path(paths$models, "fit_tier3_perop_comp"),
    file_refit = "on_change"
  )

  ate_perop <- compute_gcomp_ate(fit_perop, df_disc, "treatment", "binary",
                                  lower_is_better = TRUE)
  tier3_results$perop_comp <- list(
    label = "Peroperative complications",
    ate_draws = ate_perop,
    ate_summary = summarize_ate(ate_perop$ate),
    gating = perop_gate
  )
} else {
  cat("  Descriptive only (rate < 5%)\n")
  tier3_results$perop_comp <- list(
    label = "Peroperative complications",
    descriptive_only = TRUE,
    gating = perop_gate
  )
}

# --- Patient-Reported Complications (3-month) ---
cat("\n--- Patient-Reported Complications (3-month) ---\n")

if (ptcomp_gate$gating != "Descriptive only") {
  cat(sprintf("  EPV: %.1f — using horseshoe priors for adaptive regularization\n",
              ptcomp_gate$epv))

  fit_ptcomp <- brm(
    bf(as.formula(paste("pt_comp_any_3m ~ treatment +", cov_string))),
    data = df_disc,
    family = bernoulli(),
    prior = priors_binary,
    chains = mcmc_settings$chains,
    iter = mcmc_settings$iter,
    warmup = mcmc_settings$warmup,
    cores = min(mcmc_settings$chains, n_cores),
    seed = mcmc_settings$seed,
    control = list(adapt_delta = 0.99, max_treedepth = 14),
    file = file.path(paths$models, "fit_tier3_pt_comp"),
    file_refit = "on_change"
  )

  ate_ptcomp <- compute_gcomp_ate(fit_ptcomp, df_disc, "treatment", "binary",
                                   lower_is_better = TRUE)
  tier3_results$pt_comp <- list(
    label = "Patient-reported complications (3m)",
    ate_draws = ate_ptcomp,
    ate_summary = summarize_ate(ate_ptcomp$ate),
    gating = ptcomp_gate
  )
} else {
  cat("  Descriptive only (rate < 5%)\n")
  tier3_results$pt_comp <- list(
    label = "Patient-reported complications (3m)",
    descriptive_only = TRUE,
    gating = ptcomp_gate
  )
}

# --- Reoperations ---
cat("\n--- Reoperations During Hospital Stay ---\n")

if (reop_gate$gating != "Descriptive only") {
  fit_reop <- brm(
    bf(as.formula(paste("reop_during_stay ~ treatment +", cov_string))),
    data = df_disc,
    family = bernoulli(),
    prior = priors_binary,
    chains = mcmc_settings$chains,
    iter = mcmc_settings$iter,
    warmup = mcmc_settings$warmup,
    cores = min(mcmc_settings$chains, n_cores),
    seed = mcmc_settings$seed,
    control = list(adapt_delta = 0.99, max_treedepth = 14),
    file = file.path(paths$models, "fit_tier3_reop"),
    file_refit = "on_change"
  )

  ate_reop <- compute_gcomp_ate(fit_reop, df_disc, "treatment", "binary",
                                 lower_is_better = TRUE)
  tier3_results$reop <- list(
    label = "Reoperation during stay",
    ate_draws = ate_reop,
    ate_summary = summarize_ate(ate_reop$ate),
    gating = reop_gate
  )
} else {
  cat("  Descriptive only (rate < 5%)\n")
  tier3_results$reop <- list(
    label = "Reoperation during stay",
    descriptive_only = TRUE,
    gating = reop_gate
  )
}

# =============================================================================
# COMPLICATION COMPONENTS (descriptive, SAP Section 6.4)
# =============================================================================

cat("\n=== Individual Complication Categories (Descriptive) ===\n")

comp_vars <- c("perop_dural_tear", "perop_nerve_injury", "perop_wrong_level",
               "perop_bleeding", "perop_other",
               "pt_comp_wound_superficial", "pt_comp_wound_deep",
               "pt_comp_uti", "pt_comp_pneumonia", "pt_comp_dvt",
               "pt_comp_pe", "pt_comp_bleeding")

comp_descriptive <- map_dfr(comp_vars, function(v) {
  by_group <- df_disc %>%
    group_by(treatment) %>%
    dplyr::summarise(
      n_events = sum(.data[[v]] == 1, na.rm = TRUE),
      n_obs = sum(!is.na(.data[[v]])),
      rate = n_events / n_obs,
      .groups = "drop"
    )

  tibble(
    complication = v,
    ELD_events = by_group$n_events[by_group$treatment == "ELD"],
    ELD_n = by_group$n_obs[by_group$treatment == "ELD"],
    ELD_rate = by_group$rate[by_group$treatment == "ELD"],
    MSD_events = by_group$n_events[by_group$treatment == "MSD"],
    MSD_n = by_group$n_obs[by_group$treatment == "MSD"],
    MSD_rate = by_group$rate[by_group$treatment == "MSD"]
  )
})

print(comp_descriptive)
write.csv(comp_descriptive, file.path(paths$tables, "complication_components_disc.csv"),
          row.names = FALSE)

# =============================================================================
# SUMMARY TABLE (Table 4)
# =============================================================================

cat("\n=== Tier 3 Summary ===\n")
tier3_table <- map_dfr(tier3_results, function(r) {
  if (isTRUE(r$descriptive_only)) {
    tibble(
      Outcome = r$label,
      ATE = NA_real_,
      CrI_lo = NA_real_,
      CrI_hi = NA_real_,
      P_Superiority = NA_real_,
      Gating = r$gating$gating,
      Note = "Descriptive only (rate < 5%)"
    )
  } else {
    tibble(
      Outcome = r$label,
      ATE = r$ate_summary$mean,
      CrI_lo = r$ate_summary$cri_lo,
      CrI_hi = r$ate_summary$cri_hi,
      P_Superiority = r$ate_summary$p_superiority,
      Gating = ifelse(is.null(r$gating$gating), "Continuous", r$gating$gating),
      Note = ifelse(gate_open, "Confirmatory (gate open)", "Descriptive (gate closed)")
    )
  }
})

print(tier3_table)
write.csv(tier3_table, file.path(paths$tables, "table4_tier3_superiority.csv"),
          row.names = FALSE)

saveRDS(tier3_results, file.path(paths$output, "tier3_results.rds"))

cat("\nTier 3 analysis complete.\n")
