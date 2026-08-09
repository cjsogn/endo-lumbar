# =============================================================================
# ENDO-LUMBAR: 20 Calendar-time sensitivity for the perioperative outcomes
#
# Calendar time was pre-specified in the analysis plan as a restricted cubic
# spline interacted with treatment and was omitted from the implemented models
# (Deviation 3). The primary outcome already has a calendar-time sensitivity
# analysis; the three Tier-3 superiority outcomes did not. Surgery year is the
# largest between-arm imbalance in the cohort (standardised mean difference
# 0.53), and day-surgery pathways expanded over the study period, so this is
# the confounder those three estimates are most exposed to.
#
# This script does two things:
#   (a) refits day surgery, postoperative length of stay and patient-reported
#       complications with calendar time added to the adjustment set, using a
#       restricted cubic spline with 3 knots as the analysis plan specified;
#   (b) refits the primary outcome with the same spline form, so that the
#       calendar-time sensitivity analysis reported for the primary endpoint
#       matches the pre-specified specification and not only the linear
#       approximation used in script 11.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

library(splines)

out_dir <- file.path(paths$results, "supplement")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds")) %>% standardize_covs()

# Restricted cubic spline with 3 knots = natural cubic spline basis on 2 df.
# Centred and scaled first so the basis is numerically well behaved.
df$calendar_time <- as.numeric(df$surgery_date - as.Date("2022-01-01"))
df$calendar_time_z <- as.numeric(scale(df$calendar_time))
ct_basis <- ns(df$calendar_time_z, df = 2)
df$ct1 <- ct_basis[, 1]
df$ct2 <- ct_basis[, 2]

cov_string_ct <- paste(cov_string, "+ ct1 + ct2")

cat("=== Calendar-time sensitivity for Tier 3 and the primary outcome ===\n")
cat(sprintf("Surgery year by arm:\n"))
print(round(prop.table(table(df$treatment, format(df$surgery_date, "%Y")), 1), 3))

res <- list()

# -----------------------------------------------------------------------------
# (a1) Day surgery
# -----------------------------------------------------------------------------
cat("\n--- Day surgery, calendar time added ---\n")
fit_ds <- brm(
  bf(as.formula(paste("day_surgery ~ treatment +", cov_string_ct))),
  data = df, family = bernoulli(), prior = priors_binary_standard,
  chains = mcmc_settings$chains, iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup, cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = file.path(out_dir, "fit_tier3_day_surgery_caltime"), file_refit = "on_change"
)
a <- compute_gcomp_ate(fit_ds, df, "treatment", "binary", lower_is_better = FALSE)
s <- summarize_ate(a$ate)
cv <- check_convergence(fit_ds)
res$day_surgery <- tibble(
  Outcome = "Day surgery rate", Effect = s$mean, CrI_lo = s$cri_lo, CrI_hi = s$cri_hi,
  P_Superiority = mean(a$ate > 0), Rhat = cv$rhat_max, Divergent = cv$n_divergent
)
cat(sprintf("  RD %.3f [%.3f, %.3f], P(sup) %.4f\n", s$mean, s$cri_lo, s$cri_hi,
            mean(a$ate > 0)))

# -----------------------------------------------------------------------------
# (a2) Postoperative length of stay (ordinal)
# -----------------------------------------------------------------------------
cat("\n--- Postoperative length of stay, calendar time added ---\n")
df_los <- df %>%
  filter(!is.na(los_postop)) %>%
  mutate(los_ordinal = ordered(pmin(los_postop, 3L)))

fit_los <- brm(
  bf(as.formula(paste("los_ordinal ~ treatment +", cov_string_ct))),
  data = df_los, family = cumulative("logit"),
  prior = c(prior(normal(0, 2), class = "b"), prior(normal(0, 4), class = "Intercept")),
  chains = mcmc_settings$chains, iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup, cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  file = file.path(out_dir, "fit_tier3_los_caltime"), file_refit = "on_change"
)
log_or <- -as_draws_df(fit_los)$b_treatmentELD   # positive = shorter stay after ELD
cv <- check_convergence(fit_los)
res$los <- tibble(
  Outcome = "Postoperative length of stay (log OR)", Effect = mean(log_or),
  CrI_lo = unname(quantile(log_or, 0.025)), CrI_hi = unname(quantile(log_or, 0.975)),
  P_Superiority = mean(log_or > 0), Rhat = cv$rhat_max, Divergent = cv$n_divergent
)
cat(sprintf("  log OR %.3f [%.3f, %.3f] (OR %.2f [%.2f, %.2f]), P(sup) %.4f\n",
            mean(log_or), quantile(log_or, 0.025), quantile(log_or, 0.975),
            exp(mean(log_or)), exp(quantile(log_or, 0.025)), exp(quantile(log_or, 0.975)),
            mean(log_or > 0)))

# -----------------------------------------------------------------------------
# (a3) Patient-reported complications at 3 months
# -----------------------------------------------------------------------------
cat("\n--- Patient-reported complications, calendar time added ---\n")
fit_pc <- brm(
  bf(as.formula(paste("pt_comp_any_3m ~ treatment +", cov_string_ct))),
  data = df, family = bernoulli(), prior = priors_binary_standard,
  chains = mcmc_settings$chains, iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup, cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  file = file.path(out_dir, "fit_tier3_ptcomp_caltime"), file_refit = "on_change"
)
a <- compute_gcomp_ate(fit_pc, df, "treatment", "binary", lower_is_better = TRUE)
s <- summarize_ate(a$ate)
cv <- check_convergence(fit_pc)
res$pt_comp <- tibble(
  Outcome = "Patient-reported complications 3 months", Effect = s$mean,
  CrI_lo = s$cri_lo, CrI_hi = s$cri_hi, P_Superiority = mean(a$ate > 0),
  Rhat = cv$rhat_max, Divergent = cv$n_divergent
)
cat(sprintf("  RD %.3f [%.3f, %.3f], P(sup) %.4f\n", s$mean, s$cri_lo, s$cri_hi,
            mean(a$ate > 0)))

tier3_ct <- bind_rows(res)
write.csv(tier3_ct, file.path(paths$tables, "table_tier3_calendar_time.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# (b) Primary outcome with the pre-specified spline form
# -----------------------------------------------------------------------------
cat("\n--- Primary outcome, treatment x calendar-time spline ---\n")
df$odi_3m_zib <- transform_for_zib(df$odi_3m, upper = 100)
df_cc <- df %>% filter(!is.na(odi_3m))

cov_string_ct_int <- paste(cov_string_ct, "+ treatment:ct1 + treatment:ct2")

fit_prim_ct <- brm(
  bf(as.formula(paste("odi_3m_zib ~ treatment +", cov_string_ct_int)),
     zi ~ treatment + odi_baseline_z),
  data = df_cc, family = zero_inflated_beta(), prior = priors_zib,
  chains = mcmc_settings$chains, iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup, cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta,
                 max_treedepth = mcmc_settings$max_treedepth),
  file = file.path(out_dir, "fit_primary_caltime_spline"), file_refit = "on_change"
)
a <- compute_gcomp_ate(fit_prim_ct, newdata = df, outcome_type = "continuous",
                       lower_is_better = TRUE, scale_factor = 100)
s <- summarize_ate(a$ate, ni_margin = ni_margins$odi)
cv <- check_convergence(fit_prim_ct)
cat(sprintf("  ATE %.3f [%.3f, %.3f], P(NI) %.4f (Rhat %.4f, div %d)\n",
            s$mean, s$cri_lo, s$cri_hi, s$p_ni, cv$rhat_max, cv$n_divergent))

prim_ct <- tibble(
  Analysis = "Treatment x calendar-time spline (3 knots)",
  ATE = s$mean, CrI_lo = s$cri_lo, CrI_hi = s$cri_hi, P_NI = s$p_ni
)
write.csv(prim_ct, file.path(paths$tables, "table_primary_calendar_spline.csv"),
          row.names = FALSE)

cat("\nCalendar-time sensitivity complete.\n")
print(as.data.frame(tier3_ct), row.names = FALSE, digits = 4)
