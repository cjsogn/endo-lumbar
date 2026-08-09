# =============================================================================
# ENDO-LUMBAR: 22 Operating time with the two erroneous registrations removed
#
# Two records in the analysed cohort carry an operating time of 580 minutes or
# more, one in each arm. They are separated from the rest of the distribution by
# a wide gap (the next highest value is 269 minutes), and both are internally
# contradictory: the 585-minute endoscopic case and the 580-minute microsurgical
# case are both recorded as day surgery, the endoscopic one with a single
# postoperative night and the microsurgical one with none. A procedure lasting
# nearly ten hours is not followed by same-day discharge, so both are treated as
# erroneous registrations and set to missing.
#
# The learning-curve script already applied this rule; the descriptive means,
# the corridor comparison and the adjusted Tier-4 model did not, so the two
# records were in some reported numbers and out of others. This script applies
# the rule once, everywhere, and writes the values the manuscript reads.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

OP_TIME_MAX <- 500     # minutes; above this a registration is treated as erroneous

out_dir <- file.path(paths$results, "supplement")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df_raw <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds")) %>% standardize_covs()

excluded <- df_raw %>%
  filter(operating_time >= OP_TIME_MAX) %>%
  select(treatment, operating_time, los_postop, day_surgery)
cat("=== Records excluded as erroneous operating-time registrations ===\n")
print(as.data.frame(excluded), row.names = FALSE)

df <- df_raw %>%
  mutate(operating_time = ifelse(operating_time >= OP_TIME_MAX, NA_real_, operating_time))

# -----------------------------------------------------------------------------
# (a) Descriptive summary by arm
# -----------------------------------------------------------------------------
cat("\n=== Operating time by arm, erroneous records removed ===\n")
desc <- df %>%
  group_by(treatment) %>%
  summarise(n = sum(!is.na(operating_time)),
            mean = mean(operating_time, na.rm = TRUE),
            sd = sd(operating_time, na.rm = TRUE),
            median = median(operating_time, na.rm = TRUE),
            q25 = quantile(operating_time, 0.25, na.rm = TRUE),
            q75 = quantile(operating_time, 0.75, na.rm = TRUE),
            max = max(operating_time, na.rm = TRUE), .groups = "drop")
print(as.data.frame(desc), row.names = FALSE, digits = 4)

# -----------------------------------------------------------------------------
# (b) Corridor comparison within the endoscopic arm
# -----------------------------------------------------------------------------
cat("\n=== Corridor comparison ===\n")
corr <- df %>%
  filter(treatment == "ELD") %>%
  mutate(corridor = c("Interlaminar", NA, "Transforaminal")[haven::zap_labels(OpTilgangV3)]) %>%
  filter(!is.na(corridor)) %>%
  group_by(corridor) %>%
  summarise(n = n(), n_time = sum(!is.na(operating_time)),
            mean = mean(operating_time, na.rm = TRUE),
            median = median(operating_time, na.rm = TRUE), .groups = "drop")
print(as.data.frame(corr), row.names = FALSE, digits = 4)

# -----------------------------------------------------------------------------
# (c) Learning curve, both cohorts
# -----------------------------------------------------------------------------
cat("\n=== Learning curve ===\n")
cohort_ids <- df$ForlopsID
lc_all <- readRDS(file.path(paths$data_clean, "df_all.rds")) %>%
  filter(ForlopsID %in% cohort_ids, treatment == "ELD") %>%
  arrange(surgery_date) %>%
  mutate(case_no = row_number(),
         operating_time = ifelse(operating_time >= OP_TIME_MAX, NA_real_, operating_time),
         pure_disc = !(stenosis_central == 1 | stenosis_lateral == 1 |
                       stenosis_foraminal == 1))
lc <- lc_all %>% filter(pure_disc) %>% mutate(quartile = ntile(case_no, 4))

rho_of <- function(d, y) {
  ct <- suppressWarnings(cor.test(d$case_no, d[[y]], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(!is.na(d[[y]])))
}
r_pure <- rho_of(lc, "operating_time")
r_all  <- rho_of(lc_all, "operating_time")
r_odi  <- rho_of(lc, "odi_3m")

qtab <- lc %>% group_by(quartile) %>%
  summarise(n = sum(!is.na(operating_time)),
            mean = mean(operating_time, na.rm = TRUE),
            median = median(operating_time, na.rm = TRUE), .groups = "drop")
print(as.data.frame(qtab), row.names = FALSE, digits = 4)
cat(sprintf("Pure disc (n=%d): rho %.3f, p %.3f\n", r_pure[["n"]], r_pure[["rho"]], r_pure[["p"]]))
cat(sprintf("All ELD  (n=%d): rho %.3f, p %.3f\n", r_all[["n"]], r_all[["rho"]], r_all[["p"]]))
cat(sprintf("ODI 3m   (n=%d): rho %.3f, p %.3f\n", r_odi[["n"]], r_odi[["rho"]], r_odi[["p"]]))

lc_out <- data.frame(
  analysis = c("Pure disc herniation subset", "All endoscopic cases",
               "Pure disc, disability at 3 months"),
  n = c(r_pure[["n"]], r_all[["n"]], r_odi[["n"]]),
  rho = c(r_pure[["rho"]], r_all[["rho"]], r_odi[["rho"]]),
  p = c(r_pure[["p"]], r_all[["p"]], r_odi[["p"]]))
write.csv(lc_out, file.path(paths$tables, "table_learning_curve_robustness.csv"),
          row.names = FALSE)
write.csv(qtab, file.path(paths$tables, "table_learning_curve_quartiles_clean.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# (d) Adjusted difference, refitted without the two records
# -----------------------------------------------------------------------------
cat("\n=== Adjusted difference in operating time ===\n")
fit <- brm(
  bf(as.formula(paste("operating_time | mi() ~ treatment +", cov_string))),
  data = df, family = gaussian(),
  prior = c(prior(normal(0, 30), class = "b", coef = "treatmentELD"),
            prior(normal(0, 2), class = "b"),
            prior(student_t(3, 0, 30), class = "sigma"),
            prior(normal(60, 30), class = "Intercept")),
  chains = mcmc_settings$chains, iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup, cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed, control = list(adapt_delta = 0.95),
  file = file.path(out_dir, "fit_tier4_optime_clean"), file_refit = "on_change"
)
ate <- compute_gcomp_ate(fit, df, "treatment", "continuous", lower_is_better = TRUE)
s <- summarize_ate(ate$ate)
cv <- check_convergence(fit)
cat(sprintf("Adjusted difference (MSD - ELD): %.1f min (95%% CrI %.1f to %.1f), Rhat %.4f, div %d\n",
            s$mean, s$cri_lo, s$cri_hi, cv$rhat_max, cv$n_divergent))

summ <- desc %>%
  mutate(adjusted_diff_msd_minus_eld = s$mean,
         adjusted_lo = s$cri_lo, adjusted_hi = s$cri_hi,
         n_excluded = nrow(excluded))
write.csv(summ, file.path(paths$tables, "table_operating_time_clean.csv"), row.names = FALSE)
write.csv(corr, file.path(paths$tables, "table_corridor_operating_time_clean.csv"),
          row.names = FALSE)

cat("\nOperating-time cleaning complete.\n")
