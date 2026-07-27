# =============================================================================
# ENDO-LUMBAR: 19 Supplementary analyses added for the Neurospine submission
#
#   (a) Surgical level added to the adjustment set. Level was pre-specified in
#       the analysis plan but omitted from the implemented covariate set.
#   (b) Multiplicity adjustment (Holm and Bonferroni) recomputed on the current
#       posterior probabilities. The earlier table was built from a superseded
#       run and no longer matched Table 2.
#   (c) Loss-to-follow-up comparison at 3 months, tabulated for the supplement.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

out_dir <- file.path(paths$results, "supplement")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds")) %>% standardize_covs()

# Operated level was pre-specified as a covariate but was not carried into the
# analysis dataset, so it is merged back in from the full extract by ForlopsID.
lvl <- readRDS(file.path(paths$data_clean, "df_all.rds")) %>%
  select(ForlopsID, level_L34, level_L45, level_L5S1)
df <- left_join(df, lvl, by = "ForlopsID")
stopifnot(!any(is.na(df$level_L45)))

# =============================================================================
# (a) Surgical level sensitivity analysis
# =============================================================================

cat("=== Surgical level sensitivity analysis ===\n")
cat("Level distribution by arm:\n")
print(df %>% group_by(treatment) %>%
        summarise(L3_L4 = mean(level_L34), L4_L5 = mean(level_L45),
                  L5_S1 = mean(level_L5S1), .groups = "drop"))

df$odi_3m_zib <- transform_for_zib(df$odi_3m, upper = 100)
df_full <- df
df_cc <- df %>% filter(!is.na(odi_3m))

cov_string_level <- paste(cov_string, "+ level_L34 + level_L45 + level_L5S1")

fit_level <- brm(
  bf(as.formula(paste("odi_3m_zib ~ treatment +", cov_string_level)),
     zi ~ treatment + odi_baseline_z),
  data = df_cc, family = zero_inflated_beta(), prior = priors_zib,
  chains = mcmc_settings$chains, iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup, cores = min(mcmc_settings$chains, n_cores),
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta,
                 max_treedepth = mcmc_settings$max_treedepth),
  file = file.path(out_dir, "fit_primary_plus_level"), file_refit = "on_change"
)

conv <- check_convergence(fit_level)
ate <- compute_gcomp_ate(fit_level, newdata = df_full, outcome_type = "continuous",
                         lower_is_better = TRUE, scale_factor = 100)
s <- summarize_ate(ate$ate, ni_margin = ni_margins$odi)

cat(sprintf("\nPrimary + surgical level: ATE %.3f [%.3f, %.3f], P(NI) %.4f (Rhat %.4f, div %d)\n",
            s$mean, s$cri_lo, s$cri_hi, s$p_ni, conv$rhat_max, conv$n_divergent))

lvl_tab <- tibble(
  Analysis = c("Primary (27 covariates)", "Primary + surgical level (30 covariates)"),
  ATE = c(0.308573408972219, s$mean),
  CrI_lo = c(-2.8518479775768, s$cri_lo),
  CrI_hi = c(3.53525198407138, s$cri_hi),
  P_NI = c(1, s$p_ni)
)
write.csv(lvl_tab, file.path(paths$tables, "table_level_sensitivity.csv"),
          row.names = FALSE)

# =============================================================================
# (b) Multiplicity adjustment, recomputed on the current posteriors
# =============================================================================

cat("\n=== Multiplicity adjustment (recomputed) ===\n")

t2 <- read.csv(file.path(paths$tables, "table2_primary_results.csv"), check.names = FALSE)
t3 <- read.csv(file.path(paths$tables, "table3_tier2_effectiveness.csv"), check.names = FALSE)

mult <- bind_rows(
  tibble(Outcome = "ODI 3 months (primary)", P_NI = as.numeric(t2$`P(NI)`)),
  tibble(Outcome = t3$Outcome, P_NI = t3$P_NI)
)

# Treat 1 - P(NI) as the evidence against non-inferiority and adjust that.
mult <- mult %>%
  mutate(p_one_sided = pmax(1 - P_NI, 1 / 4000),
         p_holm = p.adjust(p_one_sided, method = "holm"),
         p_bonf = p.adjust(p_one_sided, method = "bonferroni"),
         PNI_holm = 1 - p_holm,
         PNI_bonf = 1 - p_bonf,
         pass_unadjusted = P_NI > 0.95,
         pass_holm = PNI_holm > 0.95,
         pass_bonferroni = PNI_bonf > 0.95)

write.csv(mult, file.path(paths$tables, "table_multiplicity.csv"), row.names = FALSE)
print(as.data.frame(mult[, c("Outcome", "P_NI", "PNI_holm", "PNI_bonf",
                             "pass_unadjusted", "pass_holm", "pass_bonferroni")]),
      row.names = FALSE, digits = 4)
cat(sprintf("\nPass unadjusted: %d/%d | Holm: %d/%d | Bonferroni: %d/%d\n",
            sum(mult$pass_unadjusted), nrow(mult),
            sum(mult$pass_holm), nrow(mult),
            sum(mult$pass_bonferroni), nrow(mult)))

# =============================================================================
# (c) Loss to follow-up at 3 months
# =============================================================================

cat("\n=== Loss to follow-up at 3 months ===\n")
df$ltfu_3m <- as.integer(is.na(df$odi_3m))

vars <- c("age", "bmi", "odi_baseline", "eq5d_baseline",
          "nrs_back_baseline", "nrs_leg_baseline")
ltfu <- map_dfr(vars, function(v) {
  a <- df[[v]][df$ltfu_3m == 0]; b <- df[[v]][df$ltfu_3m == 1]
  tt <- t.test(a, b)
  tibble(Variable = v,
         Followed = sprintf("%.2f (%.2f)", mean(a, na.rm = TRUE), sd(a, na.rm = TRUE)),
         Lost = sprintf("%.2f (%.2f)", mean(b, na.rm = TRUE), sd(b, na.rm = TRUE)),
         p = sprintf("%.3f", tt$p.value))
})
sex_t <- chisq.test(table(df$sex, df$ltfu_3m))
tx_t  <- chisq.test(table(df$treatment, df$ltfu_3m))
ltfu <- bind_rows(ltfu,
  tibble(Variable = "Female sex, %",
         Followed = sprintf("%.1f", 100 * mean(df$sex[df$ltfu_3m == 0] == "Female")),
         Lost = sprintf("%.1f", 100 * mean(df$sex[df$ltfu_3m == 1] == "Female")),
         p = sprintf("%.3f", sex_t$p.value)),
  tibble(Variable = "Endoscopic treatment, %",
         Followed = sprintf("%.1f", 100 * mean(df$treatment[df$ltfu_3m == 0] == "ELD")),
         Lost = sprintf("%.1f", 100 * mean(df$treatment[df$ltfu_3m == 1] == "ELD")),
         p = sprintf("%.3f", tx_t$p.value)))

write.csv(ltfu, file.path(paths$tables, "table_ltfu_3m.csv"), row.names = FALSE)
print(as.data.frame(ltfu), row.names = FALSE)
cat(sprintf("\nFollowed %d, lost %d\n", sum(df$ltfu_3m == 0), sum(df$ltfu_3m == 1)))

cat("\nSupplementary analyses complete.\n")
