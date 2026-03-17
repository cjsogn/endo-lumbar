# =============================================================================
# ENDO-LUMBAR: 17 ELD Approach Comparison (Transforaminal vs Interlaminar)
# Exploratory analysis comparing interlaminar vs transforaminal endoscopic
# approaches within ELD patients
# =============================================================================

source("/Users/cjsogn/ENDO_LUMBAR_COMPLETE_CASE/scripts/00_config_cc.R")

df <- readRDS(file.path(paths$data_clean, "df_all.rds"))

cat("=== ELD Approach Comparison: Interlaminar vs Transforaminal ===\n")

# --- Subset to ELD patients with Midline or Wiltse approach ------------------
# Data uses "Midline" (interlaminar) and "Wiltse" (transforaminal) internally.
# Keep internal labels for model fitting (cached brms models use approachWiltse),
# then relabel for all display output.

df_eld <- df %>%
  filter(treatment == "ELD", approach %in% c("Midline", "Wiltse")) %>%
  mutate(approach = factor(approach, levels = c("Midline", "Wiltse")))

# Display label mapping
approach_display <- c("Midline" = "Interlaminar", "Wiltse" = "Transforaminal")

n_interlam <- sum(df_eld$approach == "Midline")
n_transfora <- sum(df_eld$approach == "Wiltse")

cat(sprintf("ELD patients: %d (Interlaminar=%d, Transforaminal=%d)\n",
            nrow(df_eld), n_interlam, n_transfora))

# =============================================================================
# 1. DESCRIPTIVE TABLE: Baseline characteristics by approach
# =============================================================================

cat("\n--- Baseline characteristics by approach ---\n")

baseline_vars <- c("age", "sex", "bmi", "odi_baseline", "nrs_back_baseline",
                    "nrs_leg_baseline", "eq5d_baseline", "prior_surgery_any",
                    "depression_anxiety", "smoking", "motor_deficit",
                    "symptom_duration_leg", "asa_cat", "multilevel",
                    "day_surgery", "operating_time")

# Relabel for table display
df_eld_display <- df_eld %>%
  mutate(approach = recode(approach, "Midline" = "Interlaminar", "Wiltse" = "Transforaminal"))

tbl1 <- CreateTableOne(
  vars = baseline_vars,
  strata = "approach",
  data = df_eld_display,
  test = FALSE
)
tbl1_print <- print(tbl1, printToggle = FALSE, noSpaces = TRUE, showAllLevels = TRUE)
cat("\n")
print(tbl1)

write.csv(tbl1_print,
          file.path(paths$tables, "table_eld_approach_baseline.csv"))

# =============================================================================
# 2. UNADJUSTED OUTCOME COMPARISON
# =============================================================================

cat("\n--- Unadjusted outcomes by approach ---\n")

unadj_summary <- df_eld_display %>%
  group_by(approach) %>%
  dplyr::summarise(
    n = n(),
    odi_3m_n       = sum(!is.na(odi_3m)),
    odi_3m_mean    = mean(odi_3m, na.rm = TRUE),
    odi_3m_sd      = sd(odi_3m, na.rm = TRUE),
    nrs_back_3m_n  = sum(!is.na(nrs_back_3m)),
    nrs_back_3m_mean = mean(nrs_back_3m, na.rm = TRUE),
    nrs_back_3m_sd = sd(nrs_back_3m, na.rm = TRUE),
    nrs_leg_3m_n   = sum(!is.na(nrs_leg_3m)),
    nrs_leg_3m_mean = mean(nrs_leg_3m, na.rm = TRUE),
    nrs_leg_3m_sd  = sd(nrs_leg_3m, na.rm = TRUE),
    eq5d_3m_n      = sum(!is.na(eq5d_3m)),
    eq5d_3m_mean   = mean(eq5d_3m, na.rm = TRUE),
    eq5d_3m_sd     = sd(eq5d_3m, na.rm = TRUE),
    op_time_mean   = mean(operating_time, na.rm = TRUE),
    op_time_sd     = sd(operating_time, na.rm = TRUE),
    day_surg_pct   = 100 * mean(day_surgery == 1, na.rm = TRUE),
    perop_comp_pct = 100 * mean(perop_comp_any == 1, na.rm = TRUE),
    pt_comp_3m_pct = 100 * mean(pt_comp_any_3m == 1, na.rm = TRUE),
    .groups = "drop"
  )
print(unadj_summary)

write.csv(unadj_summary,
          file.path(paths$tables, "table_eld_approach_unadjusted.csv"),
          row.names = FALSE)

# =============================================================================
# 3. ADJUSTED BAYESIAN MODELS
# =============================================================================

# With n=29 in the transforaminal group, use a restricted set of ~6 covariates
# to avoid overfitting. Wider treatment prior: Normal(0, 15).

cat("\n--- Adjusted Bayesian models ---\n")

# Prepare analysis data (keep internal Midline/Wiltse for brms coefficient names)
df_model <- df_eld %>%
  mutate(
    approach_num = as.numeric(approach == "Wiltse"),
    sex_num = as.numeric(sex)
  )

# --- Priors for small-sample comparison ---
approach_priors <- c(
  set_prior("normal(0, 15)", class = "b", coef = "approachWiltse"),
  set_prior("normal(0, 2)", class = "b"),
  set_prior("student_t(3, 0, 15)", class = "sigma")
)

approach_priors_binary <- c(
  set_prior("normal(0, 1.5)", class = "b", coef = "approachWiltse"),
  set_prior("normal(0, 2)", class = "b")
)

# --- Model formulas ---
formula_continuous <- function(outcome) {
  bf(as.formula(paste(outcome, "~ approach + age + sex + bmi + odi_baseline + nrs_leg_baseline + prior_surgery_any")))
}

# --- Fit continuous outcome models ---
continuous_outcomes <- list(
  list(name = "odi_3m", label = "ODI 3m", lower_better = TRUE),
  list(name = "nrs_back_3m", label = "NRS Back 3m", lower_better = TRUE),
  list(name = "nrs_leg_3m", label = "NRS Leg 3m", lower_better = TRUE),
  list(name = "eq5d_3m", label = "EQ-5D 3m", lower_better = FALSE),
  list(name = "operating_time", label = "Operating Time", lower_better = TRUE)
)

results_list <- list()

for (oc in continuous_outcomes) {
  cat(sprintf("\nFitting model: %s\n", oc$label))

  # Filter to non-missing outcome
  df_fit <- df_model %>% filter(!is.na(.data[[oc$name]]))
  cat(sprintf("  N=%d (Interlaminar=%d, Transforaminal=%d)\n",
              nrow(df_fit),
              sum(df_fit$approach == "Midline"),
              sum(df_fit$approach == "Wiltse")))

  fit <- brm(
    formula = formula_continuous(oc$name),
    data = df_fit,
    prior = approach_priors,
    chains = mcmc_settings$chains,
    iter = mcmc_settings$iter,
    warmup = mcmc_settings$warmup,
    cores = n_cores,
    seed = mcmc_settings$seed,
    control = list(adapt_delta = mcmc_settings$adapt_delta,
                   max_treedepth = mcmc_settings$max_treedepth),
    file = file.path(paths$models, paste0("approach_", oc$name)),
    silent = 2, refresh = 0
  )

  # Extract approach effect (brms labels it approachWiltse)
  post <- as_draws_df(fit)
  approach_effect <- post$b_approachWiltse

  if (oc$lower_better) {
    delta <- -approach_effect
  } else {
    delta <- approach_effect
  }

  res <- tibble(
    outcome = oc$label,
    n_interlaminar = sum(df_fit$approach == "Midline"),
    n_transforaminal = sum(df_fit$approach == "Wiltse"),
    transforaminal_effect_mean = mean(approach_effect),
    transforaminal_effect_median = median(approach_effect),
    cri_lo = quantile(approach_effect, 0.025),
    cri_hi = quantile(approach_effect, 0.975),
    p_transforaminal_better = if (oc$lower_better) mean(approach_effect < 0) else mean(approach_effect > 0),
    rhat_max = max(brms::rhat(fit), na.rm = TRUE)
  )
  results_list[[oc$name]] <- res

  cat(sprintf("  Transforaminal effect: %.2f [%.2f, %.2f], P(TF better)=%.3f, Rhat=%.4f\n",
              res$transforaminal_effect_mean, res$cri_lo, res$cri_hi,
              res$p_transforaminal_better, res$rhat_max))
}

# --- Fit binary outcome model: day surgery ---
cat("\nFitting model: Day Surgery (binary)\n")

df_fit_ds <- df_model %>% filter(!is.na(day_surgery))

fit_ds <- brm(
  formula = bf(day_surgery ~ approach + age + sex + bmi + odi_baseline +
                 nrs_leg_baseline + prior_surgery_any,
               family = bernoulli()),
  data = df_fit_ds,
  prior = approach_priors_binary,
  chains = mcmc_settings$chains,
  iter = mcmc_settings$iter,
  warmup = mcmc_settings$warmup,
  cores = n_cores,
  seed = mcmc_settings$seed,
  control = list(adapt_delta = mcmc_settings$adapt_delta,
                 max_treedepth = mcmc_settings$max_treedepth),
  file = file.path(paths$models, "approach_day_surgery"),
  silent = 2, refresh = 0
)

post_ds <- as_draws_df(fit_ds)
ds_effect <- post_ds$b_approachWiltse

res_ds <- tibble(
  outcome = "Day Surgery",
  n_interlaminar = sum(df_fit_ds$approach == "Midline"),
  n_transforaminal = sum(df_fit_ds$approach == "Wiltse"),
  transforaminal_effect_mean = mean(ds_effect),
  transforaminal_effect_median = median(ds_effect),
  cri_lo = quantile(ds_effect, 0.025),
  cri_hi = quantile(ds_effect, 0.975),
  p_transforaminal_better = mean(ds_effect > 0),
  rhat_max = max(brms::rhat(fit_ds), na.rm = TRUE)
)
results_list[["day_surgery"]] <- res_ds

cat(sprintf("  Transforaminal effect (log-odds): %.2f [%.2f, %.2f], P(TF better)=%.3f\n",
            res_ds$transforaminal_effect_mean, res_ds$cri_lo, res_ds$cri_hi,
            res_ds$p_transforaminal_better))

# --- Combine results ---
results_df <- bind_rows(results_list)
print(results_df)

write.csv(results_df,
          file.path(paths$tables, "table_eld_approach_adjusted.csv"),
          row.names = FALSE)

# =============================================================================
# 4. SUMMARY FIGURE: Grouped bar chart
# =============================================================================

cat("\n--- Creating approach comparison figure ---\n")

plot_data <- df_eld_display %>%
  select(approach, odi_3m, nrs_back_3m, nrs_leg_3m, eq5d_3m, operating_time) %>%
  pivot_longer(-approach, names_to = "outcome", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(approach, outcome) %>%
  dplyr::summarise(
    mean = mean(value),
    sd = sd(value),
    n = n(),
    se = sd / sqrt(n),
    ci_lo = mean - 1.96 * se,
    ci_hi = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(
    outcome_label = case_when(
      outcome == "odi_3m" ~ "ODI\n3 months",
      outcome == "nrs_back_3m" ~ "NRS Back\n3 months",
      outcome == "nrs_leg_3m" ~ "NRS Leg\n3 months",
      outcome == "eq5d_3m" ~ "EQ-5D\n3 months",
      outcome == "operating_time" ~ "Op. Time\n(min)"
    ),
    outcome_label = factor(outcome_label,
                           levels = c("ODI\n3 months", "NRS Back\n3 months",
                                      "NRS Leg\n3 months", "EQ-5D\n3 months",
                                      "Op. Time\n(min)"))
  )

approach_colors <- c("Interlaminar" = "#4393C3", "Transforaminal" = "#D6604D")

p_approach <- ggplot(plot_data, aes(x = approach, y = mean, fill = approach)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.15, linewidth = 0.5, color = "grey20") +
  facet_wrap(~ outcome_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = approach_colors) +
  labs(
    title = "ELD Approach Comparison: Interlaminar vs Transforaminal",
    subtitle = sprintf("Unadjusted means with 95%% CIs (Interlaminar n=%d, Transforaminal n=%d)",
                       n_interlam, n_transfora),
    x = NULL,
    y = "Mean (95% CI)",
    fill = "Approach"
  ) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 11)
  )

save_fig(p_approach, "eld_approach_comparison.png", width = 10, height = 5)
cat("  Saved: eld_approach_comparison.png\n")

# =============================================================================
# 5. NOTE ON CONFOUNDING WITH LEARNING CURVE
# =============================================================================

cat("\n--- Learning curve confounding check ---\n")

df_eld_ordered <- df_eld %>% arrange(surgery_date) %>% mutate(case_num = row_number())

early <- df_eld_ordered %>% filter(case_num <= nrow(df_eld_ordered)/2)
late <- df_eld_ordered %>% filter(case_num > nrow(df_eld_ordered)/2)

cat(sprintf("  Early cases (1-%d): Interlaminar=%d, Transforaminal=%d\n",
            nrow(early),
            sum(early$approach == "Midline"),
            sum(early$approach == "Wiltse")))
cat(sprintf("  Late cases (%d-%d): Interlaminar=%d, Transforaminal=%d\n",
            nrow(early) + 1, nrow(df_eld_ordered),
            sum(late$approach == "Midline"),
            sum(late$approach == "Wiltse")))
cat("  Note: Transforaminal approach may be more concentrated in early cases,\n")
cat("  partially confounding approach with learning curve effects.\n")

# =============================================================================
# DONE
# =============================================================================

cat("\n=== ELD Approach Comparison Complete ===\n")
cat("  Tables saved to:", paths$tables, "\n")
cat("  Figure saved to:", paths$figures, "\n")
