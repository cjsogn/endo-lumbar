# =============================================================================
# ENDO-LUMBAR: 16 Supplementary Frequentist Analysis (TMLE)
# Targeted Maximum Likelihood Estimation with SuperLearner
# =============================================================================
#
# TMLE with SuperLearner ensemble, MICE pooling, and IPCW for 12m outcomes
# Sign convention: positive ATE = ELD superior (matches Bayesian)

# =============================================================================
# SECTION 1: SETUP
# =============================================================================

cat("=============================================================================\n")
cat("ENDO-LUMBAR: Supplementary Frequentist Analysis (TMLE)\n")
cat("=============================================================================\n\n")

# Source config for covariates, NI margins, paths
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

# Load TMLE-specific packages
suppressPackageStartupMessages({
  library(SuperLearner)
  library(tmle)
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(mice)
  library(sandwich)
  library(lmtest)
  library(parallel)
})

set.seed(20260204L)

# Load data
cat("Loading data...\n")
df_full <- readRDS(file.path(paths$data_clean, "df_disc_imp.rds"))

# 12-month eligibility is defined by surgery date alone (see script 17). The
# earlier rule also admitted late-operated patients who happened to respond,
# which conditions the analysis set on the outcome being observed.
f12_datebased <- file.path(paths$data_clean, "df_disc_12m_eligible_datebased.rds")
df_12m <- if (file.exists(f12_datebased)) {
  readRDS(f12_datebased)
} else {
  readRDS(file.path(paths$data_clean, "df_disc_12m_eligible.rds"))
}

cat(sprintf("  Full sample: n=%d (ELD=%d, MSD=%d)\n",
            nrow(df_full), sum(df_full$treatment == "ELD"), sum(df_full$treatment == "MSD")))
cat(sprintf("  12m eligible: n=%d (ELD=%d, MSD=%d)\n",
            nrow(df_12m), sum(df_12m$treatment == "ELD"), sum(df_12m$treatment == "MSD")))

# Load Bayesian results for comparison
bayesian_primary <- read.csv(file.path(paths$tables, "table2_primary_results.csv"))
bayesian_t2      <- read.csv(file.path(paths$tables, "table3_tier2_effectiveness.csv"))

# Output paths
out_tables  <- paths$tables
out_figures <- paths$figures
out_results <- paths$results
for (d in c(out_tables, out_figures, out_results)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# =============================================================================
# SECTION 2: SUPERLEARNER LIBRARY
# =============================================================================

cat("\nConfiguring SuperLearner library...\n")

# Full ensemble for primary analysis
sl_lib <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.xgboost")
sl_cvcontrol <- list(V = 10L)

# Simpler library for MICE-pooled analysis (runtime consideration)
sl_lib_mice <- c("SL.glm", "SL.glmnet")
sl_cvcontrol_mice <- list(V = 5L)

cat(sprintf("  Primary learners: %s (V=%d)\n", paste(sl_lib, collapse = ", "), sl_cvcontrol$V))
cat(sprintf("  MICE learners: %s (V=%d)\n", paste(sl_lib_mice, collapse = ", "), sl_cvcontrol_mice$V))

# =============================================================================
# SECTION 3: COVARIATE MATRIX PREPARATION
# =============================================================================

cat("\nPreparing covariate matrices...\n")

#' Build design matrix from ps_covariates for SuperLearner compatibility
#' Converts factors to dummy variables, removes intercept
#' @param dat data frame containing all covariates
#' @return numeric matrix suitable for SuperLearner/tmle
build_X_matrix <- function(dat) {
  X_df <- dat[, ps_covariates, drop = FALSE]

  # Convert integer columns that should be numeric
  for (col in names(X_df)) {
    if (is.integer(X_df[[col]])) {
      X_df[[col]] <- as.numeric(X_df[[col]])
    }
  }

  # model.matrix handles factor -> dummy conversion
  # Remove intercept column
  mm <- model.matrix(~ . - 1, data = X_df)

  # Replace any remaining NAs with column medians (safety net)
  for (j in seq_len(ncol(mm))) {
    na_idx <- is.na(mm[, j])
    if (any(na_idx)) {
      mm[na_idx, j] <- median(mm[, j], na.rm = TRUE)
    }
  }

  as.data.frame(mm)
}

X_full <- build_X_matrix(df_full)
X_12m  <- build_X_matrix(df_12m)

cat(sprintf("  Full sample: %d observations, %d covariates\n", nrow(X_full), ncol(X_full)))
cat(sprintf("  12m eligible: %d observations, %d covariates\n", nrow(X_12m), ncol(X_12m)))

# =============================================================================
# SECTION 4: TMLE ESTIMATOR FUNCTION
# =============================================================================

#' Run TMLE for a single outcome
#'
#' @param Y outcome vector (may contain NAs; complete cases used)
#' @param A treatment vector (0=MSD, 1=ELD)
#' @param W covariate data frame
#' @param family "gaussian" for continuous, "binomial" for binary
#' @param lower_is_better if TRUE, positive ATE means ELD has lower (better) value
#' @param ni_margin NI margin (positive value)
#' @param sl_lib SuperLearner library
#' @param cvcontrol CV control list
#' @return list with ATE, SE, CI, p-value, NI conclusion
run_tmle <- function(Y, A, W, family = "gaussian",
                     lower_is_better = TRUE, ni_margin = NULL,
                     sl_lib = c("SL.glm", "SL.glmnet", "SL.ranger", "SL.xgboost"),
                     cvcontrol = list(V = 10L),
                     use_cvQinit = TRUE) {
  # cvQinit = TRUE is the tmle package default and is required for the
  # influence-curve variance to be valid when the SuperLearner library contains
  # data-adaptive learners (ranger, xgboost). With cvQinit = FALSE the initial
  # outcome regression is evaluated in-sample, the residuals entering the
  # influence curve are shrunk by overfitting, and the standard errors are
  # anti-conservative. This is still standard TMLE, not CV-TMLE: only the
  # initial Q is cross-fitted, while the targeting step and the variance
  # formula are unchanged.

  # Complete cases only
  cc <- complete.cases(Y)
  Y_cc <- Y[cc]
  A_cc <- A[cc]
  W_cc <- W[cc, , drop = FALSE]
  n_cc <- length(Y_cc)

  if (n_cc < 50) {
    warning(sprintf("Only %d complete cases. Results may be unstable.", n_cc))
  }

  # For binary outcomes, ensure Y is 0/1
  if (family == "binomial") {
    Y_cc <- as.numeric(Y_cc)
    if (!all(Y_cc %in% c(0, 1))) {
      stop("Binary outcome contains values other than 0 and 1")
    }
  }

  # Run TMLE
  # tmle computes E[Y|A=1] - E[Y|A=0] = E[Y|ELD] - E[Y|MSD]
  tmle_fit <- tryCatch({
    tmle::tmle(
      Y = Y_cc,
      A = A_cc,
      W = W_cc,
      family = family,
      Q.SL.library = sl_lib,
      g.SL.library = sl_lib,
      cvQinit = use_cvQinit,
      V.Q = cvcontrol$V,
      V.g = cvcontrol$V
    )
  }, error = function(e) {
    # Fallback to simpler library if SuperLearner fails
    warning(sprintf("Full SL library failed: %s\nFalling back to SL.glm + SL.glmnet", e$message))
    tmle::tmle(
      Y = Y_cc,
      A = A_cc,
      W = W_cc,
      family = family,
      Q.SL.library = c("SL.glm", "SL.glmnet"),
      g.SL.library = c("SL.glm", "SL.glmnet"),
      cvQinit = use_cvQinit,
      V.Q = cvcontrol$V,
      V.g = cvcontrol$V
    )
  })

  # Extract results
  if (family == "gaussian") {
    est <- tmle_fit$estimates$ATE
  } else {
    est <- tmle_fit$estimates$ATE
  }

  raw_ate <- est$psi
  raw_se  <- sqrt(est$var.psi)
  raw_ci  <- est$CI

  # Apply sign convention: positive = ELD superior
  # tmle computes E[Y|A=1] - E[Y|A=0] = ELD - MSD
  # For lower-is-better: we want MSD - ELD, so negate
  # For higher-is-better: ELD - MSD is already correct
  if (lower_is_better) {
    ate <- -raw_ate
    ci_lo <- -raw_ci[2]
    ci_hi <- -raw_ci[1]
  } else {
    ate <- raw_ate
    ci_lo <- raw_ci[1]
    ci_hi <- raw_ci[2]
  }
  se <- raw_se  # SE is symmetric

  # NI test: reject H0 (inferiority) if CI_lower > -margin
  ni_conclusion <- NA_character_
  p_ni <- NA_real_
  if (!is.null(ni_margin)) {
    # One-sided test: H0: ATE <= -margin vs H1: ATE > -margin
    # Test statistic: (ATE - (-margin)) / SE
    z_ni <- (ate - (-ni_margin)) / se
    p_ni <- 1 - pnorm(z_ni)  # one-sided p-value (small = reject H0)
    ni_conclusion <- ifelse(ci_lo > -ni_margin, "NI demonstrated", "NI not demonstrated")
  }

  list(
    ate = ate,
    se = se,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    p_ni = p_ni,
    ni_conclusion = ni_conclusion,
    ni_margin = ni_margin,
    n = n_cc,
    n_total = length(Y),
    family = family,
    tmle_fit = tmle_fit
  )
}

# =============================================================================
# SECTION 5: PROPENSITY SCORE DIAGNOSTICS
# =============================================================================

cat("\n=== Propensity Score Diagnostics ===\n")

# Fit SuperLearner PS model on full sample for diagnostics
A_full <- df_full$treatment_num
ps_sl <- SuperLearner(
  Y = A_full,
  X = X_full,
  family = binomial(),
  SL.library = sl_lib,
  cvControl = list(V = 10L)
)

ps_pred <- ps_sl$SL.predict[, 1]

cat(sprintf("  PS range: [%.4f, %.4f]\n", min(ps_pred), max(ps_pred)))
cat(sprintf("  PS < 0.01: %d observations\n", sum(ps_pred < 0.01)))
cat(sprintf("  PS > 0.99: %d observations\n", sum(ps_pred > 0.99)))
cat(sprintf("  PS < 0.05: %d observations\n", sum(ps_pred < 0.05)))
cat(sprintf("  PS > 0.95: %d observations\n", sum(ps_pred > 0.95)))

# SuperLearner coefficient weights
cat("\nSuperLearner PS model weights:\n")
print(round(coef(ps_sl), 4))

# --- PS Overlap Plot ---
cat("\nGenerating PS overlap plot...\n")
ps_df <- data.frame(
  ps = ps_pred,
  treatment = df_full$treatment
)

p_ps <- ggplot(ps_df, aes(x = ps, fill = treatment)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_rug(aes(color = treatment), alpha = 0.3, sides = "b") +
  scale_fill_manual(values = tx_fills, name = "Treatment") +
  scale_color_manual(values = tx_colors, name = "Treatment") +
  labs(
    title = "Propensity Score Overlap",
    subtitle = "SuperLearner ensemble (GLM, LASSO, Random Forest, XGBoost)",
    x = "Propensity Score (probability of ELD)",
    y = "Density"
  ) +
  theme_pub +
  theme(legend.position = "top")

ggsave(file.path(out_figures, "fig_ps_overlap_tmle.png"),
       p_ps, width = 7, height = 5, dpi = 300, bg = "white")
cat("  Saved: fig_ps_overlap_tmle.png\n")

# =============================================================================
# SECTION 6: ALL 17 OUTCOMES - COMPLETE-CASE TMLE
# =============================================================================

cat("\n=== Complete-Case TMLE for All 17 Outcomes ===\n")

# Define all outcomes with their specifications
outcome_specs <- data.frame(
  outcome = c(
    "odi_3m", "odi_12m",
    "nrs_back_3m", "nrs_back_12m",
    "nrs_leg_3m", "nrs_leg_12m",
    "eq5d_3m", "eq5d_12m",
    "responder_3m",
    "rtw_3m", "rtw_12m",
    "analgesic_3m", "analgesic_12m",
    "satisfied_3m", "satisfied_12m",
    "gpe_success_3m", "gpe_success_12m"
  ),
  label = c(
    "ODI 3 months", "ODI 12 months",
    "NRS back pain 3 months", "NRS back pain 12 months",
    "NRS leg pain 3 months", "NRS leg pain 12 months",
    "EQ-5D 3 months", "EQ-5D 12 months",
    "Responder 3 months",
    "Return to work 3 months", "Return to work 12 months",
    "Analgesic use 3 months", "Analgesic use 12 months",
    "Satisfaction 3 months", "Satisfaction 12 months",
    "GPE success 3 months", "GPE success 12 months"
  ),
  family = c(
    "gaussian", "gaussian",
    "gaussian", "gaussian",
    "gaussian", "gaussian",
    "gaussian", "gaussian",
    "binomial",
    "binomial", "binomial",
    "binomial", "binomial",
    "binomial", "binomial",
    "binomial", "binomial"
  ),
  lower_is_better = c(
    TRUE, TRUE,
    TRUE, TRUE,
    TRUE, TRUE,
    FALSE, FALSE,
    FALSE,
    FALSE, FALSE,
    TRUE, TRUE,
    FALSE, FALSE,
    FALSE, FALSE
  ),
  ni_margin = c(
    7, 7,
    1.0, 1.0,
    1.0, 1.0,
    0.05, 0.05,
    0.10,
    0.10, 0.10,
    0.10, 0.10,
    0.10, 0.10,
    0.10, 0.10
  ),
  timepoint = c(
    "3m", "12m",
    "3m", "12m",
    "3m", "12m",
    "3m", "12m",
    "3m",
    "3m", "12m",
    "3m", "12m",
    "3m", "12m",
    "3m", "12m"
  ),
  stringsAsFactors = FALSE
)

# Run TMLE for each outcome sequentially
tmle_results <- list()

for (i in seq_len(nrow(outcome_specs))) {
  spec <- outcome_specs[i, ]
  cat(sprintf("\n  [%d/%d] %s (%s)...\n", i, nrow(outcome_specs), spec$label, spec$family))

  # Select appropriate dataset
  if (spec$timepoint == "12m") {
    dat <- df_12m
    W <- X_12m
  } else {
    dat <- df_full
    W <- X_full
  }

  Y <- dat[[spec$outcome]]
  A <- dat$treatment_num

  # Run TMLE
  res <- tryCatch(
    run_tmle(
      Y = Y, A = A, W = W,
      family = spec$family,
      lower_is_better = spec$lower_is_better,
      ni_margin = spec$ni_margin,
      sl_lib = sl_lib,
      cvcontrol = sl_cvcontrol
    ),
    error = function(e) {
      cat(sprintf("    ERROR: %s\n", e$message))
      list(
        ate = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
        p_ni = NA_real_, ni_conclusion = "Error", ni_margin = spec$ni_margin,
        n = sum(!is.na(Y)), n_total = length(Y), family = spec$family, tmle_fit = NULL
      )
    }
  )

  cat(sprintf("    N=%d, ATE=%.3f, 95%% CI [%.3f, %.3f], NI: %s\n",
              res$n, res$ate, res$ci_lo, res$ci_hi, res$ni_conclusion))

  tmle_results[[spec$outcome]] <- c(res, list(
    outcome = spec$outcome,
    label = spec$label,
    lower_is_better = spec$lower_is_better,
    timepoint = spec$timepoint
  ))
}

# Compile complete-case results table
cc_table <- do.call(rbind, lapply(names(tmle_results), function(nm) {
  r <- tmle_results[[nm]]
  data.frame(
    Outcome = r$label,
    N = r$n,
    ATE = round(r$ate, 4),
    SE = round(r$se, 4),
    CI_lo = round(r$ci_lo, 4),
    CI_hi = round(r$ci_hi, 4),
    NI_Margin = r$ni_margin,
    NI_Conclusion = r$ni_conclusion,
    stringsAsFactors = FALSE
  )
}))

write.csv(cc_table, file.path(out_tables, "table_tmle_results.csv"), row.names = FALSE)
cat("\nSaved: table_tmle_results.csv\n")

# =============================================================================
# SECTION 7: MICE POOLING (RUBIN'S RULES)
# =============================================================================

cat("\n=== MICE Imputation + TMLE Pooling (Rubin's Rules) ===\n")

m_imp <- 20L
cat(sprintf("  Creating %d imputations for outcome missingness...\n", m_imp))

#' Run MICE for outcome imputation, then pool TMLE across imputations
#' @param dat dataset
#' @param W_matrix covariate matrix
#' @param outcome_col outcome variable name
#' @param family "gaussian" or "binomial"
#' @param lower_is_better sign convention
#' @param ni_margin NI margin
#' @param m number of imputations
#' @return list with pooled ATE, SE, CI, NI conclusion
run_mice_tmle <- function(dat, W_matrix, outcome_col, family,
                          lower_is_better, ni_margin, m = 20L) {

  Y <- dat[[outcome_col]]
  A <- dat$treatment_num
  n_miss <- sum(is.na(Y))

  if (n_miss == 0) {
    # No missing data; run single TMLE
    cat(sprintf("      No missing outcomes. Running single TMLE.\n"))
    res <- run_tmle(Y, A, W_matrix, family = family,
                    lower_is_better = lower_is_better,
                    ni_margin = ni_margin,
                    sl_lib = sl_lib_mice, cvcontrol = sl_cvcontrol_mice)
    return(list(
      ate_pooled = res$ate,
      se_pooled = res$se,
      ci_lo = res$ci_lo,
      ci_hi = res$ci_hi,
      ni_conclusion = res$ni_conclusion,
      m_imputations = 1,
      n = res$n,
      ates = res$ate,
      ses = res$se
    ))
  }

  cat(sprintf("      %d/%d outcomes missing (%.1f%%). Running MICE...\n",
              n_miss, length(Y), 100 * n_miss / length(Y)))

  # Build imputation dataset: outcome + treatment + covariates
  imp_data <- data.frame(Y = Y, A = A, W_matrix)

  # Set up MICE
  # Use predictive mean matching for continuous, logistic for binary
  meth <- rep("", ncol(imp_data))
  names(meth) <- names(imp_data)
  if (family == "gaussian") {
    meth["Y"] <- "pmm"  # predictive mean matching
  } else {
    meth["Y"] <- "logreg"  # logistic regression
  }
  # All other variables are complete (covariates already imputed, treatment always observed)

  # Run MICE quietly
  mice_obj <- tryCatch(
    mice::mice(imp_data, m = m, method = meth, maxit = 10,
               printFlag = FALSE, seed = 20260204L),
    error = function(e) {
      warning(sprintf("MICE failed: %s", e$message))
      return(NULL)
    }
  )

  if (is.null(mice_obj)) {
    # Fallback to complete-case
    cat("      MICE failed. Falling back to complete-case.\n")
    res <- run_tmle(Y, A, W_matrix, family = family,
                    lower_is_better = lower_is_better,
                    ni_margin = ni_margin,
                    sl_lib = sl_lib_mice, cvcontrol = sl_cvcontrol_mice)
    return(list(
      ate_pooled = res$ate, se_pooled = res$se,
      ci_lo = res$ci_lo, ci_hi = res$ci_hi,
      ni_conclusion = res$ni_conclusion,
      m_imputations = 0, n = res$n,
      ates = res$ate, ses = res$se
    ))
  }

  # Run TMLE on each imputed dataset and collect ATEs + SEs
  ates <- numeric(m)
  ses  <- numeric(m)

  for (k in seq_len(m)) {
    imp_k <- mice::complete(mice_obj, k)
    Y_k <- imp_k$Y
    A_k <- imp_k$A
    W_k <- imp_k[, -(1:2), drop = FALSE]  # everything except Y and A

    res_k <- tryCatch(
      run_tmle(Y_k, A_k, W_k, family = family,
               lower_is_better = lower_is_better,
               ni_margin = ni_margin,
               sl_lib = sl_lib_mice, cvcontrol = sl_cvcontrol_mice),
      error = function(e) {
        warning(sprintf("TMLE failed on imputation %d: %s", k, e$message))
        list(ate = NA_real_, se = NA_real_)
      }
    )

    ates[k] <- res_k$ate
    ses[k]  <- res_k$se
  }

  # Pool via Rubin's rules
  valid <- !is.na(ates) & !is.na(ses)
  m_valid <- sum(valid)

  if (m_valid < 2) {
    warning("Fewer than 2 valid imputations. Cannot pool.")
    return(list(
      ate_pooled = NA_real_, se_pooled = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_,
      ni_conclusion = "Pooling failed", m_imputations = m_valid,
      n = length(Y), ates = ates, ses = ses
    ))
  }

  ates_v <- ates[valid]
  ses_v  <- ses[valid]

  ate_pooled <- mean(ates_v)
  W_var <- mean(ses_v^2)                          # within-imputation variance
  B_var <- var(ates_v)                             # between-imputation variance
  T_var <- W_var + (1 + 1 / m_valid) * B_var      # total variance
  se_pooled <- sqrt(T_var)

  # Barnard-Rubin adjusted degrees of freedom
  r_m <- (1 + 1 / m_valid) * B_var / W_var        # relative increase in variance
  df_old <- (m_valid - 1) * (1 + 1 / r_m)^2       # original Rubin df

  # Barnard-Rubin adjustment for small samples
  n_obs <- sum(!is.na(Y))
  df_obs <- (n_obs - 1) * (1 + W_var / ((1 + 1 / m_valid) * B_var))  # observed-data df
  if (is.finite(df_old) && is.finite(df_obs) && df_old > 0 && df_obs > 0) {
    df_adj <- (df_old * df_obs) / (df_old + df_obs)
  } else {
    df_adj <- Inf  # use normal approximation
  }

  # CI using t-distribution with adjusted df
  if (is.finite(df_adj) && df_adj > 0) {
    t_crit <- qt(0.975, df_adj)
  } else {
    t_crit <- qnorm(0.975)
  }
  ci_lo <- ate_pooled - t_crit * se_pooled
  ci_hi <- ate_pooled + t_crit * se_pooled

  # NI test
  ni_conclusion <- ifelse(ci_lo > -ni_margin, "NI demonstrated", "NI not demonstrated")

  cat(sprintf("      Pooled: ATE=%.3f, SE=%.3f, 95%% CI [%.3f, %.3f], NI: %s\n",
              ate_pooled, se_pooled, ci_lo, ci_hi, ni_conclusion))

  list(
    ate_pooled = ate_pooled,
    se_pooled = se_pooled,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    ni_conclusion = ni_conclusion,
    m_imputations = m_valid,
    n = length(Y),
    df = df_adj,
    W_var = W_var,
    B_var = B_var,
    T_var = T_var,
    ates = ates,
    ses = ses
  )
}

# Run MICE-pooled TMLE for all 17 outcomes
mice_results <- list()

for (i in seq_len(nrow(outcome_specs))) {
  spec <- outcome_specs[i, ]
  cat(sprintf("\n  [%d/%d] %s (MICE + TMLE)...\n", i, nrow(outcome_specs), spec$label))

  # Select appropriate dataset
  if (spec$timepoint == "12m") {
    dat <- df_12m
    W <- X_12m
  } else {
    dat <- df_full
    W <- X_full
  }

  res <- tryCatch(
    run_mice_tmle(
      dat = dat, W_matrix = W,
      outcome_col = spec$outcome,
      family = spec$family,
      lower_is_better = spec$lower_is_better,
      ni_margin = spec$ni_margin,
      m = m_imp
    ),
    error = function(e) {
      cat(sprintf("    ERROR: %s\n", e$message))
      list(
        ate_pooled = NA_real_, se_pooled = NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_,
        ni_conclusion = "Error", m_imputations = 0,
        n = nrow(if (spec$timepoint == "12m") df_12m else df_full),
        ates = NA, ses = NA
      )
    }
  )

  mice_results[[spec$outcome]] <- c(res, list(
    outcome = spec$outcome,
    label = spec$label,
    lower_is_better = spec$lower_is_better,
    timepoint = spec$timepoint
  ))
}

# Compile MICE-pooled results table
mice_table <- do.call(rbind, lapply(names(mice_results), function(nm) {
  r <- mice_results[[nm]]
  data.frame(
    Outcome = r$label,
    ATE_pooled = round(r$ate_pooled, 4),
    SE_pooled = round(r$se_pooled, 4),
    CI_lo = round(r$ci_lo, 4),
    CI_hi = round(r$ci_hi, 4),
    NI_Conclusion = r$ni_conclusion,
    m_imputations = r$m_imputations,
    stringsAsFactors = FALSE
  )
}))

write.csv(mice_table, file.path(out_tables, "table_tmle_mice_pooled.csv"), row.names = FALSE)
cat("\nSaved: table_tmle_mice_pooled.csv\n")

# =============================================================================
# SECTION 8: IPCW FOR 12-MONTH OUTCOMES
# =============================================================================

cat("\n=== IPCW-Weighted TMLE for 12-Month Outcomes ===\n")

# IPCW model: P(observed at 12m | X, A)

R_12m <- as.integer(!is.na(df_12m$odi_12m))
A_12m <- df_12m$treatment_num

cat(sprintf("  12m observation rate: %.1f%% (ELD: %.1f%%, MSD: %.1f%%)\n",
            100 * mean(R_12m),
            100 * mean(R_12m[A_12m == 1]),
            100 * mean(R_12m[A_12m == 0])))

# Fit IPCW model (simpler library for speed)
ipcw_sl <- SuperLearner(
  Y = R_12m,
  X = data.frame(A = A_12m, X_12m),
  family = binomial(),
  SL.library = sl_lib_mice,
  cvControl = list(V = 5L)
)

p_obs <- ipcw_sl$SL.predict[, 1]

# Stabilized weights: P(R=1|A) / P(R=1|X,A)
p_obs_marginal <- tapply(R_12m, A_12m, mean)
p_obs_marg <- ifelse(A_12m == 1, p_obs_marginal["1"], p_obs_marginal["0"])
ipcw_weights <- p_obs_marg / pmax(p_obs, 0.01)  # floor at 0.01 to prevent extreme weights

# Trim at 99th percentile
trim_99 <- quantile(ipcw_weights, 0.99)
ipcw_weights <- pmin(ipcw_weights, trim_99)

cat(sprintf("  IPCW weight range: [%.2f, %.2f] (trimmed at %.2f)\n",
            min(ipcw_weights), max(ipcw_weights), trim_99))
cat(sprintf("  IPCW SL weights: %s\n",
            paste(names(coef(ipcw_sl)), round(coef(ipcw_sl), 3), sep = "=", collapse = ", ")))

# 12m outcomes for IPCW analysis
outcomes_12m <- outcome_specs[outcome_specs$timepoint == "12m", ]

ipcw_results <- list()

for (i in seq_len(nrow(outcomes_12m))) {
  spec <- outcomes_12m[i, ]
  cat(sprintf("\n  [%d/%d] %s (IPCW-TMLE)...\n", i, nrow(outcomes_12m), spec$label))

  Y <- df_12m[[spec$outcome]]
  A <- df_12m$treatment_num

  # Use tmle's obsWeights parameter for IPCW-weighted estimation
  # Restrict to observed cases and apply IPCW weights
  obs_idx <- !is.na(Y)
  Y_obs <- Y[obs_idx]
  A_obs <- A[obs_idx]
  W_obs <- X_12m[obs_idx, , drop = FALSE]
  wt_obs <- ipcw_weights[obs_idx]

  # Run IPCW-weighted TMLE using obsWeights parameter
  raw_fit <- tryCatch({
    tmle::tmle(
      Y = Y_obs,
      A = A_obs,
      W = W_obs,
      family = spec$family,
      Q.SL.library = sl_lib_mice,
      g.SL.library = sl_lib_mice,
      cvQinit = TRUE,
      V.Q = 5L,
      V.g = 5L,
      obsWeights = wt_obs
    )
  }, error = function(e) {
    # Fallback to simpler library
    cat(sprintf("    SL failed, falling back: %s\n", e$message))
    tryCatch(
      tmle::tmle(
        Y = Y_obs, A = A_obs, W = W_obs,
        family = spec$family,
        Q.SL.library = c("SL.glm"),
        g.SL.library = c("SL.glm"),
        cvQinit = TRUE,
        V.Q = 5L, V.g = 5L,
        obsWeights = wt_obs
      ),
      error = function(e2) NULL
    )
  })

  if (is.null(raw_fit)) {
    cat(sprintf("    ERROR: IPCW-TMLE failed completely\n"))
    ipcw_results[[spec$outcome]] <- list(
      outcome = spec$outcome, label = spec$label,
      ate = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      ni_conclusion = "Error", n_obs = sum(obs_idx), n_eligible = nrow(df_12m)
    )
    next
  }

  est <- raw_fit$estimates$ATE
  raw_ate <- est$psi
  raw_se  <- sqrt(est$var.psi)
  raw_ci  <- est$CI

  # Apply sign convention
  if (spec$lower_is_better) {
    ate_ipcw <- -raw_ate
    ci_lo_ipcw <- -raw_ci[2]
    ci_hi_ipcw <- -raw_ci[1]
  } else {
    ate_ipcw <- raw_ate
    ci_lo_ipcw <- raw_ci[1]
    ci_hi_ipcw <- raw_ci[2]
  }

  ni_ipcw <- ifelse(ci_lo_ipcw > -spec$ni_margin, "NI demonstrated", "NI not demonstrated")

  cat(sprintf("    N_obs=%d, ATE=%.3f, SE=%.3f, 95%% CI [%.3f, %.3f], NI: %s\n",
              sum(obs_idx), ate_ipcw, raw_se, ci_lo_ipcw, ci_hi_ipcw, ni_ipcw))

  ipcw_results[[spec$outcome]] <- list(
    outcome = spec$outcome,
    label = spec$label,
    ate = ate_ipcw,
    se = raw_se,
    ci_lo = ci_lo_ipcw,
    ci_hi = ci_hi_ipcw,
    ni_conclusion = ni_ipcw,
    n_obs = sum(obs_idx),
    n_eligible = nrow(df_12m)
  )
}

# =============================================================================
# SECTION 8b: TMLE FOR PERIOPERATIVE (TIER 3) OUTCOMES
# =============================================================================
# Prior-free confirmation of the gated superiority findings. Length of stay is
# handled as the mean number of postoperative nights rather than as an ordinal
# outcome, since TMLE targets a difference in means.

cat("\n=== TMLE for Perioperative (Tier 3) Outcomes ===\n")

tier3_specs <- list(
  list(outcome = "day_surgery",  label = "Day surgery",
       family = "binomial", lower_is_better = FALSE, ni_margin = 0.10),
  list(outcome = "los_postop",   label = "LOS postop (nights)",
       family = "gaussian", lower_is_better = TRUE,  ni_margin = NULL),
  list(outcome = "pt_comp_any_3m", label = "Patient-reported complications 3m",
       family = "binomial", lower_is_better = TRUE,  ni_margin = 0.10)
)

A_full <- as.integer(df_full$treatment == "ELD")
tier3_tmle <- list()

for (i in seq_along(tier3_specs)) {
  spec <- tier3_specs[[i]]
  if (!spec$outcome %in% names(df_full)) {
    cat(sprintf("  [%d/%d] %s: column not found, skipped\n",
                i, length(tier3_specs), spec$label))
    next
  }
  cat(sprintf("  [%d/%d] %s...\n", i, length(tier3_specs), spec$label))

  res <- tryCatch(
    run_tmle(as.numeric(df_full[[spec$outcome]]), A_full, X_full,
             family = spec$family,
             lower_is_better = spec$lower_is_better,
             ni_margin = spec$ni_margin,
             sl_lib = sl_lib, cvcontrol = sl_cvcontrol),
    error = function(e) {
      cat(sprintf("    ERROR: %s\n", e$message)); NULL
    }
  )
  if (is.null(res)) next

  cat(sprintf("    N=%d, ATE=%.4f, SE=%.4f, 95%% CI [%.4f, %.4f]\n",
              res$n, res$ate, res$se, res$ci_lo, res$ci_hi))

  tier3_tmle[[spec$outcome]] <- tibble(
    Outcome = spec$label, N = res$n, ATE = res$ate, SE = res$se,
    CI_lo = res$ci_lo, CI_hi = res$ci_hi, Family = spec$family
  )
}

if (length(tier3_tmle) > 0) {
  tier3_tmle_tab <- bind_rows(tier3_tmle)
  write.csv(tier3_tmle_tab, file.path(out_tables, "table_tmle_tier3.csv"),
            row.names = FALSE)
  cat("  Saved: table_tmle_tier3.csv\n")
}

# =============================================================================
# SECTION 9: OUTPUT TABLES
# =============================================================================

cat("\n=== Generating Output Tables ===\n")

# --- Table: Bayesian vs TMLE comparison ---
# Build comparison from Bayesian table3 and TMLE complete-case results

# Parse Bayesian results into a lookup
bayesian_lookup <- bayesian_t2 %>%
  select(Outcome, ATE, CrI_lo, CrI_hi, P_NI, NI_Conclusion)

# Add primary result
# Read brms primary result from table2
brms_primary <- read.csv(file.path(out_tables, "table2_primary_results.csv"),
                         stringsAsFactors = FALSE)
bayesian_primary_row <- data.frame(
  Outcome = "ODI 3 months",
  ATE = as.numeric(brms_primary$ATE[1]),
  CrI_lo = as.numeric(gsub("\\[|\\]", "", strsplit(brms_primary$X95..CrI[1], ",")[[1]][1])),
  CrI_hi = as.numeric(gsub("\\[|\\]", "", strsplit(brms_primary$X95..CrI[1], ",")[[1]][2])),
  P_NI = 1.000,
  NI_Conclusion = "NI demonstrated",
  stringsAsFactors = FALSE
)
bayesian_lookup <- rbind(bayesian_primary_row, bayesian_lookup)

# Normalised key for label matching. The Bayesian table carries parenthetical
# definitions (for example "Responder 3 months (>=30% or >=10pt)") that the
# earlier substring match failed on, which left that row blank in the
# comparison table.
norm_label <- function(x) {
  x <- sub("\\s*\\(.*\\)\\s*$", "", x)
  tolower(trimws(gsub("\\s+", " ", x)))
}
bayesian_lookup$key <- norm_label(bayesian_lookup$Outcome)

# Create comparison table
comparison_rows <- list()
for (i in seq_len(nrow(outcome_specs))) {
  spec <- outcome_specs[i, ]
  tmle_r <- tmle_results[[spec$outcome]]

  # Find matching Bayesian result
  bay_match <- bayesian_lookup[bayesian_lookup$key == norm_label(spec$label), ]

  bay_ate <- if (nrow(bay_match) > 0) bay_match$ATE[1] else NA_real_
  bay_lo  <- if (nrow(bay_match) > 0) bay_match$CrI_lo[1] else NA_real_
  bay_hi  <- if (nrow(bay_match) > 0) bay_match$CrI_hi[1] else NA_real_
  bay_pni <- if (nrow(bay_match) > 0) bay_match$P_NI[1] else NA_real_

  # Concordance: both demonstrate NI or both don't
  tmle_ni <- tmle_r$ni_conclusion
  bay_ni  <- if (nrow(bay_match) > 0) bay_match$NI_Conclusion[1] else NA_character_
  concordance <- ifelse(!is.na(tmle_ni) & !is.na(bay_ni),
                        ifelse(grepl("demonstrated", tmle_ni) == grepl("demonstrated", bay_ni),
                               "Concordant", "Discordant"),
                        NA_character_)

  comparison_rows[[spec$outcome]] <- data.frame(
    Outcome = spec$label,
    Bayesian_ATE = round(bay_ate, 3),
    Bayesian_CrI_lo = round(bay_lo, 3),
    Bayesian_CrI_hi = round(bay_hi, 3),
    Bayesian_P_NI = round(bay_pni, 3),
    TMLE_ATE = round(tmle_r$ate, 3),
    TMLE_CI_lo = round(tmle_r$ci_lo, 3),
    TMLE_CI_hi = round(tmle_r$ci_hi, 3),
    TMLE_NI = tmle_ni,
    Concordance = concordance,
    stringsAsFactors = FALSE
  )
}

comparison_table <- do.call(rbind, comparison_rows)
write.csv(comparison_table, file.path(out_tables, "table_bayesian_vs_tmle.csv"), row.names = FALSE)
cat("  Saved: table_bayesian_vs_tmle.csv\n")

# --- IPCW table ---
ipcw_table <- do.call(rbind, lapply(names(ipcw_results), function(nm) {
  r <- ipcw_results[[nm]]
  data.frame(
    Outcome = r$label,
    N_observed = r$n_obs,
    N_eligible = r$n_eligible,
    ATE = round(r$ate, 4),
    SE = round(r$se, 4),
    CI_lo = round(r$ci_lo, 4),
    CI_hi = round(r$ci_hi, 4),
    NI_Conclusion = r$ni_conclusion,
    stringsAsFactors = FALSE
  )
}))
write.csv(ipcw_table, file.path(out_tables, "table_tmle_ipcw_12m.csv"), row.names = FALSE)
cat("  Saved: table_tmle_ipcw_12m.csv\n")

# =============================================================================
# SECTION 10: FIGURES
# =============================================================================

cat("\n=== Generating Figures ===\n")

# --- Figure: Bayesian vs TMLE Forest Plot ---
cat("  Creating Bayesian vs TMLE forest plot...\n")

# Prepare forest plot data
forest_data <- list()

for (i in seq_len(nrow(outcome_specs))) {
  spec <- outcome_specs[i, ]
  tmle_r <- tmle_results[[spec$outcome]]

  # Bayesian match
  bay_match <- bayesian_lookup[bayesian_lookup$Outcome == spec$label, ]
  if (nrow(bay_match) == 0) {
    bay_match <- bayesian_lookup[grepl(gsub(" \\d+m.*", "", spec$label),
                                       bayesian_lookup$Outcome, ignore.case = TRUE) &
                                 grepl(spec$timepoint, bayesian_lookup$Outcome), ]
  }

  if (nrow(bay_match) > 0) {
    forest_data[[length(forest_data) + 1]] <- data.frame(
      outcome = spec$label,
      method = "Bayesian",
      ate = bay_match$ATE[1],
      ci_lo = bay_match$CrI_lo[1],
      ci_hi = bay_match$CrI_hi[1],
      family = spec$family,
      ni_margin = spec$ni_margin,
      stringsAsFactors = FALSE
    )
  }

  forest_data[[length(forest_data) + 1]] <- data.frame(
    outcome = spec$label,
    method = "TMLE",
    ate = tmle_r$ate,
    ci_lo = tmle_r$ci_lo,
    ci_hi = tmle_r$ci_hi,
    family = spec$family,
    ni_margin = spec$ni_margin,
    stringsAsFactors = FALSE
  )
}

forest_df <- do.call(rbind, forest_data)

# Order outcomes
outcome_order <- rev(outcome_specs$label)
forest_df$outcome <- factor(forest_df$outcome, levels = outcome_order)
forest_df$method <- factor(forest_df$method, levels = c("Bayesian", "TMLE"))

# Split into continuous and binary panels
forest_cont <- forest_df[forest_df$family == "gaussian", ]
forest_bin  <- forest_df[forest_df$family == "binomial", ]

# --- Panel A: Continuous outcomes ---
p_cont <- ggplot(forest_cont, aes(x = ate, y = outcome, color = method, shape = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi),
                  position = position_dodge(width = 0.5),
                  linewidth = 0.5, size = 0.4) +
  scale_color_manual(values = c("Bayesian" = "#2166AC", "TMLE" = "#D95F02"), name = "Method") +
  scale_shape_manual(values = c("Bayesian" = 16, "TMLE" = 17), name = "Method") +
  labs(
    title = "A. Continuous Outcomes",
    x = "Average Treatment Effect (positive = ELD superior)",
    y = NULL
  ) +
  theme_pub +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 11, face = "bold")
  )

# --- Panel B: Binary outcomes ---
p_bin <- ggplot(forest_bin, aes(x = ate, y = outcome, color = method, shape = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = -0.10, linetype = "dotted", color = "red", linewidth = 0.5) +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi),
                  position = position_dodge(width = 0.5),
                  linewidth = 0.5, size = 0.4) +
  scale_color_manual(values = c("Bayesian" = "#2166AC", "TMLE" = "#D95F02"), name = "Method") +
  scale_shape_manual(values = c("Bayesian" = 16, "TMLE" = 17), name = "Method") +
  labs(
    title = "B. Binary Outcomes (Risk Difference)",
    x = "Average Treatment Effect (positive = ELD superior)",
    y = NULL
  ) +
  annotate("text", x = -0.10, y = 0.5, label = "NI margin", color = "red",
           size = 3, hjust = 1.1, fontface = "italic") +
  theme_pub +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 11, face = "bold")
  )

# Combine panels
p_forest <- p_cont / p_bin +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Bayesian vs TMLE: Treatment Effect Estimates",
    subtitle = "Disc herniation population (ELD vs MSD)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

ggsave(file.path(out_figures, "fig_bayesian_vs_tmle_forest.png"),
       p_forest, width = 10, height = 12, dpi = 300, bg = "white")
cat("  Saved: fig_bayesian_vs_tmle_forest.png\n")

# =============================================================================
# SECTION 11: SAVE ALL RESULTS
# =============================================================================

cat("\n=== Saving Complete Results ===\n")

tmle_all_results <- list(
  tmle_complete_case = tmle_results,
  tmle_mice_pooled = mice_results,
  tmle_ipcw_12m = ipcw_results,
  outcome_specs = outcome_specs,
  ps_model = list(
    sl_weights = coef(ps_sl),
    ps_range = range(ps_pred),
    n_extreme_low = sum(ps_pred < 0.01),
    n_extreme_high = sum(ps_pred > 0.99)
  ),
  ipcw_model = list(
    sl_weights = coef(ipcw_sl),
    weight_range = range(ipcw_weights),
    weight_trim = trim_99,
    obs_rate_overall = mean(R_12m),
    obs_rate_eld = mean(R_12m[A_12m == 1]),
    obs_rate_msd = mean(R_12m[A_12m == 0])
  ),
  settings = list(
    sl_library_primary = sl_lib,
    sl_library_mice = sl_lib_mice,
    cv_folds_primary = sl_cvcontrol$V,
    cv_folds_mice = sl_cvcontrol_mice$V,
    mice_m = m_imp,
    seed = 20260204L,
    n_full = nrow(df_full),
    n_12m = nrow(df_12m)
  )
)

saveRDS(tmle_all_results, file.path(out_results, "tmle_results.rds"))
cat("  Saved: tmle_results.rds\n")

# Also copy tables and figures to the old output paths for compatibility.
# Skip when the source and destination are the same file.
copy_if_different <- function(from, to_dir) {
  to <- file.path(to_dir, basename(from))
  if (normalizePath(from, mustWork = FALSE) != normalizePath(to, mustWork = FALSE)) {
    file.copy(from, to, overwrite = TRUE)
  }
}
for (f in list.files(out_tables, pattern = "tmle", full.names = TRUE)) {
  copy_if_different(f, paths$tables)
}
for (f in list.files(out_figures, pattern = "tmle|ps_overlap", full.names = TRUE)) {
  copy_if_different(f, paths$figures)
}

cat("\n=============================================================================\n")
cat("TMLE analysis complete.\n")
cat("=============================================================================\n")
cat(sprintf("\nOutput files:\n"))
cat(sprintf("  Tables:  %s\n", out_tables))
cat(sprintf("    - table_tmle_results.csv (complete-case)\n"))
cat(sprintf("    - table_tmle_mice_pooled.csv (MICE-pooled, m=%d)\n", m_imp))
cat(sprintf("    - table_bayesian_vs_tmle.csv (comparison)\n"))
cat(sprintf("    - table_tmle_ipcw_12m.csv (IPCW sensitivity)\n"))
cat(sprintf("  Figures: %s\n", out_figures))
cat(sprintf("    - fig_ps_overlap_tmle.png\n"))
cat(sprintf("    - fig_bayesian_vs_tmle_forest.png\n"))
cat(sprintf("  Results: %s\n", out_results))
cat(sprintf("    - tmle_results.rds\n"))
