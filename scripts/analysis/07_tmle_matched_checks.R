# Purpose: Run observed-case, outcome-imputed and censoring-weighted TMLE.
# Inputs: Prepared three- and twelve-month cohorts and shared covariate definitions.
# Outputs: TMLE estimates, pooled imputation results, fits and learner/weight checks.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
suppressPackageStartupMessages({
  library(SuperLearner)
  library(tmle)
  library(mice)
})
options(warn = 1) # Keep diagnostic warnings in each endpoint log.
PS <- sub("_z$", "", COVS)
# These wrappers set CPU threads for the native forest and boosting learners.
# Learner hyperparameters are otherwise unchanged.
SL.ranger.concurrent <- function(...) SuperLearner::SL.ranger(...,
  num.threads = as.integer(Sys.getenv("ENDO_CHAIN_CORES", "1"))
)
predict.SL.ranger.concurrent <- SuperLearner:::predict.SL.ranger
SL.xgboost.concurrent <- function(...) SuperLearner::SL.xgboost(...,
  nthread = as.integer(Sys.getenv("ENDO_CHAIN_CORES", "1"))
)
predict.SL.xgboost.concurrent <- SuperLearner:::predict.SL.xgboost
full_lib <- c("SL.glm", "SL.glmnet", "SL.ranger.concurrent", "SL.xgboost.concurrent")
small_lib <- c("SL.glm", "SL.glmnet")
effectiveness <- c(
  "odi_3m", "odi_12m", "nrs_back_3m", "nrs_back_12m",
  "nrs_leg_3m", "nrs_leg_12m", "eq5d_3m", "eq5d_12m", "responder_3m",
  "rtw_3m", "rtw_12m", "analgesic_3m", "analgesic_12m", "satisfied_3m",
  "satisfied_12m", "gpe_success_3m", "gpe_success_12m"
)
continuous <- c(
  "odi_3m", "odi_12m", "nrs_back_3m", "nrs_back_12m",
  "nrs_leg_3m", "nrs_leg_12m", "eq5d_3m", "eq5d_12m", "los_postop"
)
lower_outcomes <- c(
  setdiff(continuous, c("eq5d_3m", "eq5d_12m")), "analgesic_3m",
  "analgesic_12m", "pt_comp_any_3m"
)
ids <- c(
  paste0(c(effectiveness, "day_surgery", "los_postop", "pt_comp_any_3m"), "_calendar"),
  paste0(c(effectiveness, "day_surgery"), "_mi"), paste0(grep("12m$", effectiveness, value = TRUE), "_ipcw"),
  "day_surgery_original"
)
stopifnot(length(ids) == 47L, !anyDuplicated(ids))
extract_result <- function(fit, lower) {
  s <- fit$estimates$ATE
  sign <- if (lower) -1 else 1
  data.frame(
    mean = sign * s$psi, se = sqrt(s$var.psi),
    lower = if (lower) -s$CI[2] else s$CI[1], upper = if (lower) -s$CI[1] else s$CI[2]
  )
}
run_check <- function(id) {
  con <- file(file.path(ROOT, "10_logs", paste0("tmle_", id, ".log")), "wt")
  sink(con)
  sink(con, type = "message")
  on.exit(
    {
      sink(type = "message")
      sink()
      close(con)
    },
    add = TRUE
  )
  set.seed(SEED)
  y <- sub("_(calendar|original|mi|ipcw)$", "", id)
  d <- readRDS(file.path(ROOT, "02_data/derived", if (grepl("12m", y))
    "cohort_12m_date_eligible.rds" else "cohort_3m.rds"))
  vars <- c(PS, if (!grepl("original$", id)) c("ct1", "ct2"))
  W <- as.data.frame(model.matrix(~ . - 1, data = d[vars]))
  stopifnot(nrow(W) == nrow(d), !anyNA(W))
  Y <- as.numeric(haven::zap_labels(d[[y]]))
  A <- as.integer(d$treatment == "ELD")
  cc <- !is.na(Y)
  family <- if (y %in% continuous) "gaussian" else "binomial"
  lower <- y %in% lower_outcomes
  stopifnot(any(cc), all(is.finite(Y[cc])))
  if (family == "binomial") stopifnot(all(Y[cc] %in% 0:1))
  if (grepl("_mi$", id)) {
    impdat <- data.frame(Y = if (family == "binomial") factor(Y, levels = 0:1) else Y, A = A, W)
    meth <- setNames(rep("", ncol(impdat)), names(impdat))
    meth["Y"] <- if (family == "binomial") "logreg" else "pmm"
    imp <- mice(impdat, m = 20, method = meth, maxit = 10, printFlag = FALSE, seed = SEED)
    saveRDS(imp, file.path(ROOT, "02_data/derived", paste0(id, "_mice_m20.rds")))
    all <- vector("list", imp$m)
    for (k in seq_len(imp$m)) {
      dk <- mice::complete(imp, k)
      yy <- if (is.factor(dk$Y)) as.numeric(as.character(dk$Y)) else dk$Y
      stopifnot(!anyNA(yy))
      fit <- tmle::tmle(
        Y = yy, A = dk$A, W = dk[, -c(1, 2), drop = FALSE], family = family,
        Q.SL.library = small_lib, g.SL.library = small_lib, cvQinit = TRUE, V.Q = 5, V.g = 5
      )
      stopifnot(
        length(fit$Qinit$SL.library) == length(small_lib),
        length(fit$g$SL.library) == length(small_lib)
      )
      all[[k]] <- extract_result(fit, lower)
    }
    est <- do.call(rbind, all)
    stopifnot(nrow(est) == 20, all(is.finite(est$mean)), all(is.finite(est$se)))
    pooled <- mice::pool.scalar(est$mean, est$se^2, n = nrow(d), k = ncol(W) + 2)
    saveRDS(list(estimates = est, pool = pooled), file.path(ROOT, "04_results", paste0("tmle_", id, ".rds")))
    s <- data.frame(
      mean = pooled$qbar, se = sqrt(pooled$t), lower = pooled$qbar - qt(.975, pooled$df) * sqrt(pooled$t),
      upper = pooled$qbar + qt(.975, pooled$df) * sqrt(pooled$t), imputations = pooled$m
    )
    write_csv(est, paste0("08_qa/tmle_", id, "_all_imputations.csv"))
  } else {
    weights <- rep(1, sum(cc))
    lib <- full_lib
    V <- 10
    if (grepl("_ipcw$", id)) {
      # Use ODI response as the shared visit-response proxy for twelve-month endpoints.
      # Item-specific missingness is not separately modelled by this sensitivity.
      R <- as.integer(!is.na(d$odi_12m))
      ps <- SuperLearner(Y = R, X = data.frame(A = A, W), family = binomial(), SL.library = small_lib, cvControl = list(V = 5))
      prob <- as.vector(ps$SL.predict)
      marg <- tapply(R, A, mean)
      wt <- ifelse(A == 1, marg["1"], marg["0"]) / pmax(prob, .01)
      cutoff <- unname(quantile(wt, .99))
      wt <- pmin(wt, cutoff)
      weights <- wt[cc]
      lib <- small_lib
      V <- 5
      saveRDS(
        list(probability = prob, weights = wt, observed = cc, response_proxy = R, treatment = A, trim = cutoff, learner = ps),
        file.path(ROOT, "08_qa", paste0("tmle_", id, "_weights.rds"))
      )
      write_csv(
        data.frame(
          min_probability = min(prob), max_probability = max(prob),
          max_weight = max(weights), trim = cutoff, effective_n = sum(weights)^2 / sum(weights^2)
        ),
        paste0("04_results/tmle_", id, "_weight_summary.csv")
      )
      if (y == "odi_12m") {
        file.copy(file.path(ROOT, "08_qa", paste0("tmle_", id, "_weights.rds")),
          file.path(ROOT, "08_qa/tmle_odi12_ipcw_weights.rds"),
          overwrite = TRUE
        )
        file.copy(file.path(ROOT, paste0("04_results/tmle_", id, "_weight_summary.csv")),
          file.path(ROOT, "04_results/tmle_odi12_ipcw_weight_summary.csv"),
          overwrite = TRUE
        )
      }
    }
    fit <- tmle::tmle(
      Y = Y[cc], A = A[cc], W = W[cc, , drop = FALSE], family = family,
      Q.SL.library = lib, g.SL.library = lib, cvQinit = TRUE, V.Q = V, V.g = V, obsWeights = weights
    )
    # No fallback library or failed-imputation exclusion is permitted.
    stopifnot(length(fit$Qinit$SL.library) == length(lib), length(fit$g$SL.library) == length(lib))
    saveRDS(list(
      fit = fit, outcome = y, observed = cc, ids = d$ForlopsID[cc], covariates = vars,
      learners = lib, folds = V, seed = SEED
    ), file.path(ROOT, "03_models", paste0("tmle_", id, ".rds")))
    s <- extract_result(fit, lower)
    write_csv(
      data.frame(
        component = c(rep("Q", length(fit$Qinit$coef)), rep("g", length(fit$g$coef))),
        learner = c(names(fit$Qinit$coef), names(fit$g$coef)), weight = c(fit$Qinit$coef, fit$g$coef)
      ),
      paste0("08_qa/tmle_", id, "_ensemble_weights.csv")
    )
  }
  s$id <- id
  s$observed <- sum(cc)
  s$eligible <- nrow(d)
  s$outcome <- y
  s$family <- family
  s$lower_better <- lower
  write_csv(s, paste0("04_results/tmle_", id, "_summary.csv"))
  print(s)
}
if (Sys.getenv("ENDO_LOAD_ONLY", "") != "1") {
  if (Sys.getenv("ENDO_MODEL_WORKER", "") == "1") {
    for (id in commandArgs(trailingOnly = TRUE)) run_check(id)
  } else {
    write_csv(data.frame(id = ids), "08_qa/tmle_checks_manifest.csv")
    Sys.setenv(ENDO_MANIFEST = "tmle_checks_manifest.csv", ENDO_SCRIPT = "07_tmle_matched_checks.R", ENDO_BATCH = "tmle")
    s <- system2("python3", shQuote(file.path(CODE_DIR, "04_model_scheduler.py")))
    if (s != 0) stop("TMLE check failure. Required learners and all imputations must succeed.")
  }
}
