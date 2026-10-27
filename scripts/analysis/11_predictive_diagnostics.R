# Purpose: Assess predictive distributions, influence and prior/posterior dispersion.
# Inputs: Combined Bayesian endpoint list and fitted outcome models.
# Outputs: Numerical predictive summaries, PSIS-LOO and dispersion diagnostics.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
tab <- read.csv(file.path(ROOT, "04_results/calendar_all_bayesian_results.csv"))
stats <- function(y) c(
  mean = mean(y), sd = sd(y), median = median(y), q10 = unname(quantile(y, .1)),
  q90 = unname(quantile(y, .9)), zero_fraction = mean(y == 0)
)
pp <- llsum <- dispersion <- list()
for (id in tab$id) {
  fit <- readRDS(file.path(ROOT, "03_models", paste0(id, "_calendar.rds")))
  # Compare posterior dispersion of mean-model slopes with their stated normal
  # priors. This describes prior influence, not a pass/fail test of fit.
  slopes <- fixef(fit)
  pr <- as.data.frame(prior_summary(fit))
  pr <- pr[pr$class == "b" & pr$dpar == "" & pr$nlpar == "", ]
  generic <- pr$prior[pr$coef == ""][1]
  terms <- rownames(slopes)
  terms <- terms[!startsWith(terms, "Intercept") & !startsWith(terms, "zi_")]
  dispersion[[id]] <- do.call(rbind, lapply(terms, function(term) {
    p <- pr$prior[pr$coef == term]
    if (!length(p) || !nzchar(p[1])) p <- generic
    stopifnot(grepl("^normal\\(", p[1]))
    sd0 <- as.numeric(sub(".*,[[:space:]]*([^)]*)\\)", "\\1", p[1]))
    data.frame(
      id = id, coefficient = term, prior_sd = sd0,
      posterior_sd = slopes[term, "Est.Error"],
      posterior_to_prior_sd = slopes[term, "Est.Error"] / sd0
    )
  }))
  response <- all.vars(fit$formula$formula)[1]
  y <- fit$data[[response]]
  is_observed <- !is.na(y)
  sf <- if (id %in% c("odi_3m", "odi_12m")) 100 else if (grepl("nrs_", id)) 10 else 1
  if (is.factor(y)) y <- as.integer(y)
  set.seed(SEED)
  pred <- posterior_predict(fit, ndraws = 400)
  stopifnot(ncol(pred) == length(y))
  pred <- pred[, is_observed, drop = FALSE] * sf
  y <- y[is_observed] * sf
  ob <- stats(y)
  sim <- t(apply(pred, 1, stats))
  pp[[id]] <- data.frame(
    id = id, statistic = names(ob), observed = ob,
    predictive_median = apply(sim, 2, median), predictive_lower = apply(sim, 2, quantile, .025),
    predictive_upper = apply(sim, 2, quantile, .975), p_predictive_ge_observed = colMeans(sweep(sim, 2, ob, ">="))
  )
  saveRDS(list(observed = y, predicted = pred), file.path(ROOT, "08_qa", paste0(id, "_ppc.rds")))
  # Pointwise PSIS diagnostics use only actually observed outcome likelihoods.
  # No comparison of ELPD between differently transformed outcome families.
  ll <- log_lik(fit)
  stopifnot(ncol(ll) == length(is_observed))
  ll <- ll[, is_observed, drop = FALSE]
  nch <- posterior::nchains(posterior::as_draws_array(fit))
  chain <- rep(seq_len(nch), each = nrow(ll) / nch)
  reff <- loo::relative_eff(exp(ll), chain_id = chain, cores = 1)
  lo <- loo::loo(ll, r_eff = reff, cores = 1)
  saveRDS(lo, file.path(ROOT, "08_qa", paste0(id, "_loo.rds")))
  pk <- loo::pareto_k_values(lo)
  llsum[[id]] <- data.frame(
    id = id, n_observed = length(y), max_pareto_k = max(pk),
    n_k_gt_07 = sum(pk > .7), n_k_gt_1 = sum(pk > 1), elpd_loo = lo$estimates["elpd_loo", "Estimate"]
  )
  cat(id, "PPC and PSIS completed\n")
  rm(fit, pred, ll, lo)
  gc(verbose = FALSE)
}
write_csv(do.call(rbind, pp), "08_qa/calendar_posterior_predictive_statistics.csv")
write_csv(do.call(rbind, llsum), "08_qa/calendar_psis_loo_diagnostics.csv")
write_csv(do.call(rbind, dispersion), "08_qa/calendar_prior_posterior_dispersion.csv")
