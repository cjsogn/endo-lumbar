# Purpose: Fit secondary models with calendar terms and treatment interactions.
# Inputs: Endpoint definitions from 03_model_jobs.R and prepared cohorts.
# Outputs: Models, standardised effects, predictive summaries and diagnostics.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "03_model_jobs.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  if (any(!args %in% names(jobs))) stop("Unknown model identifier.")
  jobs <- jobs[args]
}
stopifnot(length(jobs) > 0)
prior_predictive <- function(b) {
  # Direct independent draws from the exact proper priors. The intercept is
  # centered on the fitting design, matching brms. This avoids MCMC on priors.
  set.seed(SEED)
  n_draw <- 4000L
  sd_mu <- switch(b$job$family,
    zib = .5,
    binary = .5,
    eq = 2,
    negative = 2,
    optime = 2
  )
  tx_sd <- switch(b$job$family,
    zib = 1,
    binary = 1,
    eq = .5,
    negative = 2,
    optime = 30
  )
  ic <- switch(b$job$family,
    zib = c(0, 3),
    binary = c(0, 5),
    eq = c(.5, .5),
    negative = c(2, 2),
    optime = c(60, 30)
  )
  X <- model.matrix(as.formula(paste("~", RHS)), b$data)[, -1, drop = FALSE]
  X <- sweep(X, 2, colMeans(X))
  B <- matrix(rnorm(n_draw * ncol(X), 0, sd_mu), n_draw, ncol(X), dimnames = list(NULL, colnames(X)))
  B[, "treatmentELD"] <- rnorm(n_draw, 0, tx_sd)
  if (b$job$family == "eq") B[, "odi_baseline_z"] <- rnorm(n_draw, 0, .5)
  eta <- B %*% t(X) + rnorm(n_draw, ic[1], ic[2])
  if (b$job$family == "zib") {
    Z <- model.matrix(~ treatment + baseline_zi, b$data)[, -1, drop = FALSE]
    Z <- sweep(Z, 2, colMeans(Z))
    zi <- plogis(matrix(rnorm(n_draw * ncol(Z)), n_draw, ncol(Z)) %*% t(Z) + rnorm(n_draw, 0, 1.5))
    mu <- plogis(eta)
    phi <- rgamma(n_draw, 2, .1)
    pred <- matrix(rbeta(length(mu), as.vector(mu * phi), as.vector((1 - mu) * phi)), n_draw)
    pred[matrix(runif(length(zi)), n_draw) < zi] <- 0
    pred <- pred * b$job$scale
  } else if (b$job$family == "binary") {
    pred <- matrix(rbinom(length(eta), 1, plogis(eta)), n_draw)
  } else {
    sigma_scale <- switch(b$job$family,
      eq = .3,
      negative = 3,
      optime = 30
    )
    sigma <- abs(rt(n_draw, 3)) * sigma_scale
    pred <- eta + matrix(rnorm(length(eta)), n_draw) * sigma
  }
  stats <- data.frame(
    mean = rowMeans(pred), sd = apply(pred, 1, sd),
    zero_fraction = rowMeans(pred == 0)
  )
  saveRDS(stats, file.path(ROOT, "08_qa", paste0(b$job$id, "_prior_predictive.rds")))
  write_csv(
    data.frame(statistic = names(stats), t(apply(stats, 2, quantile, c(.025, .5, .975)))),
    paste0("08_qa/", b$job$id, "_prior_predictive_summary.csv")
  )
  invisible(TRUE)
}

# Compile unique Stan programs sequentially. Sampling below uses up to fourteen
# cores, but concurrent C++ compilers would create unnecessary memory pressure.
if (Sys.getenv("ENDO_MODEL_WORKER", "") != "1") {
  seen <- character()
  manifest <- list()
  for (id in names(jobs)) {
    b <- build_job(jobs[[id]])
    prior_predictive(b)
    code <- make_stancode(b$formula, data = b$data, family = b$family, prior = b$prior)
    stan_path <- cmdstanr::write_stan_file(code)
    if (!stan_path %in% seen) {
      cat("Compiling", id, "\n")
      flush.console()
      cmdstanr::cmdstan_model(stan_path, quiet = TRUE)
      seen <- c(seen, stan_path)
    }
    manifest[[id]] <- data.frame(
      id = id, family = b$job$family, observed = b$observed,
      fit_rows = nrow(b$data), target_rows = nrow(b$target), stan_file = stan_path,
      formula = paste(deparse(b$formula$formula), collapse = " "), seed = SEED
    )
    capture.output(list(formula = b$formula, prior = b$prior),
      file = file.path(ROOT, "08_qa", paste0(id, "_new_model_specification.txt"))
    )
  }
  write_csv(do.call(rbind, manifest), "08_qa/calendar_new_models_manifest.csv")
}

run_job <- function(id, cores) {
  log <- file.path(ROOT, "10_logs", paste0("model_", id, ".log"))
  con <- file(log, "wt")
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
  b <- build_job(jobs[[id]])
  dest <- file.path(ROOT, "03_models", paste0(id, "_calendar"))
  # Reuse only a completed model with matching code/data/priors.
  fit <- brm(b$formula,
    data = b$data, family = b$family, prior = b$prior,
    chains = 4, iter = 2000, warmup = 1000, cores = cores, seed = SEED,
    control = list(adapt_delta = .95, max_treedepth = 12), backend = "cmdstanr",
    file = dest, file_refit = "on_change", refresh = 500
  )
  dg <- diagnostics(fit, id, 12)
  write_csv(dg, paste0("08_qa/", id, "_calendar_diagnostics_initial.csv"))
  if (!dg$pass) {
    cat("Diagnostics require additional sampling and stricter integration controls.\n")
    saveRDS(fit, paste0(dest, "_initial.rds"))
    fit <- update(fit,
      iter = 4000, warmup = 2000, cores = cores, seed = SEED,
      control = list(adapt_delta = .99, max_treedepth = 14), file = NULL, refresh = 1000
    )
    saveRDS(fit, paste0(dest, ".rds"))
    dg <- diagnostics(fit, id, 14)
  }
  write_csv(dg, paste0("08_qa/", id, "_calendar_diagnostics.csv"))
  if (!dg$pass) stop(paste("NOT READY: unresolved diagnostics for", id))
  draws <- gcomp(fit, b$target, b$job$lower, b$job$scale)
  saveRDS(draws, file.path(ROOT, "04_results", paste0(id, "_calendar_draws.rds")))
  s <- cbind(
    id = id, observed = b$observed, target_n = nrow(b$target),
    summarize_effect(draws$delta, b$job$margin)
  )
  write_csv(s, paste0("04_results/", id, "_calendar_summary.csv"))
  set.seed(SEED)
  p <- posterior_predict(fit, ndraws = 400) * b$job$scale
  obs <- b$data$y * b$job$scale
  # Gaussian latent-outcome rows are excluded from observed-data PPC summaries.
  p <- p[, !is.na(obs), drop = FALSE]
  obs <- obs[!is.na(obs)]
  saveRDS(list(observed = obs, predicted = p), file.path(ROOT, "08_qa", paste0(id, "_ppc.rds")))
  capture.output(sessionInfo(), file = file.path(ROOT, "10_logs", paste0(id, "_session.txt")))
  print(s)
  list(id = id, pass = TRUE)
}
if (Sys.getenv("ENDO_MODEL_WORKER", "") == "1") {
  cores <- as.integer(Sys.getenv("ENDO_CHAIN_CORES"))
  for (id in names(jobs)) run_job(id, cores)
} else {
  status <- system2("python3", shQuote(file.path(CODE_DIR, "04_model_scheduler.py")))
  if (status != 0) stop("At least one model worker failed. Read its endpoint log.")
  cat("Calendar model batch finished.\n")
}
