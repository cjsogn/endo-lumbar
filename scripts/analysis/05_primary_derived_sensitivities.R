# Purpose: Calculate E-values and outcome-shift sensitivities from primary draws.
# Inputs: Prepared three-month cohort and primary ODI posterior contrasts.
# Outputs: E-values, pattern-mixture shifts and tipping-point summaries in 04_results.

source(file.path(Sys.getenv("ENDO_CODE_DIR"), "00_config.R"))
d <- readRDS(file.path(ROOT, "02_data/derived/cohort_3m.rds"))
x <- readRDS(file.path(ROOT, "04_results/odi_3m_calendar_draws.rds"))$delta
stopifnot(length(x) > 0)
sd_y <- sd(d$odi_3m, na.rm = TRUE)
ev <- function(distance) {
  rr <- exp(.91 * pmax(distance, 0) / sd_y)
  rr + sqrt(rr * (rr - 1))
}
tab <- do.call(rbind, lapply(list(calendar = x), function(z) {
  lo <- unname(quantile(z, .025))
  hi <- unname(quantile(z, .975))
  data.frame(
    sd_observed_odi = sd_y, mean = mean(z), lower = lo, upper = hi,
    evalue_null_point = ev(abs(mean(z))),
    evalue_null_interval = if (lo <= 0 && hi >= 0) 1 else ev(min(abs(c(lo, hi)))),
    evalue_margin_point = ev(mean(z) + 7), evalue_margin_interval = ev(lo + 7)
  )
}))
write_csv(data.frame(model = rownames(tab), tab), "04_results/evalues_calendar.csv")
n_eld <- sum(is.na(d$odi_3m) & d$treatment == "ELD")
n_msd <- sum(is.na(d$odi_3m) & d$treatment == "MSD")
N <- nrow(d)
# Apply factual-arm shifts to conditional means without bounds clipping.
# This sensitivity calculation does not model observed ODI outside 0 to 100.
pm <- do.call(rbind, lapply(seq(-8, 8, 2), function(delta)
  cbind(shift = delta, summarize_effect(x + delta * (n_msd - n_eld) / N, 7))))
tp <- do.call(rbind, lapply(seq(0, 25, .5), function(delta)
  cbind(adversarial_shift = delta, summarize_effect(x - delta * (n_msd + n_eld) / N, 7))))
write_csv(pm, "04_results/pattern_mixture_calendar.csv")
write_csv(tp, "04_results/tipping_point_calendar.csv")
tipped <- tp$adversarial_shift[tp$p_ni <= .95]
write_csv(data.frame(
  eld_missing = n_eld, msd_missing = n_msd, target_n = N,
  criterion = "PNI <= 0.95 (registered declaration requires > 0.95)",
  first_grid_value = if (length(tipped)) min(tipped) else NA_real_,
  grid_max = 25, grid_step = .5
), "04_results/tipping_point_calendar_summary.csv")
