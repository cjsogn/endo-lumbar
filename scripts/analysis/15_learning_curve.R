# =============================================================================
# ENDO-LUMBAR: 15 Learning Curve Analysis (ELD Only)
# Operating time and clinical outcomes by cumulative case experience
#
# Figures produced:
#   1. learning_curve_eld.png         - All ELD (n=124): Op time + ODI 3m
#   2. learning_curve_pure_disc.png   - Pure disc only (no stenosis, n=90)
#   3. learning_curve_calendar.png    - Calendar time: covariate + interaction
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("scripts", "analysis", "00_config.R"))

df <- readRDS(file.path(paths$data_clean, "df_all.rds"))

cat("=== Learning Curve Analysis (ELD) ===\n")

# --- Prepare ALL ELD data ordered by surgery date ----------------------------
df_eld_all <- df %>%
  filter(treatment == "ELD") %>%
  arrange(surgery_date) %>%
  mutate(
    case_num = row_number(),
    approach = recode_factor(approach,
                            "Midline" = "Interlaminar", "Wiltse" = "Transforaminal",
                            "Other" = "Other",
                            .ordered = FALSE),
    has_stenosis = (stenosis_central == 1 | stenosis_lateral == 1 | stenosis_foraminal == 1),
    pure_disc = !has_stenosis
  )

# Exclude implausible operating time outlier (>=500 min)
n_op_outlier <- sum(df_eld_all$operating_time >= 500, na.rm = TRUE)
df_eld_all$operating_time[df_eld_all$operating_time >= 500] <- NA
cat(sprintf("Operating time outliers set to NA: %d (values >= 500 min)\n", n_op_outlier))

n_eld <- nrow(df_eld_all)
n_pure <- sum(df_eld_all$pure_disc)
n_stenosis <- sum(df_eld_all$has_stenosis)

cat(sprintf("Total ELD: %d\n", n_eld))
cat(sprintf("  Pure disc herniation: %d\n", n_pure))
cat(sprintf("  Concomitant stenosis: %d\n", n_stenosis))
cat(sprintf("Date range: %s to %s\n",
            min(df_eld_all$surgery_date), max(df_eld_all$surgery_date)))

# MSD reference values (all MSD patients)
msd_op_time_mean <- mean(df$operating_time[df$treatment == "MSD"], na.rm = TRUE)
msd_odi_3m_mean  <- mean(df$odi_3m[df$treatment == "MSD"], na.rm = TRUE)

cat(sprintf("MSD reference - Op time: %.1f min, ODI 3m: %.1f\n",
            msd_op_time_mean, msd_odi_3m_mean))

# --- Helper: assign quartiles ------------------------------------------------
assign_quartiles <- function(d) {
  n <- nrow(d)
  d %>% mutate(quartile = cut(case_num,
                               breaks = c(0, ceiling(n/4), ceiling(n/2),
                                          ceiling(3*n/4), n),
                               labels = c("Q1", "Q2", "Q3", "Q4"),
                               include.lowest = TRUE))
}

# --- Helper: compute quartile stats ------------------------------------------
compute_quartile_stats <- function(d) {
  d %>%
    group_by(quartile) %>%
    dplyr::summarise(
      n = n(),
      case_range = paste0(min(case_num), "-", max(case_num)),
      date_range = paste0(min(surgery_date), " to ", max(surgery_date)),
      op_time_n    = sum(!is.na(operating_time)),
      op_time_mean = mean(operating_time, na.rm = TRUE),
      op_time_sd   = sd(operating_time, na.rm = TRUE),
      op_time_median = median(operating_time, na.rm = TRUE),
      odi_3m_n    = sum(!is.na(odi_3m)),
      odi_3m_mean = mean(odi_3m, na.rm = TRUE),
      odi_3m_sd   = sd(odi_3m, na.rm = TRUE),
      nrs_leg_3m_n    = sum(!is.na(nrs_leg_3m)),
      nrs_leg_3m_mean = mean(nrs_leg_3m, na.rm = TRUE),
      nrs_leg_3m_sd   = sd(nrs_leg_3m, na.rm = TRUE),
      n_interlaminar = sum(approach == "Interlaminar", na.rm = TRUE),
      n_transforaminal  = sum(approach == "Transforaminal", na.rm = TRUE),
      pct_transforaminal = 100 * n_transforaminal / n,
      n_stenosis = sum(has_stenosis, na.rm = TRUE),
      .groups = "drop"
    )
}

# --- Helper: build operating time panel (uses ALL cases) ---------------------
build_optime_panel <- function(d, q_stats, n_total, subtitle_extra = "") {
  q_labels <- q_stats %>%
    mutate(
      x_pos = (as.numeric(gsub(".*-(\\d+)", "\\1", case_range)) +
                 as.numeric(gsub("(\\d+)-.*", "\\1", case_range))) / 2,
      label = sprintf("%.0f min", op_time_mean)
    )
  approach_colors <- c("Interlaminar" = "#4393C3", "Transforaminal" = "#D6604D", "Other" = "grey50")

  ggplot(d, aes(x = case_num, y = operating_time)) +
    geom_hline(yintercept = msd_op_time_mean,
               linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
    annotate("text", x = max(d$case_num) * 0.98, y = msd_op_time_mean + 4,
             label = sprintf("MSD mean (%.0f min)", msd_op_time_mean),
             hjust = 1, size = 3, color = "#B2182B", fontface = "italic") +
    annotate("rect",
             xmin = c(0.5, ceiling(n_total/4)+0.5, ceiling(n_total/2)+0.5, ceiling(3*n_total/4)+0.5),
             xmax = c(ceiling(n_total/4)+0.5, ceiling(n_total/2)+0.5, ceiling(3*n_total/4)+0.5, n_total+0.5),
             ymin = -Inf, ymax = Inf,
             fill = rep(c("grey95", "white"), 2), alpha = 0.5) +
    geom_smooth(method = "loess", span = 0.5, se = TRUE,
                color = "#2166AC", fill = "#92C5DE", alpha = 0.3, linewidth = 1) +
    geom_point(aes(color = approach), size = 2, alpha = 0.7) +
    scale_color_manual(values = approach_colors, name = "Approach") +
    geom_text(data = q_labels,
              aes(x = x_pos, y = max(d$operating_time, na.rm = TRUE) * 0.95, label = label),
              inherit.aes = FALSE, size = 3, fontface = "bold", color = "grey30") +
    geom_text(data = q_labels,
              aes(x = x_pos, y = max(d$operating_time, na.rm = TRUE) * 0.90,
                  label = as.character(quartile)),
              inherit.aes = FALSE, size = 2.8, color = "grey50") +
    scale_x_continuous(breaks = seq(0, n_total, by = 20),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = "Cumulative ELD case number", y = "Operating time (minutes)",
         title = paste0("Operating Time (all cases, n=", nrow(d), ")", subtitle_extra)) +
    theme(
      legend.position = c(0.85, 0.75),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )
}

# --- Helper: build ODI panel (observed cases only) ---------------------------
build_odi_panel <- function(d, n_total, cor_result, subtitle_extra = "") {
  d_odi <- d %>% filter(!is.na(odi_3m))

  ggplot(d_odi, aes(x = case_num, y = odi_3m)) +
    geom_hline(yintercept = msd_odi_3m_mean,
               linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
    annotate("text", x = max(d$case_num, na.rm = TRUE) * 0.98,
             y = msd_odi_3m_mean + 2,
             label = sprintf("MSD mean (%.1f)", msd_odi_3m_mean),
             hjust = 1, size = 3, color = "#B2182B", fontface = "italic") +
    annotate("rect",
             xmin = c(0.5, ceiling(n_total/4)+0.5, ceiling(n_total/2)+0.5, ceiling(3*n_total/4)+0.5),
             xmax = c(ceiling(n_total/4)+0.5, ceiling(n_total/2)+0.5, ceiling(3*n_total/4)+0.5, n_total+0.5),
             ymin = -Inf, ymax = Inf,
             fill = rep(c("grey95", "white"), 2), alpha = 0.5) +
    geom_smooth(method = "loess", span = 0.6, se = TRUE,
                color = "#2166AC", fill = "#92C5DE", alpha = 0.3, linewidth = 1) +
    geom_point(color = "#2166AC", size = 2, alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, n_total, by = 20),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = "Cumulative ELD case number", y = "ODI score at 3 months",
         title = paste0("ODI at 3 Months (n=", nrow(d_odi), " with follow-up)", subtitle_extra)) +
    annotate("text", x = min(d_odi$case_num) + 2,
             y = max(d_odi$odi_3m, na.rm = TRUE) * 0.95,
             label = sprintf("Spearman \u03c1 = %.3f, p = %.3f",
                             cor_result$estimate, cor_result$p.value),
             hjust = 0, size = 3, color = "grey30")
}

# #############################################################################
# FIGURE 1: ALL ELD PATIENTS (n=124)
# #############################################################################

cat("\n===== Figure 1: All ELD patients (n=124) =====\n")

df_eld_all <- assign_quartiles(df_eld_all)
q_stats_all <- compute_quartile_stats(df_eld_all)

cat("Quartile statistics (all ELD):\n")
print(q_stats_all)
write.csv(q_stats_all, file.path(paths$tables, "table_learning_curve_quartiles.csv"),
          row.names = FALSE)

# Spearman correlations
cor_optime_all <- cor.test(df_eld_all$case_num, df_eld_all$operating_time,
                           method = "spearman", exact = FALSE)
df_odi_all <- df_eld_all %>% filter(!is.na(odi_3m))
cor_odi_all <- cor.test(df_odi_all$case_num, df_odi_all$odi_3m,
                        method = "spearman", exact = FALSE)

cat(sprintf("\nSpearman (all ELD):\n  Op time: rho=%.3f, p=%.4f\n  ODI 3m:  rho=%.3f, p=%.4f\n",
            cor_optime_all$estimate, cor_optime_all$p.value,
            cor_odi_all$estimate, cor_odi_all$p.value))

p1_optime <- build_optime_panel(df_eld_all, q_stats_all, n_eld)
p1_odi    <- build_odi_panel(df_eld_all, n_eld, cor_odi_all)

p1_combined <- p1_optime + p1_odi +
  plot_annotation(
    title = "Learning Curve Analysis: All ELD Cases",
    subtitle = sprintf("n = %d consecutive ELD cases ordered by surgery date", n_eld),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, color = "grey30"))
  )

save_fig(p1_combined, "learning_curve_eld.png", width = 12, height = 5.5)
cat("  Saved: learning_curve_eld.png\n")

# #############################################################################
# FIGURE 2: PURE DISC HERNIATION ONLY (no concomitant stenosis)
# #############################################################################

cat("\n===== Figure 2: Pure disc herniation ELD (n=90) =====\n")

df_pure <- df_eld_all %>%
  filter(pure_disc) %>%
  mutate(case_num = row_number())  # Re-number sequentially

n_pure_eld <- nrow(df_pure)
df_pure <- assign_quartiles(df_pure)
q_stats_pure <- compute_quartile_stats(df_pure)

cat("Quartile statistics (pure disc):\n")
print(q_stats_pure)
write.csv(q_stats_pure, file.path(paths$tables, "table_learning_curve_quartiles_pure_disc.csv"),
          row.names = FALSE)

# MSD pure disc reference
df_msd_pure <- df %>%
  filter(treatment == "MSD",
         !(stenosis_central == 1 | stenosis_lateral == 1 | stenosis_foraminal == 1))
msd_pure_optime <- mean(df_msd_pure$operating_time, na.rm = TRUE)
msd_pure_odi3m  <- mean(df_msd_pure$odi_3m, na.rm = TRUE)
cat(sprintf("MSD pure disc reference - Op time: %.1f min, ODI 3m: %.1f\n",
            msd_pure_optime, msd_pure_odi3m))

cor_optime_pure <- cor.test(df_pure$case_num, df_pure$operating_time,
                            method = "spearman", exact = FALSE)
df_odi_pure <- df_pure %>% filter(!is.na(odi_3m))
cor_odi_pure <- cor.test(df_odi_pure$case_num, df_odi_pure$odi_3m,
                         method = "spearman", exact = FALSE)

cat(sprintf("Spearman (pure disc):\n  Op time: rho=%.3f, p=%.4f\n  ODI 3m:  rho=%.3f, p=%.4f\n",
            cor_optime_pure$estimate, cor_optime_pure$p.value,
            cor_odi_pure$estimate, cor_odi_pure$p.value))

p2_optime <- build_optime_panel(df_pure, q_stats_pure, n_pure_eld)
p2_odi    <- build_odi_panel(df_pure, n_pure_eld, cor_odi_pure)

p2_combined <- p2_optime + p2_odi +
  plot_annotation(
    title = "Learning Curve Analysis: Pure Disc Herniation ELD Cases",
    subtitle = sprintf("n = %d ELD cases without concomitant stenosis (excluded %d with stenosis)",
                       n_pure_eld, n_stenosis),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, color = "grey30"))
  )

save_fig(p2_combined, "learning_curve_pure_disc.png", width = 12, height = 5.5)
cat("  Saved: learning_curve_pure_disc.png\n")

# #############################################################################
# FIGURE 3: CALENDAR TIME (covariate and interaction models)
# #############################################################################

cat("\n===== Figure 3: Calendar time learning curve =====\n")

# Use surgery_date as numeric days for plotting; calendar_time for model
df_eld_all <- df_eld_all %>%
  mutate(surgery_date_num = as.numeric(surgery_date),
         cal_time_centered = calendar_time - mean(calendar_time))

# --- Panel A: Operating time vs calendar date (all cases) --------------------

# Fit calendar time models for operating time
# Model 1: calendar time as linear covariate
lm_caltime_lin <- lm(operating_time ~ calendar_time, data = df_eld_all)
cat("\nLinear calendar time model (operating time):\n")
cat(sprintf("  Coefficient: %.3f min/day (SE=%.3f), p=%.4f\n",
            coef(lm_caltime_lin)[2], summary(lm_caltime_lin)$coefficients[2,2],
            summary(lm_caltime_lin)$coefficients[2,4]))

# Model 2: calendar time with natural cubic spline (NCS, 3 df)
library(splines)
lm_caltime_rcs <- lm(operating_time ~ ns(calendar_time, df = 3), data = df_eld_all)
cat(sprintf("  NCS model R-sq: %.3f (vs linear R-sq: %.3f)\n",
            summary(lm_caltime_rcs)$r.squared,
            summary(lm_caltime_lin)$r.squared))

# Generate prediction grid
pred_grid <- tibble(
  calendar_time = seq(min(df_eld_all$calendar_time), max(df_eld_all$calendar_time), length.out = 200),
  surgery_date = as.Date(seq(min(df_eld_all$surgery_date_num),
                              max(df_eld_all$surgery_date_num), length.out = 200),
                          origin = "1970-01-01")
)

# Linear predictions
pred_lin <- predict(lm_caltime_lin, newdata = pred_grid, interval = "confidence")
pred_grid$fit_lin <- pred_lin[, "fit"]
pred_grid$lwr_lin <- pred_lin[, "lwr"]
pred_grid$upr_lin <- pred_lin[, "upr"]

# RCS predictions
pred_rcs <- predict(lm_caltime_rcs, newdata = pred_grid, interval = "confidence")
pred_grid$fit_rcs <- pred_rcs[, "fit"]
pred_grid$lwr_rcs <- pred_rcs[, "lwr"]
pred_grid$upr_rcs <- pred_rcs[, "upr"]

approach_colors <- c("Interlaminar" = "#4393C3", "Transforaminal" = "#D6604D", "Other" = "grey50")

p3a <- ggplot(df_eld_all, aes(x = surgery_date, y = operating_time)) +
  geom_hline(yintercept = msd_op_time_mean,
             linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
  annotate("text", x = max(df_eld_all$surgery_date) - 30, y = msd_op_time_mean + 4,
           label = sprintf("MSD mean (%.0f min)", msd_op_time_mean),
           hjust = 1, size = 3, color = "#B2182B", fontface = "italic") +
  # RCS fit (flexible)
  geom_ribbon(data = pred_grid, aes(x = surgery_date, ymin = lwr_rcs, ymax = upr_rcs),
              inherit.aes = FALSE, fill = "#2166AC", alpha = 0.15) +
  geom_line(data = pred_grid, aes(x = surgery_date, y = fit_rcs),
            inherit.aes = FALSE, color = "#2166AC", linewidth = 1) +
  # Linear fit
  geom_line(data = pred_grid, aes(x = surgery_date, y = fit_lin),
            inherit.aes = FALSE, color = "#D6604D", linewidth = 0.8, linetype = "dashed") +
  # Points
  geom_point(aes(color = approach), size = 2, alpha = 0.7) +
  scale_color_manual(values = approach_colors, name = "Approach") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "4 months") +
  labs(x = "Surgery date", y = "Operating time (minutes)",
       title = sprintf("Operating Time vs Calendar Date (all cases, n=%d)", n_eld)) +
  annotate("text", x = min(df_eld_all$surgery_date) + 30,
           y = max(df_eld_all$operating_time, na.rm = TRUE) * 0.92,
           label = sprintf("Linear: %.2f min/day, p=%.3f\nNCS (3 df): R\u00b2=%.3f",
                           coef(lm_caltime_lin)[2],
                           summary(lm_caltime_lin)$coefficients[2,4],
                           summary(lm_caltime_rcs)$r.squared),
           hjust = 0, size = 3, color = "grey30") +
  theme(
    legend.position = c(0.88, 0.78),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# --- Panel B: ODI 3m vs calendar date (observed cases) -----------------------

df_eld_odi_cal <- df_eld_all %>% filter(!is.na(odi_3m))

lm_odi_lin <- lm(odi_3m ~ calendar_time, data = df_eld_odi_cal)
lm_odi_rcs <- lm(odi_3m ~ ns(calendar_time, df = 3), data = df_eld_odi_cal)

cat(sprintf("\nLinear calendar time model (ODI 3m):\n  Coefficient: %.3f/day (SE=%.3f), p=%.4f\n",
            coef(lm_odi_lin)[2], summary(lm_odi_lin)$coefficients[2,2],
            summary(lm_odi_lin)$coefficients[2,4]))
cat(sprintf("  NCS model R-sq: %.3f (vs linear R-sq: %.3f)\n",
            summary(lm_odi_rcs)$r.squared, summary(lm_odi_lin)$r.squared))

pred_grid_odi <- tibble(
  calendar_time = seq(min(df_eld_odi_cal$calendar_time),
                      max(df_eld_odi_cal$calendar_time), length.out = 200),
  surgery_date = as.Date(seq(min(df_eld_odi_cal$surgery_date_num),
                              max(df_eld_odi_cal$surgery_date_num), length.out = 200),
                          origin = "1970-01-01")
)

pred_odi_lin <- predict(lm_odi_lin, newdata = pred_grid_odi, interval = "confidence")
pred_grid_odi$fit_lin <- pred_odi_lin[, "fit"]
pred_grid_odi$lwr_lin <- pred_odi_lin[, "lwr"]
pred_grid_odi$upr_lin <- pred_odi_lin[, "upr"]

pred_odi_rcs <- predict(lm_odi_rcs, newdata = pred_grid_odi, interval = "confidence")
pred_grid_odi$fit_rcs <- pred_odi_rcs[, "fit"]
pred_grid_odi$lwr_rcs <- pred_odi_rcs[, "lwr"]
pred_grid_odi$upr_rcs <- pred_odi_rcs[, "upr"]

p3b <- ggplot(df_eld_odi_cal, aes(x = surgery_date, y = odi_3m)) +
  geom_hline(yintercept = msd_odi_3m_mean,
             linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
  annotate("text", x = max(df_eld_odi_cal$surgery_date) - 30,
           y = msd_odi_3m_mean + 2,
           label = sprintf("MSD mean (%.1f)", msd_odi_3m_mean),
           hjust = 1, size = 3, color = "#B2182B", fontface = "italic") +
  # RCS fit
  geom_ribbon(data = pred_grid_odi, aes(x = surgery_date, ymin = lwr_rcs, ymax = upr_rcs),
              inherit.aes = FALSE, fill = "#2166AC", alpha = 0.15) +
  geom_line(data = pred_grid_odi, aes(x = surgery_date, y = fit_rcs),
            inherit.aes = FALSE, color = "#2166AC", linewidth = 1) +
  # Linear fit
  geom_line(data = pred_grid_odi, aes(x = surgery_date, y = fit_lin),
            inherit.aes = FALSE, color = "#D6604D", linewidth = 0.8, linetype = "dashed") +
  # Points
  geom_point(color = "#2166AC", size = 2, alpha = 0.5) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "4 months") +
  labs(x = "Surgery date", y = "ODI score at 3 months",
       title = sprintf("ODI at 3 Months vs Calendar Date (n=%d with follow-up)",
                       nrow(df_eld_odi_cal))) +
  annotate("text", x = min(df_eld_odi_cal$surgery_date) + 30,
           y = max(df_eld_odi_cal$odi_3m, na.rm = TRUE) * 0.92,
           label = sprintf("Linear: %.3f/day, p=%.3f\nNCS (3 df): R\u00b2=%.3f",
                           coef(lm_odi_lin)[2],
                           summary(lm_odi_lin)$coefficients[2,4],
                           summary(lm_odi_rcs)$r.squared),
           hjust = 0, size = 3, color = "grey30") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Legend for model fits
p3a <- p3a +
  annotate("segment", x = min(df_eld_all$surgery_date) + 15,
           xend = min(df_eld_all$surgery_date) + 70,
           y = max(df_eld_all$operating_time, na.rm = TRUE) * 0.82,
           yend = max(df_eld_all$operating_time, na.rm = TRUE) * 0.82,
           color = "#2166AC", linewidth = 1) +
  annotate("text", x = min(df_eld_all$surgery_date) + 75,
           y = max(df_eld_all$operating_time, na.rm = TRUE) * 0.82,
           label = "NCS (3 df)", hjust = 0, size = 2.8, color = "#2166AC") +
  annotate("segment", x = min(df_eld_all$surgery_date) + 15,
           xend = min(df_eld_all$surgery_date) + 70,
           y = max(df_eld_all$operating_time, na.rm = TRUE) * 0.77,
           yend = max(df_eld_all$operating_time, na.rm = TRUE) * 0.77,
           color = "#D6604D", linewidth = 0.8, linetype = "dashed") +
  annotate("text", x = min(df_eld_all$surgery_date) + 75,
           y = max(df_eld_all$operating_time, na.rm = TRUE) * 0.77,
           label = "Linear", hjust = 0, size = 2.8, color = "#D6604D")

p3_combined <- p3a + p3b +
  plot_annotation(
    title = "Learning Curve by Calendar Time: Linear Covariate vs Natural Cubic Spline",
    subtitle = sprintf("ELD cases (n = %d), calendar time modeled as linear trend and NCS with 3 df",
                       n_eld),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, color = "grey30"))
  )

save_fig(p3_combined, "learning_curve_calendar.png", width = 13, height = 5.5)
cat("  Saved: learning_curve_calendar.png\n")

# #############################################################################
# SUMMARY TABLE: All correlations and comparisons
# #############################################################################

cat("\n===== Summary =====\n")

summary_tbl <- tibble(
  analysis = c("All ELD", "All ELD", "Pure disc", "Pure disc",
               "Calendar (linear)", "Calendar (linear)"),
  outcome = c("Operating time", "ODI 3m", "Operating time", "ODI 3m",
              "Operating time", "ODI 3m"),
  n = c(n_eld, nrow(df_odi_all), n_pure_eld, nrow(df_odi_pure),
        n_eld, nrow(df_eld_odi_cal)),
  spearman_rho = c(cor_optime_all$estimate, cor_odi_all$estimate,
                   cor_optime_pure$estimate, cor_odi_pure$estimate,
                   NA, NA),
  spearman_p = c(cor_optime_all$p.value, cor_odi_all$p.value,
                 cor_optime_pure$p.value, cor_odi_pure$p.value,
                 NA, NA),
  linear_coef = c(NA, NA, NA, NA,
                  coef(lm_caltime_lin)[2], coef(lm_odi_lin)[2]),
  linear_p = c(NA, NA, NA, NA,
               summary(lm_caltime_lin)$coefficients[2,4],
               summary(lm_odi_lin)$coefficients[2,4]),
  q1_mean = c(q_stats_all$op_time_mean[1], q_stats_all$odi_3m_mean[1],
              q_stats_pure$op_time_mean[1], q_stats_pure$odi_3m_mean[1],
              NA, NA),
  q4_mean = c(q_stats_all$op_time_mean[4], q_stats_all$odi_3m_mean[4],
              q_stats_pure$op_time_mean[4], q_stats_pure$odi_3m_mean[4],
              NA, NA)
)

print(summary_tbl)
write.csv(summary_tbl, file.path(paths$tables, "table_learning_curve_summary.csv"),
          row.names = FALSE)

# =============================================================================
# DONE
# =============================================================================

cat("\n=== Learning Curve Analysis Complete ===\n")
cat("Figures saved:\n")
cat("  1. learning_curve_eld.png        (all ELD, n=124)\n")
cat("  2. learning_curve_pure_disc.png  (pure disc, n=90)\n")
cat("  3. learning_curve_calendar.png   (calendar time models)\n")
cat("Tables saved:\n")
cat("  - table_learning_curve_quartiles.csv\n")
cat("  - table_learning_curve_quartiles_pure_disc.csv\n")
cat("  - table_learning_curve_summary.csv\n")
