## 06_plots.R
## All figures for MSc thesis — HR/HRV and everyday listening
## Output: figures/
## -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggeffects)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)

dir.create("figures", showWarnings = FALSE)

# =============================================================
# LOAD DATA AND MODELS
# =============================================================

data <- readRDS("ema_physio_aligned.rds") %>%
  filter(between(hour(completed_ts), 7, 23)) %>%
  filter(!is.na(hr_norm), !is.na(hrv_norm)) %>%
  mutate(
    participant_id    = factor(participant_id),
    listening_context = factor(listening_context_label,
                               levels = c("Baseline Morning",
                                          "Baseline Evening",
                                          "Conversation",
                                          "Focused Listening",
                                          "Passive Listening")),
    emotional_state_INT        = as.numeric(emotional_state),
    fatigue_cognitive_load_INT = as.numeric(fatigue_cognitive_load),
    acoustic_perception_INT    = as.numeric(acoustic_perception),
    listening_effort_INT       = as.numeric(listening_effort),
    time_of_day = case_when(
      hour(completed_ts) < 12 ~ "morning",
      hour(completed_ts) < 18 ~ "afternoon",
      TRUE                    ~ "evening"
    ),
    time_of_day = factor(time_of_day,
                         levels = c("morning", "afternoon", "evening")),
    hour_cont = hour(completed_ts) + minute(completed_ts) / 60
  )

data_nb <- data %>%
  filter(listening_context %in%
           c("Conversation", "Focused Listening", "Passive Listening")) %>%
  droplevels()

physio_raw <- readRDS("physio_processed.rds")

model_rq1_hr         <- readRDS("model_rq1_hr.rds")
model_rq1_hrv        <- readRDS("model_rq1_hrv.rds")
model_rq1_hr_ext     <- readRDS("model_rq1_hr_ext.rds")
model_rq1_hrv_ext    <- readRDS("model_rq1_hrv_ext.rds")
model_rq2_hr_cat     <- readRDS("model_rq2_hr_cat.rds")
model_rq2_hrv_cat    <- readRDS("model_rq2_hrv_cat.rds")
model_rq2_hrv_spline <- readRDS("model_rq2_hrv_spline.rds")
model_baseline_hr    <- readRDS("model_baseline_hr.rds")
model_baseline_hrv   <- readRDS("model_baseline_hrv.rds")
model_rq3_hr         <- readRDS("model_rq3_hr.rds")
model_rq3_hrv        <- readRDS("model_rq3_hrv.rds")
model_rq3_hr_add     <- readRDS("model_rq3_hr_add.rds")
model_rq3_hrv_add    <- readRDS("model_rq3_hrv_add.rds")

# =============================================================
# THEME AND COLOURS
# =============================================================

context_colors <- c(
  "Baseline Morning"  = "#5F5E5A",
  "Baseline Evening"  = "#3C3489",
  "Conversation"      = "#185FA5",
  "Focused Listening" = "#0F6E56",
  "Passive Listening" = "#993C1D"
)

theme_sofia <- theme_minimal(base_size = 16) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(size = 13, color = "grey40"),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 13),
    legend.position  = "top",
    legend.title     = element_text(size = 13),
    legend.text      = element_text(size = 12),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 13)
  )

theme_set(theme_sofia)

# =============================================================
# HELPER: save figure
# =============================================================

save_fig <- function(plot, filename, width = 12, height = 5) {
  ggsave(paste0("figures/", filename),
         plot, width = width, height = height, dpi = 150)
  message("Saved: figures/", filename)
}

# =============================================================
# FIG 7: Raw HR and HRV distribution per participant
# (Figure 7 in thesis)
# =============================================================

p_raw_hr <- physio_raw %>%
  rename(participant_id = ID) %>%
  filter(!is.na(MeanHR)) %>%
  ggplot(aes(x = reorder(participant_id, MeanHR,
                          FUN = median), y = MeanHR,
             fill = participant_id)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7,
               show.legend = FALSE) +
  labs(title = "A — Raw HR distribution per participant",
       x = "Participant", y = "HR (bpm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_raw_hrv <- physio_raw %>%
  rename(participant_id = ID) %>%
  filter(!is.na(niHR)) %>%
  ggplot(aes(x = reorder(participant_id, niHR,
                          FUN = median), y = niHR,
             fill = participant_id)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7,
               show.legend = FALSE) +
  labs(title = "B — Raw HRV distribution per participant",
       x = "Participant", y = "RMSSD proxy (ms)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

save_fig(p_raw_hr + p_raw_hrv,
         "fig07_raw_hr_hrv_per_participant.png",
         width = 14, height = 6)

# =============================================================
# FIG 8: HR and HRV distribution at EMA moments
# =============================================================

p_distr_hr <- ggplot(data, aes(x = hr_mean)) +
  geom_histogram(bins = 30, fill = "#185FA5",
                 color = "white", alpha = 0.85) +
  labs(title = "A — HR distribution at EMA moments",
       x = "HR (bpm)", y = "Count")

p_distr_hrv <- ggplot(data, aes(x = rmssd_mean)) +
  geom_histogram(bins = 30, fill = "#0F6E56",
                 color = "white", alpha = 0.85) +
  labs(title = "B — HRV distribution at EMA moments",
       x = "RMSSD proxy (ms)", y = "Count")

save_fig(p_distr_hr + p_distr_hrv,
         "fig08_hr_hrv_distribution_ema.png")

# =============================================================
# FIG 9: EMA overview — responses per context and time of day
# =============================================================

p_ema_count <- data %>%
  count(listening_context_label) %>%
  ggplot(aes(x = reorder(listening_context_label, -n),
             y = n, fill = listening_context_label)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(values = context_colors) +
  labs(title = "A — EMA responses per listening context",
       x = NULL, y = "Number of responses") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_ema_time <- data %>%
  mutate(hour = hour(completed_ts)) %>%
  count(hour, listening_context_label) %>%
  group_by(hour) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = hour, y = pct,
             fill = listening_context_label)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = context_colors, name = NULL) +
  scale_x_continuous(breaks = seq(7, 23, by = 2),
                     labels = paste0(seq(7, 23, by = 2), ":00")) +
  labs(title = "B — Context distribution across the day",
       x = "Hour of day", y = "% of responses") +
  theme(legend.position = "bottom")

save_fig(p_ema_count + p_ema_time,
         "fig09_ema_overview.png",
         width = 14, height = 6)

# =============================================================
# FIG 10: Normalised HR and HRV across contexts (boxplot)
# =============================================================

p_ctx_hr <- data %>%
  ggplot(aes(x = listening_context, y = hr_norm,
             fill = listening_context)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(values = context_colors,
                    guide = "none") +
  labs(title = "A — Normalised HR across listening contexts",
       x = NULL, y = "ΔHR (bpm)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_ctx_hrv <- data %>%
  ggplot(aes(x = listening_context, y = hrv_norm,
             fill = listening_context)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(values = context_colors,
                    guide = "none") +
  labs(title = "B — Normalised HRV across listening contexts",
       x = NULL, y = "ΔHRV (ms)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_fig(p_ctx_hr + p_ctx_hrv,
         "fig10_context_hr_hrv.png")

# =============================================================
# FIG 11: Heatmap — mean normalised HR and HRV per participant
# =============================================================

heatmap_data <- data %>%
  group_by(participant_id, listening_context) %>%
  summarise(hr_mean_norm  = mean(hr_norm,  na.rm = TRUE),
            hrv_mean_norm = mean(hrv_norm, na.rm = TRUE),
            .groups = "drop")

p_heat_hr <- heatmap_data %>%
  ggplot(aes(x = listening_context,
             y = participant_id,
             fill = hr_mean_norm)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low  = "#2166AC",
                       mid  = "white",
                       high = "#B2182B",
                       midpoint = 0,
                       name = "ΔHR (bpm)") +
  labs(title = "A — Mean normalised HR",
       x = NULL, y = "Participant") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid  = element_blank())

p_heat_hrv <- heatmap_data %>%
  ggplot(aes(x = listening_context,
             y = participant_id,
             fill = hrv_mean_norm)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low  = "#2166AC",
                       mid  = "white",
                       high = "#B2182B",
                       midpoint = 0,
                       name = "ΔHRV (ms)") +
  labs(title = "B — Mean normalised HRV",
       x = NULL, y = "Participant") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid  = element_blank())

save_fig(p_heat_hr + p_heat_hrv,
         "fig11_heatmap_participant_context.png",
         width = 14, height = 7)

# =============================================================
# FIG 12: Circadian variation of HR and HRV (median ± IQR)
# =============================================================

circadian_data <- data %>%
  mutate(hour = hour(completed_ts)) %>%
  group_by(hour) %>%
  summarise(
    hr_med  = median(hr_norm,  na.rm = TRUE),
    hr_q1   = quantile(hr_norm,  0.25, na.rm = TRUE),
    hr_q3   = quantile(hr_norm,  0.75, na.rm = TRUE),
    hrv_med = median(hrv_norm, na.rm = TRUE),
    hrv_q1  = quantile(hrv_norm, 0.25, na.rm = TRUE),
    hrv_q3  = quantile(hrv_norm, 0.75, na.rm = TRUE),
    n       = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 5)

p_circ_hr <- circadian_data %>%
  ggplot(aes(x = hour)) +
  geom_ribbon(aes(ymin = hr_q1, ymax = hr_q3),
              fill = "#185FA5", alpha = 0.15) +
  geom_line(aes(y = hr_med), color = "#185FA5",
            linewidth = 1.2) +
  geom_point(aes(y = hr_med, size = n),
             color = "#185FA5", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_x_continuous(breaks = seq(7, 22, by = 2),
                     labels = paste0(seq(7, 22, by = 2), ":00")) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(title = "A — Circadian variation of HR",
       x = NULL, y = "ΔHR (bpm)")

p_circ_hrv <- circadian_data %>%
  ggplot(aes(x = hour)) +
  geom_ribbon(aes(ymin = hrv_q1, ymax = hrv_q3),
              fill = "#0F6E56", alpha = 0.15) +
  geom_line(aes(y = hrv_med), color = "#0F6E56",
            linewidth = 1.2) +
  geom_point(aes(y = hrv_med, size = n),
             color = "#0F6E56", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_x_continuous(breaks = seq(7, 22, by = 2),
                     labels = paste0(seq(7, 22, by = 2), ":00")) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(title = "B — Circadian variation of HRV",
       x = "Hour of day", y = "ΔHRV (ms)")

save_fig(p_circ_hr / p_circ_hrv,
         "fig12_circadian_hr_hrv.png",
         width = 10, height = 8)

# =============================================================
# FIG 13: Morning vs Evening baseline — group level
# =============================================================

baseline_group <- data %>%
  filter(listening_context %in%
           c("Baseline Morning", "Baseline Evening")) %>%
  group_by(listening_context) %>%
  summarise(
    hr_mean  = mean(hr_norm,  na.rm = TRUE),
    hr_se    = sd(hr_norm,  na.rm = TRUE) / sqrt(n()),
    hrv_mean = mean(hrv_norm, na.rm = TRUE),
    hrv_se   = sd(hrv_norm, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

p_base_hr <- baseline_group %>%
  ggplot(aes(x = listening_context, y = hr_mean,
             fill = listening_context)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = hr_mean - hr_se,
                    ymax = hr_mean + hr_se),
                width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("Baseline Morning" = "#5F5E5A",
               "Baseline Evening" = "#3C3489"),
    guide  = "none"
  ) +
  labs(title = "A — HR", x = NULL, y = "ΔHR (bpm)")

p_base_hrv <- baseline_group %>%
  ggplot(aes(x = listening_context, y = hrv_mean,
             fill = listening_context)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = hrv_mean - hrv_se,
                    ymax = hrv_mean + hrv_se),
                width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("Baseline Morning" = "#5F5E5A",
               "Baseline Evening" = "#3C3489"),
    guide  = "none"
  ) +
  labs(title = "B — HRV", x = NULL, y = "ΔHRV (ms)")

save_fig(p_base_hr + p_base_hrv,
         "fig13_morning_evening_group.png",
         width = 8, height = 5)

# =============================================================
# FIG 14: Morning vs Evening — per participant (bar chart)
# =============================================================

baseline_ind <- data %>%
  filter(listening_context %in%
           c("Baseline Morning", "Baseline Evening")) %>%
  group_by(participant_id, listening_context) %>%
  summarise(hr_mean  = mean(hr_norm,  na.rm = TRUE),
            hrv_mean = mean(hrv_norm, na.rm = TRUE),
            .groups  = "drop")

p_ind_hr <- baseline_ind %>%
  ggplot(aes(x = participant_id, y = hr_mean,
             fill = listening_context)) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("Baseline Morning" = "#5F5E5A",
               "Baseline Evening" = "#3C3489"),
    name   = NULL
  ) +
  labs(title = "HR — Morning vs Evening per participant",
       x = "Participant", y = "ΔHR (bpm)") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 8),
        legend.position = "top")

p_ind_hrv <- baseline_ind %>%
  ggplot(aes(x = participant_id, y = hrv_mean,
             fill = listening_context)) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("Baseline Morning" = "#5F5E5A",
               "Baseline Evening" = "#3C3489"),
    name   = NULL
  ) +
  labs(title = "HRV — Morning vs Evening per participant",
       x = "Participant", y = "ΔHRV (ms)") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 8),
        legend.position = "none")

save_fig(p_ind_hr / p_ind_hrv,
         "fig14_morning_evening_individual.png",
         width = 12, height = 10)

# =============================================================
# FIG 15: Emotional valence × HR and HRV (line plot ± SE)
# =============================================================

valence_line <- data %>%
  filter(!is.na(emotional_state_INT)) %>%
  group_by(emotional_state_INT) %>%
  summarise(
    hr_mean  = mean(hr_norm,  na.rm = TRUE),
    hr_se    = sd(hr_norm,  na.rm = TRUE) / sqrt(n()),
    hrv_mean = mean(hrv_norm, na.rm = TRUE),
    hrv_se   = sd(hrv_norm, na.rm = TRUE) / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  )

p_val_hr <- valence_line %>%
  ggplot(aes(x = emotional_state_INT, y = hr_mean)) +
  geom_ribbon(aes(ymin = hr_mean - hr_se,
                  ymax = hr_mean + hr_se),
              fill = "#185FA5", alpha = 0.2) +
  geom_line(color = "#185FA5", linewidth = 1.2) +
  geom_point(aes(size = n), color = "#185FA5") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_x_continuous(breaks = 0:5,
                     labels = c("0\nVery\nFrustrated",
                                "1","2","3","4",
                                "5\nVery\nGood")) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(title = "A — HR across emotional valence levels",
       x = "Emotional state", y = "ΔHR (bpm)")

p_val_hrv <- valence_line %>%
  ggplot(aes(x = emotional_state_INT, y = hrv_mean)) +
  geom_ribbon(aes(ymin = hrv_mean - hrv_se,
                  ymax = hrv_mean + hrv_se),
              fill = "#0F6E56", alpha = 0.2) +
  geom_line(color = "#0F6E56", linewidth = 1.2) +
  geom_point(aes(size = n), color = "#0F6E56") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_x_continuous(breaks = 0:5,
                     labels = c("0\nVery\nFrustrated",
                                "1","2","3","4",
                                "5\nVery\nGood")) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(title = "B — HRV across emotional valence levels",
       x = "Emotional state", y = "ΔHRV (ms)")

save_fig(p_val_hr + p_val_hrv,
         "fig15_valence_line.png")

# =============================================================
# FIG 16: Emotional valence × context (boxplot)
# =============================================================

valence_ctx <- data %>%
  filter(!is.na(emotional_state_INT)) %>%
  mutate(valence_group = case_when(
    emotional_state_INT <= 1 ~ "Negative\n(0-1)",
    emotional_state_INT <= 3 ~ "Neutral\n(2-3)",
    TRUE                     ~ "Positive\n(4-5)"
  ),
  valence_group = factor(valence_group,
                         levels = c("Negative\n(0-1)",
                                    "Neutral\n(2-3)",
                                    "Positive\n(4-5)")))

p_val_ctx_hr <- valence_ctx %>%
  ggplot(aes(x = listening_context,
             y = hr_norm, fill = valence_group)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(
    values = c("Negative\n(0-1)" = "#E76F51",
               "Neutral\n(2-3)"  = "#E9C46A",
               "Positive\n(4-5)" = "#2A9D8F"),
    name = "Emotional valence"
  ) +
  labs(title = "A — HR by context and emotional valence",
       x = NULL, y = "ΔHR (bpm)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p_val_ctx_hrv <- valence_ctx %>%
  ggplot(aes(x = listening_context,
             y = hrv_norm, fill = valence_group)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(
    values = c("Negative\n(0-1)" = "#E76F51",
               "Neutral\n(2-3)"  = "#E9C46A",
               "Positive\n(4-5)" = "#2A9D8F"),
    name = "Emotional valence"
  ) +
  labs(title = "B — HRV by context and emotional valence",
       x = NULL, y = "ΔHRV (ms)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

save_fig(
  (p_val_ctx_hr + p_val_ctx_hrv) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  "fig16_valence_context_boxplot.png",
  width = 14, height = 6
)

# =============================================================
# FIG 17: Within-person correlations (valence × HR and HRV)
# =============================================================

within_corr <- data %>%
  filter(!is.na(emotional_state_INT)) %>%
  group_by(participant_id) %>%
  summarise(
    r_hr  = cor(emotional_state_INT, hr_norm,
                use = "complete.obs"),
    r_hrv = cor(emotional_state_INT, hrv_norm,
                use = "complete.obs"),
    .groups = "drop"
  ) %>%
  arrange(r_hr) %>%
  mutate(participant_id = factor(participant_id,
                                 levels = participant_id))

p_corr_hr <- within_corr %>%
  ggplot(aes(x = r_hr, y = participant_id,
             fill = r_hr > 0)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("TRUE" = "#2A9D8F", "FALSE" = "#E76F51"),
    guide  = "none"
  ) +
  labs(title = "A — Valence vs HR",
       x = "Within-person Pearson r", y = "Participant")

p_corr_hrv <- within_corr %>%
  arrange(r_hrv) %>%
  mutate(participant_id = factor(participant_id,
                                 levels = participant_id)) %>%
  ggplot(aes(x = r_hrv, y = participant_id,
             fill = r_hrv > 0)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("TRUE" = "#2A9D8F", "FALSE" = "#E76F51"),
    guide  = "none"
  ) +
  labs(title = "B — Valence vs HRV",
       x = "Within-person Pearson r", y = "Participant")

save_fig(p_corr_hr + p_corr_hrv,
         "fig17_within_person_correlations.png",
         width = 12, height = 7)

# =============================================================
# FIG 18: Cognitive fatigue × context (boxplot)
# =============================================================

fatigue_data <- data %>%
  filter(listening_context %in%
           c("Conversation", "Focused Listening",
             "Passive Listening"),
         !is.na(fatigue_cognitive_load_INT)) %>%
  mutate(
    fatigue_level = case_when(
      fatigue_cognitive_load_INT <= 1 ~ "Low (0-1)",
      fatigue_cognitive_load_INT <= 3 ~ "Medium (2-3)",
      TRUE                            ~ "High (4-5)"
    ),
    fatigue_level = factor(fatigue_level,
                           levels = c("Low (0-1)",
                                      "Medium (2-3)",
                                      "High (4-5)"))
  )

p_fat_hr <- fatigue_data %>%
  ggplot(aes(x = listening_context, y = hr_norm,
             fill = fatigue_level)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(
    values = c("Low (0-1)"    = "#2A9D8F",
               "Medium (2-3)" = "#E9C46A",
               "High (4-5)"   = "#E76F51"),
    name = "Cognitive fatigue"
  ) +
  labs(title = "A — HR by context and cognitive fatigue",
       x = NULL, y = "ΔHR (bpm)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p_fat_hrv <- fatigue_data %>%
  ggplot(aes(x = listening_context, y = hrv_norm,
             fill = fatigue_level)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(
    values = c("Low (0-1)"    = "#2A9D8F",
               "Medium (2-3)" = "#E9C46A",
               "High (4-5)"   = "#E76F51"),
    name = "Cognitive fatigue"
  ) +
  labs(title = "B — HRV by context and cognitive fatigue",
       x = NULL, y = "ΔHRV (ms)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

save_fig(
  (p_fat_hr + p_fat_hrv) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  "fig18_fatigue_context.png",
  width = 14, height = 6
)

# =============================================================
# FIG 32: EMA compliance per participant
# =============================================================

p_compliance <- data %>%
  count(participant_id) %>%
  ggplot(aes(x = reorder(participant_id, -n),
             y = n, fill = n)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
  scale_fill_gradient(low = "#E9C46A", high = "#E76F51") +
  geom_hline(yintercept = mean(data %>%
                                 count(participant_id) %>%
                                 pull(n)),
             linetype = "dashed", color = "gray40") +
  labs(title = "Number of valid EMA responses per participant",
       subtitle = "Dashed line = group mean",
       x = "Participant", y = "N responses") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 9))

save_fig(p_compliance,
         "fig32_ema_compliance.png",
         width = 10, height = 5)

# =============================================================
# FIG 33: EMA response density heatmap (day × hour)
# =============================================================

p_density <- data %>%
  mutate(
    day  = as.Date(completed_ts),
    hour = hour(completed_ts)
  ) %>%
  count(day, hour) %>%
  ggplot(aes(x = hour, y = day, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#E9C46A", high = "#E76F51",
                      name = "N EMA") +
  scale_x_continuous(breaks = seq(7, 22, by = 2),
                     labels = paste0(seq(7, 22, by = 2),
                                     ":00")) +
  labs(title = "EMA response density across days and hours",
       x = "Hour of day", y = "Date") +
  theme(panel.grid = element_blank())

save_fig(p_density,
         "fig33_ema_density_heatmap.png",
         width = 10, height = 7)

# =============================================================
# FIG 34: Acoustic perception × context (boxplot) — APPENDIX
# =============================================================

acoustic_data <- data %>%
  filter(!is.na(acoustic_perception_INT)) %>%
  mutate(
    acoustic_group = case_when(
      acoustic_perception_INT <= 1 ~ "Very Quiet\n(0-1)",
      acoustic_perception_INT <= 3 ~ "Moderate\n(2-3)",
      TRUE                         ~ "Very Noisy\n(4-5)"
    ),
    acoustic_group = factor(acoustic_group,
                            levels = c("Very Quiet\n(0-1)",
                                       "Moderate\n(2-3)",
                                       "Very Noisy\n(4-5)"))
  )

p_ac_hr <- acoustic_data %>%
  ggplot(aes(x = listening_context, y = hr_norm,
             fill = acoustic_group)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(
    values = c("Very Quiet\n(0-1)" = "#2A9D8F",
               "Moderate\n(2-3)"   = "#E9C46A",
               "Very Noisy\n(4-5)" = "#E76F51"),
    name = "Acoustic perception"
  ) +
  labs(title = "A — HR by context and acoustic perception",
       x = NULL, y = "ΔHR (bpm)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p_ac_hrv <- acoustic_data %>%
  ggplot(aes(x = listening_context, y = hrv_norm,
             fill = acoustic_group)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(
    values = c("Very Quiet\n(0-1)" = "#2A9D8F",
               "Moderate\n(2-3)"   = "#E9C46A",
               "Very Noisy\n(4-5)" = "#E76F51"),
    name = "Acoustic perception"
  ) +
  labs(title = "B — HRV by context and acoustic perception",
       x = NULL, y = "ΔHRV (ms)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

save_fig(
  (p_ac_hr + p_ac_hrv) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  "fig34_acoustic_context.png",
  width = 14, height = 6
)

# =============================================================
# FIG 35: Listening effort × context (boxplot) — APPENDIX
# =============================================================

effort_data <- data %>%
  filter(!is.na(listening_effort_INT)) %>%
  mutate(
    effort_group = case_when(
      listening_effort_INT <= 1 ~ "Low\n(0-1)",
      listening_effort_INT <= 3 ~ "Medium\n(2-3)",
      TRUE                      ~ "High\n(4-5)"
    ),
    effort_group = factor(effort_group,
                          levels = c("Low\n(0-1)",
                                     "Medium\n(2-3)",
                                     "High\n(4-5)"))
  )

p_eff_hr <- effort_data %>%
  ggplot(aes(x = listening_context, y = hr_norm,
             fill = effort_group)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(
    values = c("Low\n(0-1)"    = "#2A9D8F",
               "Medium\n(2-3)" = "#E9C46A",
               "High\n(4-5)"   = "#E76F51"),
    name = "Listening effort"
  ) +
  labs(title = "A — HR by context and listening effort",
       x = NULL, y = "ΔHR (bpm)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p_eff_hrv <- effort_data %>%
  ggplot(aes(x = listening_context, y = hrv_norm,
             fill = effort_group)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_fill_manual(
    values = c("Low\n(0-1)"    = "#2A9D8F",
               "Medium\n(2-3)" = "#E9C46A",
               "High\n(4-5)"   = "#E76F51"),
    name = "Listening effort"
  ) +
  labs(title = "B — HRV by context and listening effort",
       x = NULL, y = "ΔHRV (ms)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

save_fig(
  (p_eff_hr + p_eff_hrv) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  "fig35_effort_context.png",
  width = 14, height = 6
)

# =============================================================
# FIG 36: Cumulative exposure × HR and HRV at evening baseline
# APPENDIX
# =============================================================

cumulative_data <- data %>%
  filter(listening_context_label == "Baseline Evening",
         !is.na(cumulative_listening_exposure)) %>%
  mutate(cumexp = as.numeric(cumulative_listening_exposure))

p_cum_hr <- cumulative_data %>%
  ggplot(aes(x = factor(round(cumexp)), y = hr_norm)) +
  geom_boxplot(fill = "#5F5E5A", alpha = 0.7,
               outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  labs(title = "A — HR at Evening Baseline by cumulative exposure",
       x = "Cumulative listening exposure (0-5)",
       y = "ΔHR (bpm)")

p_cum_hrv <- cumulative_data %>%
  ggplot(aes(x = factor(round(cumexp)), y = hrv_norm)) +
  geom_boxplot(fill = "#3C3489", alpha = 0.7,
               outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  labs(title = "B — HRV at Evening Baseline by cumulative exposure",
       x = "Cumulative listening exposure (0-5)",
       y = "ΔHRV (ms)")

save_fig(p_cum_hr + p_cum_hrv,
         "fig36_cumulative_exposure.png",
         width = 12, height = 5)

# =============================================================
# FIG 37: Emmeans RQ1 — context (with Tukey arrows)
# APPENDIX
# =============================================================

emm_rq1_hr  <- emmeans(model_rq1_hr,  ~ listening_context)
emm_rq1_hrv <- emmeans(model_rq1_hrv, ~ listening_context)

p_emm_rq1_hr <- plot(emm_rq1_hr, comparisons = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50") +
  labs(title = "A — EMM HR by listening context",
       x = "ΔHR (bpm)", y = "Listening context")

p_emm_rq1_hrv <- plot(emm_rq1_hrv, comparisons = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50") +
  labs(title = "B — EMM HRV by listening context",
       x = "ΔHRV (ms)", y = "Listening context")

save_fig(p_emm_rq1_hr / p_emm_rq1_hrv,
         "fig37_emmeans_rq1.png",
         width = 8, height = 8)

# =============================================================
# FIG 38: Coefficient plot — extended model RQ1
# APPENDIX
# =============================================================

tidy_hr <- tidy(model_rq1_hr_ext,
                conf.int = TRUE,
                effects  = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = recode(term,
      "listening_contextFocused Listening" = "Focused Listening",
      "listening_contextPassive Listening" = "Passive Listening",
      "acoustic_perception_INT"            = "Acoustic perception",
      "fatigue_cognitive_load_INT"         = "Cognitive fatigue",
      "listening_effort_INT"               = "Listening effort",
      "emotional_state_INT"                = "Emotional valence"
    ),
    significant = p.value < 0.05,
    signal = "HR"
  )

tidy_hrv <- tidy(model_rq1_hrv_ext,
                 conf.int = TRUE,
                 effects  = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = recode(term,
      "listening_contextFocused Listening" = "Focused Listening",
      "listening_contextPassive Listening" = "Passive Listening",
      "acoustic_perception_INT"            = "Acoustic perception",
      "fatigue_cognitive_load_INT"         = "Cognitive fatigue",
      "listening_effort_INT"               = "Listening effort",
      "emotional_state_INT"                = "Emotional valence"
    ),
    significant = p.value < 0.05,
    signal = "HRV"
  )

coef_data <- bind_rows(tidy_hr, tidy_hrv)

p_coef <- coef_data %>%
  ggplot(aes(x = estimate,
             y = reorder(term, estimate),
             color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_color_manual(
    values = c("TRUE"  = "#E76F51",
               "FALSE" = "#264653"),
    labels = c("TRUE"  = "p < .05",
               "FALSE" = "p ≥ .05"),
    name   = NULL
  ) +
  facet_wrap(~ signal, scales = "free_x") +
  labs(title = "Fixed effect coefficients — RQ1 extended model",
       subtitle = "Non-baseline contexts only (N = 379). Reference: Conversation.",
       x = "β with 95% CI", y = NULL)

save_fig(p_coef,
         "fig38_coef_plot_rq1_ext.png",
         width = 12, height = 6)

# =============================================================
# FIG 39: Emmeans RQ2 — time of day (with Tukey arrows)
# APPENDIX
# =============================================================

emm_rq2_hr  <- emmeans(model_rq2_hr_cat,  ~ time_of_day)
emm_rq2_hrv <- emmeans(model_rq2_hrv_cat, ~ time_of_day)

p_emm_rq2_hr <- plot(emm_rq2_hr, comparisons = TRUE) +
  scale_y_discrete(limits = c("evening", "afternoon", "morning")) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50") +
  labs(title = "A — HR by time of day",
       x = "ΔHR (bpm)", y = "Time of day")

p_emm_rq2_hrv <- plot(emm_rq2_hrv, comparisons = TRUE) +
  scale_y_discrete(limits = c("evening", "afternoon", "morning")) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50") +
  labs(title = "B — HRV by time of day",
       x = "ΔHRV (ms)", y = "Time of day")

save_fig(p_emm_rq2_hr / p_emm_rq2_hrv,
         "fig39_emmeans_rq2.png",
         width = 8, height = 8)

# =============================================================
# FIG 40: Emmeans RQ3 — emotional state (with Tukey arrows)
# APPENDIX
# =============================================================

data_cat <- data %>%
  filter(!is.na(emotional_state_INT)) %>%
  mutate(
    valence_cat = case_when(
      emotional_state_INT <= 1 ~ "Negative\n(0-1)",
      emotional_state_INT <= 3 ~ "Neutral\n(2-3)",
      TRUE                     ~ "Positive\n(4-5)"
    ),
    valence_cat = factor(valence_cat,
                         levels = c("Negative\n(0-1)",
                                    "Neutral\n(2-3)",
                                    "Positive\n(4-5)"))
  )

model_rq3_hr_cat <- lmer(
  hr_norm ~ valence_cat + (1 | participant_id),
  data = data_cat, REML = FALSE
)
model_rq3_hrv_cat <- lmer(
  hrv_norm ~ valence_cat + (1 | participant_id),
  data = data_cat, REML = FALSE
)

emm_rq3_hr  <- emmeans(model_rq3_hr_cat,  ~ valence_cat)
emm_rq3_hrv <- emmeans(model_rq3_hrv_cat, ~ valence_cat)

p_emm_rq3_hr <- plot(emm_rq3_hr, comparisons = TRUE) +
  scale_y_discrete(limits = c("Positive\n(4-5)",
                               "Neutral\n(2-3)",
                               "Negative\n(0-1)")) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50") +
  labs(title = "A — HR by emotional state",
       x = "ΔHR (bpm)", y = "Emotional state")

p_emm_rq3_hrv <- plot(emm_rq3_hrv, comparisons = TRUE) +
  scale_y_discrete(limits = c("Positive\n(4-5)",
                               "Neutral\n(2-3)",
                               "Negative\n(0-1)")) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50") +
  labs(title = "B — HRV by emotional state",
       x = "ΔHRV (ms)", y = "Emotional state")

save_fig(p_emm_rq3_hr / p_emm_rq3_hrv,
         "fig40_emmeans_rq3.png",
         width = 8, height = 8)

# =============================================================
# Q-Q PLOTS — model diagnostics (APPENDIX)
# =============================================================

qq_plot <- function(model, title) {
  res <- residuals(model)
  ggplot(data.frame(res = res), aes(sample = res)) +
    stat_qq(color = "#264653", alpha = 0.5, size = 1) +
    stat_qq_line(color = "#E76F51", linewidth = 0.8) +
    labs(title = title,
         x = "Theoretical quantiles",
         y = "Sample quantiles") +
    theme(panel.grid.minor = element_blank())
}

fig_qq <- (
  qq_plot(model_rq1_hr,      "RQ1 — HR") |
  qq_plot(model_rq1_hrv,     "RQ1 — HRV")
) / (
  qq_plot(model_rq2_hr_cat,  "RQ2 — HR") |
  qq_plot(model_rq2_hrv_cat, "RQ2 — HRV")
) / (
  qq_plot(model_rq3_hr_add,  "RQ3 — HR") |
  qq_plot(model_rq3_hrv_add, "RQ3 — HRV")
) +
  plot_annotation(
    title = "Q-Q plots of model residuals",
    subtitle = "Points should follow the red line if residuals are normally distributed"
  )

save_fig(fig_qq,
         "fig_qq_plots.png",
         width = 10, height = 12)

# =============================================================
# HRV CIRCADIAN SPLINE — model-based (APPENDIX)
# =============================================================

pred_hrv_spline <- ggpredict(model_rq2_hrv_spline,
                             terms = "hour_cont [all]")

p_spline <- plot(pred_hrv_spline) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_x_continuous(
    breaks = seq(7, 23, by = 2),
    labels = paste0(seq(7, 23, by = 2), ":00")
  ) +
  labs(
    title = "HRV circadian pattern — natural cubic spline (model-based)",
    subtitle = "Shaded area = 95% CI",
    x = "Hour of day",
    y = "Predicted ΔHRV (ms)"
  )

save_fig(p_spline,
         "fig_circadian_spline_hrv.png",
         width = 8, height = 5)

message("All figures saved to figures/")
