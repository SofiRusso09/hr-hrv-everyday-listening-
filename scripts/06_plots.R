## 06_plots.R
## Descriptive and model-based figures
## Output: figures/
## -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggeffects)
library(patchwork)
library(lme4)
library(lmerTest)

dir.create("figures", showWarnings = FALSE)

data                 <- readRDS("ema_physio_aligned.rds")
model_rq1_hr         <- readRDS("model_rq1_hr.rds")
model_rq1_hrv        <- readRDS("model_rq1_hrv.rds")
model_rq1_hr_ext     <- readRDS("model_rq1_hr_ext.rds")
model_rq2_hr_cat     <- readRDS("model_rq2_hr_cat.rds")
model_rq2_hrv_cat    <- readRDS("model_rq2_hrv_cat.rds")
model_rq2_hrv_spline <- readRDS("model_rq2_hrv_spline.rds")
model_rq3_hr         <- readRDS("model_rq3_hr.rds")
model_rq3_hrv        <- readRDS("model_rq3_hrv.rds")
model_rq3_hr_add     <- readRDS("model_rq3_hr_add.rds")
model_baseline_hrv   <- readRDS("model_baseline_hrv.rds")

# --- Colour palette ---
context_colors <- c(
  "Baseline Morning"  = "#E9C46A",
  "Conversation"      = "#E76F51",
  "Focused Listening" = "#264653",
  "Passive Listening" = "#2A9D8F",
  "Baseline Evening"  = "#F4A261"
)

theme_set(theme_minimal(base_size = 13))

# =============================================================
# Fig 1: Raw HR and HRV distribution per participant
# =============================================================
p_raw_hr <- data %>%
  group_by(participant_id) %>%
  summarise(hr_med = median(hr_mean, na.rm = TRUE),
            hr_q1  = quantile(hr_mean, 0.25, na.rm = TRUE),
            hr_q3  = quantile(hr_mean, 0.75, na.rm = TRUE)) %>%
  ggplot(aes(x = reorder(participant_id, hr_med), y = hr_med)) +
  geom_point(color = "#264653", size = 2) +
  geom_errorbar(aes(ymin = hr_q1, ymax = hr_q3), width = 0.3, color = "#264653") +
  labs(title = "Raw HR per participant", x = "Participant", y = "HR (bpm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_raw_hrv <- data %>%
  group_by(participant_id) %>%
  summarise(hrv_med = median(hrv_proxy_rmssd, na.rm = TRUE),
            hrv_q1  = quantile(hrv_proxy_rmssd, 0.25, na.rm = TRUE),
            hrv_q3  = quantile(hrv_proxy_rmssd, 0.75, na.rm = TRUE)) %>%
  ggplot(aes(x = reorder(participant_id, hrv_med), y = hrv_med)) +
  geom_point(color = "#E76F51", size = 2) +
  geom_errorbar(aes(ymin = hrv_q1, ymax = hrv_q3), width = 0.3, color = "#E76F51") +
  labs(title = "Raw HRV per participant", x = "Participant", y = "RMSSD proxy (ms)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("figures/fig_raw_per_participant.png",
       p_raw_hr + p_raw_hrv, width = 14, height = 5, dpi = 150)

# =============================================================
# Fig 2: Normalised HR and HRV across contexts (boxplot)
# =============================================================
p_hr <- data %>%
  filter(!is.na(listening_context_label)) %>%
  ggplot(aes(x = listening_context_label, y = hr_norm, fill = listening_context_label)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = context_colors, guide = "none") +
  labs(title = "Normalised HR across listening contexts",
       x = NULL, y = "HR deviation from morning baseline (bpm)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_hrv <- data %>%
  filter(!is.na(listening_context_label)) %>%
  ggplot(aes(x = listening_context_label, y = hrv_norm, fill = listening_context_label)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = context_colors, guide = "none") +
  labs(title = "Normalised HRV across listening contexts",
       x = NULL, y = "HRV deviation from morning baseline (ms)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("figures/fig_context_hr_hrv.png",
       p_hr + p_hrv, width = 12, height = 5, dpi = 150)

# =============================================================
# Fig 3: Heatmap — mean normalised HR per participant × context
# =============================================================
heatmap_data <- data %>%
  filter(!is.na(listening_context_label)) %>%
  group_by(participant_id, listening_context_label) %>%
  summarise(hr_mean_norm  = mean(hr_norm,  na.rm = TRUE),
            hrv_mean_norm = mean(hrv_norm, na.rm = TRUE),
            .groups = "drop")

p_heat_hr <- heatmap_data %>%
  ggplot(aes(x = listening_context_label, y = participant_id, fill = hr_mean_norm)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2A9D8F", mid = "white", high = "#E76F51",
                       midpoint = 0, name = "ΔHR (bpm)") +
  labs(title = "Mean normalised HR — participant × context",
       x = NULL, y = "Participant") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_heat_hrv <- heatmap_data %>%
  ggplot(aes(x = listening_context_label, y = participant_id, fill = hrv_mean_norm)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#E76F51", mid = "white", high = "#264653",
                       midpoint = 0, name = "ΔHRV (ms)") +
  labs(title = "Mean normalised HRV — participant × context",
       x = NULL, y = "Participant") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("figures/fig_heatmap_context.png",
       p_heat_hr + p_heat_hrv, width = 14, height = 7, dpi = 150)

# =============================================================
# Fig 4: Circadian plot — HR and HRV median ± IQR per hour
# =============================================================
circadian_data <- data %>%
  filter(!is.na(hr_norm), !is.na(hrv_norm)) %>%
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
  )

p_circ_hr <- circadian_data %>%
  ggplot(aes(x = hour)) +
  geom_ribbon(aes(ymin = hr_q1, ymax = hr_q3), fill = "#264653", alpha = 0.2) +
  geom_line(aes(y = hr_med), color = "#264653", linewidth = 1) +
  geom_point(aes(y = hr_med, size = n), color = "#264653") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = seq(7, 23, by = 2),
                     labels = paste0(seq(7, 23, by = 2), ":00")) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(title = "A — Circadian variation of HR",
       x = "Hour of day", y = "ΔHR (bpm)")

p_circ_hrv <- circadian_data %>%
  ggplot(aes(x = hour)) +
  geom_ribbon(aes(ymin = hrv_q1, ymax = hrv_q3), fill = "#E76F51", alpha = 0.2) +
  geom_line(aes(y = hrv_med), color = "#E76F51", linewidth = 1) +
  geom_point(aes(y = hrv_med, size = n), color = "#E76F51") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = seq(7, 23, by = 2),
                     labels = paste0(seq(7, 23, by = 2), ":00")) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(title = "B — Circadian variation of HRV",
       x = "Hour of day", y = "ΔHRV (ms)")

ggsave("figures/fig_circadian_hr_hrv.png",
       p_circ_hr + p_circ_hrv, width = 12, height = 5, dpi = 150)

# =============================================================
# Fig 5: Morning vs Evening baseline per participant
# =============================================================
baseline_data <- data %>%
  filter(listening_context_label %in% c("Baseline Morning", "Baseline Evening")) %>%
  group_by(participant_id, listening_context_label) %>%
  summarise(hr_mean_norm  = mean(hr_norm,  na.rm = TRUE),
            hrv_mean_norm = mean(hrv_norm, na.rm = TRUE),
            .groups = "drop")

p_bm_hr <- baseline_data %>%
  ggplot(aes(x = participant_id, y = hr_mean_norm,
             fill = listening_context_label)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Baseline Morning" = "#E9C46A",
                               "Baseline Evening" = "#F4A261"),
                    name = NULL) +
  labs(title = "HR — Morning vs Evening Baseline per participant",
       x = "Participant", y = "ΔHR (bpm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "top")

p_bm_hrv <- baseline_data %>%
  ggplot(aes(x = participant_id, y = hrv_mean_norm,
             fill = listening_context_label)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Baseline Morning" = "#E9C46A",
                               "Baseline Evening" = "#F4A261"),
                    name = NULL) +
  labs(title = "HRV — Morning vs Evening Baseline per participant",
       x = "Participant", y = "ΔHRV (ms)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "top")

ggsave("figures/fig_morning_evening_per_participant.png",
       p_bm_hr / p_bm_hrv, width = 12, height = 10, dpi = 150)

# =============================================================
# Fig 6: Emotional valence × HR and HRV (line plot)
# =============================================================
valence_data <- data %>%
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

p_val_hr <- valence_data %>%
  ggplot(aes(x = emotional_state_INT, y = hr_mean)) +
  geom_ribbon(aes(ymin = hr_mean - hr_se, ymax = hr_mean + hr_se),
              fill = "#264653", alpha = 0.2) +
  geom_line(color = "#264653", linewidth = 1) +
  geom_point(aes(size = n), color = "#264653") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = 0:5,
                     labels = c("0\nVery\nFrustrated","1","2","3","4","5\nVery\nGood")) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(title = "A — HR across emotional valence levels",
       x = "Emotional state", y = "ΔHR (bpm)")

p_val_hrv <- valence_data %>%
  ggplot(aes(x = emotional_state_INT, y = hrv_mean)) +
  geom_ribbon(aes(ymin = hrv_mean - hrv_se, ymax = hrv_mean + hrv_se),
              fill = "#E76F51", alpha = 0.2) +
  geom_line(color = "#E76F51", linewidth = 1) +
  geom_point(aes(size = n), color = "#E76F51") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = 0:5,
                     labels = c("0\nVery\nFrustrated","1","2","3","4","5\nVery\nGood")) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(title = "B — HRV across emotional valence levels",
       x = "Emotional state", y = "ΔHRV (ms)")

ggsave("figures/fig_valence_hr_hrv.png",
       p_val_hr + p_val_hrv, width = 12, height = 5, dpi = 150)

# =============================================================
# Fig 7: Within-person correlations (valence × HR and HRV)
# =============================================================
within_person_corr <- data %>%
  filter(!is.na(emotional_state_INT)) %>%
  group_by(participant_id) %>%
  summarise(
    r_hr  = cor(emotional_state_INT, hr_norm,  use = "complete.obs"),
    r_hrv = cor(emotional_state_INT, hrv_norm, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  pivot_longer(c(r_hr, r_hrv), names_to = "signal", values_to = "r") %>%
  mutate(signal = recode(signal, r_hr = "HR", r_hrv = "HRV"))

p_corr <- within_person_corr %>%
  ggplot(aes(x = reorder(participant_id, r), y = r, fill = r > 0)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ signal) +
  scale_fill_manual(values = c("TRUE" = "#2A9D8F", "FALSE" = "#E76F51"),
                    guide = "none") +
  labs(title = "Within-person correlations: emotional valence × HR / HRV",
       x = "Participant", y = "Pearson r") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("figures/fig_within_person_corr.png",
       p_corr, width = 10, height = 5, dpi = 150)

# =============================================================
# Fig 8: Cognitive fatigue × context (boxplot)
# =============================================================
fatigue_data <- data %>%
  filter(listening_context_label %in%
           c("Conversation", "Focused Listening", "Passive Listening"),
         !is.na(fatigue_cognitive_load_INT)) %>%
  mutate(fatigue_level = case_when(
    fatigue_cognitive_load_INT <= 1 ~ "Low (0-1)",
    fatigue_cognitive_load_INT <= 3 ~ "Medium (2-3)",
    TRUE                            ~ "High (4-5)"
  ),
  fatigue_level = factor(fatigue_level,
                         levels = c("Low (0-1)", "Medium (2-3)", "High (4-5)")))

p_fat_hr <- fatigue_data %>%
  ggplot(aes(x = listening_context_label, y = hr_norm, fill = fatigue_level)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Low (0-1)"    = "#2A9D8F",
                               "Medium (2-3)" = "#E9C46A",
                               "High (4-5)"   = "#E76F51"),
                    name = "Fatigue level") +
  labs(title = "A — HR by context and cognitive fatigue",
       x = NULL, y = "ΔHR (bpm)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p_fat_hrv <- fatigue_data %>%
  ggplot(aes(x = listening_context_label, y = hrv_norm, fill = fatigue_level)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Low (0-1)"    = "#2A9D8F",
                               "Medium (2-3)" = "#E9C46A",
                               "High (4-5)"   = "#E76F51"),
                    name = "Fatigue level") +
  labs(title = "B — HRV by context and cognitive fatigue",
       x = NULL, y = "ΔHRV (ms)") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave("figures/fig_fatigue_context.png",
       p_fat_hr + p_fat_hrv, width = 12, height = 5, dpi = 150)

# =============================================================
# Fig 9: HRV circadian spline (model-based)
# =============================================================
pred_hrv_spline <- ggpredict(model_rq2_hrv_spline, terms = "hour_cont [all]")

p_circ_spline <- plot(pred_hrv_spline) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = seq(7, 23, by = 2),
                     labels = paste0(seq(7, 23, by = 2), ":00")) +
  labs(title = "Circadian variation of HRV (natural cubic spline — model-based)",
       x = "Hour of day", y = "Predicted ΔHRV (ms)")

ggsave("figures/fig_circadian_spline.png",
       p_circ_spline, width = 8, height = 5, dpi = 150)

message("All figures saved to figures/")
