## 06_plots.R
## Descriptive and model-based figures
## Output: figures/
## -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggeffects)
library(patchwork)
library(lme4)

dir.create("figures", showWarnings = FALSE)

data          <- readRDS("ema_physio_aligned.rds")
model_rq1_hr  <- readRDS("model_rq1_hr.rds")
model_rq1_hrv <- readRDS("model_rq1_hrv.rds")
model_rq2_hrv_cat    <- readRDS("model_rq2_hrv_cat.rds")
model_rq2_hrv_spline <- readRDS("model_rq2_hrv_spline.rds")
model_rq3_hr  <- readRDS("model_rq3_hr.rds")

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
# Figure: HR and HRV across listening contexts (boxplot)
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

ggsave("figures/fig_context_hr_hrv.png", p_hr + p_hrv, width = 12, height = 5, dpi = 150)

# =============================================================
# Figure: Circadian variation of HRV (spline)
# =============================================================
pred_hrv_spline <- ggpredict(model_rq2_hrv_spline, terms = "hour_cont [all]")

p_circ <- plot(pred_hrv_spline) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = seq(7, 23, by = 2),
                     labels = paste0(seq(7, 23, by = 2), ":00")) +
  labs(title = "Circadian variation of HRV (natural cubic spline)",
       x = "Hour of day", y = "Predicted HRV deviation (ms)")

ggsave("figures/fig_circadian_hrv.png", p_circ, width = 8, height = 5, dpi = 150)

# =============================================================
# Figure: Within-person correlations (valence × HR and HRV)
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
  scale_fill_manual(values = c("TRUE" = "#2A9D8F", "FALSE" = "#E76F51"), guide = "none") +
  labs(title = "Within-person correlations: emotional valence × HR / HRV",
       x = "Participant", y = "Pearson r") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("figures/fig_within_person_corr.png", p_corr, width = 10, height = 5, dpi = 150)

message("All figures saved to figures/")
