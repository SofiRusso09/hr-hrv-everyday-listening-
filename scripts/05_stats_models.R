## 05_stats_models.R
## Linear Mixed Models for RQ1, RQ2, RQ3
## All models include a random intercept per participant
## Outcomes: hr_norm (bpm), hrv_norm (ms) — normalised to morning baseline
## -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)
library(splines)
library(performance)
library(car)

# --- Load data ---
raw <- read_rds("ema_physio_aligned.rds") %>%
  filter(between(hour(completed_ts), 7, 23)) %>%
  filter(!is.na(hr_norm), !is.na(hrv_norm))

data <- raw %>%
  mutate(
    participant_id    = factor(participant_id),
    listening_context = factor(listening_context_label,
                               levels = c("Baseline Evening",
                                          "Baseline Morning",
                                          "Conversation",
                                          "Focused Listening",
                                          "Passive Listening")),
    time_of_day = case_when(
      hour(completed_ts) < 12 ~ "morning",
      hour(completed_ts) < 18 ~ "afternoon",
      TRUE                    ~ "evening"
    ),
    time_of_day = factor(time_of_day, levels = c("morning", "afternoon", "evening")),
    hour_cont   = hour(completed_ts) + minute(completed_ts) / 60
  )

# Non-baseline subset (for extended RQ1 model)
data_nb <- data %>%
  filter(listening_context %in% c("Conversation", "Focused Listening", "Passive Listening")) %>%
  droplevels()

# Baseline-only subset (for Morning vs Evening comparison)
data_baseline <- data %>%
  filter(listening_context %in% c("Baseline Morning", "Baseline Evening")) %>%
  droplevels()

# =============================================================
# ICC — null models
# =============================================================
model_null_hr  <- lmer(hr_norm  ~ 1 + (1 | participant_id), data = data, REML = TRUE)
model_null_hrv <- lmer(hrv_norm ~ 1 + (1 | participant_id), data = data, REML = TRUE)

message("ICC HR:  ", round(icc(model_null_hr)$ICC_adjusted,  3))
message("ICC HRV: ", round(icc(model_null_hrv)$ICC_adjusted, 3))

# =============================================================
# RQ1 — Listening context
# =============================================================
model_rq1_hr  <- lmer(hr_norm  ~ listening_context + (1 | participant_id), data = data, REML = FALSE)
model_rq1_hrv <- lmer(hrv_norm ~ listening_context + (1 | participant_id), data = data, REML = FALSE)

# F-tests
print(anova(model_rq1_hr,  type = 3))
print(anova(model_rq1_hrv, type = 3))

# Post-hoc Tukey
emm_hr  <- emmeans(model_rq1_hr,  ~ listening_context)
emm_hrv <- emmeans(model_rq1_hrv, ~ listening_context)
print(pairs(emm_hr,  adjust = "tukey"))
print(pairs(emm_hrv, adjust = "tukey"))

# Effect sizes
print(eff_size(emm_hr,  sigma = sigma(model_rq1_hr),  edf = df.residual(model_rq1_hr)))
print(eff_size(emm_hrv, sigma = sigma(model_rq1_hrv), edf = df.residual(model_rq1_hrv)))

# Extended model (non-baseline, EMA covariates)
model_rq1_hr_ext <- lmer(
  hr_norm ~ listening_context +
    acoustic_perception_INT + fatigue_cognitive_load_INT +
    listening_effort_INT + emotional_state_INT +
    (1 | participant_id),
  data = data_nb, REML = FALSE
)

print(anova(model_rq1_hr_ext, type = 3))
print(vif(model_rq1_hr_ext))

# =============================================================
# RQ2 — Time of day
# =============================================================
model_rq2_hr_cat  <- lmer(hr_norm  ~ time_of_day + (1 | participant_id), data = data, REML = FALSE)
model_rq2_hrv_cat <- lmer(hrv_norm ~ time_of_day + (1 | participant_id), data = data, REML = FALSE)

print(anova(model_rq2_hr_cat,  type = 3))
print(anova(model_rq2_hrv_cat, type = 3))

# Natural cubic spline (robustness check)
model_rq2_hr_spline  <- lmer(hr_norm  ~ ns(hour_cont, 3) + (1 | participant_id), data = data, REML = FALSE)
model_rq2_hrv_spline <- lmer(hrv_norm ~ ns(hour_cont, 3) + (1 | participant_id), data = data, REML = FALSE)

print(anova(model_rq2_hr_spline,  lmer(hr_norm  ~ 1 + (1 | participant_id), data = data, REML = FALSE)))
print(anova(model_rq2_hrv_spline, lmer(hrv_norm ~ 1 + (1 | participant_id), data = data, REML = FALSE)))

# Morning vs Evening baseline (formal test)
model_baseline_hr  <- lmer(hr_norm  ~ listening_context + (1 | participant_id), data = data_baseline, REML = FALSE)
model_baseline_hrv <- lmer(hrv_norm ~ listening_context + (1 | participant_id), data = data_baseline, REML = FALSE)

print(summary(model_baseline_hr))
print(summary(model_baseline_hrv))

# =============================================================
# RQ3 — Emotional valence
# =============================================================
model_rq3_hr  <- lmer(hr_norm  ~ emotional_state_INT + (1 | participant_id), data = data, REML = FALSE)
model_rq3_hrv <- lmer(hrv_norm ~ emotional_state_INT + (1 | participant_id), data = data, REML = FALSE)

# Additive model (valence + fatigue)
model_rq3_hr_add  <- lmer(hr_norm  ~ emotional_state_INT + fatigue_cognitive_load_INT + (1 | participant_id), data = data, REML = FALSE)
model_rq3_hrv_add <- lmer(hrv_norm ~ emotional_state_INT + fatigue_cognitive_load_INT + (1 | participant_id), data = data, REML = FALSE)

print(anova(model_rq3_hr_add,  type = 3))
print(anova(model_rq3_hrv_add, type = 3))



# =============================================================
# Save models
# =============================================================
saveRDS(model_rq1_hr,         "model_rq1_hr.rds")
saveRDS(model_rq1_hrv,        "model_rq1_hrv.rds")
saveRDS(model_rq1_hr_ext,     "model_rq1_hr_ext.rds")
saveRDS(model_rq2_hr_cat,     "model_rq2_hr_cat.rds")
saveRDS(model_rq2_hrv_cat,    "model_rq2_hrv_cat.rds")
saveRDS(model_rq2_hr_spline,  "model_rq2_hr_spline.rds")
saveRDS(model_rq2_hrv_spline, "model_rq2_hrv_spline.rds")
saveRDS(model_rq3_hr,         "model_rq3_hr.rds")
saveRDS(model_rq3_hrv,        "model_rq3_hrv.rds")
saveRDS(model_rq3_hr_add,     "model_rq3_hr_add.rds")
saveRDS(model_rq3_hrv_add,    "model_rq3_hrv_add.rds")


message("All models saved.")
