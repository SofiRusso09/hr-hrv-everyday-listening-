## 03_clean_physiology.R
## Apply physiological range filters and movement artefact removal
## Output: physio_processed.rds
## -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)

# --- Settings ---
HR_MIN        <- 40     # bpm
HR_MAX        <- 120    # bpm
HRV_MIN       <- 5      # ms
HRV_MAX       <- 200    # ms
STEPS_THRESH  <- 0.95   # percentile above which steps are considered high movement
HOUR_START    <- 7
HOUR_END      <- 23

physio_raw <- readRDS("physio_raw.rds")

physio_processed <- physio_raw %>%
  filter(
    between(hr_mean, HR_MIN, HR_MAX),
    between(hrv_proxy_rmssd, HRV_MIN, HRV_MAX),
    between(hour(window_start), HOUR_START, HOUR_END)
  ) %>%
  group_by(participant_id) %>%
  mutate(
    steps_threshold = quantile(steps, STEPS_THRESH, na.rm = TRUE),
    high_movement   = !is.na(steps) & steps > steps_threshold
  ) %>%
  ungroup() %>%
  filter(!high_movement) %>%
  select(-steps_threshold, -high_movement)

saveRDS(physio_processed, "physio_processed.rds")
message("Saved: physio_processed.rds — ", nrow(physio_processed), " rows remaining")
