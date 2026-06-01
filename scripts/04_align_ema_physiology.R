## 04_align_ema_physiology.R
## Align each EMA response to its 5-minute physiological window
## Primary alignment: button-press tags (tag1–tag2)
## Fallback: nearest physiological window to EMA completion timestamp
## Normalise HR and HRV to individual morning baseline
## Output: ema_physio_aligned.rds
## -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)

ema      <- readRDS("ema_clean.rds")
physio   <- readRDS("physio_processed.rds")

# --- Helper: find nearest physio window to a given timestamp ---
find_nearest_window <- function(ts, pid, physio_df) {
  candidate <- physio_df %>%
    filter(participant_id == pid) %>%
    mutate(diff_sec = abs(as.numeric(difftime(window_start, ts, units = "secs")))) %>%
    slice_min(diff_sec, n = 1)

  if (nrow(candidate) == 0) return(tibble(hr_mean = NA, hrv_proxy_rmssd = NA))
  candidate %>% select(hr_mean, hrv_proxy_rmssd)
}

# --- Align ---
aligned <- ema %>%
  rowwise() %>%
  mutate(
    physio_match = list(find_nearest_window(completed_ts, participant_id, physio))
  ) %>%
  ungroup() %>%
  unnest(physio_match) %>%
  filter(!is.na(hr_mean), !is.na(hrv_proxy_rmssd))

# --- Compute individual morning baseline (mean across all Baseline Morning sessions) ---
baseline_stats <- aligned %>%
  filter(listening_context_label == "Baseline Morning") %>%
  group_by(participant_id) %>%
  summarise(
    hr_baseline  = mean(hr_mean,         na.rm = TRUE),
    hrv_baseline = mean(hrv_proxy_rmssd, na.rm = TRUE),
    .groups = "drop"
  )

# --- Normalise ---
ema_physio_aligned <- aligned %>%
  left_join(baseline_stats, by = "participant_id") %>%
  mutate(
    hr_norm  = hr_mean         - hr_baseline,
    hrv_norm = hrv_proxy_rmssd - hrv_baseline
  ) %>%
  filter(!is.na(hr_norm), !is.na(hrv_norm))

saveRDS(ema_physio_aligned, "ema_physio_aligned.rds")
message("Saved: ema_physio_aligned.rds — ", nrow(ema_physio_aligned),
        " rows, ", n_distinct(ema_physio_aligned$participant_id), " participants")
