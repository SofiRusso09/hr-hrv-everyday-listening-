## 04_align_ema_physio.R
library(tidyverse)
library(lubridate)
library(janitor)

# ---------------------------------------------------------
# 1. Load Data
# ---------------------------------------------------------
ema    <- read_rds("ema_clean.rds")
physio <- read_rds("physio_processed.rds")
tags   <- read_rds("tags_clean.rds")

# CRITICAL: convert EMA timestamps to CET to match physio/tags timezone
ema <- ema %>%
  mutate(completed_ts = lubridate::force_tz(completed_ts, "CET"))

# Rename columns
physio <- physio %>%
  dplyr::rename(
    participant_id = ID,
    timestamp      = TimeR,
    hr             = MeanHR,
    rmssd          = niHR
  )

tags <- tags %>%
  dplyr::rename(participant_id = ID)

# Convert to plain data.frames to avoid plyr/dplyr conflicts
physio <- as.data.frame(physio)
tags   <- as.data.frame(tags)

# ---------------------------------------------------------
# 2. Function: tag1-tag2 window + nearest point fallback
#    ALL filters use R base (no dplyr) to avoid plyr conflict
# ---------------------------------------------------------
# Protocol:
#   - Participant presses button (tag1)
#   - Waits ~5 min
#   - Presses again (tag2)
#   - Fills in EMA immediately after
#
# Constraints:
#   - tag1 and tag2 must be on the same day as completed_ts
#   - tag2 must be within 30 min before completed_ts
#   - gap tag1-tag2 must be <= 30 min
#   - If no valid pair → nearest physiological point

extract_physio_window <- function(ema_row, physio_df, tags_df) {
  
  pid      <- ema_row$participant_id
  ema_time <- ema_row$completed_ts  # already in CET
  ema_date <- as.Date(ema_time)
  
  # Step 1: get all tags for this participant on the same day
  # also restrict to tags before or within 2 min after completed_ts
  pid_mask     <- tags_df$participant_id == pid
  all_pid_tags <- tags_df[pid_mask, ]
  
  if (nrow(all_pid_tags) == 0) {
    tags_same_day <- all_pid_tags
  } else {
    day_mask      <- as.Date(all_pid_tags$Timestamps) == ema_date &
      all_pid_tags$Timestamps <= ema_time + 120
    tags_same_day <- all_pid_tags[day_mask, ]
    tags_same_day <- tags_same_day[order(tags_same_day$Timestamps), ]
  }
  
  # Step 2: find tag2 candidates (within 30 min before completed_ts)
  if (nrow(tags_same_day) > 0) {
    gap_mins  <- as.numeric(difftime(ema_time, tags_same_day$Timestamps,
                                     units = "mins"))
    # positive gap_mins = tag is before ema_time
    cand_mask      <- gap_mins >= -2 & gap_mins <= 30
    tags_tag2_cand <- tags_same_day[cand_mask, ]
  } else {
    tags_tag2_cand <- tags_same_day[0, ]
  }
  
  # Step 3: use tag1-tag2 window if valid pair found
  if (nrow(tags_tag2_cand) >= 1) {
    
    tag2 <- tags_tag2_cand$Timestamps[nrow(tags_tag2_cand)]
    
    # tag1 = tag immediately before tag2 on the same day
    before_mask      <- tags_same_day$Timestamps < tag2
    tags_before_tag2 <- tags_same_day[before_mask, ]
    
    if (nrow(tags_before_tag2) >= 1) {
      
      tag1         <- tags_before_tag2$Timestamps[nrow(tags_before_tag2)]
      gap_tags_min <- as.numeric(difftime(tag2, tag1, units = "mins"))
      
      if (gap_tags_min <= 30) {
        
        # Extract physio window between tag1 and tag2 (R base)
        win_mask      <- physio_df$participant_id == pid &
          physio_df$timestamp >= tag1 &
          physio_df$timestamp <= tag2
        physio_window <- physio_df[win_mask, ]
        
        if (nrow(physio_window) > 0) {
          return(data.frame(
            tag1                = tag1,
            tag2                = tag2,
            window_duration_sec = as.numeric(difftime(tag2, tag1,
                                                      units = "secs")),
            hr_mean    = mean(physio_window$hr,    na.rm = TRUE),
            hr_sd      = sd(physio_window$hr,      na.rm = TRUE),
            rmssd_mean = mean(physio_window$rmssd, na.rm = TRUE),
            eda_mean   = mean(physio_window$eda,   na.rm = TRUE),
            eda_sd     = sd(physio_window$eda,     na.rm = TRUE),
            steps_sum  = sum(physio_window$steps,  na.rm = TRUE),
            n_samples  = nrow(physio_window),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Step 4: fallback — nearest physiological timestamp (R base)
  pid_phys_mask <- physio_df$participant_id == pid
  physio_sub    <- physio_df[pid_phys_mask, ]
  
  if (nrow(physio_sub) == 0) {
    return(data.frame(
      tag1                = as.POSIXct(NA),
      tag2                = as.POSIXct(NA),
      window_duration_sec = NA_real_,
      hr_mean    = NA_real_, hr_sd      = NA_real_,
      rmssd_mean = NA_real_,
      eda_mean   = NA_real_, eda_sd     = NA_real_,
      steps_sum  = NA_real_, n_samples  = 0L,
      stringsAsFactors = FALSE
    ))
  }
  
  nearest_idx <- which.min(abs(physio_sub$timestamp - ema_time))
  nearest_row <- physio_sub[nearest_idx, ]
  
  data.frame(
    tag1                = as.POSIXct(NA),
    tag2                = as.POSIXct(NA),
    window_duration_sec = NA_real_,
    hr_mean    = nearest_row$hr,
    hr_sd      = NA_real_,
    rmssd_mean = nearest_row$rmssd,
    eda_mean   = nearest_row$eda,
    eda_sd     = NA_real_,
    steps_sum  = nearest_row$steps,
    n_samples  = 1L,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------
# 3. Apply to all EMA rows
# ---------------------------------------------------------
physio_features_list <- purrr::map(
  seq_len(nrow(ema)),
  function(i) extract_physio_window(ema[i, ], physio, tags)
)

aligned <- ema %>%
  bind_cols(bind_rows(physio_features_list))

message("Alignment complete: ", nrow(aligned), " rows, ",
        n_distinct(aligned$participant_id), " participants.")
message("  → tag1-tag2 windows: ",
        sum(!is.na(aligned$tag1)), " / ", nrow(aligned))
message("  → fallback nearest:  ",
        sum(is.na(aligned$tag1) & aligned$n_samples == 1), " / ", nrow(aligned))
message("  → no physio data:    ",
        sum(aligned$n_samples == 0), " / ", nrow(aligned))

# ---------------------------------------------------------
# 3b. Time gap tag2 → completed_ts
# ---------------------------------------------------------
aligned <- aligned %>%
  mutate(
    tag2_physio_gap_sec = case_when(
      !is.na(tag2) ~ as.numeric(difftime(completed_ts, tag2, units = "secs")),
      TRUE         ~ NA_real_
    ),
    tag2_same_session = !is.na(tag2_physio_gap_sec) &
      tag2_physio_gap_sec >= -120 &
      tag2_physio_gap_sec <= 1800
  )

valid_gaps <- aligned$tag2_physio_gap_sec[
  !is.na(aligned$tag2_physio_gap_sec) & aligned$tag2_same_session]
if (length(valid_gaps) > 0) {
  message("Time gap tag2 → completed_ts (valid sessions):")
  message("  Median: ", round(median(valid_gaps)), " sec (",
          round(median(valid_gaps) / 60, 1), " min)")
  message("  Max:    ", round(max(valid_gaps)), " sec (",
          round(max(valid_gaps) / 60, 1), " min)")
}

# ---------------------------------------------------------
# 4. Normalisation  — morning baseline per participant
# ---------------------------------------------------------
baseline_morning <- aligned %>%
  filter(listening_context_label == "Baseline Morning") %>%
  dplyr::group_by(participant_id) %>%
  dplyr::summarise(
    hr_baseline_morning  = mean(hr_mean,    na.rm = TRUE),
    hrv_baseline_morning = mean(rmssd_mean, na.rm = TRUE),
    .groups = "drop"
  )

aligned <- aligned %>%
  left_join(baseline_morning, by = "participant_id") %>%
  mutate(
    hr_norm  = hr_mean    - hr_baseline_morning,
    hrv_norm = rmssd_mean - hrv_baseline_morning
  )


# ---------------------------------------------------------
# 5. Save
# ---------------------------------------------------------
write_rds(aligned, "ema_physio_aligned.rds")
write.csv(aligned, "ema_physio_aligned.csv", row.names = FALSE)

message("Saved: ema_physio_aligned.rds and ema_physio_aligned.csv")
