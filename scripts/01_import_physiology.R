## 01_import_physiology.R
## Import raw PPG data, compute HR and proxy RMSSD over 5-minute windows
## Output: physio_raw.rds
## -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(RHRV)

# --- Settings ---
WINDOW_MIN  <- 5        # HRV window length in minutes
HR_MIN      <- 40       # Minimum plausible HR (bpm)
HR_MAX      <- 120      # Maximum plausible HR (bpm)
DATA_DIR    <- "data/raw_physiology"

# --- Helper: compute proxy RMSSD from interpolated HR signal ---
compute_proxy_rmssd <- function(x) {
  if (length(x) < 3) return(NA_real_)
  sqrt(mean(diff(diff(x))^2, na.rm = TRUE))
}

# --- Helper: process one participant file ---
process_participant <- function(file_path, participant_id) {

  hrv <- CreateHRVData()
  hrv <- SetVerbose(hrv, FALSE)
  hrv <- LoadBeatWFDB(hrv, file_path)   # adapt to your file format
  hrv <- BuildNIHR(hrv)
  hrv <- FilterNIHR(hrv)
  hrv <- InterpolateNIHR(hrv, freqhr = 4)

  # Split into 5-minute windows
  total_sec <- length(hrv$HR) / 4
  n_windows <- floor(total_sec / (WINDOW_MIN * 60))

  windows <- map_dfr(seq_len(n_windows), function(w) {
    start_idx <- (w - 1) * WINDOW_MIN * 60 * 4 + 1
    end_idx   <- w * WINDOW_MIN * 60 * 4
    segment   <- hrv$HR[start_idx:end_idx]

    tibble(
      participant_id  = participant_id,
      window_start    = hrv$datetime[start_idx],
      hr_mean         = mean(segment, na.rm = TRUE),
      hrv_proxy_rmssd = compute_proxy_rmssd(segment)
    )
  })

  windows
}

# --- Main ---
files <- list.files(DATA_DIR, pattern = "\\.dat$", full.names = TRUE)

physio_raw <- map_dfr(files, function(f) {
  pid <- tools::file_path_sans_ext(basename(f))
  message("Processing: ", pid)
  process_participant(f, pid)
}) %>%
  filter(
    between(hr_mean, HR_MIN, HR_MAX),
    !is.na(hrv_proxy_rmssd)
  )

saveRDS(physio_raw, "physio_raw.rds")
message("Saved: physio_raw.rds — ", nrow(physio_raw), " rows")
