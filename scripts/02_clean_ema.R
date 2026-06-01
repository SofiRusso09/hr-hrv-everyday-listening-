## 02_clean_ema.R
## Standardise EMA exports, encode listening context labels
## Output: ema_clean.rds
## -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)

DATA_DIR <- "data/ema_exports"

# --- Context label mapping ---
CONTEXT_LABELS <- c(
  "1" = "Baseline Morning",
  "2" = "Conversation",
  "3" = "Focused Listening",
  "4" = "Passive Listening",
  "5" = "Baseline Evening"
)

# --- Load and combine EMA files ---
ema_raw <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE) %>%
  map_dfr(read_csv, show_col_types = FALSE)

# --- Clean ---
ema_clean <- ema_raw %>%
  rename_with(tolower) %>%
  rename_with(~ str_replace_all(., " ", "_")) %>%
  filter(
    !is.na(listening_context),
    listening_context != ""
  ) %>%
  mutate(
    completed_ts = parse_date_time(completed_ts, orders = c("ymd HMS", "dmy HMS")),
    listening_context_label = recode(
      as.character(listening_context),
      !!!CONTEXT_LABELS
    ),
    # Numeric EMA predictors (0-5 scale)
    emotional_state_INT        = as.numeric(emotional_state),
    fatigue_cognitive_load_INT = as.numeric(fatigue_cognitive_load),
    listening_effort_INT       = as.numeric(listening_effort),
    acoustic_perception_INT    = as.numeric(acoustic_perception),
    cumulative_listening_INT   = as.numeric(cumulative_listening_exposure)
  ) %>%
  filter(!is.na(completed_ts))

saveRDS(ema_clean, "ema_clean.rds")
message("Saved: ema_clean.rds — ", nrow(ema_clean), " rows, ",
        n_distinct(ema_clean$participant_id), " participants")
