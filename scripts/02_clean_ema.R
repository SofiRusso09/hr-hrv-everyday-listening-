## 02_clean_ema.R
# Import and clean EMA survey data from individual participant CSVs
# Output: ema_clean.rds

library(tidyverse)
library(lubridate)
library(janitor)

# ---------------------------------------------------------
# 1. List EMA files
# ---------------------------------------------------------
# Expected folder structure:
#   EMA_data/
#     P001.csv
#     P003.csv
#     ...

ema_folder <- "EMA_data/"

ema_files <- list.files(
  path     = ema_folder,
  pattern  = "^P[0-9]{3}\\.csv$",
  full.names = TRUE
)

message("Found ", length(ema_files), " EMA files.")

# ---------------------------------------------------------
# 2. Cleaning function for a single participant file
# ---------------------------------------------------------

clean_single_ema <- function(file_path) {
  
  # Read everything as character to avoid type conflicts
  df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
    clean_names()
  
  # Extract participant ID from filename
  participant_id <- str_extract(basename(file_path), "P[0-9]{3}")
  
  # Column name aliases (clean_names() renames some columns):
  #   ACOUSTIC_PERCEPTION_   -> acoustic_perception_
  #   LISTENING_EFFORT_      -> listening_effort_
  #   FATIGUE/COGNITIVE_LOAD -> fatigue_cognitive_load
  col_map <- list(
    listening_context             = c("listening_context"),
    acoustic_perception           = c("acoustic_perception_", "acoustic_perception"),
    listening_effort              = c("listening_effort_", "listening_effort"),
    fatigue_cognitive_load        = c("fatigue_cognitive_load"),
    cumulative_listening_exposure = c("cumulative_listening_exposure"),
    emotional_state               = c("emotional_state"),
    listening_intent              = c("listening_intent")
  )
  
  # Resolve aliases: pick first existing column for each variable
  for (v in names(col_map)) {
    existing <- intersect(col_map[[v]], names(df))
    if (length(existing) > 0) {
      df[[v]] <- df[[existing[1]]]
    } else {
      df[[v]] <- NA_character_
    }
  }
  
  # Keep only required columns
  df <- df %>%
    select(
      completed_ts,
      listening_context,
      acoustic_perception,
      listening_effort,
      fatigue_cognitive_load,
      cumulative_listening_exposure,
      emotional_state,
      listening_intent
    )
  
  # Parse timestamp (format: dd/mm/yyyy HH:MM)
  df <- df %>%
    filter(!completed_ts %in% c("<no-response>", "<not-shown>", "", NA)) %>%
    filter(!is.na(completed_ts)) %>%
    mutate(completed_ts = dmy_hm(completed_ts, quiet = TRUE)) %>%
    filter(!is.na(completed_ts))
  
  # Replace survey placeholders with NA before numeric conversion
  placeholder_to_na <- function(x) {
    if_else(x %in% c("<no-response>", "<not-shown>", ""), NA_character_, x)
  }
  
  df <- df %>%
    mutate(across(
      c(listening_context, acoustic_perception, listening_effort,
        fatigue_cognitive_load, cumulative_listening_exposure,
        emotional_state, listening_intent),
      placeholder_to_na
    ))
  
  # Convert to numeric
  df <- df %>%
    mutate(across(
      c(listening_context, acoustic_perception, listening_effort,
        fatigue_cognitive_load, cumulative_listening_exposure,
        emotional_state, listening_intent),
      ~ suppressWarnings(as.numeric(.x))
    ))
  
  # Drop rows where all measures are NA
  df <- df %>%
    filter(if_any(-completed_ts, ~ !is.na(.x)))
  
  # Add context label
  df <- df %>%
    mutate(
      listening_context_label = case_when(
        listening_context == 1 ~ "Conversation",
        listening_context == 2 ~ "Focused Listening",
        listening_context == 3 ~ "Passive Listening",
        listening_context == 4 ~ "Baseline Morning",
        listening_context == 5 ~ "Baseline Evening",
        TRUE                   ~ NA_character_
      )
    ) %>%
    relocate(listening_context_label, .after = listening_context)
  
  # Add participant ID
  df <- df %>%
    mutate(participant_id = participant_id) %>%
    relocate(participant_id)
  
  return(df)
}

# ---------------------------------------------------------
# 3. Apply to all files and combine
# ---------------------------------------------------------

ema_clean <- map_dfr(ema_files, clean_single_ema)

message("EMA cleaning complete. ",
        nrow(ema_clean), " rows across ",
        n_distinct(ema_clean$participant_id), " participants.")

# ---------------------------------------------------------
# 4. Save
# ---------------------------------------------------------

write_rds(ema_clean, "ema_clean.rds")
write.csv(ema_clean, "ema_clean..csv", row.names = FALSE)
message("Saved: ema_clean.rds")


