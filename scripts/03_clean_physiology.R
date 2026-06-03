## 03_clean_physio.R
# Clean and filter the aggregated physiological data (df.pp) produced by 01_import_physio.R
# Output: physio_processed.rds, tags_clean.rds

library(tidyverse)
library(lubridate)

# ---------------------------------------------------------
# 1. Load raw physiological data
# ---------------------------------------------------------

df.pp   <- readRDS("physio_raw.rds")
df.tags <- readRDS("tags_raw.rds")

# ---------------------------------------------------------
# 2. Filter to daytime hours (7:00 - 23:00)
# ---------------------------------------------------------

df.pp <- df.pp %>%
  filter(hour(TimeR) >= 7,
         hour(TimeR) <= 23)

# ---------------------------------------------------------
# 3. Filter HR outside physiological range (40-160 bpm)
# ---------------------------------------------------------

df.pp <- df.pp %>%
  mutate(MeanHR = if_else(MeanHR < 40 | MeanHR > 160, NA_real_, MeanHR))

# ---------------------------------------------------------
# 4. Filter RMSSD outside plausible range (0-200 ms)
# ---------------------------------------------------------

df.pp <- df.pp %>%
  mutate(niHR = if_else(niHR < 0 | niHR > 200, NA_real_, niHR))

# ---------------------------------------------------------
# 5. Filter EDA outside physiological range (0.01-40 µS)
#    Typical range for wrist-worn sensor (e.g. Empatica E4)
#    Verify with summary(df.pp$eda) before applying
# ---------------------------------------------------------

df.pp <- df.pp %>%
  mutate(eda = if_else(eda < 0.01 | eda > 40, NA_real_, eda))

# ---------------------------------------------------------
# 6. Invalidate HR and RMSSD during intense movement
#    Steps > 95th percentile likely indicate motion artifacts
# ---------------------------------------------------------

steps_threshold <- quantile(df.pp$steps, 0.95, na.rm = TRUE)

df.pp <- df.pp %>%
  mutate(
    MeanHR = if_else(!is.na(steps) & steps > steps_threshold, NA_real_, MeanHR),
    niHR   = if_else(!is.na(steps) & steps > steps_threshold, NA_real_, niHR)
  )

# ---------------------------------------------------------
# 7. Remove windows with fewer than 2 valid physiological measures
# ---------------------------------------------------------

df.pp <- df.pp %>%
  mutate(
    n_valid = (!is.na(MeanHR)) + (!is.na(niHR)) + (!is.na(eda)) + (!is.na(steps))
  ) %>%
  filter(n_valid >= 2) %>%
  select(-n_valid)

# ---------------------------------------------------------
# 8. Clean tags: keep only daytime, sort by participant and time
# ---------------------------------------------------------

df.tags <- df.tags %>%
  filter(hour(Timestamps) >= 7,
         hour(Timestamps) <= 23) %>%
  arrange(ID, Timestamps)

message("Physio cleaning complete. ",
        nrow(df.pp), " 5-min windows retained across ",
        length(unique(df.pp$ID)), " participants.")

# ---------------------------------------------------------
# 9. Save
# ---------------------------------------------------------

saveRDS(df.pp,   "physio_processed.rds")
saveRDS(df.tags, "tags_clean.rds")
write.csv(df.pp,   "physio_processed.csv", row.names = FALSE)
write.csv(df.tags, "tags_clean.csv",        row.names = FALSE)
message("Saved CSV: physio_processed.csv, tags_clean.csv")
message("Saved: physio_processed.rds, tags_clean.rds")

