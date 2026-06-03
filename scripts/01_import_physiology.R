## 01_import_physio.R
## Import raw PPG data, compute HR and proxy RMSSD over 5-minute windows
# Output: physio_raw.rds, tags_raw.rds

library(tidyverse)
library(lubridate)
library(plyr)
library(RHRV)

# ---------------------------------------------------------
# 1. Load raw CSVs for all participants
# ---------------------------------------------------------
# Expected folder structure:
#   processed_raw_data/
#     P001/
#       systolic_peaks.csv
#       steps.csv
#       tags.csv
#       eda.csv
#     P003/
#       ...

IDs <- dir("processed_raw_data/")

df.qpeaks <- df.steps <- df.tags <- df.eda <- NULL

for (i in IDs) {
  
  base <- paste0("processed_raw_data/", i, "/")
  
  if (file.exists(paste0(base, "systolic_peaks.csv"))) {
    tmp <- read.csv(paste0(base, "systolic_peaks.csv"))
    tmp$ID <- i
    df.qpeaks <- rbind(df.qpeaks, tmp)
  }
  
  if (file.exists(paste0(base, "steps.csv"))) {
    tmp <- read.csv(paste0(base, "steps.csv"))
    tmp$ID <- i
    df.steps <- rbind(df.steps, tmp)
  }
  
  if (file.exists(paste0(base, "tags.csv"))) {
    tmp <- read.csv(paste0(base, "tags.csv"))
    if (nrow(tmp) > 0) {
      tmp$ID <- i
      df.tags <- rbind(df.tags, tmp)
    }
  }
  
  if (file.exists(paste0(base, "eda.csv"))) {
    tmp <- read.csv(paste0(base, "eda.csv"))
    tmp$ID <- i
    df.eda <- rbind(df.eda, tmp)
  }
}

message("Raw CSVs loaded for ", length(IDs), " participants.")

# ---------------------------------------------------------
# 2. Convert timestamps and compute IBI
# ---------------------------------------------------------

df.qpeaks$systolic_peak_timestamp <- df.qpeaks$systolic_peak_timestamp / 1000000

# IBI in ms, filter physiologically implausible values
df.qpeaks$IBI <- c(0, diff(df.qpeaks$systolic_peak_timestamp))
df.qpeaks$IBI[df.qpeaks$IBI <= 400]  <- NA
df.qpeaks$IBI[df.qpeaks$IBI >= 1200] <- NA
df.qpeaks$dIBI <- c(0, diff(df.qpeaks$IBI))

# Convert all timestamps to POSIXct
df.qpeaks$Timestamps <- as.POSIXct(as.numeric(df.qpeaks$systolic_peak_timestamp) / 1000,
                                   origin = "1970-01-01", tz = "CET")
df.steps$Timestamps  <- as.POSIXct(as.numeric(df.steps$unix_timestamp) / 1000000,
                                   origin = "1970-01-01", tz = "CET")
df.tags$Timestamps   <- as.POSIXct(as.numeric(df.tags$tags_timestamp) / 1000000,
                                   origin = "1970-01-01", tz = "CET")
df.eda$Timestamps    <- as.POSIXct(as.numeric(df.eda$unix_timestamp) / 1000000,
                                   origin = "1970-01-01", tz = "CET")

# ---------------------------------------------------------
# 3. Compute HRV and aggregate in 5-minute windows
# ---------------------------------------------------------

IDs    <- unique(df.qpeaks$ID)
df.pp  <- NULL

for (i in IDs) {
  
  Dates <- unique(date(df.qpeaks$Timestamps[df.qpeaks$ID == i]))
  
  for (ii in Dates) {
    
    df.qpeaks.run <- df.qpeaks[date(df.qpeaks$Timestamps) == ii & df.qpeaks$ID == i, ]
    df.steps.run  <- df.steps[date(df.steps$Timestamps)   == ii & df.steps$ID  == i, ]
    df.eda.run    <- df.eda[date(df.eda$Timestamps)        == ii & df.eda$ID    == i, ]
    
    if (nrow(df.qpeaks.run) > 5) {
      
      df.qpeaks.run$beatPos <- as.numeric(df.qpeaks.run$Timestamps) -
        min(as.numeric(df.qpeaks.run$Timestamps))
      
      MD <- tryCatch({
        md <- RHRV::CreateHRVData()
        md <- RHRV::LoadBeatVector(md,
                                   df.qpeaks.run$beatPos,
                                   scale    = 1,
                                   datetime = format(df.qpeaks.run$Timestamps[1],
                                                     format = "%d/%m/%Y %H:%M:%S"))
        md <- RHRV::BuildNIHR(HRVData = md)
        md <- RHRV::InterpolateNIHR(HRVData = md)
        md
      }, error = function(e) {
        message("RHRV failed for ", i, " on ", ii, ": ", e$message)
        NULL
      })
      
      if (is.null(MD)) next
      
      df.run.run          <- data.frame(Time = MD$Beat$Time + MD$datetime,
                                        niHR = MD$Beat$niHR)
      df.run.run$TimeR    <- round_date(df.run.run$Time, "5 min")
      df.steps.run$TimeR  <- round_date(df.steps.run$Timestamps, "5 min")
      df.eda.run$TimeR    <- round_date(df.eda.run$Timestamps, "5 min")
      
      # Skip if no data after RHRV processing
      if (nrow(df.run.run) == 0) next
      
      # RMSSD proxy and mean HR per 5-min window
      # Need at least 3 points per window for diff(diff(x)) to work
      d.niHR <- tryCatch({
        tmp        <- aggregate(niHR ~ TimeR,
                                FUN  = function(x) {
                                  if (length(x) < 3) return(NA_real_)
                                  sqrt(mean(diff(diff(x))^2))
                                },
                                data = df.run.run)
        tmp$MeanHR <- aggregate(niHR ~ TimeR,
                                FUN  = function(x) mean(x, na.rm = TRUE),
                                data = df.run.run)$niHR
        tmp
      }, error = function(e) {
        message("Skipping aggregation for ", i, " on ", ii, ": ", e$message)
        NULL
      })
      
      if (is.null(d.niHR) || nrow(d.niHR) == 0) next
      
      # Aggregate steps and EDA with safety checks
      d <- d.niHR
      
      if (nrow(df.steps.run) > 5 && any(!is.na(df.steps.run$steps))) {
        d.steps <- tryCatch(
          aggregate(steps ~ TimeR, FUN = function(x) mean(x, na.rm = TRUE),
                    data = df.steps.run[!is.na(df.steps.run$steps), ]),
          error = function(e) NULL
        )
        if (!is.null(d.steps) && nrow(d.steps) > 0)
          d <- merge.data.frame(d, d.steps, by = "TimeR", all = TRUE)
      }
      
      if (nrow(df.eda.run) > 0 && any(!is.na(df.eda.run$eda))) {
        d.eda <- tryCatch(
          aggregate(eda ~ TimeR, FUN = function(x) mean(x, na.rm = TRUE),
                    data = df.eda.run[!is.na(df.eda.run$eda), ]),
          error = function(e) NULL
        )
        if (!is.null(d.eda) && nrow(d.eda) > 0)
          d <- merge.data.frame(d, d.eda, by = "TimeR", all = TRUE)
      }
      
      d$ID  <- i
      df.pp <- rbind.fill(df.pp, d)
    }
  }
}

df.pp$Date <- date(df.pp$TimeR)

# Filter to daytime only
df.pp <- df.pp[hour(df.pp$TimeR) > 7, ]

message("HRV processing complete. ", nrow(df.pp), " 5-min windows across ",
        length(unique(df.pp$ID)), " participants.")

# ---------------------------------------------------------
# 4. Save
# ---------------------------------------------------------

saveRDS(df.pp,   "physio_raw.rds")
saveRDS(df.tags, "tags_raw.rds")

write.csv(df.pp,   "physio_raw.csv",   row.names = FALSE)
write.csv(df.tags, "tags_raw.csv",     row.names = FALSE)
message("Saved: physio_raw.rds, tags_raw.rds")
