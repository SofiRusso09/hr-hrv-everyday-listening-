# HR and HRV as Physiological Indicators for Everyday Listening

**MSc Thesis — Sofia Russo**  
DTU (Technical University of Denmark)   
In collaboration with Eriksholm Research Center (Oticon)

---

## Overview

This repository contains the R analysis pipeline for the MSc thesis:

> *"Heart Rate and Heart Rate Variability as Physiological Indicators for Everyday Listening"*

The study investigates how HR and HRV vary across daily listening contexts (focused listening, conversation, passive listening, and resting baseline) using a combination of continuous wearable-based cardiac monitoring and Ecological Momentary Assessment (EMA).

**Study design:** 22 young normal-hearing adults, 7-day field study  
**Device:** Empatica EmbracePlus (PPG-based wearable)  
**EMA app:** SEMA3  
**N observations:** 690 valid EMA–physiology pairs

---

## Pipeline Structure

The analysis is organised into 6 modular scripts, each producing an `.rds` output used by the next stage.

```
01_import_physiology.R   →  physio_raw.rds
02_clean_ema.R           →  ema_clean.rds
03_clean_physiology.R    →  physio_processed.rds
04_align_ema_physiology.R →  ema_physio_aligned.rds
05_stats_models.R        →  model_*.rds
06_plots.R               →  figures/
```

---

## Scripts

| Script | Description |
|--------|-------------|
| `01_import_physiology.R` | Imports raw PPG data, extracts systolic peaks, computes HR and proxy RMSSD over 5-minute windows using the RHRV package |
| `02_clean_ema.R` | Standardises EMA column names, removes placeholder responses, encodes listening context labels |
| `03_clean_physiology.R` | Applies physiological range filters (HR: 40–120 bpm), removes high-movement artefacts, filters to 07:00–23:00 |
| `04_align_ema_physiology.R` | Aligns each EMA response to its corresponding physiological window using button-press tags; applies nearest-point fallback when tags are missing; normalises HR and HRV to individual morning baseline |
| `05_stats_models.R` | Fits Linear Mixed Models (LMMs) with random intercept per participant for RQ1 (listening context), RQ2 (time of day), and RQ3 (emotional valence); computes ICC, post-hoc Tukey contrasts, effect sizes, and R² |
| `06_plots.R` | Generates all descriptive and model-based figures |

---

## Dependencies

```r
install.packages(c(
  "tidyverse",
  "lubridate",
  "lme4",
  "lmerTest",
  "emmeans",
  "effectsize",
  "splines",
  "ggeffects",
  "performance",
  "car",
  "RHRV"
))
```

R version: 4.5.3

---

## Key Results

- Listening context significantly affects both HR (p = 0.014) and HRV (p = 0.005)
- Focused Listening associated with lower HR vs Passive Listening (β = −4.705 bpm, p = 0.007, d = 0.41)
- HRV shows significant circadian variation (p = 0.013); HR does not (p = 0.985)
- Emotional valence not significant at group level; emerges in additive model controlling for fatigue (p = 0.044)
- High inter-individual variability dominates the signal (ICC: HR = 0.154, HRV = 0.129)

---

## Data

Raw data are not included in this repository due to participant privacy.  
The pipeline expects the following input structure:

```
data/
├── raw_physiology/     # Per-participant PPG files from Empatica EmbracePlus
└── ema_exports/        # SEMA3 EMA export files (.csv)
```

---

## License

MIT License — free to use and adapt with attribution.

---


