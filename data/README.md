# Data Directory

Raw data for coral wound healing experiments (2022 and 2025).

## Experiment 1 (2022) — AUTHORITATIVE

### `airbrush_dremel.csv`
- **Source of record:** regenerated from
  `Wound_Type_exp1_2022-may16-2026-1846.xlsx` (**2026-05-16 independent re-score**;
  raw `.xlsx` lives in `data/`). Supersedes the 2026-04-21 re-score (migrated
  from the standalone 2022 wound-type repository).
- **Species:** *Acropora*, *Pocillopora*, *Porites*
- **N:** 22 fragments (paired split-colony), 154 rows, days 0,3,8,13,18,23,28
- **Treatments:** `airbrush` (tissue only), `dremel` (= scrape; tissue + skeleton)
- **Outcomes:** `regenerated` (yes/no — polyps in wound center; manuscript
  headline), `healed` (no/incomplete/yes — coenosarc/closure), `debris` (yes/no)
- **NA counts (after NA-backtrack re-score):** regenerated 0, healed 3, debris 4.
  Usable `regenerated` observations: 154 (**36 yes / 118 no**, 0 NA).
- **2026-05-16 re-score scope:** **21 cells changed across all three outcomes**
  (including `regenerated` calls in *Acropora* × Dremel and NA-backtrack
  corrections applied via t-1/t+1 neighbour rule across `healed`, `debris`,
  and `regenerated`).
- **Metadata sheet:** the `.xlsx` metadata sheet now documents **13 columns**
  (was 9), including the re-score definition note.
- Dictionary: `airbrush_dremel_metadata.csv`; sample map: `sample_key.csv`
- Prior vintages: `archive_exp1_2022/` (provenance only, not included in this snapshot). The immediately
  superseded 2026-04-21 data + metadata + sample_key are archived at
  `archive_exp1_2022/superseded_pre-may16-2026/`.

## Experiment 2 (2025) — Detailed Characterization

### `Porites microscope characterization - complete.xlsx`
- **Species:** *Porites* spp.
- **N colonies:** 15 (5 per treatment)
- **N observations:** 206
- **Time range:** Day 0 to Day 63
- **Treatments:** Airbrush, Dremel, Airbrush-Dremel
- **Outcomes:** 7 (coenosarc coverage, healed, polyps, algal plug, pink, RFP, yellow aggregations)

## Data Structure

### Experiment 2 (2025) Excel File

The Excel file contains:

**Sheet: "data"**
- `species` — Species name
- `date` — Observation date
- `day` — Days post-wounding (numeric)
- `treatment` — Treatment type (airbrush, dremel, air_drem)
- `id` — Coral colony ID (unique per colony)
- **Outcomes** (yes/no/incomplete):
  - `coenosarc_coverage` — Tissue regrowth over wound
  - `healed` — Complete healing status
  - `polyp_in_center_of_wound` — Polyp regeneration in wound center
  - `algal_plug` — Algal fouling presence
  - `pink` — Pink tissue (inflammation marker)
  - `rfp` — RFP expression (Porites only)
  - `yellow_aggregations` — Yellow granule clusters (Porites only)
- **Metadata:**
  - `photo_notes` — Notes on imaging
  - `microscope_notes` — Observation notes

**Sheet: "READ.ME"** (optional)
- Metadata and data dictionary

### Experiment 1 (2022) CSV File

`airbrush_dremel.csv` columns:
- `coral_id` — Unique coral identifier (e.g., "11a", "11b" for paired siblings)
- `parent_id` — Parent colony ID (extracted from coral_id)
- `species` — Species name (acropora, pocillopora, porites)
- `treatment` — Treatment type (airbrush, dremel)
- `day` — Days post-wounding
- `healed` — Healing status (yes, no, incomplete)
- `regenerated` — Polyps in wound center (yes/no — manuscript headline)
- `debris` — Debris present (yes/no)

## Treatments

1. **Airbrush** — Deep tissue removal using compressed air (no skeletal damage)
2. **Dremel** — Skeletal abrasion using rotary tool (retains basal tissue)
3. **Airbrush-Dremel (Combo)** — Airbrush followed by Dremel (both tissue removal + skeletal damage)

## Collection Details

- **Location:** Mo'orea, French Polynesia (2-4m depth)
- **Collection:** Colonies gently chiseled from pavement
- **Preparation:** Bases leveled with band saw, glued to labeled tiles
- **Acclimation:** 1 day in flow-through seawater (27-28°C, ~35 PSU)

## Imaging Protocol

- **Equipment:** SZX-10 Olympus stereomicroscope with fluorescence filters
- **Filters:** GFP (470nm), RFP (545nm)
- **Timing:** Day 0, daily D1-D3, then weekly until healed or study end
- **Settings:** Gain=0, exposure=100-300ms, LED intensity=50%

## Data Quality

- All wounds imaged at standardized magnifications
- Two magnifications: overview + high-resolution detail
- Corals imaged in seawater, immediately returned to tables
- Consistent imaging conditions across all timepoints

## Usage

Load data in R:

```r
library(readxl)
library(janitor)
library(tidyverse)

# Experiment 2 (2025) - Detailed characterization
porites <- read_excel("data/Porites microscope characterization - complete.xlsx",
                      sheet = "data") %>%
  clean_names()

# Experiment 1 (2022) - Authoritative re-score
exp1 <- read_csv("data/airbrush_dremel.csv")
```

---

**See main [README.md](../README.md) for analysis workflow**
