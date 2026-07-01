# Coral wound-type healing & regeneration — analysis code and data

Analysis code and data for **Seifert, Brzezinski, Osenberg & Stier (2026)**,
*"Deep tissue removal in wounds facilitates algal colonization and inhibits
healing and regeneration in tropical corals."* See `CITATION.cff` for full
author, ORCID, and citation metadata.

## Overview

Two field experiments at Mo'orea, French Polynesia test how **wound type** —
tissue-only removal (*airbrush*) vs. tissue + skeletal damage (*scrape*, called
`dremel` in the data) — affects coral wound healing and regeneration.

- **Experiment 1 (3 genera: *Acropora*, *Pocillopora*, *Porites*)** — 22 paired
  split-colony fragments, 154 observations over 28 days. Outcomes: regeneration
  (polyps in wound center), healing (coenosarc closure), debris/algal coverage.
- **Experiment 2 (*Porites* spp.)** — 206 observations, 15 colonies, 3 treatments
  (airbrush / scrape / airbrush+scrape) over 63 days; up to ~10 microscope-scored
  wound-response metrics.
- **Experiments 3 & 4** are qualitative (single-polyp regeneration; histology)
  with **no dataset or code** — nothing to reproduce here; reported descriptively
  in the manuscript.

All reported statistics are **regenerated** by the pipeline; this README does
not hardcode them. After `make all`, read the numbers from:

- Exp 1: `output/text/paper_results_summary.md`, `output/tables/`
- Exp 2: `output/exp2_tables/PUBLICATION_STATISTICS_TABLE.csv`
- Figures: `output/figures/` (Exp 1), `output/exp2_figures_main/` (Exp 2)

## Methods (condensed)

Binary wound-response outcomes are modeled with **GLMMs**. Experiment 1 uses
`lme4` (`outcome ~ treatment + species + (1 | parent_id/coral_id)`) with a
**Firth penalized-likelihood** refit (`brglm2`) as a separation-robust
sensitivity, plus an ordinal cumulative-link mixed model for 3-level healing.
Experiment 2 uses `glmmTMB` (`outcome ~ treatment * ns(day, 3) + (1 | id)`);
outcomes exhibiting complete separation are fitted with a **penalized binomial
GLMM** (weakly-informative Normal prior on fixed effects; Gelman et al. 2008)
on the real 0/1 data. Experiment 2 outcomes were specified **a priori**, so
unadjusted *p*-values are primary; a family-wise correction table
(`output/exp2_tables/MULTIPLE_TESTING_CORRECTION.csv`) is a defensive
sensitivity. DHARMa residual diagnostics are seeded for determinism.

**Honest limitation:** the DHARMa dispersion test is flagged on most penalized
*Porites* models, but dispersion is not separately identifiable for Bernoulli
outcomes, so it is reported descriptively rather than used as an acceptance
gate; an observation-level random effect is not estimable at these per-cell
sample sizes. Inference rests on the penalized GLMM with this limitation stated.

## Reproduce

```bash
Rscript scripts/install_dependencies.R   # one-time
make all                                 # full pipeline (~1 min)
make verify                              # integrity self-check
```

Requires **R ≥ 4.3** (tested on R 4.5.2). Exact package versions used in the
reference run are recorded in `output/text/sessionInfo.txt`. Run all commands
from the repository root.

## Figures reproduced by this code

This repository regenerates the manuscript's **data figures**. After `make all`:

| Manuscript figure | Regenerated file | Script |
|---|---|---|
| Fig 2B — regeneration × species × wound type, with algal-prevalence overlay | `output/figures/regenerated_composition_with_debris_overlay_*.pdf` / `.png` | `airbrush_dremel_10_15_2025.R` |
| Fig 3 data panel — *Porites* wound-response outcomes over time by treatment | `output/exp2_figures_main/figure2_porites_all_outcomes.pdf` / `.png` | `exp2_02_create_figures.R` |
| Data supplement — *Porites* early-timepoint dynamics | `output/exp2_figures_supplement/figureS2_porites_early_dynamics.pdf` / `.png` | `exp2_05_create_early_timepoint_figure.R` |

Supporting data figures (regeneration trajectories, outcome time-series panels,
ordinal-healing variant) are also written to `output/figures/`. All remaining
manuscript figures — anatomy diagrams, in-situ wound photographs,
photomicrographs, RFP/native-GFP imaging, and histology — are **images or
Illustrator composites with no code source**, and are therefore not part of
this repository by design.

## Layout

```
data/      authoritative inputs + column dictionaries (data/README.md)
scripts/   analysis pipeline (Exp 1 + Exp 2), verify_pipeline.R, Makefile driver
output/    regenerated tables, figures, and diagnostics
```

## Terminology

`dremel` (data/code) = `scrape` (manuscript) — tissue + skeletal damage.
`airbrush` = tissue-only removal.

## License

Code and data released under the MIT License (see `LICENSE`).
