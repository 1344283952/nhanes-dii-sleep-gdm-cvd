# DII × Sleep × GDM history → Cardiovascular Mortality (NHANES 2007-2018)

Reproducible analytic pipeline for the manuscript:

> **Pro-inflammatory diet, sleep disturbance, and cardiovascular mortality in US women with prior gestational diabetes: NHANES 2007-2018**
> Submitted to *Cardiovascular Diabetology*

## What this repository contains

This is the full reproducible R pipeline for our analysis of pro-inflammatory diet (Dietary Inflammatory Index, DII), sleep disturbance, and long-term mortality in **10,931 ever-pregnant US women** from six NHANES cycles (2007-2018), linked to the 2019 NCHS Public-Use Linked Mortality File.

- **Primary exposures**: DII (Shivappa 2014, 27 NHANES-available components), sleep disturbance composite (0-6 score)
- **Primary outcomes**: all-cause and cardiovascular mortality
- **Modifier**: prior gestational diabetes mellitus (GDM; RHQ162)
- **Methods**: survey-weighted Cox + logistic regression (R `survey`), CMAverse causal mediation (B=1000 bootstrap), restricted cubic splines, IPTW, E-values, subgroup interactions, multiple sensitivity analyses

## Key findings

- DII × all-cause mortality (M2 HR per 1-unit = 1.09, 95% CI 1.04-1.14; Q4 vs Q1 = 1.58, *P*-trend = 0.003)
- Sleep disorder × stroke (M2 OR = 1.99, 1.55-2.56)
- DII × education *P*~interaction~ = 0.002 (socioeconomic effect modification)
- Sleep does *not* mediate DII → mortality (PM = 0.6%, *P* = 0.41; B = 1000)
- GDM history does **not** modify DII or sleep effects (*P*~interaction~ = 0.90)

## Repository layout

```
scripts/                R analytic pipeline (run 00 → 16)
  00_install_packages.R    one-time install
  01_download_data.R       download NHANES + mortality
  02_merge_data.R          merge across 6 cycles
  03_clean_data.R          DII + sleep + GDM coding
  04_survey_design.R       svydesign setup
  05_table1.R              Table 1 baseline
  06_cox_main.R            Cox proportional hazards
  07_logistic_cvd.R        logistic regression for 6 CVD outcomes
  08_subgroup.R            stratified + interaction
  09_rcs.R                 restricted cubic splines
  11_mediation.R           CMAverse causal mediation
  12_sensitivity.R         S1-S6 sensitivity analyses
  13_predictive.R          C-index + E-value
  14_consort.R             CONSORT flow chart
  15_dag.R                 directed acyclic graph
  16_forest.R              forest plots
  run_all.R                end-to-end orchestrator
  _explore_N.R             pre-analysis sample-size exploration
  _make_delivery_v2.R      build submission zip
  _render_readme.R         pandoc md→docx batch
data/processed/         intermediate .RData files (start here to skip 01-03)
  nhanes_raw_merged.RData    raw merged across 6 cycles
  nhanes_final.RData         analytic sample N=10,931
  nhanes_design.RData        svydesign object
output/
  tables/               primary + supplementary tables (CSV + XLSX)
  figures/              CONSORT, DAG, RCS panels, forest plots
task.md                 project specification
references.bib          36 BibTeX entries (cited in manuscript)
```

## How to reproduce

Software: **R 4.6** + **Rtools 4.6** (for compiling packages); **pandoc** (for documentation rendering).

```r
# 1. Install packages (one-time, ~30-60 min including CMAverse from GitHub)
Rscript scripts/00_install_packages.R

# 2. Download raw NHANES + mortality data (~15 min)
Rscript scripts/01_download_data.R

# 3. Run the full pipeline (~20-30 min)
Rscript scripts/run_all.R
```

Or, to skip the 30-minute data download/merge step, **load the saved analytic sample directly**:

```r
load("data/processed/nhanes_final.RData")  # → nhanes_final, 10,931 × 1,008
load("data/processed/nhanes_design.RData") # → nhanes_design, svydesign object
# then run any scripts from 05_table1.R onward
```

## Reproducibility notes

- The full pipeline (data download → final figures) is deterministic; no randomized splits other than CMAverse bootstrap iterations.
- CMAverse uses B = 1000 bootstrap iterations with a fixed seed of `set.seed(20260513)` in `scripts/11_mediation.R` (reproducible to the 4th decimal).
- All survey-weighted analyses use the six-cycle dietary weight (WTDRD1 / 6) with `survey.lonely.psu = "adjust"`.
- The mediation analysis is unweighted (CMAverse limitation; documented in manuscript).

## Data

- **NHANES 2007-2018**: public domain. CDC/NCHS. https://wwwn.cdc.gov/nchs/nhanes/
- **NCHS Linked Mortality File**: public-use 2019 release. https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

This repository does **not** include the raw `.XPT` NHANES files. `scripts/01_download_data.R` will download them into `data/raw/` (which is `.gitignore`-d).

## Reproducibility

This repository contains the full analytic pipeline as submitted to *Cardiovascular Diabetology*. All scripts run end-to-end on a clean R 4.6 installation (after `00_install_packages.R`). The mediation analysis uses a fixed seed (`set.seed(20260513)`) for bootstrap reproducibility.

## Manuscript

See `manuscript.docx` (in the parent submission package).

## License

Code: MIT License (see `LICENSE`).
Data: NHANES is in the public domain (US government data).

## Citation

If you use this code, please cite:

> Li J, Sun X, Zhang J, Wang Z, Zhai L, Yu L. Pro-inflammatory diet, sleep disturbance, and cardiovascular mortality in ever-pregnant US women with and without prior gestational diabetes: a NHANES 2007–2018 cohort study. *Cardiovascular Diabetology*. (Submitted 2026; volume / pages / DOI to be assigned upon acceptance.)

## Contact

**Corresponding author**: Ling Yu, Department of Pharmacy, The Second Hospital of Jilin University, Changchun, Jilin Province, China.
(ORCID and institutional email are listed in the submission cover letter and will be confirmed at publication.)

**Co-authors**:
- Jie Li (first author) — Department of Obstetrics and Gynecology, The Second Hospital of Jilin University
- Xiubo Sun, Jing Zhang, Lijie Zhai — Department of Pharmacy, The Second Hospital of Jilin University
- Zhendong Wang — Beijing Union University

## AI tool disclosure

Per COPE 2025 guidance:
- AI-assisted tools were used for code generation and language editing of the Methods and Results sections.
- The Introduction and Discussion sections were drafted by the authors without AI-assisted text generation.
- All scientific decisions, interpretations, and conclusions are the responsibility of the authors.
