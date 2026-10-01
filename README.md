# heart-capillary-sex-disease

Do capillary endothelial cells from female and male **human** hearts differ in function-related gene programs **in heart disease**?
This is paper 2 of a two-part project. Paper 1 (healthy hearts) is in
[heart-capillary-sex-differences](https://github.com/mladelus/heart-capillary-sex-differences).

The analysis plan (`ANALYSIS_PLAN.md`) was committed **before** any expression data were analysed; the commit date is its time-stamp.

## Data
Open data only, with no application or registration: the CZ CELLxGENE Census (release 2025-11-08), adult ventricular donors
with dilated cardiomyopathy, arrhythmogenic cardiomyopathy, myocardial infarction or myocarditis.

## Analyses
| Script | What it does |
|---|---|
| `R/00_setup.R` | packages, settings, gene programs, shared functions |
| `R/01_inventory.R` | eligible donors and disease strata (metadata only) |
| `R/02_download.R` | downloads counts for included donors |
| `R/02b_download_mural.R` | downloads mural cells (pericytes + smooth muscle), one rule for all datasets |
| `R/03_ec_subtypes_pseudobulk.R` | endothelial subtypes; donor-level pseudobulks |
| `R/04_primary_sex_within_disease.R` | **primary**: female - male within disease, pooled; equivalence ±0.8 SD; permutations |
| `R/05_sex_by_disease.R` | sex × DCM interaction (same consortium) and cross-study comparison with paper 1 |
| `R/06_pericytes_xy.R` | mural cells (pericytes) and X–Y paralog dosage in disease |
| `R/07_figures.R` | figures |
| `R/08_where_sex_differences.R` | Stage 3 (exploratory): discovery in one cohort, replication in the other, by cell type; Hallmark pathways |
| `R/09_hormone_receptors.R` | Stage 3: sex-hormone receptor map; estrogen/androgen-response scores in capillaries |
| `R/10_capillary_states_variability.R` | Stage 3: capillary cell states; donor-to-donor variability |
| `R/11_combined_and_figures.R` | Health + disease combined estimates; stage 3 figures |
| `R/12_immune_and_ambient.R` | Stage 3: immune-cell programs and composition; estrogen response by cell type; ambient-RNA check |
| `R/13_cardiomyocyte_pathways_network.R` | Stage 3: pathways and STRING protein network for cardiomyocyte sex-biased genes; endothelin genes in capillaries |

## How to run
1. Install R ≥ 4.3 and the packages listed in `R/00_setup.R` (the script installs missing ones).
2. Run paper 1 first; its `pb.rds`, `aim1.rds` and `aim4.rds` must be in `~/Documents/heart_capillary_project` (or set `HCAP_DATA`).
3. In R:
```r
setwd("~/Desktop/heart-capillary-sex-disease")
source("run_all.R")
```
Outputs (tables, figures, `.rds`) are written to `~/Documents/heart_capillary_disease` (or set `HCAP2_DATA`).

## Licence
MIT (see `LICENSE`). Author: Maria Adelus (ORCID 0000-0002-9676-9214).
