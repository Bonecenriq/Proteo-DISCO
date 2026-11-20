# Proteo-DISCO

Analysis code for a proteomic characterization of DISCO (Distance of Covariance), using SomaScan aptamer proteomics in the SOMMA (Study of Muscle, Mobility and Aging) cohort. Includes data preparation, entropy scoring, phenotype associations, organ-specific decomposition, pathway enrichment, and hospitalization risk modeling.

**Status:** code accompanying a manuscript currently under review.

---

## Contents

1. [Overview](#1-overview)
2. [System requirements](#2-system-requirements)
3. [Installation guide](#3-installation-guide)
4. [Demo](#4-demo)
5. [Instructions for use](#5-instructions-for-use)
6. [Reproduction instructions](#6-reproduction-instructions)
7. [Methods notes](#7-methods-notes)
8. [License, authors, and citation](#8-license-authors-and-citation)

---

## 1. Overview

This repository contains the complete analysis pipeline used to:

1. Build the SOMMA Parent and SOMMA Jr analytic datasets (demographics, anthropometrics, physical function, cardiorespiratory fitness, mitochondrial respiration, frailty, medication use, and medical history) from the study baseline files.
2. Pre-process the SomaScan aptamer data (aptamer QC filtering, regularized Box-Cox transformation, standardization to the SOMMA Parent mean/SD), and derive the outcome variables: GDF-15 and cystatin C as circulating biomarkers, plus non-elective hospitalization events and follow-up time.
3. Compute proteomic entropy for each SOMMA Parent participant, using SOMMA Jr as the reference population. The primary measure is age-weighted DISCO; Mahalanobis-distance measures and sex-specific variants are also computed.
4. Merge the Box-Cox transformed aptamer expression data with the DISCO score (weighted, unweighted, and organ-specific variants), and run differential expression analysis (DEA) of the proteome against continuous DISCO, both unadjusted (`limma`) and covariate-adjusted (age/sex/site).
5. Build the baseline characteristics table stratified by DISCO quartile.
6. Test associations between DISCO and physical performance, cognitive & frailty, and energetic capacity & efficiency (mitochondrial function) phenotypes. Sensitivity analyses adjust for GDF-15, cystatin C, and medication count.
7. Repeat the phenotype analyses using an organ-specific decomposition of DISCO across nine organs, including organ-specific DEA and Hallmark GSEA.
8. Generate the volcano plot and GSEA bubble plot summarizing the proteome-wide and pathway-level signature of DISCO.
9. Model non-elective hospitalization risk by DISCO quartile using Cox proportional-hazards models, with whole-proteome and organ-specific hazard-ratio heatmaps, Kaplan–Meier curves, and sensitivity analyses.
10. Compare weighted vs. unweighted DISCO scoring.

The pipeline is organized as a set of sequentially-run, numbered scripts (see [Pipeline structure](#53-pipeline-structure)). Each script prints sanity-check summaries to the console and writes tables/figures to `results/`.

---

## 2. System requirements

### 2.1 Operating systems

The code is platform-independent R (plus one SAS program) and runs on Linux, macOS, and Windows wherever R ≥ 4.0 and the packages below can be installed.

### 2.2 Software dependencies

- **R** ≥ 4.0
- **SAS** 9.4 (only for `code/01_make_data.sas`)
- **Bioconductor** (installed automatically by `install_packages.R`)

| Purpose | Packages | Source |
|---|---|---|
| Data handling | `tidyverse`, `data.table`, `assertthat`, `haven`, `plyr`, `readxl`, `openxlsx`, `writexl`, `DescTools`, `moments`, `future`, `future.apply` | CRAN |
| SomaScan I/O and annotation | `SomaDataIO` | CRAN |
| | `SomaScan.db`, `AnnotationDbi` | Bioconductor |
| Statistics | `MASS`, `glmnet`, `MuMIn`, `emmeans`, `survival`, `survminer`, `timeROC`, `sandwich`, `lmtest`, `FSA` | CRAN |
| Differential expression and enrichment | `limma`, `fgsea`, `GSVA`, `clusterProfiler` | Bioconductor |
| | `msigdbr` | CRAN |
| Visualization | `ggpubr`, `ggrepel`, `ggtext`, `cowplot`, `patchwork`, `scales`, `corrplot`, `ggcorrplot`, `pheatmap`, `RColorBrewer` | CRAN |
| | `ComplexHeatmap` | Bioconductor |

`ggtext`, `writexl`, and `FSA` are optional: without `ggtext` the subscript axis labels fall back to plain text, and without `FSA` the Table 1 script uses pairwise Wilcoxon post-hoc tests instead of Dunn's test.

### 2.3 Non-standard hardware

**No non-standard hardware is required.** Because DISCO works on all pairwise protein correlations, memory use grows with the square of the number of aptamers; for a panel of several thousand aptamers (as in `03_calc_entropy.R`), at least 16 GB of RAM is recommended.

---

## 3. Installation guide

```bash
git clone https://github.com/Bonecenriq/Proteo-DISCO.git
cd Proteo-DISCO
Rscript install_packages.R
```

`install_packages.R` installs any missing CRAN and Bioconductor packages and reports any that failed. `code/99_analysis_setup.R` also stops with a list of missing packages if a required package is not installed.

---

## 4. Demo

A demo dataset is **not provided**. The analyses require individual-level SOMMA data, which cannot be redistributed (see [Data availability](#51-data-availability)), and the results depend on the specific SOMMA cohort structure (SOMMA Jr as the reference population, SOMMA Parent as the target). Researchers with authorized SOMMA data access can run the full pipeline as described in Sections [5](#5-instructions-for-use) and [6](#6-reproduction-instructions).

---

## 5. Instructions for use

### 5.1 Data availability

Individual-level SOMMA data are **not** included in this repository. You need authorized access to the SOMMA Parent and SOMMA Jr data to run the full pipeline. All local paths in the code are placeholders (`/path/to/...`).

Raw inputs required (not distributed): SOMMA Parent/Jr baseline SAS datasets and formats, `somalogic_subset_ids.csv`, the SomaScan ADAT file, `SO000_Annotations.xlsx` (sample and analyte annotations), and the hospitalization and death endpoint files.

### 5.2 Configuration

1. Set the local data paths at the top of `code/01_make_data.sas` (`PATH_BASELINE`, `PATH_BASEJR`, `PATH_BASEJR2`, `PATH_FORMATS`, `PATH_DATA`).
2. Set `PROJECT_DIR` in `code/99_analysis_setup.R`.

Repository layout:

```
Proteo-DISCO/
├── README.md
├── LICENSE
├── .gitignore            # keeps data and results out of version control
├── install_packages.R    # installs all R dependencies
└── code/
    ├── 01_make_data.sas
    ├── 02_preprocess_data.R
    ├── 03_calc_entropy.R
    ├── 04_analysis_dataset.R
    ├── 05_table1_baseline.R
    ├── 06_phenotype_associations.R
    ├── 07_organ_entropy.R
    ├── 08_proteins_pathways.R
    ├── 09_survival_analysis.R
    ├── 99_analysis_setup.R
    └── 99_entropy_methods.R
```

Create the following local folders next to `code/` before running the pipeline. They are listed in `.gitignore` and are never committed:

```
data/      # SOMMA input data (SAS datasets, ADAT, annotations, endpoint files)
output/    # entropy measures written by 03_calc_entropy.R
results/   # tables and figures written by scripts 04–09
figures/
```

### 5.3 Pipeline structure

```
SOMMA baseline files ──► 01_make_data.sas ──► 02_preprocess_data.R ──► 03_calc_entropy.R
                                                                              │
                                                                              ▼
                                                                    04_analysis_dataset.R
                                                                              │
     ┌──────────────────┬──────────────────────┬─────────────────────┬────────┴─────────────┐
     ▼                  ▼                      ▼                     ▼                      ▼
05_table1_baseline  06_phenotype_assoc ──► 07_organ_entropy ──► 09_survival_analysis  08_proteins_pathways
```

| Script | Section | Main outputs |
|---|---|---|
| `01_make_data.sas` | Analytic dataset construction (Parent & Jr; full, men, women) | `prot_subset`, `protjr_subset` (.sas7bdat) |
| `02_preprocess_data.R` | Phenotype/ADAT merge, hospitalization & death endpoints, aptamer QC, log2 / Box-Cox / inverse-normal transformations | `02_<transform>_transformed_data_with_pheno.csv` |
| `03_calc_entropy.R` | Entropy scoring vs. the SOMMA Jr reference: DISCO, Mahalanobis distance (PCA / shrinkage), quartiles, sex-specific variants | `05_entropy_measures.csv` |
| `04_analysis_dataset.R` | Parent analysis frame (n = 838), aptamer → gene symbol mapping (`SomaScan.db`), unadjusted continuous-DISCO DEA | `DEA_continuous_BOX_DISCO.csv`, `checkpoint_section2.RData` |
| `05_table1_baseline.R` | Baseline characteristics overall and by DISCO quartile, with BH-adjusted post-hoc tests | `Table1_by_DISCO_quartile.xlsx/.csv` |
| `06_phenotype_associations.R` | DISCO quartile vs. phenotypes and circulating biomarkers; GDF-15 / cystatin C sensitivity; unweighted DISCO | `Fig2_v9.png`, `FigS_sensitivity_2panel.png`, `SuppFig_unweighted_DISCO.*` |
| `07_organ_entropy.R` | Organ-specific DISCO vs. phenotypes, by quartile and per 1 SD | `Fig_organ_dot_*_BOX.*`, `Table_Fig3b_organ_outcome_continuous.csv` |
| `08_proteins_pathways.R` | Adjusted DEA, Hallmark GSEA, volcano + GSEA bubble plot, organ × protein / organ × pathway heatmaps | `DEA_lm_BOX_DISCO_adj.csv`, `GSEA_Hallmark_DISCO.csv`, `Fig_volcano_GSEA_combined.*`, `Supp_Fig_organ_combined.*` |
| `09_survival_analysis.R` | KM curves and HR heatmaps; GDF-15 / cystatin C and medication-count sensitivity; reverse-KM follow-up; weighted vs. unweighted DISCO | `Fig_hosp_KM_heatmap_BOX.*`, `Supp_Fig_hosp_HR_sensitivity_BOX.*`, `FigS_meds_combined.*`, `Table_hosp_HR_*.csv` |
| `99_analysis_setup.R` | Paths, packages, global variables, outcome definitions, `boxcox_transform()`, `invnorm_transform()` | n/a |
| `99_entropy_methods.R` | `DISCO()` and `DISCO_fast()` | n/a |

### 5.4 Running the pipeline on your data

```bash
sas code/01_make_data.sas
Rscript code/02_preprocess_data.R
Rscript code/03_calc_entropy.R
```

```r
# from the project root, in a single R session
source("code/04_analysis_dataset.R")
source("code/05_table1_baseline.R")
source("code/06_phenotype_associations.R")
source("code/07_organ_entropy.R")
source("code/08_proteins_pathways.R")
source("code/09_survival_analysis.R")
```

Scripts 04–09 share **in-memory objects** (`meta_df`, `m`, `prot_df`, `analysis_df`, plus plotting helpers defined in 06), so run them in order within a single R session.

**Using DISCO on other datasets.** `DISCO_fast()` in `code/99_entropy_methods.R` can be used on its own: pass a reference data frame, a target data frame, the feature column names, and a named vector of per-feature age associations.

---

## 6. Reproduction instructions

All quantitative results in the manuscript are produced by running the scripts in this order, with access to the SOMMA data:

`01_make_data.sas` → `02_preprocess_data.R` → `03_calc_entropy.R` → (unweighted and organ-specific DISCO variants; see table below) → `04_analysis_dataset.R` → `05_table1_baseline.R` → `06_phenotype_associations.R` → `07_organ_entropy.R` → `08_proteins_pathways.R` → `09_survival_analysis.R`

Commands for each step are in Section [5.4](#54-running-the-pipeline-on-your-data), and the script table in Section [5.3](#53-pipeline-structure) lists which script writes each table and figure file.

The unweighted and organ-specific DISCO variants are produced by variants of the entropy pipeline and must exist before scripts 06–09 are run:

| File | Used by | Content |
|---|---|---|
| `09b_entropy_measures_organ_box.csv` | 07, 09 | Organ-specific DISCO (continuous and quartile) for nine organs |
| `05z_entropy_measures_unweighted.csv` | 06, 09 | Age-unweighted DISCO |
| `Table_organ_protein_beta.csv` | 08 | Organ-specific DISCO × protein β |
| `GSEA_organ_Hallmark.csv` | 08 | Organ-specific Hallmark GSEA results |

---

## 7. Methods notes

**DISCO.** For each SOMMA Parent participant, DISCO measures how much adding that individual's aptamer profile to the SOMMA Jr reference changes the reference Pearson correlation network. Each protein pair is weighted by the absolute product of the two proteins' age associations (estimated in SOMMA Parent):

DISCO = log<sub>10</sub>[ n<sub>ref</sub><sup>2</sup> · Σ<sub>i<j</sub> w<sub>ij</sub> (R<sub>pert,ij</sub> − R<sub>ref,ij</sub>)<sup>2</sup> ],  with w<sub>ij</sub> ∝ |a<sub>i</sub>·a<sub>j</sub>|

`DISCO_fast()` gives the same result as `DISCO()` using a rank-one covariance update, which makes it practical for the full SomaScan panel. Unweighted DISCO uses equal weights across protein pairs.

**Models.**
- Phenotypes: linear regression, β in outcome SD units vs. Q1, adjusted for age, sex, clinic site, body weight, and height.
- Proteome: z-scored aptamer regressed on 1-SD DISCO with `limma`, adjusted for age, sex, and site.
- Pathways: pre-ranked GSEA on the 50 MSigDB Hallmark sets, ranked by the gene-level moderated t-statistic, BH FDR.
- Hospitalization: Cox PH model (Efron ties, robust SE) with the same adjustment set as the phenotype models.
- Sensitivity models additionally adjust for GDF-15, cystatin C, or log medication count.

---

## 8. License, authors, and citation

This code is released under the [MIT License](LICENSE), an Open Source Initiative-approved license. The license covers the analysis code only; SOMMA study data are not distributed and remain subject to the study's data access policies.

- **Tyler Mansfield**: data preparation, preprocessing, and entropy methods (scripts 01–03, `99_*`)
- **Namki Hong**: analysis datasets, statistical analyses, and figures (scripts 04–09)

Citation will be added upon publication.
