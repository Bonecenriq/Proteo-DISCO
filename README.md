# DISCO-SomaScan-Aging

Analysis code for a proteomic characterization of **DISCO (Distance of
Covariance)**, a Mahalanobis-distance-based multi-organ dysregulation
entropy score, using SomaScan aptamer proteomics. Includes phenotype
associations, organ-specific decomposition, pathway enrichment, and
hospitalization risk modeling.

> **Status:** code accompanying a manuscript currently under review.

---

## 1. Overview

This repository contains the complete analysis pipeline used to:

1. Merge SomaScan aptamer expression data (Box-Cox transformed) with the
   DISCO (Distance of Covariance) entropy score (weighted, unweighted,
   and organ-specific variants) computed in a companion pipeline.
2. Run differential expression analysis (DEA) of the proteome against
   continuous DISCO, both unadjusted (`limma`) and covariate-adjusted
   (linear models with age/sex/site).
3. Build the baseline characteristics table (Table 1) stratified by
   DISCO quartile.
4. Test associations between DISCO and physical performance, cognitive,
   and mitochondrial-function phenotypes (Figure 2), including
   sensitivity analyses adjusting for GDF-15, cystatin C, and
   medication count.
5. Repeat the phenotype and outcome analyses using an **organ-specific**
   decomposition of DISCO (Figure 3, Supplementary Figures 3–4), including
   organ-specific DEA and Hallmark GSEA.
6. Generate the volcano plot and GSEA bubble plot summarizing the
   proteome-wide and pathway-level signature of DISCO (Figure 4).
7. Model non-elective hospitalization risk by DISCO quartile using Cox
   proportional-hazards models, with organ-specific and whole-proteome
   hazard-ratio heatmaps and Kaplan–Meier curves (Figure 5, Supplementary
   Figure 5), plus sensitivity analyses.
8. Compare weighted vs. unweighted DISCO scoring.

The script is organized as a single, sequentially-run R file with clearly
numbered/titled sections (see [Pipeline structure](#3-pipeline-structure)
below); each section prints sanity-check summaries to the console and
writes tables/figures to `results/`.

## 2. Data availability

This repository contains **analysis code only** — no individual-level
phenotype, proteomic, or outcome data are included.

- The cohort used is Study of Muscle, Mobility and Aging (SOMMA) , accessed under data use agreement. Investigators can request access via SOMMA website(https://sommaonline.ucsf.edu/)


## 3. Pipeline structure

| Section in script | Output |
|---|---|
| 1. Packages | — |
| 2. Data loading & preprocessing | `results/checkpoint_section2.RData` |
| 3. Protein annotation mapping + DEA (continuous DISCO, `limma`) | `results/DEA_continuous_BOX_DISCO.csv` |
| 4. Table 1 (baseline characteristics by DISCO quartile) | `results/Table1_baseline.csv`, `Table1_by_DISCO_quartile*.csv` |
| Figure 2 + Suppl. Fig. 1 (phenotype associations; GDF-15/CysC sensitivity) | figure files, `results/*` |
| Suppl. Fig. 3 (weighted vs. unweighted DISCO) | figure files |
| Figure 3 + Suppl. Fig. 3 (organ-specific DISCO vs. outcomes) | `results/Table_Fig3b_organ_outcome_continuous.csv` |
| Adjusted continuous-DISCO DEA (age/sex/site) | `results/DEA_lm_BOX_DISCO_adj.csv` |
| Suppl. Fig. 4 (organ-specific Hallmark pathways) | figure files |
| Figure 4 (volcano + GSEA bubble plot) | `results/GSEA_Hallmark_DISCO.csv`, figure files |
| Suppl. Fig. A/B (organ × protein heatmap; organ × pathway GSEA heatmap) | `results/GSEA_organ_Hallmark.csv`, `Table_organ_protein_beta.csv` |
| Figure 5 + Suppl. Fig. 5 (KM curves, Cox HR heatmaps, hospitalization) | `results/Table_hosp_HR_BOX.csv`, figure files |
| Weighted vs. unweighted DISCO correlation | console output, scatter plot |
| Sensitivity: + medication count (log) | `results/Table_hosp_HR_meds.csv`, figure files |

All figure and table numbers above refer to the manuscript; update this
table if numbering changes at revision.

## 4. Required inputs

Place the following files in the working directory before running the
script (update `setwd()` at the top of the script accordingly):

| File | Description |
|---|---|
| `02_box_transformed_data_with_pheno.csv` | Box-Cox-transformed SomaScan aptamer expression (`seq.*` columns) merged with phenotype/covariate data; `COHORT` column used to subset the "Parent" cohort (n = 838 in this analysis). |
| `05_entropy_measures.csv` | Weighted DISCO entropy score (`BOX_DISCO_JR`, quartile `BOX_DISCO_JR_Q`, `LOG_DISCO_JR`) per subject `ID`. |
| `05z_entropy_measures_unweighted.csv` | Unweighted DISCO score, for the weighted-vs-unweighted comparison. |
| `09b_entropy_measures_organ_box.csv` | Organ-specific DISCO scores, for Figure 3 and related supplementary analyses. |

These files are produced by the upstream DISCO-scoring pipeline (see
[Data availability](#2-data-availability)) and are **not included here**.

## 5. Requirements

R ≥ 4.2 is recommended. Package dependencies:

```r
install.packages(c(
  "tidyverse", "patchwork", "cowplot", "ggrepel", "RColorBrewer",
  "cluster", "zoo", "sandwich", "lmtest", "scales", "openxlsx"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SomaScan.db", "limma", "clusterProfiler", "msigdbr", "GSVA"
))
```

`survival` ships with base R. See [`sessionInfo.txt`](sessionInfo.txt)
for exact package versions used to generate the manuscript results
(generate this by running `writeLines(capture.output(sessionInfo()),
"sessionInfo.txt")` after sourcing the script in your own environment,
and commit the real file before publishing).

## 6. Usage

```r
# 1. Edit the working-directory path in Section 2 of the script
# 2. Place the four required input CSVs in that directory
# 3. Run section-by-section (recommended on first pass) or source() the
#    whole file. Each major figure/table section is self-contained given
#    the checkpoints saved earlier in the script.
source("somma_analysis_pipeline.R")
```

Outputs are written to `results/`, which is created automatically.

## 7. Repository structure

```
.
├── README.md
├── LICENSE
├── .gitignore
├── somma_analysis_pipeline.R   # main analysis script (rename as you like)
└── results/                    # generated on run; not tracked in git
```

## 8. License

Code is released under the [MIT License](LICENSE) unless noted
otherwise. This does **not** cover the underlying cohort data, which remain subject to their original data use agreements.
