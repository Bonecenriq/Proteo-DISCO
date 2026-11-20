# ==============================================================================
# STUDY       : SOMMA
# PROJECT     : SOMMA SomaLogic
# PROGRAM     : 04_analysis_dataset.R
#
# PURPOSE     : Assemble the SOMMA Parent analysis dataset, map SomaScan aptamers
#               to gene symbols, and run the unadjusted continuous-DISCO
#               differential expression analysis.
#
# INPUTS      : - 02_box_transformed_data_with_pheno.csv (02_preprocess_data.R)
#               - 05_entropy_measures.csv (03_calc_entropy.R)
#               - SomaScan.db aptamer annotation
#
# OUTPUTS     : - checkpoint_section2.RData
#               - DEA_continuous_BOX_DISCO.csv
#               - In-memory objects used by scripts 05-09:
#                   meta_df, meta_sorted, prot_df, prot_sorted, disco_vals,
#                   sd_disco, m, mapping_tbl, anno_tbl
#
# DATA ACCESS : Individual-level SOMMA data are not included in this
#               repository. Access to the underlying study data is required
#               to reproduce the analyses.
#
# DEPENDENCIES: - 99_figure_setup.R
#               - SomaScan.db, AnnotationDbi, limma, readr, dplyr
#
# VERSION     : 1.0
# DATE        : 2026-08-15
# AUTHOR      : Namki Hong
#
# REVISION HISTORY
# 08/15/26    Creation (split out of the combined analysis script)
#
# NOTES       : Analyses are restricted to the SOMMA Parent cohort (n = 838).
#               GDF-15 and Cystatin C are treated as outcomes and are excluded
#               from the aptamer matrix in 02_preprocess_data.R.
#
#               Run this script first; scripts 05-09 expect the objects listed
#               under OUTPUTS to be present in the session.
# ==============================================================================
source("code/99_figure_setup.R")


# =============================================================================
# 1. DATA LOADING AND PREPROCESSING
# =============================================================================


# --- 1a. Protein expression (Box-Cox transformed, Parent n=838) ---
pheno_raw <- readr::read_csv(
  file.path(PATH_DATA, "02_box_transformed_data_with_pheno.csv"),
  show_col_types = FALSE
)

pheno_parent <- pheno_raw |>
  dplyr::filter(COHORT == "Parent")
cat("Parent n =", nrow(pheno_parent), "\n")

meta_cols <- names(pheno_parent)[
  !grepl("^seq\\.", names(pheno_parent))]
meta_df   <- pheno_parent |>
  dplyr::select(all_of(meta_cols))
prot_df   <- pheno_parent |>
  dplyr::select(starts_with("seq."))

cat("Proteins:", ncol(prot_df), "\n")

# --- 1b. DISCO entropy merge ---
entropy_raw <- readr::read_csv(
  file.path(PATH_OUTPUT, "05_entropy_measures.csv"),
  show_col_types = FALSE
)

meta_df <- meta_df |>
  dplyr::left_join(
    entropy_raw |>
      dplyr::select(ID, LOG_DISCO_JR,
                    BOX_DISCO_JR,
                    BOX_DISCO_JR_Q) |>
      dplyr::distinct(ID, .keep_all = TRUE),
    by = "ID"
  )

cat("BOX_DISCO_JR range:",
    round(range(meta_df$BOX_DISCO_JR,
                na.rm=TRUE), 3), "\n")
cat("BOX_DISCO_JR_Q distribution:\n")
print(table(meta_df$BOX_DISCO_JR_Q,
            useNA = "ifany"))

# --- 1c. Sort by BOX_DISCO_JR ---
order_idx   <- order(meta_df$BOX_DISCO_JR)
meta_sorted <- meta_df[order_idx, ]
prot_sorted <- as.matrix(prot_df[order_idx, ])
disco_vals  <- meta_sorted$BOX_DISCO_JR

# sd_disco (for per-SD standardization)
sd_disco <- sd(meta_df$BOX_DISCO_JR, na.rm=TRUE)
cat("sd_disco:", round(sd_disco, 4), "\n")

# working data frame
m <- meta_df |>
  dplyr::mutate(
    Q = factor(BOX_DISCO_JR_Q, levels = 1:4,
               labels = c("Q1\n(Lowest)", "Q2",
                          "Q3", "Q4\n(Highest)"))
  )

save(meta_df, meta_sorted, prot_df, prot_sorted,
     disco_vals, sd_disco, m,
     file = file.path(PATH_RESULTS, "checkpoint_section2.RData"))
cat("Checkpoint saved.\n")


# =============================================================================
# 2. PROTEIN NAME MAPPING AND UNADJUSTED DEA (continuous DISCO)
# =============================================================================

# --- 2a. SomaScan.db mapping ---
all_keys      <- keys(SomaScan.db)
seq_numbers   <- sub("^seq\\.", "", names(prot_df))
seq_numbers   <- gsub("\\.", "-", seq_numbers)
keys_nohyphen <- gsub("-", "", all_keys)
seq_nohyphen  <- gsub("-", "", seq_numbers)
match_idx     <- match(seq_nohyphen, keys_nohyphen)

mapping_tbl <- data.frame(
  original  = names(prot_df),
  converted = all_keys[match_idx]
)
cat("Mapped:", sum(!is.na(mapping_tbl$converted)),
    "/", nrow(mapping_tbl), "\n")

anno_tbl <- AnnotationDbi::select(
  SomaScan.db,
  keys    = na.omit(unique(mapping_tbl$converted)),
  columns = c("PROBEID","SYMBOL","GENENAME"),
  keytype = "PROBEID"
) |>
  dplyr::distinct(PROBEID, .keep_all = TRUE)

# --- 2b. limma DEA (continuous BOX_DISCO_JR) ---
prot_mat_all <- t(as.matrix(prot_df))
design_cont  <- model.matrix(
  ~ BOX_DISCO_JR, data = meta_df)
fit_cont     <- lmFit(prot_mat_all,
                      design_cont) |> eBayes()
results_cont <- topTable(
  fit_cont,
  coef    = "BOX_DISCO_JR",
  number  = Inf,
  sort.by = "P"
)

results_cont_anno <- results_cont |>
  as.data.frame() |>
  dplyr::mutate(
    original  = rownames(results_cont),
    converted = mapping_tbl$converted[
      match(rownames(results_cont),
            mapping_tbl$original)]
  ) |>
  dplyr::left_join(
    anno_tbl,
    by = c("converted" = "PROBEID")
  ) |>
  dplyr::mutate(
    label       = dplyr::if_else(
      is.na(SYMBOL), original, SYMBOL),
    logFC_perSD = logFC * sd_disco,
    neg_log10_p = -log10(adj.P.Val)
  )

cat("DEA done. FDR<0.05:",
    sum(results_cont$adj.P.Val < 0.05,
        na.rm=TRUE), "\n")

readr::write_csv(
  results_cont_anno |>
    dplyr::select(original, SYMBOL, GENENAME,
                  logFC, logFC_perSD,
                  AveExpr, P.Value,
                  adj.P.Val),
  file.path(PATH_RESULTS, "DEA_continuous_BOX_DISCO.csv")
)
cat("Saved: DEA_continuous_BOX_DISCO.csv\n")
