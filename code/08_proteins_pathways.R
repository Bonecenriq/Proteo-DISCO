# ==============================================================================
# STUDY       : SOMMA
# PROJECT     : SOMMA SomaLogic
# PROGRAM     : 08_proteins_pathways.R
#
# PURPOSE     : Identify individual proteins and Hallmark pathways associated with
#               proteomic entropy, and produce Figure 4 and Supplementary Figure 4
#               together with Supplementary Tables S2 and S3.
#
# INPUTS      : - prot_df, meta_df, mapping_tbl, anno_tbl (04_analysis_dataset.R)
#               - Table_organ_protein_beta.csv (organ x protein betas; see NOTES)
#               - GSEA_organ_Hallmark.csv (organ-specific GSEA; see NOTES)
#
# OUTPUTS     : - DEA_lm_BOX_DISCO_adj.csv       (Supplementary Table S2)
#               - GSEA_Hallmark_DISCO.csv        (Supplementary Table S3)
#               - Fig_GSEA_Hallmark_DISCO.png / .pdf
#               - Fig_volcano_GSEA_combined.png / .pdf   (Figure 4)
#               - Supp_Fig_organ_protein_heatmap.png / .pdf
#               - Supp_Fig_GSEA_organ_heatmap.png / .pdf
#               - Supp_Fig_organ_combined.png / .pdf     (Supplementary Figure 4)
#
# DATA ACCESS : Individual-level SOMMA data are not included in this
#               repository. Access to the underlying study data is required
#               to reproduce the analyses.
#
# DEPENDENCIES: - 99_figure_setup.R
#               - 04_analysis_dataset.R (objects prot_df, meta_df, mapping_tbl,
#                 anno_tbl)
#               - limma, clusterProfiler, msigdbr, ggrepel, patchwork
#
# VERSION     : 1.0
# DATE        : 2026-08-15
# AUTHOR      : Namki Hong
#
# REVISION HISTORY
# 08/15/26    Creation (split out of the combined analysis script)
#
# NOTES       : The adjusted differential expression analysis models each aptamer,
#               standardized to unit SD, on 1-SD of DISCO with adjustment for age,
#               sex and clinic site. Betas are therefore standardized coefficients.
#
#               Gene set enrichment uses the 50 MSigDB Hallmark sets, pre-ranked by
#               the limma moderated t-statistic collapsed to the gene level
#               (strongest aptamer per symbol).
#
#               Supplementary Figure 4 requires Table_organ_protein_beta.csv and
#               GSEA_organ_Hallmark.csv, the organ-specific counterparts of the
#               tables produced here. Those files are generated outside this
#               repository; see README.md.
# ==============================================================================
source("code/99_figure_setup.R")


# =============================================================================
# 1. ADJUSTED CONTINUOUS-DISCO DEA (limma)  ->  Supplementary Table S2
#   Model : aptamer_z ~ DISCO_z + age + sex + site   (age/sex/site adjusted)
#   beta  : per 1-SD DISCO, in aptamer SD units (standardized coefficient)
#   Output: results_v14/DEA_lm_BOX_DISCO_adj.csv
#   Prereqs (04_analysis_dataset.R): prot_df, meta_df, mapping_tbl, anno_tbl
# =============================================================================
stopifnot(exists("prot_df"), exists("meta_df"),
          exists("mapping_tbl"), exists("anno_tbl"))

# --- model covariates (complete-case on DISCO + covariates) ---
mdat <- meta_df |>
  dplyr::mutate(SITE_f = factor(SITE)) |>
  dplyr::select(BOX_DISCO_JR, D1AGE2, EEFEMALE, SITE_f)
keep <- stats::complete.cases(
  mdat[, c("BOX_DISCO_JR", "D1AGE2", "EEFEMALE", "SITE_f")])
mdat <- mdat[keep, ]
cat("DEA (adjusted) n =", nrow(mdat), "\n")

# --- response: each aptamer standardized to unit SD (=> beta in aptamer-SD units) ---
X      <- as.matrix(prot_df)[keep, , drop = FALSE]      # participants x aptamers
sd_col <- apply(X, 2, stats::sd, na.rm = TRUE)
X      <- X[, sd_col > 0, drop = FALSE]                 # drop zero-variance aptamers
Xz     <- scale(X)                                      # z-score each aptamer (column)
prot_mat_adj <- t(Xz)                                   # aptamers x participants

# --- predictor: DISCO standardized to 1-SD ---
disco_z <- as.numeric(scale(mdat$BOX_DISCO_JR))

design_adj <- model.matrix(
  ~ disco_z + mdat$D1AGE2 + mdat$EEFEMALE + mdat$SITE_f)
colnames(design_adj)[2] <- "DISCO_z"

fit_adj <- limma::lmFit(prot_mat_adj, design_adj) |> limma::eBayes()
res_adj <- limma::topTable(fit_adj, coef = "DISCO_z",
                           number = Inf, sort.by = "P",
                           confint = TRUE)   # adds CI.L / CI.R

# --- annotate (reuse the mapping/annotation from 04_analysis_dataset.R) ---
res_adj_anno <- res_adj |>
  as.data.frame() |>
  dplyr::mutate(
    original  = rownames(res_adj),
    converted = mapping_tbl$converted[
      match(original, mapping_tbl$original)]
  ) |>
  dplyr::left_join(anno_tbl, by = c("converted" = "PROBEID")) |>
  dplyr::mutate(
    beta        = logFC,          # per 1-SD DISCO, aptamer-SD units (standardized)
    beta_lo     = CI.L,
    beta_hi     = CI.R,
    FDR         = adj.P.Val,
    neg_log10_q = -log10(adj.P.Val),
    label       = dplyr::if_else(is.na(SYMBOL), original, SYMBOL)
  ) |>
  dplyr::arrange(P.Value)

# --- save Supplementary Table S2 ---
readr::write_csv(
  res_adj_anno |>
    dplyr::select(original, SYMBOL, GENENAME, label,
                  beta, beta_lo, beta_hi,
                  t, P.Value, FDR, neg_log10_q, AveExpr),
  file.path(PATH_RESULTS, "DEA_lm_BOX_DISCO_adj.csv")
)
cat("Saved: DEA_lm_BOX_DISCO_adj.csv\n")

# --- sanity check against the manuscript thresholds ---
n_up <- sum(res_adj_anno$FDR < 0.001 & res_adj_anno$beta >  0.2, na.rm = TRUE)
n_dn <- sum(res_adj_anno$FDR < 0.001 & res_adj_anno$beta < -0.2, na.rm = TRUE)
cat(sprintf("|beta|>0.2 & FDR<0.001  ->  UP: %d | DOWN: %d\n", n_up, n_dn))

# =============================================================
# Figure 4. Volcano plot and bubble plot
# =============================================================


col_up      <- "#C0392B"
col_dn      <- "#2471A3"
col_up_dark <- "#7B241C"
col_dn_dark <- "#1A5276"
col_ns      <- "grey80"


# =============================================================
# 1. DEA data load
# =============================================================
dea_df <- readr::read_csv(
  file.path(PATH_RESULTS, "DEA_lm_BOX_DISCO_adj.csv"),
  show_col_types=FALSE
) |>
  dplyr::mutate(
    sig_group = dplyr::case_when(
      FDR < 0.001 & beta >  0.2 ~ "UP",
      FDR < 0.001 & beta < -0.2 ~ "DOWN",
      TRUE ~ "NS"
    ),
    sig_group = factor(sig_group,
                       levels=c("UP","DOWN","NS"))
  )

cat("DEA n:", nrow(dea_df), "\n")
cat("UP:", sum(dea_df$sig_group=="UP"),
    "| DOWN:", sum(dea_df$sig_group=="DOWN"), "\n")


# =============================================================
# 2. Label protein definition
# =============================================================
protein_name_map <- c(
  # UP — Known aging markers
  "B2M"      = "\u03b22-Microglobulin",
  "IGFBP4"   = "IGFBP-4",
  "EDA2R"    = "EDA2R",
  "CCL2"     = "MCP-1",
  "FSTL3"    = "FSTL-3",
  "TNFRSF1B" = "sTNFR2",
  "TNFRSF1A" = "sTNFR1",
  "CHI3L1"   = "YKL-40",
  # UP — Novel / Emerging
  "TIMP1"    = "TIMP-1",
  "VSIG4"    = "VSIG4",
  "SVEP1"    = "SVEP1",
  "EFEMP1"   = "Fibulin-3",
  "PTGDS"    = "PTGDS",
  "PYCARD"   = "ASC",
  "METRNL"   = "Meteorin-like",
  # DOWN
  "VASP"     = "VASP",
  "EIF4E"    = "eIF4E",
  "RAB7A"    = "RAB7A",
  "NAP1L1"   = "NAP1L1",
  "LYN"      = "LYN"
)

label_df <- dea_df |>
  dplyr::filter(SYMBOL %in% names(protein_name_map)) |>
  dplyr::group_by(SYMBOL) |>
  dplyr::slice_max(abs(beta), n=1,
                   with_ties=FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    prot_name = protein_name_map[SYMBOL]
  )

cat("Label proteins:", nrow(label_df), "\n")


# =============================================================
# 3. Panel a: Volcano plot
# =============================================================
p_volcano <- ggplot2::ggplot(
  dea_df,
  ggplot2::aes(x=beta, y=neg_log10_q)
) +
  # 
  ggplot2::geom_vline(xintercept=0,
                      color="grey55", linewidth=0.25) +
  ggplot2::geom_hline(
    yintercept=-log10(0.001),
    linetype="dashed",
    color="grey40", linewidth=0.35) +
  ggplot2::geom_vline(
    xintercept=c(-0.2, 0.2),
    linetype="dashed",
    color="grey40", linewidth=0.35) +
  
  # NS
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d, sig_group=="NS"),
    color=col_ns, size=0.7, alpha=0.30) +
  
  # Sig 
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d,
                            sig_group=="UP",
                            !SYMBOL %in% names(protein_name_map)),
    color=col_up, size=0.7, alpha=0.55) +
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d,
                            sig_group=="DOWN",
                            !SYMBOL %in% names(protein_name_map)),
    color=col_dn, size=0.7, alpha=0.55) +
  
  # 
  ggplot2::geom_point(
    data=dplyr::filter(label_df,
                       sig_group=="UP"),
    ggplot2::aes(x=beta, y=neg_log10_q),
    color=col_up_dark, size=0.7,
    alpha=1.0, shape=16) +
  ggplot2::geom_point(
    data=dplyr::filter(label_df,
                       sig_group=="DOWN"),
    ggplot2::aes(x=beta, y=neg_log10_q),
    color=col_dn_dark, size=0.7,
    alpha=1.0, shape=16) +
  
  # UP label
  ggrepel::geom_text_repel(
    data=dplyr::filter(label_df,
                       sig_group=="UP"),
    ggplot2::aes(
      x     = beta,
      y     = neg_log10_q,
      label = prot_name
    ),
    color              = col_up_dark,
    size               = 2.4,
    fontface           = "italic",
    direction          = "y",
    nudge_x            = 0.62 -
      dplyr::filter(label_df,
                    sig_group=="UP")$beta,
    hjust              = 0,
    segment.size       = 0.3,
    segment.alpha      = 0.5,
    segment.curvature  = -0.1,
    box.padding        = 0.15,
    max.overlaps       = 50,
    min.segment.length = 0,
    show.legend        = FALSE,
    seed               = 42
  ) +
  
  # DOWN label
  ggrepel::geom_text_repel(
    data=dplyr::filter(label_df,
                       sig_group=="DOWN"),
    ggplot2::aes(
      x     = beta,
      y     = neg_log10_q,
      label = prot_name
    ),
    color              = col_dn_dark,
    size               = 2.4,
    fontface           = "italic",
    direction          = "y",
    nudge_x            = -0.62 -
      dplyr::filter(label_df,
                    sig_group=="DOWN")$beta,
    hjust              = 1,
    segment.size       = 0.3,
    segment.alpha      = 0.5,
    segment.curvature  = 0.1,
    box.padding        = 0.15,
    max.overlaps       = 50,
    min.segment.length = 0,
    show.legend        = FALSE,
    seed               = 42
  ) +
  
  # UP/DOWN count
  ggplot2::annotate("text",
                    x= 0.22,
                    y= max(dea_df$neg_log10_q,
                           na.rm=TRUE) * 0.97,
                    label=paste0("UP: ",
                                 sum(dea_df$sig_group=="UP")),
                    color=col_up_dark, size=2.8,
                    hjust=0, fontface="bold") +
  ggplot2::annotate("text",
                    x=-0.52,
                    y= max(dea_df$neg_log10_q,
                           na.rm=TRUE) * 0.97,
                    label=paste0("DOWN: ",
                                 sum(dea_df$sig_group=="DOWN")),
                    color=col_dn_dark, size=2.8,
                    hjust=0, fontface="bold") +
  
  ggplot2::scale_x_continuous(
    name   = expression(
      beta~"(per 1-SD DISCO, protein SD units)"),
    limits = c(-0.85, 0.85),
    breaks = seq(-0.8, 0.8, 0.4),
    labels = c("-0.8","-0.4","0","0.4","0.8")
  ) +
  ggplot2::scale_y_continuous(
    name  = expression(-log[10](FDR)),
    expand= ggplot2::expansion(
      mult=c(0.02, 0.06))
  ) +
  ggplot2::labs(tag="a") +
  ggplot2::theme_bw(base_size=10) +
  ggplot2::theme(
    plot.tag         = ggplot2::element_text(
      face="bold", size=13),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      color="grey92", linewidth=0.3),
    plot.background  = ggplot2::element_rect(
      fill="white", color=NA),
    plot.margin      = ggplot2::margin(5, 2, 5, 5)
  )

# =============================================================================
# 3. HALLMARK GSEA ON CONTINUOUS DISCO  ->  Supplementary Table S3
#   Ranking : limma moderated t-statistic for the DISCO term (Table S2),
#             collapsed to gene level (strongest aptamer = max |t| per SYMBOL).
#   Sets    : 50 MSigDB Hallmark gene sets (msigdbr, Homo sapiens, "H").
#   Method  : clusterProfiler::GSEA (pre-ranked), BH FDR across Hallmark sets.
#   Output  : results_v14/GSEA_Hallmark_DISCO.csv   (all 50 pathways retained)
#   Prereqs : res_adj_anno from the Table S2 block above (columns SYMBOL, t).
# =============================================================================
stopifnot(exists("res_adj_anno"),
          all(c("SYMBOL", "t") %in% names(res_adj_anno)))
suppressMessages({library(clusterProfiler); library(msigdbr); library(stringr)})

# --- gene-level ranking metric: DISCO t-statistic, one value per gene symbol --
rank_tbl <- res_adj_anno |>
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "", !is.na(t)) |>
  dplyr::group_by(SYMBOL) |>
  dplyr::slice_max(order_by = abs(t), n = 1, with_ties = FALSE) |>  # strongest aptamer/gene
  dplyr::ungroup()

geneList <- rank_tbl$t
names(geneList) <- rank_tbl$SYMBOL
geneList <- sort(geneList, decreasing = TRUE)   # GSEA needs decreasing, unique names
cat("Ranked genes for GSEA:", length(geneList),
    " (from", nrow(res_adj_anno), "aptamers)\n")

# --- Hallmark gene sets as TERM2GENE (handle msigdbr old/new API) --------------
h_sets <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens", category   = "H"),
  error = function(e)
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"))
h_sets <- h_sets |>
  dplyr::select(gs_name, gene_symbol) |>
  dplyr::distinct()
cat("Hallmark sets:", dplyr::n_distinct(h_sets$gs_name), "\n")

# --- pre-ranked GSEA (clusterProfiler defaults; keep ALL sets for the table) ---
set.seed(42)
gsea_res <- clusterProfiler::GSEA(
  geneList      = geneList,
  TERM2GENE     = h_sets,
  pvalueCutoff  = 1,          # retain every pathway; the figure code filters at FDR<0.25
  pAdjustMethod = "BH",       # Benjamini-Hochberg
  seed          = TRUE,
  verbose       = FALSE
)

# --- format to the exact schema consumed by Figure 4b + organ-heatmap reader ---
gsea_out <- as.data.frame(gsea_res) |>
  dplyr::mutate(
    gs_label  = stringr::str_to_title(gsub("_", " ", sub("^HALLMARK_", "", ID))),
    direction = dplyr::if_else(NES > 0, "UP (high DISCO)", "DOWN (high DISCO)")
  ) |>
  dplyr::arrange(dplyr::desc(NES)) |>
  dplyr::select(ID, gs_label, direction, setSize,
                enrichmentScore, NES, pvalue, p.adjust, qvalue,
                rank, leading_edge, core_enrichment)

readr::write_csv(gsea_out, file.path(PATH_RESULTS, "GSEA_Hallmark_DISCO.csv"))
cat("Saved: GSEA_Hallmark_DISCO.csv  (", nrow(gsea_out), "Hallmark pathways )\n")

# --- sanity checks vs. manuscript / Figure 4b expectations --------------------
cat("  FDR<0.05 :", sum(gsea_out$p.adjust < 0.05), "\n")
cat("  FDR<0.25 :", sum(gsea_out$p.adjust < 0.25), "\n")
cat("  UP / DOWN (FDR<0.05): ",
    sum(gsea_out$p.adjust < 0.05 & gsea_out$NES > 0), "/",
    sum(gsea_out$p.adjust < 0.05 & gsea_out$NES < 0), "\n")
cat("  top UP  :", paste(head(gsea_out$gs_label[gsea_out$NES > 0], 3), collapse = ", "), "\n")
cat("  top DOWN:", paste(head(rev(gsea_out$gs_label[gsea_out$NES < 0]), 3), collapse = ", "), "\n")


# =============================================================
# 4. GSEA data load
# =============================================================
gsea_raw <- readr::read_csv(
  file.path(PATH_RESULTS, "GSEA_Hallmark_DISCO.csv"),
  show_col_types=FALSE
)

gsea_df <- gsea_raw |>
  dplyr::filter(p.adjust < 0.25) |>
  dplyr::mutate(
    # protein labeling
    gs_label = gsub("Dna",     "DNA",       gs_label),
    gs_label = gsub("Il2",     "IL2",       gs_label),
    gs_label = gsub("Il6",     "IL6",       gs_label),
    gs_label = gsub("Pi3k",    "PI3K",      gs_label),
    gs_label = gsub("Mtorc1",  "mTORC1",    gs_label),
    gs_label = gsub("Mtor",    "mTOR",      gs_label),
    gs_label = gsub("Myc",     "MYC",       gs_label),
    gs_label = gsub("Kras",    "KRAS",      gs_label),
    gs_label = gsub("G2m",     "G2M",       gs_label),
    gs_label = gsub("Jak Stat3","JAK-STAT3",gs_label),
    gs_label = gsub("Stat5",   "STAT5",     gs_label),
    gs_label = gsub("Akt",     "AKT",       gs_label),
    # FDR<0.05: UP/DOWN 
    pt_color = dplyr::case_when(
      p.adjust < 0.05 &
        direction == "UP (high DISCO)"   ~ "UP",
      p.adjust < 0.05 &
        direction == "DOWN (high DISCO)" ~ "DOWN",
      TRUE ~ "NS"
    ),
    alpha_pt = dplyr::if_else(
      p.adjust < 0.05, 1.0, 0.45)
  ) |>
  dplyr::arrange(dplyr::desc(NES)) |>
  dplyr::mutate(
    gs_label = factor(gs_label,
                      levels=gs_label)
  )

cat("GSEA pathways (FDR<0.25):",
    nrow(gsea_df), "\n")
cat("  FDR<0.05:",
    sum(gsea_df$p.adjust < 0.05), "\n")
cat("  FDR 0.05-0.25:",
    sum(gsea_df$p.adjust >= 0.05), "\n")


# =============================================================
# 5. Panel b: GSEA Bubble plot
# =============================================================
p_gsea <- ggplot2::ggplot(
  gsea_df,
  ggplot2::aes(
    x     = NES,
    y     = gs_label,
    size  = -log10(p.adjust),
    color = pt_color,
    alpha = alpha_pt
  )
) +
  ggplot2::geom_vline(
    xintercept=0,
    color="grey40", linewidth=0.4) +
  
  ggplot2::geom_point(shape=16) +
  
  ggplot2::scale_color_manual(
    values=c(
      "UP"   = col_up,
      "DOWN" = col_dn,
      "NS"   = "grey60"
    ),
    name  = "Direction\n(high DISCO)",
    guide = ggplot2::guide_legend(
      order=1,
      override.aes=list(
        size=4, alpha=1,
        color=c(col_dn, "grey60", col_up)))
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_size_continuous(
    name   = "FDR q-value",
    range  = c(2.5, 9),
    breaks = c(-log10(0.05),
               -log10(0.01),
               -log10(0.001),
               -log10(0.0001)),
    labels = c("<0.05","<0.01",
               "<0.001","<0.0001"),
    guide  = ggplot2::guide_legend(order=2)
  ) +
  ggplot2::scale_x_continuous(
    name  ="Normalized Enrichment Score (NES)",
    limits=c(
      min(gsea_df$NES) - 0.35,
      max(gsea_df$NES) + 0.35),
    breaks=seq(-2.5, 2.5, 0.5)
  ) +
  ggplot2::scale_y_discrete(name=NULL) +
  ggplot2::labs(
    tag     = "b",
  ) +
  ggplot2::theme_bw(base_size=10) +
  ggplot2::theme(
    plot.tag           = ggplot2::element_text(
      face="bold", size=13),
    plot.subtitle      = ggplot2::element_text(
      size=7.5, color="grey40",
      lineheight=1.2),
    plot.caption       = ggplot2::element_text(
      size=7, color="grey50",
      hjust=0, lineheight=1.2),
    axis.text.y        = ggplot2::element_text(
      size=9, color="black"),
    axis.text.x        = ggplot2::element_text(
      size=9),
    axis.title.x       = ggplot2::element_text(
      size=10),
    panel.grid.major.y = ggplot2::element_line(
      color="grey90", linewidth=0.3),
    panel.grid.major.x = ggplot2::element_line(
      color="grey90", linewidth=0.3),
    panel.grid.minor   = ggplot2::element_blank(),
    legend.position    = "right",
    legend.text        = ggplot2::element_text(
      size=8),
    legend.title       = ggplot2::element_text(
      size=8.5),
    plot.background    = ggplot2::element_rect(
      fill="white", color=NA),
    plot.margin        = ggplot2::margin(5, 5, 5, 2)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(order=1),
    size  = ggplot2::guide_legend(order=2)
  )

# GSEA single save
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_GSEA_Hallmark_DISCO.pdf"),
  plot=p_gsea, width=200, height=180,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_GSEA_Hallmark_DISCO.png"),
  plot=p_gsea, width=200, height=180,
  units="mm", dpi=300)
cat("Saved: Fig_GSEA_Hallmark_DISCO\n")


# =============================================================
# 6. Combined figure (Panel a + b)
# =============================================================
fig_combined <- p_volcano + p_gsea +
  patchwork::plot_layout(
    widths=c(1.1, 1)
  ) +
  patchwork::plot_annotation(
    theme=ggplot2::theme(
      plot.title      = ggplot2::element_text(
        face="bold", size=12,
        margin=ggplot2::margin(b=3)),
      plot.subtitle   = ggplot2::element_text(
        size=7.5, color="grey40",
        margin=ggplot2::margin(b=5)),
      plot.caption    = ggplot2::element_text(
        size=6.5, color="grey50",
        hjust=0, lineheight=1.2),
      plot.background = ggplot2::element_rect(
        fill="white", color=NA),
      plot.margin     = ggplot2::margin(8,8,6,8)
    )
  )

ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_volcano_GSEA_combined.pdf"),
  plot  = fig_combined,
  width = 380, height = 175,
  units = "mm", device = cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_volcano_GSEA_combined.png"),
  plot  = fig_combined,
  width = 380, height = 175,
  units = "mm", dpi   = 300)
cat("Saved: Fig_volcano_GSEA_combined\n")



# =============================================================
# Supplementary Figure 4. Organ-specific Entropy and hallmark pathways
# =============================================================


organs <- c("blood","bone","brain","heart",
            "kidney","liver","lung",
            "lymphoid","skeletal")

organ_labels <- c(
  "blood"    = "Blood Vessel",
  "bone"     = "Bone",
  "brain"    = "Brain",
  "heart"    = "Heart",
  "kidney"   = "Kidney",
  "liver"    = "Liver",
  "lung"     = "Lung",
  "lymphoid" = "Lymphoid",
  "skeletal" = "Skeletal Muscle"
)

col_up <- "#C0392B"
col_dn <- "#2471A3"


# =============================================================
# SUPP FIG A: Organ × Protein Heatmap
# =============================================================

# ── A-1. Protein selection: Full DISCO DEA top 3% ──────────────────
dea_df <- readr::read_csv(
  file.path(PATH_RESULTS, "DEA_lm_BOX_DISCO_adj.csv"),
  show_col_types=FALSE
) |>
  dplyr::filter(!is.na(SYMBOL))

sig_up <- dea_df |>
  dplyr::filter(FDR < 0.001, beta > 0.2) |>
  dplyr::group_by(SYMBOL) |>
  dplyr::slice_max(beta, n=1,
                   with_ties=FALSE) |>
  dplyr::ungroup() |>
  dplyr::slice_max(beta, n=24) |>
  dplyr::arrange(dplyr::desc(beta))

sig_dn <- dea_df |>
  dplyr::filter(FDR < 0.001, beta < -0.2) |>
  dplyr::group_by(SYMBOL) |>
  dplyr::slice_min(beta, n=1,
                   with_ties=FALSE) |>
  dplyr::ungroup() |>
  dplyr::slice_min(beta, n=51) |>
  dplyr::arrange(beta)

target_proteins <- dplyr::bind_rows(
  sig_up, sig_dn
) |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::mutate(
    full_disco_dir = dplyr::if_else(
      beta > 0, "UP","DOWN")
  )

cat("Target proteins:",
    nrow(target_proteins), "\n")
cat("UP:", sum(target_proteins$full_disco_dir=="UP"),
    "| DOWN:",
    sum(target_proteins$full_disco_dir=="DOWN"), "\n")

# ── A-2. Organ × protein beta load  ───────────────────────────
organ_beta_df <- readr::read_csv(
  file.path(PATH_RESULTS, "Table_organ_protein_beta.csv"),
  show_col_types=FALSE
) |>
  dplyr::filter(!is.na(SYMBOL))

cat("organ_beta_df rows:",
    nrow(organ_beta_df), "\n")

# ── A-3. protein sorting: UP , DOWN  ───────────────────────
protein_order <- c(
  target_proteins |>
    dplyr::filter(full_disco_dir=="UP",
                  !is.na(SYMBOL)) |>
    dplyr::arrange(dplyr::desc(beta)) |>
    dplyr::pull(SYMBOL),
  target_proteins |>
    dplyr::filter(full_disco_dir=="DOWN",
                  !is.na(SYMBOL)) |>
    dplyr::arrange(beta) |>
    dplyr::pull(SYMBOL)
) |> unique()

# ── A-4. Full Proteome beta column ─────────────────────────────
heat_full <- target_proteins |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::mutate(
    organ    = "full",
    sig_mark = dplyr::case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
  ) |>
  dplyr::select(SYMBOL, organ,
                beta, sig_mark)

# ── A-5. Organ beta column ─────────────────────────────────────
heat_organ <- organ_beta_df |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::mutate(
    sig_mark = dplyr::case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
  ) |>
  dplyr::select(SYMBOL, organ,
                beta, sig_mark)

# ── A-6. merge + factor config ────────────────────────────────
heat_long <- dplyr::bind_rows(
  heat_full, heat_organ
) |>
  dplyr::mutate(
    SYMBOL = factor(
      SYMBOL,
      levels=rev(protein_order)
    ),
    organ_label = factor(
      dplyr::recode(organ,
                    "full"     = "Full\nProteome",
                    "blood"    = "Blood\nVessel",
                    "bone"     = "Bone",
                    "brain"    = "Brain",
                    "heart"    = "Heart",
                    "kidney"   = "Kidney",
                    "liver"    = "Liver",
                    "lung"     = "Lung",
                    "lymphoid" = "Lymphoid",
                    "skeletal" = "Skeletal\nMuscle"
      ),
      levels=c("Full\nProteome",
               "Blood\nVessel","Bone","Brain",
               "Heart","Kidney","Liver","Lung",
               "Lymphoid","Skeletal\nMuscle")
    )
  ) |>
  dplyr::filter(!is.na(SYMBOL))

# UP/DOWN margins
n_up <- sum(target_proteins$full_disco_dir=="UP",
            na.rm=TRUE)
n_dn <- sum(target_proteins$full_disco_dir=="DOWN",
            na.rm=TRUE)
hline_y <- n_dn + 0.5

# ── A-7. Plot ─────────────────────────────────────────────────
p_heat_prot <- ggplot2::ggplot(
  heat_long,
  ggplot2::aes(
    x    = organ_label,
    y    = SYMBOL,
    fill = beta
  )
) +
  ggplot2::geom_tile(
    color="white", linewidth=0.3) +
  
  # sig asterisk
  ggplot2::geom_text(
    ggplot2::aes(label=sig_mark),
    size=1.8, vjust=0.75,
    color="white"
  ) +
  
  # Full Proteome line (vertical)
  ggplot2::annotate("segment",
                    x=1.5, xend=1.5,
                    y=0.5, yend=n_up+n_dn+0.5,
                    color="grey20", linewidth=0.8
  ) +
  
  # UP/DOWN line (horizontal)
  ggplot2::geom_hline(
    yintercept=hline_y,
    color="grey20", linewidth=0.8
  ) +
  
  # UP/DOWN label
  ggplot2::annotate("text",
                    x    = 10.7,
                    y    = c(n_dn/2, n_dn + n_up/2),
                    label= c("DOWN\n(full DISCO)",
                             "UP\n(full DISCO)"),
                    color= c(col_dn, col_up),
                    size = 2.8, fontface="bold",
                    hjust= 0
  ) +
  
  ggplot2::scale_fill_gradient2(
    low      = col_dn,
    mid      = "white",
    high     = col_up,
    midpoint = 0,
    limits   = c(-0.35, 0.35),
    oob      = scales::squish,
    name     = expression(
      beta~"(per 1-SD entropy)")
  ) +
  ggplot2::scale_x_discrete(
    name  = NULL,
    expand= ggplot2::expansion(
      add=c(0.5, 1.5))
  ) +
  ggplot2::scale_y_discrete(name=NULL) +
  ggplot2::labs(
    title   = NULL,
    subtitle= NULL,
    caption = paste0(
      "* FDR<0.05, ** FDR<0.01, ",
      "*** FDR<0.001. ",
      "Color: \u03b2 per 1-SD entropy. ",
      "Color scale capped at \u00b10.35."
    )
  ) +
  ggplot2::coord_cartesian(clip="off") +
  ggplot2::theme_minimal(base_size=9) +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(
      face="bold", size=11),
    plot.subtitle   = ggplot2::element_text(
      size=7.5, color="grey40"),
    plot.caption    = ggplot2::element_text(
      size=7, color="grey50", hjust=0),
    axis.text.x     = ggplot2::element_text(
      size=8, angle=35, hjust=1,
      color="black", lineheight=0.85),
    axis.text.y     = ggplot2::element_text(
      size=6.5, color="black",
      face="italic"),
    legend.position  = "bottom",
    legend.key.width = ggplot2::unit(
      1.5,"cm"),
    legend.title     = ggplot2::element_text(
      size=8),
    legend.text      = ggplot2::element_text(
      size=7.5),
    panel.grid       = ggplot2::element_blank(),
    plot.background  = ggplot2::element_rect(
      fill="white", color=NA),
    plot.margin      = ggplot2::margin(
      8, 65, 8, 8)
  )

ggplot2::ggsave(
  file.path(PATH_RESULTS, "Supp_Fig_organ_protein_heatmap.pdf"),
  plot=p_heat_prot, width=190, height=250,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Supp_Fig_organ_protein_heatmap.png"),
  plot=p_heat_prot, width=190, height=250,
  units="mm", dpi=300)
cat("Saved: Supp_Fig_organ_protein_heatmap\n")


# =============================================================
# SUPP FIG B: Organ × Pathway GSEA Heatmap
# =============================================================

# ── B-1. data load ─────────────────────────────────────────
gsea_full <- readr::read_csv(
  file.path(PATH_RESULTS, "GSEA_Hallmark_DISCO.csv"),
  show_col_types=FALSE
)

gsea_organ <- readr::read_csv(
  file.path(PATH_RESULTS, "GSEA_organ_Hallmark.csv"),
  show_col_types=FALSE
)

# ── B-2.  pathway select ───────────────────────────────────
# Full FDR<0.05 + organ 2 or more FDR<0.05
pw_full_sig <- gsea_full |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::pull(ID)

pw_organ_sig <- gsea_organ |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::count(ID) |>
  dplyr::filter(n >= 2) |>
  dplyr::pull(ID)

show_pws <- union(pw_full_sig, pw_organ_sig)
cat("presented pathway:", length(show_pws), "\n")

# ── B-3. Nomenclature ──────────────────────────────────────
clean_label <- function(x) {
  x |>
    gsub("HALLMARK_", "", x=_) |>
    gsub("_", " ", x=_) |>
    stringr::str_to_title() |>
    gsub("Dna",     "DNA",       x=_) |>
    gsub("Il2",     "IL2",       x=_) |>
    gsub("Il6",     "IL6",       x=_) |>
    gsub("Pi3k",    "PI3K",      x=_) |>
    gsub("Mtorc1",  "mTORC1",    x=_) |>
    gsub("Mtor",    "mTOR",      x=_) |>
    gsub("Myc",     "MYC",       x=_) |>
    gsub("Kras",    "KRAS",      x=_) |>
    gsub("G2m",     "G2M",       x=_) |>
    gsub("Jak Stat3","JAK-STAT3",x=_) |>
    gsub("Stat5",   "STAT5",     x=_) |>
    gsub("Akt",     "AKT",       x=_)
}

# ── B-4. Full Proteome column ──────────────────────────────────
full_col <- gsea_full |>
  dplyr::filter(ID %in% show_pws) |>
  dplyr::mutate(
    organ       = "full",
    organ_label = "Full Proteome",
    NES_show    = dplyr::if_else(
      p.adjust < 0.05, NES, NA_real_),
    sig_mark    = dplyr::case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01  ~ "**",
      p.adjust < 0.05  ~ "*",
      TRUE             ~ ""
    ),
    gs_label    = clean_label(ID)
  )

# ── B-5. Organ column ──────────────────────────────────────────
organ_col <- gsea_organ |>
  dplyr::filter(ID %in% show_pws) |>
  dplyr::mutate(
    organ_label = organ_labels[organ],
    NES_show    = dplyr::if_else(
      p.adjust < 0.05, NES, NA_real_),
    sig_mark    = dplyr::case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01  ~ "**",
      p.adjust < 0.05  ~ "*",
      TRUE             ~ ""
    ),
    gs_label    = clean_label(ID)
  )

# ── B-6. pathway sorting: Full NES  ─────────────────────────
pw_order <- full_col |>
  dplyr::arrange(dplyr::desc(NES)) |>
  dplyr::pull(gs_label) |>
  unique()

# ── B-7. merge ───────────────────────────────────────────────
plot_df <- dplyr::bind_rows(
  full_col |>
    dplyr::select(organ, organ_label,
                  gs_label, NES,
                  NES_show, sig_mark,
                  p.adjust),
  organ_col |>
    dplyr::select(organ, organ_label,
                  gs_label, NES,
                  NES_show, sig_mark,
                  p.adjust)
) |>
  dplyr::mutate(
    gs_label = factor(gs_label,
                      levels=pw_order),
    organ_label = factor(
      organ_label,
      levels=c("Full Proteome",
               organ_labels[organs]))
  )

# ── B-8. Plot ─────────────────────────────────────────────────
p_heat_gsea <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x    = organ_label,
    y    = gs_label,
    fill = NES_show
  )
) +
  # NA tile (FDR≥0.05) → gray
  ggplot2::geom_tile(
    data=\(d) dplyr::filter(d,
                            is.na(NES_show)),
    fill="grey92", color="white",
    linewidth=0.4
  ) +
  ggplot2::geom_tile(
    data=\(d) dplyr::filter(d,
                            !is.na(NES_show)),
    color="white", linewidth=0.4
  ) +
  
  # Full Proteome line
  ggplot2::annotate("segment",
                    x=1.5, xend=1.5,
                    y=0.5, yend=length(pw_order)+0.5,
                    color="grey20", linewidth=0.8
  ) +
  
  ggplot2::scale_fill_gradient2(
    low      = col_dn,
    mid      = "white",
    high     = col_up,
    midpoint = 0,
    limits   = c(-2.6, 2.6),
    na.value = "grey92",
    name     = "NES"
  ) +
  ggplot2::scale_x_discrete(
    name  = NULL,
    expand= ggplot2::expansion(
      add=c(0.5, 0.5))
  ) +
  ggplot2::scale_y_discrete(
    name=NULL, limits=rev) +
  ggplot2::labs(
    title   = NULL,
    subtitle= NULL,
    caption = paste0(
      "Pathways shown: FDR<0.05 in ",
      "full proteome or \u22652 organs. ",
      "NES: Normalized Enrichment Score."
    )
  ) +
  ggplot2::theme_minimal(base_size=10) +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(
      face="bold", size=11),
    plot.subtitle   = ggplot2::element_text(
      size=7.5, color="grey40"),
    plot.caption    = ggplot2::element_text(
      size=7, color="grey50", hjust=0),
    axis.text.x     = ggplot2::element_text(
      size=8.5, angle=35, hjust=1,
      color="black", lineheight=0.85),
    axis.text.y     = ggplot2::element_text(
      size=9, color="black"),
    legend.position  = "bottom",
    legend.key.width = ggplot2::unit(
      1.5,"cm"),
    legend.title     = ggplot2::element_text(
      size=8.5),
    legend.text      = ggplot2::element_text(
      size=8),
    panel.grid       = ggplot2::element_blank(),
    plot.background  = ggplot2::element_rect(
      fill="white", color=NA),
    plot.margin      = ggplot2::margin(
      8, 10, 8, 8)
  )

ggplot2::ggsave(
  file.path(PATH_RESULTS, "Supp_Fig_GSEA_organ_heatmap.pdf"),
  plot=p_heat_gsea, width=220, height=180,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Supp_Fig_GSEA_organ_heatmap.png"),
  plot=p_heat_gsea, width=220, height=180,
  units="mm", dpi=300)
cat("Saved: Supp_Fig_GSEA_organ_heatmap\n")

# =============================================================
# Combined: Supp Fig A (protein) + Supp Fig B (GSEA)
# =============================================================

# tag add
p_heat_prot <- p_heat_prot +
  ggplot2::labs(tag="a")

p_heat_gsea <- p_heat_gsea +
  ggplot2::labs(tag="b")

supp_combined <- p_heat_prot + p_heat_gsea +
  patchwork::plot_layout(
    ncol   = 2,
    widths = c(1, 1.2)   
  ) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill="white", color=NA),
      plot.margin     = ggplot2::margin(8,8,6,8)
    )
  )

ggplot2::ggsave(
  file.path(PATH_RESULTS, "Supp_Fig_organ_combined.pdf"),
  plot  = supp_combined,
  width = 420, height = 260,
  units = "mm", device = cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Supp_Fig_organ_combined.png"),
  plot  = supp_combined,
  width = 420, height = 260,
  units = "mm", dpi   = 300)
cat("Saved: Supp_Fig_organ_combined\n")
