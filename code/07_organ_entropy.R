# ==============================================================================
# STUDY       : SOMMA
# PROJECT     : SOMMA SomaLogic
# PROGRAM     : 07_organ_entropy.R
#
# PURPOSE     : Estimate associations between organ-specific proteomic entropy
#               and phenotypic outcomes, and produce Figure 3 and Supplementary
#               Figure 3.
#
# INPUTS      : - meta_df (04_analysis_dataset.R)
#               - 09b_entropy_measures_organ_box.csv (organ-specific DISCO;
#                 see NOTES)
#
# OUTPUTS     : - Fig_organ_dot_Q2Q3Q4_BOX.png / .pdf   (Figure 3a)
#               - Fig_organ_dot_continuous_BOX.png / .pdf (Figure 3b)
#               - Table_Fig3b_organ_outcome_continuous.csv
#               - In-memory object analysis_df, reused by 09_survival_analysis.R
#
# DATA ACCESS : Individual-level SOMMA data are not included in this
#               repository. Access to the underlying study data is required
#               to reproduce the analyses.
#
# DEPENDENCIES: - 99_figure_setup.R
#               - 04_analysis_dataset.R (object meta_df)
#               - ggplot2, ggtext (optional), survival
#
# VERSION     : 1.0
# DATE        : 2026-08-15
# AUTHOR      : Namki Hong
#
# REVISION HISTORY
# 08/15/26    Creation (split out of the combined analysis script)
#
# NOTES       : Nine organ-specific entropy measures (blood vessel, bone, brain,
#               heart, kidney, liver, lung, lymphoid, skeletal muscle) are compared
#               with the full-proteome measure. Models are adjusted for age, sex,
#               clinic site, body weight and height.
#
#               Requires 09b_entropy_measures_organ_box.csv, an organ-panel version
#               of the 03_calc_entropy.R output produced outside this repository;
#               see README.md.
#
#               The two blocks labelled 'FIX' re-join the organ-specific continuous
#               measures and rebuild the standardized entropy variables before
#               Figure 3b is recomputed. They are retained from the original script
#               so that the published figure can be reproduced exactly.
# ==============================================================================
source("code/99_figure_setup.R")

stopifnot(exists("meta_df"))


# =============================================================================
# FIGURE 3: ORGAN SPECIFIC ENTROPY AND OUTCOMES & SUPPLEMENTARY FIGURE 3
# =============================================================================

# ── 0.  ──────────────────────────────────────────────
organs <- c("blood","bone","brain","heart",
            "kidney","liver","lung",
            "lymphoid","skeletal")

organ_order <- c(
  "Full\nProteome",
  "Blood\nVessel","Bone","Brain","Heart",
  "Kidney","Liver","Lung","Lymphoid",
  "Skeletal\nMuscle"
)

# ── CHANGE: "VO2 peak" -> "VO2peak" ───────────────

use_md  <- requireNamespace("ggtext", quietly = TRUE)
lab_vo2 <- if (use_md) "VO<sub>2</sub>peak (mL/min)" else
  "VO\u2082peak (mL/min)"
axis_y_elem <- if (use_md)
  ggtext::element_markdown(size = 8.5, color = "black") else
    ggplot2::element_text(size = 8.5, color = "black")

outcome_order_fig <- c(
  lab_vo2,                            # <- CHANGE
  "Walking speed (m/s)",
  "Leg peak power (W)",
  "Daily step count",
  "DSST score",                       # ← Digit Symbol → DSST score
  "Vigor-to-Frailty",
  "Max OxPhos (pmol/s/mg)",
  "Cost-capacity ratio"
)

outcomes_cont <- list(
  VO2    = "TTPKVO2U",
  Speed  = "NF400MPACE",
  Power  = "LEPEAKPWR2",
  Steps  = "ACSCFUAV",
  Digit  = "DSCORR",
  VF     = "FT0V2FN",
  OxPhos = "REMOXPHOS",
  CCR    = "TTSSPKVO2",
  GDF15  = "SOGDF15",
  CysC   = "SOCYSC"
)

outcome_labels_tbl <- c(
  VO2    = lab_vo2,                   # <- CHANGE
  Speed  = "Walking speed (m/s)",
  Power  = "Leg peak power (W)",
  Steps  = "Daily step count",
  Digit  = "DSST score",              # ← Digit Symbol → DSST score
  VF     = "Vigor-to-Frailty",
  OxPhos = "Max OxPhos (pmol/s/mg)",
  CCR    = "Cost-capacity ratio",
  GDF15  = "GDF-15 (log)",
  CysC   = "Cystatin-C (log)"
)

base_cov  <- "D1AGE2 + EEFEMALE + SITE_f + HWWGT + HWHGT"
col_higher <- "#C0392B"
col_lower  <- "#2471A3"
alpha_sig  <- 0.88

dir.create(PATH_RESULTS, showWarnings=FALSE)

# ── 1. BOX  ─────────────────────
organ_box_raw <- readr::read_csv(
  file.path(PATH_OUTPUT, "09b_entropy_measures_organ_box.csv"),
  show_col_types=FALSE
)

# quartile/continuous
q_vars_box    <- paste0("BOX_DISCO_JR_Q_", organs)
cont_vars_box <- paste0("BOX_DISCO_JR_", organs)

analysis_df <- meta_df |>
  dplyr::left_join(
    organ_box_raw |>
      dplyr::select(ID,
                    all_of(q_vars_box),
                    all_of(cont_vars_box)),
    by = "ID"
  ) |>
  dplyr::mutate(
    SITE_f    = factor(SITE),
    GDF15_std = (SOGDF15 - mean(SOGDF15, na.rm=TRUE)) /
      sd(SOGDF15, na.rm=TRUE),
    CysC_std  = (SOCYSC  - mean(SOCYSC,  na.rm=TRUE)) /
      sd(SOCYSC,  na.rm=TRUE),
    # Full proteome standardized
    BOX_DISCO_STD = (BOX_DISCO_JR -
                       mean(BOX_DISCO_JR, na.rm=TRUE)) /
      sd(BOX_DISCO_JR, na.rm=TRUE)
  )

# quartile factor
for (qv in q_vars_box) {
  analysis_df[[qv]] <- factor(
    analysis_df[[qv]], levels=1:4)
}

# organ continuous standardize
for (organ_nm in organs) {
  cv <- paste0("BOX_DISCO_JR_", organ_nm)
  sv <- paste0(cv, "_std")
  analysis_df[[sv]] <-
    (analysis_df[[cv]] -
       mean(analysis_df[[cv]], na.rm=TRUE)) /
    sd(analysis_df[[cv]], na.rm=TRUE)
}

cat("Analysis n =", nrow(analysis_df), "\n")

# ── 2. Helper  ────────────────────────────────────────────
# 2a. Quartile β (Q2/Q3/Q4 vs Q1)
get_beta_all_q <- function(outcome_var, q_var,
                           data, extra="",
                           unadj=FALSE) {
  d <- data |>
    dplyr::filter(!is.na(.data[[outcome_var]]),
                  !is.na(.data[[q_var]]))
  if (nrow(d) < 50) return(NULL)
  sd_out <- sd(d[[outcome_var]], na.rm=TRUE)
  fml <- if (unadj) {
    as.formula(paste0(outcome_var, " ~ ", q_var))
  } else {
    as.formula(paste0(
      outcome_var, " ~ ", q_var, " + ", base_cov,
      if (nchar(extra)>0) paste0(" + ",extra) else ""))
  }
  fit <- tryCatch(lm(fml, data=d),
                  error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  coef_tbl <- as.data.frame(summary(fit)$coefficients)
  coef_tbl$term <- rownames(coef_tbl)
  q_rows <- coef_tbl |>
    dplyr::filter(grepl(paste0(q_var,"[234]"), term)) |>
    dplyr::mutate(
      quartile = paste0("Q", sub(
        paste0(".*",q_var), "", term))
    )
  if (nrow(q_rows)==0) return(NULL)
  data.frame(
    quartile = q_rows$quartile,
    beta     = q_rows$Estimate / sd_out,
    beta_lo  = (q_rows$Estimate -
                  1.96*q_rows$`Std. Error`) / sd_out,
    beta_hi  = (q_rows$Estimate +
                  1.96*q_rows$`Std. Error`) / sd_out,
    p        = q_rows$`Pr(>|t|)`,
    n        = nrow(d)
  )
}

# 2b. Continuous β (per-1-SD)
get_beta_cont <- function(outcome_var, disco_var,
                          data) {
  d <- data |>
    dplyr::filter(!is.na(.data[[outcome_var]]),
                  !is.na(.data[[disco_var]]))
  if (nrow(d) < 50) return(NULL)
  sd_out   <- sd(d[[outcome_var]], na.rm=TRUE)
  sd_disco <- sd(d[[disco_var]],   na.rm=TRUE)
  fml <- as.formula(paste0(
    outcome_var, " ~ ", disco_var,
    " + ", base_cov))
  fit <- tryCatch(lm(fml, data=d),
                  error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  coef_tbl <- as.data.frame(summary(fit)$coefficients)
  coef_tbl$term <- rownames(coef_tbl)
  row <- coef_tbl |>
    dplyr::filter(term == disco_var)
  if (nrow(row)==0) return(NULL)
  data.frame(
    beta    = row$Estimate * sd_disco / sd_out,
    beta_lo = (row$Estimate -
                 1.96*row$`Std. Error`) *
      sd_disco / sd_out,
    beta_hi = (row$Estimate +
                 1.96*row$`Std. Error`) *
      sd_disco / sd_out,
    p       = row$`Pr(>|t|)`,
    n       = nrow(d)
  )
}

# Organ recode
recode_organ <- function(x) {
  dplyr::recode(x,
                "blood"    = "Blood\nVessel",
                "bone"     = "Bone",
                "brain"    = "Brain",
                "heart"    = "Heart",
                "kidney"   = "Kidney",
                "liver"    = "Liver",
                "lung"     = "Lung",
                "lymphoid" = "Lymphoid",
                "skeletal" = "Skeletal\nMuscle",
                "full"     = "Full\nProteome"
  )
}

# ── 3. Organ-specific quartile β  (BOX) ──────────────────
cat("Computing organ-specific quartile β (BOX)...\n")
results_baseline_list <- list()

for (organ_nm in organs) {
  q_var <- paste0("BOX_DISCO_JR_Q_", organ_nm)
  for (nm in names(outcomes_cont)) {
    res <- get_beta_all_q(
      outcomes_cont[[nm]], q_var, analysis_df)
    if (!is.null(res)) {
      results_baseline_list[[
        length(results_baseline_list)+1]] <-
        data.frame(
          organ         = organ_nm,
          outcome       = nm,
          outcome_label = outcome_labels_tbl[nm],
          model         = "Base",
          res
        )
    }
  }
}

results_baseline_df <- dplyr::bind_rows(
  results_baseline_list) |>
  dplyr::mutate(
    sig       = p < 0.05,
    direction = dplyr::if_else(beta >= 0,
                               "Higher","Lower"),
    abs_beta  = abs(beta)
  )
cat("Baseline rows:", nrow(results_baseline_df), "\n")

# ── 4. Continuous β (BOX) ───────────────────────────────
cat("Computing continuous β (BOX)...\n")
results_cont_list <- list()

# Full proteome
for (nm in names(outcomes_cont)) {
  res <- get_beta_cont(
    outcomes_cont[[nm]],
    "BOX_DISCO_STD", analysis_df)
  if (!is.null(res)) {
    results_cont_list[[
      length(results_cont_list)+1]] <-
      data.frame(organ="full", outcome=nm,
                 outcome_label=outcome_labels_tbl[nm],
                 res)
  }
}

# Organ-specific
for (organ_nm in organs) {
  sv <- paste0("BOX_DISCO_JR_", organ_nm, "_std")
  for (nm in names(outcomes_cont)) {
    res <- get_beta_cont(
      outcomes_cont[[nm]], sv, analysis_df)
    if (!is.null(res)) {
      results_cont_list[[
        length(results_cont_list)+1]] <-
        data.frame(organ=organ_nm, outcome=nm,
                   outcome_label=outcome_labels_tbl[nm],
                   res)
    }
  }
}

results_cont_df <- dplyr::bind_rows(
  results_cont_list) |>
  dplyr::mutate(
    sig       = p < 0.05,
    direction = dplyr::if_else(beta >= 0,
                               "Higher","Lower"),
    abs_beta  = abs(beta)
  )
cat("Continuous rows:", nrow(results_cont_df), "\n")

# ── 5. Full proteome quartile β (std_df2) ────────────────────
full_q_df <- std_df2 |>
  dplyr::mutate(
    organ         = "full",
    organ_label   = "Full\nProteome",
    quartile      = paste0("Q", quartile),
    p             = p_value,
    sig           = p_value < 0.05,
    direction     = dplyr::if_else(beta >= 0,
                                   "Higher","Lower"),
    abs_beta      = abs(beta),
    outcome_label = unname(outcome_labels_tbl[outcome])
  ) |>
  dplyr::select(organ, organ_label, outcome_label,
                outcome, quartile, beta, beta_lo,
                beta_hi, p, sig, direction, abs_beta)

# ── 6. Supp data prep ─────────────────────────────────
make_supp_data <- function(q) {
  dplyr::bind_rows(
    full_q_df |>
      dplyr::filter(quartile == q) |>
      dplyr::mutate(
        organ_label   = factor("Full\nProteome",
                               levels=organ_order),
        outcome_label = factor(
          as.character(outcome_label),
          levels=outcome_order_fig)
      ),
    results_baseline_df |>
      dplyr::filter(model=="Base", quartile==q) |>
      dplyr::mutate(
        organ_label = factor(
          recode_organ(organ),
          levels=organ_order),
        outcome_label = factor(
          dplyr::recode(
            as.character(outcome_label),
            "VO2 peak (mL/min)" = lab_vo2),   # <- CHANGE
          levels=outcome_order_fig)
      )
  ) |>
    dplyr::filter(!is.na(outcome_label),
                  !is.na(organ_label))
}

supp_q2 <- make_supp_data("Q2")
supp_q3 <- make_supp_data("Q3")

# Continuous supp
cont_supp <- dplyr::bind_rows(
  results_cont_df |>
    dplyr::filter(organ=="full") |>
    dplyr::mutate(
      organ_label   = "Full\nProteome",
      outcome_label = dplyr::recode(
        as.character(outcome_label),
        "VO2 peak (mL/min)" = lab_vo2)        # <- CHANGE
    ),
  results_cont_df |>
    dplyr::filter(organ != "full") |>
    dplyr::mutate(
      organ_label   = recode_organ(organ),
      outcome_label = dplyr::recode(
        as.character(outcome_label),
        "VO2 peak (mL/min)" = lab_vo2)        # <- CHANGE
    )
) |>
  dplyr::mutate(
    organ_label   = factor(organ_label,
                           levels=organ_order),
    outcome_label = factor(outcome_label,
                           levels=outcome_order_fig)
  ) |>
  dplyr::filter(!is.na(outcome_label),
                !is.na(organ_label))

cat("supp_q2:", nrow(supp_q2), "\n")
cat("supp_q3:", nrow(supp_q3), "\n")
cat("cont_supp:", nrow(cont_supp), "\n")

# ── 7. Figure 1: Q2/Q3/Q4 Dot plot ───────────────────────────
dot_df <- dplyr::bind_rows(
  full_q_df |>
    dplyr::filter(
      quartile %in% c("Q2","Q3","Q4"),
      !outcome_label %in% c("GDF-15 (log)",
                            "Cystatin-C (log)")
    ) |>
    dplyr::mutate(organ_label="Full\nProteome"),
  results_baseline_df |>
    dplyr::filter(
      model    == "Base",
      quartile %in% c("Q2","Q3","Q4"),
      !outcome_label %in% c("GDF-15 (log)",
                            "Cystatin-C (log)")
    ) |>
    dplyr::mutate(
      organ_label = recode_organ(organ),
      outcome_label = dplyr::recode(
        as.character(outcome_label),
        "VO2 peak (mL/min)" = lab_vo2)        # <- CHANGE
    )
) |>
  dplyr::mutate(
    organ_label = factor(organ_label,
                         levels=organ_order),
    outcome_label = factor(outcome_label,
                           levels=outcome_order_fig),
    quartile  = factor(quartile,
                       levels=c("Q2","Q3","Q4")),
    direction = dplyr::if_else(beta >= 0,
                               "Higher","Lower"),
    sig       = p < 0.05,
    abs_beta  = abs(beta),
    q_offset  = dplyr::case_when(
      quartile=="Q2" ~ -0.27,
      quartile=="Q3" ~  0.00,
      quartile=="Q4" ~  0.27
    ),
    organ_num = as.numeric(organ_label) + q_offset
  ) |>
  dplyr::filter(!is.na(outcome_label),
                !is.na(organ_label))

q_label_df <- dot_df |>
  dplyr::distinct(organ_label, quartile, organ_num) |>
  dplyr::mutate(q_txt=as.character(quartile)) |>
  dplyr::arrange(organ_label, quartile)

n_organs <- length(levels(dot_df$organ_label))
x_breaks <- seq_len(n_organs)
x_labels <- levels(dot_df$organ_label)

p_dot_q234 <- ggplot2::ggplot(
  dot_df,
  ggplot2::aes(x=organ_num, y=outcome_label)
) +
  ggplot2::annotate("rect",
                    xmin=0.55, xmax=1.45, ymin=0.5, ymax=8.5,
                    fill="grey93", alpha=0.7) +
  ggplot2::annotate("rect",
                    xmin=seq(2.5, n_organs-0.5, 2),
                    xmax=seq(3.5, n_organs+0.5, 2),
                    ymin=0.5, ymax=8.5,
                    fill="grey96", alpha=0.6) +
  ggplot2::annotate("segment",
                    x=1.5, xend=1.5, y=0.5, yend=8.5,
                    color="grey30", linewidth=0.6) +
  ggplot2::annotate("segment",
                    x    = as.numeric(outer(2:n_organs,
                                            c(-0.135,0.135),"+")),
                    xend = as.numeric(outer(2:n_organs,
                                            c(-0.135,0.135),"+")),
                    y=0.5, yend=8.5,
                    color="grey85", linewidth=0.2,
                    linetype="dotted") +
  ggplot2::annotate("segment",
                    x=c(0.865,1.135), xend=c(0.865,1.135),
                    y=0.5, yend=8.5,
                    color="grey75", linewidth=0.2,
                    linetype="dotted") +
  ggplot2::geom_vline(
    xintercept=seq(1.5, n_organs-0.5, 1),
    color="grey80", linewidth=0.25) +
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d, !sig),
    ggplot2::aes(size=abs_beta, color=direction),
    shape=1, stroke=0.5, alpha=0.25) +
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d,
                            sig & direction=="Higher"),
    ggplot2::aes(size=abs_beta, color=direction),
    shape=16, alpha=alpha_sig) +
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d,
                            sig & direction=="Lower"),
    ggplot2::aes(size=abs_beta, color=direction),
    shape=16, alpha=alpha_sig) +
  ggplot2::geom_text(
    data=\(d) dplyr::filter(d, sig),
    ggplot2::aes(x=organ_num, y=outcome_label,
                 label=sprintf("%.2f", beta),
                 color=direction),
    size=2.0, vjust=3.4, hjust=0.5,
    fontface="plain",
    show.legend=FALSE,
    inherit.aes=FALSE) +
  ggplot2::geom_text(
    data=q_label_df,
    ggplot2::aes(x=organ_num, y=0.22,
                 label=q_txt),
    size=1.8, color="grey50", hjust=0.5,
    show.legend=FALSE, inherit.aes=FALSE) +
  ggplot2::annotate("text",
                    x=n_organs+0.6, y=c(2.0,5.5,7.5),
                    label=c("Mitochondrial\nfunction",
                            "Cognitive &\nfrailty",
                            "Physical\nperformance"),
                    size=2.3, hjust=0,
                    color=c("#BA7517","#534AB7","#0F6E56"),
                    fontface="bold", lineheight=0.85) +
  ggplot2::scale_color_manual(
    values=c("Higher"=col_higher,
             "Lower" =col_lower),
    name  =expression(beta~"direction (vs Q1)"),
    labels=c("Higher","Lower"),
    guide =ggplot2::guide_legend(
      order=1,
      override.aes=list(
        shape=16, size=4, alpha=0.9))
  ) +
  ggplot2::scale_size_continuous(
    name  =expression("|"*beta*"| (SD units)"),
    range =c(0.8, 7.0),
    breaks=c(0.1, 0.2, 0.3, 0.5),
    labels=c("0.1","0.2","0.3","0.5")) +
  ggplot2::scale_x_continuous(
    breaks=x_breaks, labels=x_labels,
    limits=c(0.4, n_organs+0.5),
    expand=ggplot2::expansion(mult=c(0,0)),
    name=paste0(
      "Organ-specific proteomic entropy\n",
      "(Q2 / Q3 / Q4 vs Q1)")
  ) +
  ggplot2::scale_y_discrete(name=NULL, limits=rev) +
  ggplot2::labs(
  ) +
  ggplot2::coord_cartesian(clip="off") +
  ggplot2::theme_classic(base_size=10) +
  ggplot2::theme(
    axis.line        = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(
      size=7.5, color="black", angle=35,
      hjust=1, lineheight=0.85),
    axis.text.y      = axis_y_elem,            # <- CHANGE
    axis.title.x     = ggplot2::element_text(
      size=9, margin=ggplot2::margin(t=10)),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    legend.text      = ggplot2::element_text(size=8),
    legend.title     = ggplot2::element_text(size=8.5),
    legend.key.size  = ggplot2::unit(0.4,"cm"),
    panel.grid       = ggplot2::element_blank(),
    plot.background  = ggplot2::element_rect(
      fill="white", color=NA),
    plot.caption     = ggplot2::element_text(
      size=6.5, color="grey40", hjust=0,
      margin=ggplot2::margin(t=5)),
    plot.margin      = ggplot2::margin(8,60,10,5)) +
  ggplot2::guides(
    color=ggplot2::guide_legend(order=1),
    size =ggplot2::guide_legend(order=2))

ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_organ_dot_Q2Q3Q4_BOX.pdf"),
  plot=p_dot_q234, width=220, height=140,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_organ_dot_Q2Q3Q4_BOX.png"),
  plot=p_dot_q234, width=220, height=140,
  units="mm", dpi=300)
cat("Saved: Fig_organ_dot_Q2Q3Q4_BOX\n")

# ── 8. Figure 3: Continuous Dot plot ─────────────────────────
dot_cont_df <- results_cont_df |>
  dplyr::filter(
    !outcome_label %in% c("GDF-15 (log)",
                          "Cystatin-C (log)")
  ) |>
  dplyr::mutate(
    organ_label   = factor(recode_organ(organ),
                           levels=organ_order),
    outcome_label = factor(
      as.character(outcome_label),
      levels=outcome_order_fig),
    organ_num     = as.numeric(organ_label)
  ) |>
  dplyr::filter(!is.na(outcome_label),
                !is.na(organ_label))

# check
cat("Continuous outcome levels:\n")
print(levels(dot_cont_df$outcome_label))

p_dot_cont <- ggplot2::ggplot(
  dot_cont_df,
  ggplot2::aes(x=organ_num, y=outcome_label)
) +
  ggplot2::annotate("rect",
                    xmin=0.55, xmax=1.45, ymin=0.5, ymax=8.5,
                    fill="grey93", alpha=0.7) +
  ggplot2::annotate("rect",
                    xmin=seq(2.5, n_organs-0.5, 2),
                    xmax=seq(3.5, n_organs+0.5, 2),
                    ymin=0.5, ymax=8.5,
                    fill="grey96", alpha=0.6) +
  ggplot2::annotate("segment",
                    x=1.5, xend=1.5, y=0.5, yend=8.5,
                    color="grey30", linewidth=0.6) +
  ggplot2::geom_vline(
    xintercept=seq(1.5, n_organs-0.5, 1),
    color="grey80", linewidth=0.25) +
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d, !sig),
    ggplot2::aes(size=abs_beta, color=direction),
    shape=1, stroke=0.5, alpha=0.25) +
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d,
                            sig & direction=="Higher"),
    ggplot2::aes(size=abs_beta, color=direction),
    shape=16, alpha=alpha_sig) +
  ggplot2::geom_point(
    data=\(d) dplyr::filter(d,
                            sig & direction=="Lower"),
    ggplot2::aes(size=abs_beta, color=direction),
    shape=16, alpha=alpha_sig) +
  ggplot2::geom_text(
    data=\(d) dplyr::filter(d, sig),
    ggplot2::aes(x=organ_num, y=outcome_label,
                 label=sprintf("%.2f", beta),
                 color=direction),
    size=2.0, vjust=3.4, hjust=0.5,
    show.legend=FALSE, inherit.aes=FALSE) +
  ggplot2::annotate("text",
                    x=n_organs+0.6, y=c(2.0,5.5,7.5),
                    label=c("Mitochondrial\nfunction",
                            "Cognitive &\nfrailty",
                            "Physical\nperformance"),
                    size=2.3, hjust=0,
                    color=c("#BA7517","#534AB7","#0F6E56"),
                    fontface="bold", lineheight=0.85) +
  ggplot2::scale_color_manual(
    values=c("Higher"=col_higher,
             "Lower" =col_lower),
    name  =expression(beta~"direction"),
    labels=c("Higher","Lower"),
    guide =ggplot2::guide_legend(
      order=1,
      override.aes=list(
        shape=16, size=4, alpha=0.9))
  ) +
  ggplot2::scale_size_continuous(
    name  =expression("|"*beta*"| (SD units)"),
    range =c(0.8, 7.0),
    breaks=c(0.05, 0.1, 0.2, 0.3),
    labels=c("0.05","0.1","0.2","0.3")) +
  ggplot2::scale_x_continuous(
    breaks=x_breaks, labels=x_labels,
    limits=c(0.4, n_organs+0.5),
    expand=ggplot2::expansion(mult=c(0,0)),
    name=paste0(
      "Organ-specific proteomic entropy\n",
      "(per 1-SD continuous)")
  ) +
  ggplot2::scale_y_discrete(    name   = NULL,
                                limits = rev(outcome_order_fig)
  ) +
  ggplot2::labs(
    caption=paste0()
  ) +
  ggplot2::coord_cartesian(clip="off") +
  ggplot2::theme_classic(base_size=10) +
  ggplot2::theme(
    axis.line        = ggplot2::element_blank(),
    axis.ticks       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(
      size=7.5, color="black", angle=35,
      hjust=1, lineheight=0.85),
    axis.text.y      = axis_y_elem,            # <- CHANGE
    axis.title.x     = ggplot2::element_text(
      size=9, margin=ggplot2::margin(t=10)),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    legend.text      = ggplot2::element_text(size=8),
    legend.title     = ggplot2::element_text(size=8.5),
    legend.key.size  = ggplot2::unit(0.4,"cm"),
    panel.grid       = ggplot2::element_blank(),
    plot.background  = ggplot2::element_rect(
      fill="white", color=NA),
    plot.caption     = ggplot2::element_text(
      size=6.5, color="grey40", hjust=0,
      margin=ggplot2::margin(t=5)),
    plot.margin      = ggplot2::margin(8,60,10,5)) +
  ggplot2::guides(
    color=ggplot2::guide_legend(order=1),
    size =ggplot2::guide_legend(order=2))

ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_organ_dot_continuous_BOX.pdf"),
  plot=p_dot_cont, width=200, height=130,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_organ_dot_continuous_BOX.png"),
  plot=p_dot_cont, width=200, height=130,
  units="mm", dpi=300)
cat("Saved: Fig_organ_dot_continuous_BOX\n")

cat("\n=== complete ===\n")
cat("Exported files:\n")
cat("  Fig_organ_dot_Q2Q3Q4_BOX.png\n")
cat("  Fig_organ_dot_continuous_BOX.png\n")

# =============================================================================
# FIX: organ-specific continuous DISCO variable- analysis_df re-join
# =============================================================================
# 1) Continuous DISCO omit recheck
cont_vars_box <- paste0("BOX_DISCO_JR_", organs)
missing_cont <- setdiff(cont_vars_box, names(analysis_df))
cat("count omitted continuous DISCO var:\n")
print(missing_cont)

# 2) organ_box_raw to analysis df join
if (length(missing_cont) > 0) {
  analysis_df <- analysis_df |>
    dplyr::left_join(
      organ_box_raw |>
        dplyr::select(ID, dplyr::all_of(missing_cont)),
      by = "ID"
    )
  cat("\n=== analysis_df var check ===\n")
  print(grep("BOX_DISCO_JR_[a-z]+$", names(analysis_df),
             value = TRUE))
}

# =============================================================================
# FIX 2: BOX_DISCO_STD (full proteome standardized) var recreate
# =============================================================================
# Full proteome standardized var recreate
analysis_df$BOX_DISCO_STD <-
  (analysis_df$BOX_DISCO_JR -
     mean(analysis_df$BOX_DISCO_JR, na.rm = TRUE)) /
  sd(analysis_df$BOX_DISCO_JR, na.rm = TRUE)

cat("BOX_DISCO_STD recreate\n")
cat("Mean:", round(mean(analysis_df$BOX_DISCO_STD, na.rm=TRUE), 4), "\n")
cat("SD:",   round(sd(analysis_df$BOX_DISCO_STD,   na.rm=TRUE), 4), "\n")
cat("N non-missing:", sum(!is.na(analysis_df$BOX_DISCO_STD)), "\n")

# Organ-specific standardized recreate
for (organ_nm in organs) {
  cv <- paste0("BOX_DISCO_JR_", organ_nm)
  sv <- paste0(cv, "_std")
  if (cv %in% names(analysis_df)) {
    analysis_df[[sv]] <-
      (analysis_df[[cv]] -
         mean(analysis_df[[cv]], na.rm = TRUE)) /
      sd(analysis_df[[cv]], na.rm = TRUE)
  } else {
    cat("WARNING:", cv, "still missing\n")
  }
}

cat("\n=== Standardized var list ===\n")
print(grep("_std$|BOX_DISCO_STD", names(analysis_df),
           value = TRUE))

# =============================================================================
#  Section (B) Figure 3b recalc
# =============================================================================
cat("\nRecomputing continuous β...\n")
results_cont_list <- list()

# Full proteome
for (nm in names(outcomes_cont)) {
  res <- get_beta_cont(outcomes_cont[[nm]],
                       "BOX_DISCO_STD", analysis_df)
  if (!is.null(res)) {
    results_cont_list[[length(results_cont_list)+1]] <-
      data.frame(organ = "full", outcome = nm,
                 outcome_label = outcome_labels_tbl[nm], res)
  }
}

# Organ-specific
for (organ_nm in organs) {
  sv <- paste0("BOX_DISCO_JR_", organ_nm, "_std")
  for (nm in names(outcomes_cont)) {
    res <- get_beta_cont(outcomes_cont[[nm]], sv, analysis_df)
    if (!is.null(res)) {
      results_cont_list[[length(results_cont_list)+1]] <-
        data.frame(organ = organ_nm, outcome = nm,
                   outcome_label = outcome_labels_tbl[nm], res)
    }
  }
}

results_cont_df <- dplyr::bind_rows(results_cont_list) |>
  dplyr::mutate(
    sig       = p < 0.05,
    direction = dplyr::if_else(beta >= 0, "Higher", "Lower"),
    abs_beta  = abs(beta)
  )

cat("Continuous rows:", nrow(results_cont_df), "\n")
cat("Expected: 10 organs × 10 outcomes = 100\n")

# Export Figure 3b
fig3b_export <- results_cont_df |>
  dplyr::filter(!outcome %in% c("GDF15", "CysC")) |>
  dplyr::mutate(
    organ_label = ifelse(organ == "full",
                         "Full Proteome",
                         recode_organ(organ)),
    organ_label = gsub("\n", " ", organ_label),
    across(c(beta, beta_lo, beta_hi), ~ round(.x, 3)),
    p_value = format.pval(p, digits = 3, eps = 0.001)
  ) |>
  dplyr::select(organ_label, outcome, outcome_label,
                beta, beta_lo, beta_hi, p_value, n)

write.csv(fig3b_export,
          file.path(PATH_RESULTS, "Table_Fig3b_organ_outcome_continuous.csv"),
          row.names = FALSE)
cat("Saved: Table_Fig3b_organ_outcome_continuous.csv\n")
