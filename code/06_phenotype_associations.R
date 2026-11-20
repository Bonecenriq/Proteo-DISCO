# ==============================================================================
# STUDY       : SOMMA
# PROJECT     : SOMMA SomaLogic
# PROGRAM     : 06_phenotype_associations.R
#
# PURPOSE     : Estimate associations between whole-proteome entropy (DISCO) and
#               phenotypic outcomes, and produce Figure 2 and Supplementary
#               Figures 1 and 2.
#
# INPUTS      : - m, meta_df (04_analysis_dataset.R)
#               - 05z_entropy_measures_unweighted.csv (age-unweighted DISCO;
#                 see NOTES)
#
# OUTPUTS     : - Fig2_v9.png                    (Figure 2: scatter, forest, biomarkers)
#               - Fig2_panel_b_base.png
#               - Fig2_panel_b_GDF15.png
#               - Fig2_panel_b_CysC.png
#               - Fig2_sensitivity_3panel.png
#               - FigS_sensitivity_2panel.png    (Supplementary Figure 1)
#               - SuppFig_unweighted_DISCO.png / .pdf   (Supplementary Figure 2)
#               - SuppFig_unweighted_DISCO_stdbeta.csv
#               - In-memory objects reused by 09_survival_analysis.R:
#                   outcome_labels, outcome_order, domain_colors, q_colors_main,
#                   q_shapes_v5, domain_label_df, make_forest_b(), save_fig()
#
# DATA ACCESS : Individual-level SOMMA data are not included in this
#               repository. Access to the underlying study data is required
#               to reproduce the analyses.
#
# DEPENDENCIES: - 99_figure_setup.R
#               - 04_analysis_dataset.R (objects m, meta_df)
#               - ggplot2, cowplot, ggtext (optional, for subscript axis labels)
#
# VERSION     : 1.0
# DATE        : 2026-08-15
# AUTHOR      : Namki Hong
#
# REVISION HISTORY
# 08/15/26    Creation (split out of the combined analysis script)
#
# NOTES       : All models are linear regressions adjusted for age, sex, clinic
#               site, body weight and height. Betas are expressed in SD units of
#               each outcome, with quartile 1 as the reference.
#
#               Supplementary Figure 2 requires 05z_entropy_measures_unweighted.csv,
#               an age-unweighted DISCO version of the 03_calc_entropy.R output.
#               That file is produced outside this repository; see README.md.
#
#               get_std_beta() is redefined in the Supplementary Figure 2 block so
#               that the quartile variable can be parameterised. The original
#               definition is restored by re-running the Figure 2 block.
# ==============================================================================
source("code/99_figure_setup.R")

stopifnot(exists("m"), exists("meta_df"))


# =============================================================================
# FIGURE 2: PHENOTYPE ASSOCIATION & SUPPLEMENTARY FIGURE 1   [v9 - full]
#   1) "VO2 peak" -> "VO2peak" (no space; ggtext subscript)
#   2) "Mitochondrial function" -> "Energetic capacity and efficiency"
#   3) Panel b: "P for trend < 0.001 for all outcomes" note inside
#      the green (Physical performance) band
#   4) FIX: explicit scale_y_discrete() (avoids "Discrete value
#      supplied to a continuous scale")
# =============================================================================

# --- 1a. get_std_beta function ---
get_std_beta <- function(outcome_var, data,
                         extra_str = "") {
  sd_out <- sd(data[[outcome_var]], na.rm=TRUE)
  data2  <- data |>
    dplyr::filter(
      !is.na(.data[[outcome_var]])) |>
    dplyr::mutate(
      Q_fac     = factor(BOX_DISCO_JR_Q,
                         levels = 1:4),
      GDF15_std = (SOGDF15 -
                     mean(SOGDF15,
                          na.rm=TRUE)) /
        sd(SOGDF15, na.rm=TRUE),
      CysC_std  = (SOCYSC -
                     mean(SOCYSC,
                          na.rm=TRUE)) /
        sd(SOCYSC, na.rm=TRUE)
    )
  base_cov <- paste0(
    "D1AGE2 + EEFEMALE + SITE + ",
    "HWWGT + HWHGT")
  fml <- as.formula(paste0(
    outcome_var, " ~ Q_fac + ", base_cov,
    if (nchar(extra_str) > 0)
      paste0(" + ", extra_str) else ""
  ))
  fit      <- lm(fml, data = data2)
  coef_tbl <- as.data.frame(
    summary(fit)$coefficients)
  coef_tbl$term <- rownames(coef_tbl)
  q_rows <- coef_tbl |>
    dplyr::filter(grepl("Q_fac[234]", term))
  data.frame(
    quartile = c(2, 3, 4),
    beta     = q_rows$Estimate / sd_out,
    beta_lo  = (q_rows$Estimate -
                  1.96*q_rows$`Std. Error`) /
      sd_out,
    beta_hi  = (q_rows$Estimate +
                  1.96*q_rows$`Std. Error`) /
      sd_out,
    p_value  = q_rows$`Pr(>|t|)`,
    sd_out   = sd_out
  )
}

# --- 1b. P-for-trend (quartile as ordinal 1-4) ---
get_p_trend <- function(outcome_var, data,
                        extra_str = "") {
  data2 <- data |>
    dplyr::filter(!is.na(.data[[outcome_var]]),
                  !is.na(BOX_DISCO_JR_Q)) |>
    dplyr::mutate(
      Q_num     = as.numeric(BOX_DISCO_JR_Q),
      GDF15_std = (SOGDF15 - mean(SOGDF15, na.rm=TRUE)) /
        sd(SOGDF15, na.rm=TRUE),
      CysC_std  = (SOCYSC - mean(SOCYSC, na.rm=TRUE)) /
        sd(SOCYSC, na.rm=TRUE)
    )
  fml <- as.formula(paste0(
    outcome_var, " ~ Q_num + D1AGE2 + EEFEMALE + ",
    "SITE + HWWGT + HWHGT",
    if (nchar(extra_str) > 0)
      paste0(" + ", extra_str) else ""))
  cf <- summary(lm(fml, data = data2))$coefficients
  cf["Q_num", "Pr(>|t|)"]
}

# --- Outcome lists ---
outcomes_all <- list(
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
outcomes_main <- outcomes_all[
  !names(outcomes_all) %in% c("GDF15","CysC")]

# ── CHANGE 1: "VO2 peak" -> "VO2peak" ────────────────────────

use_md <- requireNamespace("ggtext", quietly = TRUE)

lab_vo2  <- if (use_md) "VO<sub>2</sub>peak (mL/min)" else
  "VO\u2082peak (mL/min)"
lab_dsst <- if (use_md)
  "DSST score (0-133)<br>(Lower = poorer cognition)" else
    "DSST score (0-133)\n(Lower = poorer cognition)"
lab_vf   <- if (use_md)
  "Vigor-to-Frailty (0-12)<br>(Higher = more frail)" else
    "Vigor-to-Frailty (0-12)\n(Higher = more frail)"

# y element (markdown element_markdown)
axis_y_elem <- if (use_md)
  ggtext::element_markdown(size = 8.5, lineheight = 0.95) else
    element_text(size = 8.5)

outcome_labels <- c(
  VO2    = lab_vo2,
  Speed  = "Walking speed (m/s)",
  Power  = "Leg peak power (W)",
  Steps  = "Daily step count",
  Digit  = lab_dsst,
  VF     = lab_vf,
  OxPhos = "Max OxPhos (pmol/s/mg)",
  CCR    = "Cost-capacity ratio (%)",
  GDF15  = "GDF-15 (log)",
  CysC   = "Cystatin-C (log)"
)
outcome_order <- c(
  lab_vo2,                            # <- CHANGE 1
  "Walking speed (m/s)",
  "Leg peak power (W)",
  "Daily step count",
  lab_dsst,
  lab_vf,
  "Max OxPhos (pmol/s/mg)",
  "Cost-capacity ratio (%)"
)

# ── CHANGE 2: domain name ───────────────────────────────────
#   "Mitochondrial function" -> "Energetic capacity and efficiency"
DOMAIN_ENERGY <- "Energetic capacity and efficiency"

# --- Compute std_df2 ---
std_df2 <- lapply(names(outcomes_all),
                  function(nm) {
                    df <- get_std_beta(outcomes_all[[nm]], m)
                    df$outcome <- nm
                    df
                  }) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    quartile_label = factor(
      paste0("Q", quartile),
      levels = c("Q4","Q3","Q2")),
    domain = dplyr::case_when(
      outcome %in% c("VO2","Speed",
                     "Power","Steps") ~
        "Physical performance",
      outcome %in% c("Digit","VF") ~
        "Cognitive & frailty",
      outcome %in% c("OxPhos","CCR") ~
        DOMAIN_ENERGY,                      # <- CHANGE 2
      TRUE ~ "Biomarker"),
    # build as character first, then convert to factor
    outcome_label = outcome_labels[outcome],  # character
    outcome_label = factor(
      outcome_label,
      levels = c(
        rev(outcome_order),     # the 8 outcomes used in panel b
        "GDF-15 (log)",
        "Cystatin-C (log)"
      )
    )
  )

# check
cat("=== Std beta Q4 vs Q1 ===\n")
std_df2 |>
  dplyr::filter(quartile == 4) |>
  dplyr::select(outcome, outcome_label,
                beta, beta_lo, beta_hi, p_value) |>
  dplyr::mutate(across(c(beta,beta_lo,beta_hi),
                       round, 3),
                p_value = round(p_value, 4)) |>
  as.data.frame() |>
  print()

# --- Color/shape settings ---
domain_colors <- c(
  "Physical performance"              = "#0F6E56",
  "Cognitive & frailty"               = "#534AB7",
  "Energetic capacity and efficiency" = "#BA7517"   # <- CHANGE 2
)

q_colors_main <- c("Q2"="#74C476","Q3"="#FD8D3C","Q4"="#D73027")
q_fills_bio    <- c("Q2"="#74C476","Q3"="#FD8D3C","Q4"="#D73027")
q_shapes_v5  <- c("Q2"=21, "Q3"=24, "Q4"=23)
q_colors_scatter <- c(
  "1" = "#4575B4",  # blue   (Q1 lowest)
  "2" = "#74C476",  # green  (Q2)
  "3" = "#FD8D3C",  # orange (Q3)
  "4" = "#D73027"   # red    (Q4 highest)
)

# ── CHANGE 2: in-panel domain label ─────────────────────────
domain_label_df <- data.frame(
  y     = c(1.5, 3.5, 6.5),
  label = c("Energetic capacity\n& efficiency",  # <- CHANGE 2
            "Cognitive &\nfrailty",
            "Physical\nperformance"),
  color = unname(domain_colors[c(
    DOMAIN_ENERGY,
    "Cognitive & frailty",
    "Physical performance")])
)

std_main <- std_df2 |>
  dplyr::filter(domain != "Biomarker") |>
  dplyr::mutate(
    quartile_label = factor(
      paste0("Q", quartile),
      levels = c("Q4","Q3","Q2"))
  )

# ── std_bio ─────────────────────────────────────────────────
std_bio <- std_df2 |>
  dplyr::filter(domain == "Biomarker") |>
  dplyr::mutate(
    outcome_label = factor(
      as.character(outcome_label),
      levels = c("Cystatin-C (log)",
                 "GDF-15 (log)")
    ),
    quartile_label = factor(
      paste0("Q", quartile),
      levels = c("Q4","Q3","Q2"))
  )

cat("\n=== std_bio check ===\n")
std_bio |>
  dplyr::select(outcome, outcome_label,
                quartile, beta) |>
  as.data.frame() |>
  print()

# ── ADD: P-for-trend  (8 outcome) ────────────
ptrend_tbl <- data.frame(
  outcome = names(outcomes_main),
  p_trend = vapply(
    outcomes_main,
    function(v) get_p_trend(v, m),
    numeric(1))
)
cat("\n=== P for trend (Q1->Q4, adjusted) ===\n")
print(transform(ptrend_tbl,
                p_trend = format.pval(p_trend,
                                      digits = 2,
                                      eps = 1e-4)),
      row.names = FALSE)
cat("All < 0.001? ",
    all(ptrend_tbl$p_trend < 0.001), "\n")

shared_legend_b <- scale_shape_manual(
  values = q_shapes_v5,
  breaks = c("Q2","Q3","Q4"),
  labels = c("Q2","Q3","Q4 (highest entropy)"),
  name   = "Quartile vs. Q1 (reference)"
)
shared_fill_b <- scale_fill_manual(
  values = q_colors_main,
  breaks = c("Q2","Q3","Q4"),
  labels = c("Q2","Q3","Q4 (highest entropy)"),
  name   = "Quartile vs. Q1 (reference)"
)

# --- Panel a: Scatter ---
p_a <- ggplot(
  m |> dplyr::filter(!is.na(BOX_DISCO_JR_Q)),
  aes(x=D1AGE2, y=BOX_DISCO_JR,
      color=factor(BOX_DISCO_JR_Q))
) +
  geom_point(size=0.85, alpha=0.5) +
  geom_smooth(aes(group=1), method="lm",
              color="grey20", linewidth=0.9,
              se=TRUE, fill="grey70",
              alpha=0.2) +
  scale_color_manual(
    values = q_colors_scatter,
    labels = c("Q1 (lowest)","Q2",
               "Q3","Q4 (highest)"),
    name   = "Proteomic entropy\n(DISCO) quartile",
    guide  = guide_legend(
      nrow=2,
      override.aes=list(size=2.2))
  ) +
  labs(title="a", x="Age (years)",
       y="Proteomic entropy (DISCO)") +
  theme_bw(base_size=9) +
  theme(
    plot.title      = element_text(
      face="bold", size=10, hjust=0),
    legend.position = "bottom",
    legend.text     = element_text(size=7.5),
    legend.title    = element_text(
      size=8, face="bold"),
    legend.key.size = unit(3.5,"mm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color="grey92", linewidth=0.3),
    plot.margin = margin(3,5,3,3,"mm")
  )

# --- Panel b: Forest (physical / cognitive / energetic) ---
p_b <- ggplot(
  std_main,
  aes(x=beta, y=outcome_label,
      color=domain, shape=quartile_label,
      fill=quartile_label)
) +
  annotate("rect", xmin=-Inf, xmax=Inf,
           ymin=4.5, ymax=8.5,
           fill="#E1F5EE", alpha=0.5) +
  annotate("rect", xmin=-Inf, xmax=Inf,
           ymin=2.5, ymax=4.5,
           fill="#EEEDFE", alpha=0.5) +
  annotate("rect", xmin=-Inf, xmax=Inf,
           ymin=0.5, ymax=2.5,
           fill="#FAEEDA", alpha=0.5) +
  geom_vline(xintercept=0, color="grey30",
             linewidth=0.5, linetype="dashed") +
  annotate("text", x=0.02, y=Inf,
           label="Q1 (ref)", size=2.3,
           color="grey40", hjust=0, vjust=1.5,
           fontface="italic") +
  # ── P-for-trend  ────────
annotate("text", x=0.88, y=8.28,
         label="P for trend < 0.001 for all outcomes",
         size=2.4, color="grey25",
         hjust=1, vjust=0.5,
         fontface="italic") +
  geom_errorbarh(
    aes(xmin=beta_lo, xmax=beta_hi,
        color=domain),
    height=0, linewidth=0.55,
    position=position_dodge(width=0.65)
  ) +
  geom_point(
    aes(color=domain), size=2.3, stroke=0.5,
    position=position_dodge(width=0.65)
  ) +
  geom_text(
    data=domain_label_df,
    aes(x=0.88, y=y, label=label,
        color=I(color)),
    hjust=1, size=2.2, lineheight=0.85,
    fontface="bold", inherit.aes=FALSE
  ) +
  scale_color_manual(
    values=domain_colors, guide="none") +
  shared_legend_b + shared_fill_b +
  guides(
    shape=guide_legend(
      nrow=1,
      override.aes=list(size=2.5)),
    fill=guide_legend(nrow=1)
  ) +
  scale_x_continuous(
    breaks=seq(-0.8,0.8,0.2),
    limits=c(-0.9,0.9),
    labels=function(x) sprintf("%.1f",x)
  ) +
  # ── FIX
scale_y_discrete() +
  labs(title="b",
       x="Adjusted difference vs. Q1 (outcome SD units)",
       y=NULL) +
  coord_cartesian(clip="off") +
  theme_bw(base_size=9) +
  theme(
    plot.title         = element_text(
      face="bold", size=10, hjust=0),
    axis.text.y        = axis_y_elem,
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size=7.5),
    legend.title       = element_text(
      size=8, face="bold"),
    legend.key.size    = unit(4,"mm"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      color="grey88", linewidth=0.3),
    plot.margin        = margin(3,5,3,0,"mm")
  )

# --- Panel c: Biomarker ---
p_c <- ggplot(
  std_bio,
  aes(x=beta, y=outcome_label,
      color=outcome_label,
      shape=quartile_label,
      fill=quartile_label)
) +
  annotate("rect", xmin=-Inf, xmax=Inf,
           ymin=-Inf, ymax=Inf,
           fill="#FAECE7", alpha=0.4) +
  geom_vline(xintercept=0, color="grey30",
             linewidth=0.5, linetype="dashed") +
  annotate("text", x=0.02, y=Inf,
           label="Q1 (ref)", size=2.3,
           color="grey40", hjust=0, vjust=1.5,
           fontface="italic") +
  geom_errorbarh(
    aes(xmin=beta_lo, xmax=beta_hi),
    height=0, linewidth=0.55,
    position=position_dodge(width=0.65)
  ) +
  geom_point(
    size=2.3, stroke=0.5,
    position=position_dodge(width=0.65)
  ) +
  scale_color_manual(
    values=c("GDF-15 (log)"="#993C1D",
             "Cystatin-C (log)"="#BA7517"),
    guide="none"
  ) +
  scale_shape_manual(values=q_shapes_v5,
                     guide="none") +
  scale_fill_manual(values=q_fills_bio,
                    guide="none") +
  scale_x_continuous(
    breaks=seq(0,1.4,0.2),
    limits=c(-0.15,1.45),
    labels=function(x) sprintf("%.1f",x)
  ) +
  labs(title="c  Circulating biomarkers",
       x="Adjusted difference vs. Q1 (outcome SD units)",
       y=NULL) +
  theme_bw(base_size=9) +
  theme(
    plot.title         = element_text(
      face="bold", size=9, hjust=0),
    axis.text.y        = axis_y_elem,
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      color="grey88", linewidth=0.3),
    plot.margin        = margin(3,5,3,3,"mm")
  )

# --- Combine Fig2 ---
legend_b <- cowplot::get_legend(
  p_b + theme(legend.position="bottom",
              legend.direction="horizontal"))
legend_a <- cowplot::get_legend(
  p_a + theme(legend.position="bottom"))
left_top <- cowplot::plot_grid(
  p_a + theme(legend.position="none"),
  legend_a, ncol=1, rel_heights=c(1,0.15))
left_col <- cowplot::plot_grid(
  left_top,
  p_c + theme(legend.position="none"),
  ncol=1, rel_heights=c(1.7,1))
right_col <- cowplot::plot_grid(
  p_b + theme(legend.position="none"),
  legend_b, ncol=1, rel_heights=c(1,0.08))
fig2_body <- cowplot::plot_grid(
  left_col, right_col,
  ncol=2, rel_widths=c(0.32,0.68))
fig2_v9 <- cowplot::plot_grid(
  fig2_body,
  cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("Linear regression adjusted ",
             "for age, sex, clinic site, ",
             "body weight, and height. ",
             "Error bars = 95% CI. ",
             "\u03b2 = mean difference vs. ",
             "Q1 in SD units of each outcome."),
      x=0.01, y=0.5, hjust=0, vjust=0.5,
      size=6.5, color="grey40",
      fontface="plain"),
  ncol=1, rel_heights=c(1,0.06))

png(file.path(PATH_RESULTS, "Fig2_v9.png"),
    width=220/25.4, height=150/25.4,
    units="in", res=300)
print(fig2_v9)
dev.off()
cat("Saved: results_v14/Fig2_v9.png\n")

# --- Sensitivity: + GDF-15 / + CysC ---
make_std_df <- function(extra_str="",
                        model_label="Base") {
  lapply(names(outcomes_main), function(nm) {
    df <- get_std_beta(
      outcomes_main[[nm]], m, extra_str)
    df$outcome       <- nm
    df$model         <- model_label
    df$outcome_label <- outcome_labels[nm]
    df
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      quartile_label = factor(
        paste0("Q", quartile),
        levels=c("Q4","Q3","Q2")),
      domain = dplyr::case_when(
        outcome %in% c("VO2","Speed",
                       "Power","Steps") ~
          "Physical performance",
        outcome %in% c("Digit","VF") ~
          "Cognitive & frailty",
        TRUE ~ DOMAIN_ENERGY),          # <- CHANGE 2
      outcome_label = factor(
        outcome_label,
        levels=rev(outcome_order))
    )
}

make_forest_b <- function(df, title_lbl,
                          sub_txt=NULL) {
  ggplot(df,
         aes(x=beta, y=outcome_label,
             color=domain,
             shape=quartile_label,
             fill=quartile_label)) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=4.5, ymax=8.5,
             fill="#E1F5EE", alpha=0.5) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=2.5, ymax=4.5,
             fill="#EEEDFE", alpha=0.5) +
    annotate("rect", xmin=-Inf, xmax=Inf,
             ymin=0.5, ymax=2.5,
             fill="#FAEEDA", alpha=0.5) +
    geom_vline(xintercept=0,
               color="grey30", linewidth=0.5,
               linetype="dashed") +
    annotate("text", x=0.02, y=Inf,
             label="Q1 (ref)", size=2.3,
             color="grey40", hjust=0,
             vjust=1.5, fontface="italic") +
    geom_errorbarh(
      aes(xmin=beta_lo, xmax=beta_hi,
          color=domain),
      height=0, linewidth=0.55,
      position=position_dodge(width=0.65)
    ) +
    geom_point(
      aes(color=domain), size=2.3,
      stroke=0.5,
      position=position_dodge(width=0.65)
    ) +
    geom_text(
      data=domain_label_df,
      aes(x=0.88, y=y, label=label,
          color=I(color)),
      hjust=1, size=2.2, lineheight=0.85,
      fontface="bold", inherit.aes=FALSE
    ) +
    scale_color_manual(
      values=domain_colors, guide="none") +
    scale_shape_manual(
      values=q_shapes_v5,
      breaks=c("Q2","Q3","Q4"),
      labels=c("Q2","Q3",
               "Q4 (highest entropy)"),
      name="Quartile vs. Q1 (reference)") +
    scale_fill_manual(
      values=q_colors_main,
      breaks=c("Q2","Q3","Q4"),
      labels=c("Q2","Q3",
               "Q4 (highest entropy)"),
      name="Quartile vs. Q1 (reference)") +
    guides(
      shape=guide_legend(
        nrow=1,
        override.aes=list(size=2.5)),
      fill=guide_legend(nrow=1)) +
    scale_x_continuous(
      breaks=seq(-0.8,0.8,0.2),
      limits=c(-0.9,0.9),
      labels=function(x) sprintf("%.1f",x)
    ) +
    scale_y_discrete() +                # <- FIX 
    labs(title=title_lbl, subtitle=sub_txt,
         x="Adjusted difference vs. Q1 (outcome SD units)",
         y=NULL) +
    coord_cartesian(clip="off") +
    theme_bw(base_size=9) +
    theme(
      plot.title         = element_text(
        face="bold", size=10, hjust=0),
      plot.subtitle      = element_text(
        size=7.5, color="grey40"),
      axis.text.y        = axis_y_elem,
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      legend.text        = element_text(
        size=7.5),
      legend.title       = element_text(
        size=8, face="bold"),
      legend.key.size    = unit(4,"mm"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(
        color="grey88", linewidth=0.3),
      plot.margin        = margin(
        3,5,3,0,"mm"))
}

save_fig <- function(p, fname, w=120, h=130) {
  leg <- cowplot::get_legend(
    p + theme(legend.position="bottom",
              legend.direction="horizontal"))
  combined <- cowplot::plot_grid(
    p + theme(legend.position="none"),
    leg, ncol=1, rel_heights=c(1,0.1))
  png(fname, width=w/25.4, height=h/25.4,
      units="in", res=300)
  print(combined)
  dev.off()
  cat("Saved:", fname, "\n")
}

std_v1 <- make_std_df("",          "Base")
std_v2 <- make_std_df("GDF15_std", "Base + GDF-15")
std_v3 <- make_std_df("CysC_std",  "Base + Cystatin-C")

fig_b_v1 <- make_forest_b(std_v1, "b",
                          "Adjusted: age, sex, site, height, weight")
fig_b_v2 <- make_forest_b(std_v2, "b",
                          "Adjusted: base + GDF-15 (continuous)")
fig_b_v3 <- make_forest_b(std_v3, "b",
                          "Adjusted: base + Cystatin-C (continuous)")

save_fig(fig_b_v1,
         file.path(PATH_RESULTS, "Fig2_panel_b_base.png"))
save_fig(fig_b_v2,
         file.path(PATH_RESULTS, "Fig2_panel_b_GDF15.png"))
save_fig(fig_b_v3,
         file.path(PATH_RESULTS, "Fig2_panel_b_CysC.png"))

leg_shared <- cowplot::get_legend(
  fig_b_v1 + theme(
    legend.position="bottom",
    legend.direction="horizontal"))
fig2_sensitivity <- cowplot::plot_grid(
  cowplot::plot_grid(
    fig_b_v1 + theme(
      legend.position="none",
      plot.subtitle=element_text(size=6.5)),
    fig_b_v2 + theme(
      legend.position="none",
      plot.subtitle=element_text(size=6.5)),
    fig_b_v3 + theme(
      legend.position="none",
      plot.subtitle=element_text(size=6.5)),
    ncol=3),
  leg_shared,
  ncol=1, rel_heights=c(1,0.08))
png(file.path(PATH_RESULTS, "Fig2_sensitivity_3panel.png"),
    width=330/25.4, height=130/25.4,
    units="in", res=300)
print(fig2_sensitivity)
dev.off()
cat("Saved: results_v14/Fig2_sensitivity_3panel.png\n")

# --- Supplementary Figure 1 (2-panel) ---
std_v1 <- make_std_df("",          "Base")
std_v2 <- make_std_df("GDF15_std", "Base + GDF-15")
std_v3 <- make_std_df("CysC_std",  "Base + Cystatin-C")

fig_b_v1 <- make_forest_b(std_v1, NULL,
                          "Adjusted: age, sex, site, height, weight")
fig_b_v2 <- make_forest_b(std_v2, "a",
                          "Adjusted: base + GDF-15 (continuous)")
fig_b_v3 <- make_forest_b(std_v3, "b",
                          "Adjusted: base + Cystatin-C (continuous)")

save_fig(fig_b_v2,
         file.path(PATH_RESULTS, "Fig2_panel_b_GDF15.png"))
save_fig(fig_b_v3,
         file.path(PATH_RESULTS, "Fig2_panel_b_CysC.png"))

leg_shared <- cowplot::get_legend(
  fig_b_v1 + theme(
    legend.position="bottom",
    legend.direction="horizontal"))

# GDF-15 / Cystatin-C adjusted panels
fig2_sensitivity <- cowplot::plot_grid(
  cowplot::plot_grid(
    fig_b_v2 + theme(
      legend.position="none",
      plot.subtitle=element_text(size=6.5)),
    fig_b_v3 + theme(
      legend.position="none",
      plot.subtitle=element_text(size=6.5)),
    ncol=2),
  leg_shared,
  ncol=1, rel_heights=c(1,0.08))
png(file.path(PATH_RESULTS, "FigS_sensitivity_2panel.png"),
    width=220/25.4, height=130/25.4,
    units="in", res=300)
print(fig2_sensitivity)
dev.off()
cat("Saved: results_v14/FigS_sensitivity_2panel.png\n")

# =============================================================================
# Supplementary Figure 2: Age-unweighted DISCO 
# =============================================================================


out_dir <- PATH_RESULTS

# ---- 0. Merge unweighted DISCO into the analysis data `m` -------------------
# `m` is the same analysis data frame used for Figure 2 (weighted BOX_DISCO_JR
# + all outcomes). We add the unweighted columns with a _UNW tag.
unw <- readr::read_csv(file.path(PATH_OUTPUT, "05z_entropy_measures_unweighted.csv"), show_col_types = FALSE) |>
  dplyr::distinct(ID, .keep_all = TRUE) |>
  dplyr::rename_with(~ paste0(.x, "_UNW"), .cols = -ID)
# -> BOX_DISCO_JR_UNW (continuous), BOX_DISCO_JR_Q_UNW (quartile 1-4)

stopifnot(exists("m"))
m_unw <- as.data.frame(m) |> dplyr::inner_join(unw, by = "ID")  # dedicated frame; do NOT overwrite shared m
cat("merged n =", nrow(m_unw), "\n")

Q_VAR <- "BOX_DISCO_JR_Q_UNW"   # <<< unweighted quartile drives the forest

# ---- 1. get_std_beta (quartile variable parameterized) ----------------------
get_std_beta <- function(outcome_var, data, q_var = Q_VAR, extra_str = "") {
  sd_out <- sd(data[[outcome_var]], na.rm = TRUE)
  data2  <- data |>
    dplyr::filter(!is.na(.data[[outcome_var]])) |>
    dplyr::mutate(
      Q_fac     = factor(.data[[q_var]], levels = 1:4),
      GDF15_std = (SOGDF15 - mean(SOGDF15, na.rm = TRUE)) / sd(SOGDF15, na.rm = TRUE),
      CysC_std  = (SOCYSC  - mean(SOCYSC,  na.rm = TRUE)) / sd(SOCYSC,  na.rm = TRUE)
    )
  base_cov <- "D1AGE2 + EEFEMALE + SITE + HWWGT + HWHGT"
  fml <- as.formula(paste0(
    outcome_var, " ~ Q_fac + ", base_cov,
    if (nchar(extra_str) > 0) paste0(" + ", extra_str) else ""))
  fit      <- lm(fml, data = data2)
  coef_tbl <- as.data.frame(summary(fit)$coefficients)
  coef_tbl$term <- rownames(coef_tbl)
  q_rows <- coef_tbl |> dplyr::filter(grepl("Q_fac[234]", term))
  data.frame(
    quartile = c(2, 3, 4),
    beta     = q_rows$Estimate / sd_out,
    beta_lo  = (q_rows$Estimate - 1.96 * q_rows$`Std. Error`) / sd_out,
    beta_hi  = (q_rows$Estimate + 1.96 * q_rows$`Std. Error`) / sd_out,
    p_value  = q_rows$`Pr(>|t|)`,
    sd_out   = sd_out
  )
}

# ---- 2. Outcomes / labels / ordering (identical to Figure 2) ----------------
outcomes_all <- list(
  VO2 = "TTPKVO2U", Speed = "NF400MPACE", Power = "LEPEAKPWR2", Steps = "ACSCFUAV",
  Digit = "DSCORR", VF = "FT0V2FN", OxPhos = "REMOXPHOS", CCR = "TTSSPKVO2",
  GDF15 = "SOGDF15", CysC = "SOCYSC"
)

# ── CHANGE 1: "VO2 peak" -> "VO2peak"  ─────────────────────────────

use_md <- requireNamespace("ggtext", quietly = TRUE)

lab_vo2  <- if (use_md) "VO<sub>2</sub>peak (mL/min)" else
  "VO\u2082peak (mL/min)"
lab_dsst <- if (use_md)
  "DSST score (0-133)<br>(Lower = poorer cognition)" else
    "DSST score (0-133)\n(Lower = poorer cognition)"
lab_vf   <- if (use_md)
  "Vigor-to-Frailty (0-12)<br>(Higher = more frail)" else
    "Vigor-to-Frailty (0-12)\n(Higher = more frail)"

# y element (element_markdown)
axis_y_elem <- if (use_md)
  ggtext::element_markdown(size = 8.5, lineheight = 0.95) else
    element_text(size = 8.5)

outcome_labels <- c(
  VO2 = lab_vo2, Speed = "Walking speed (m/s)",
  Power = "Leg peak power (W)", Steps = "Daily step count",
  Digit = lab_dsst,
  VF = lab_vf,
  OxPhos = "Max OxPhos (pmol/s/mg)", CCR = "Cost-capacity ratio (%)",
  GDF15 = "GDF-15 (log)", CysC = "Cystatin-C (log)"
)
outcome_order <- c(
  lab_vo2, "Walking speed (m/s)", "Leg peak power (W)",
  "Daily step count", lab_dsst,
  lab_vf,
  "Max OxPhos (pmol/s/mg)", "Cost-capacity ratio (%)"
)

# ── CHANGE 2: "Mitochondrial function" -> "Energetic capacity & efficiency" ──
DOMAIN_ENERGY <- "Energetic capacity & efficiency"

std_df2 <- lapply(names(outcomes_all), function(nm) {
  d <- get_std_beta(outcomes_all[[nm]], m_unw); d$outcome <- nm; d
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    quartile_label = factor(paste0("Q", quartile), levels = c("Q4", "Q3", "Q2")),
    domain = dplyr::case_when(
      outcome %in% c("VO2", "Speed", "Power", "Steps") ~ "Physical performance",
      outcome %in% c("Digit", "VF")                    ~ "Cognitive & frailty",
      outcome %in% c("OxPhos", "CCR")                  ~ DOMAIN_ENERGY,   # <- CHANGE 2
      TRUE                                             ~ "Biomarker"),
    outcome_label = outcome_labels[outcome],
    outcome_label = factor(outcome_label,
                           levels = c(rev(outcome_order), "GDF-15 (log)", "Cystatin-C (log)"))
  )

cat("=== Unweighted-DISCO Std beta, Q4 vs Q1 ===\n")
std_df2 |> dplyr::filter(quartile == 4) |>
  dplyr::select(outcome, beta, beta_lo, beta_hi, p_value) |>
  dplyr::mutate(dplyr::across(c(beta, beta_lo, beta_hi), \(x) round(x, 3)),
                p_value = round(p_value, 4)) |>
  as.data.frame() |> print()

# ---- 3. Aesthetics (identical to Figure 2) ----------------------------------
domain_colors <- c("Physical performance"            = "#0F6E56",
                   "Cognitive & frailty"             = "#534AB7",
                   "Energetic capacity & efficiency" = "#BA7517")  # <- CHANGE 2
q_colors_main <- c("Q2" = "#74C476", "Q3" = "#FD8D3C", "Q4" = "#D73027")
q_shapes_v5   <- c("Q2" = 21, "Q3" = 24, "Q4" = 23)
q_colors_scatter <- c("1" = "#4575B4", "2" = "#74C476", "3" = "#FD8D3C", "4" = "#D73027")

domain_label_df <- data.frame(
  y = c(1.5, 3.5, 6.5),
  label = c("Energetic capacity\n& efficiency",      # <- CHANGE 2
            "Cognitive &\nfrailty", "Physical\nperformance"),
  color = unname(domain_colors[c(DOMAIN_ENERGY, "Cognitive & frailty",
                                 "Physical performance")])
)

std_main <- std_df2 |>
  dplyr::filter(domain != "Biomarker") |>
  dplyr::mutate(quartile_label = factor(paste0("Q", quartile), levels = c("Q4", "Q3", "Q2")))

shared_shape_b <- scale_shape_manual(
  values = q_shapes_v5, breaks = c("Q2", "Q3", "Q4"),
  labels = c("Q2", "Q3", "Q4 (highest entropy)"),
  name = "Quartile vs. Q1 (reference)")
shared_fill_b <- scale_fill_manual(
  values = q_colors_main, breaks = c("Q2", "Q3", "Q4"),
  labels = c("Q2", "Q3", "Q4 (highest entropy)"),
  name = "Quartile vs. Q1 (reference)")

# ---- 4. Panel a: weighted vs unweighted DISCO -------------------------------
dA <- m_unw[stats::complete.cases(m_unw$BOX_DISCO_JR, m_unw$BOX_DISCO_JR_UNW), ]
r_p   <- cor(dA$BOX_DISCO_JR, dA$BOX_DISCO_JR_UNW)
rho_s <- cor(dA$BOX_DISCO_JR, dA$BOX_DISCO_JR_UNW, method = "spearman")

p_a <- ggplot(dA[!is.na(dA[[Q_VAR]]), ],
              aes(x = BOX_DISCO_JR, y = BOX_DISCO_JR_UNW,
                  color = factor(.data[[Q_VAR]]))) +
  geom_point(size = 0.85, alpha = 0.5) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x,
              color = "grey20", linewidth = 0.9, se = TRUE, fill = "grey70", alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.6,
           label = sprintf("Pearson~italic(r)==%.2f", r_p), parse = TRUE,
           size = 3, color = "grey20") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.09, vjust = 3.4,
           label = paste0("n = ", nrow(dA)), size = 2.9, color = "grey40") +
  scale_color_manual(values = q_colors_scatter,
                     labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                     name = "Unweighted DISCO\nquartile",
                     guide = guide_legend(nrow = 2, override.aes = list(size = 2.2))) +
  labs(title = "a", x = "Weighted DISCO", y = "Unweighted DISCO") +
  theme_bw(base_size = 9) +
  theme(aspect.ratio = 1,                       # <<< square plotting panel
        plot.title = element_text(face = "bold", size = 10, hjust = 0),
        legend.position = "bottom", legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8, face = "bold"),
        legend.key.size = unit(3.5, "mm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        plot.margin = margin(3, 5, 3, 3, "mm"))

# ---- 5. Panel b: forest (identical format to Figure 2b) ---------------------
p_b <- ggplot(std_main,
              aes(x = beta, y = outcome_label,
                  color = domain, shape = quartile_label, fill = quartile_label)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 4.5, ymax = 8.5, fill = "#E1F5EE", alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 4.5, fill = "#EEEDFE", alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 2.5, fill = "#FAEEDA", alpha = 0.5) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 0.02, y = Inf, label = "Q1 (ref)", size = 2.3,
           color = "grey40", hjust = 0, vjust = 1.5, fontface = "italic") +
  geom_errorbarh(aes(xmin = beta_lo, xmax = beta_hi, color = domain),
                 height = 0, linewidth = 0.55, position = position_dodge(width = 0.65)) +
  geom_point(aes(color = domain), size = 2.3, stroke = 0.5,
             position = position_dodge(width = 0.65)) +
  geom_text(data = domain_label_df, aes(x = 0.88, y = y, label = label, color = I(color)),
            hjust = 1, size = 2.2, lineheight = 0.85, fontface = "bold", inherit.aes = FALSE) +
  scale_color_manual(values = domain_colors, guide = "none") +
  shared_shape_b + shared_fill_b +
  guides(shape = guide_legend(nrow = 1, override.aes = list(size = 2.5)),
         fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), limits = c(-0.9, 0.9),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_y_discrete() +                      # <- CHANGE 3 (avoids the discrete/continuous scale error)
  labs(title = "b", x = "Adjusted difference vs. Q1 (outcome SD units)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
        axis.text.y = axis_y_elem,          # <- CHANGE 1 (markdown labels)
        legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8, face = "bold"),
        legend.key.size = unit(4, "mm"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3),
        plot.margin = margin(3, 5, 3, 0, "mm"))

# ---- 6. Combine & save ------------------------------------------------------
legend_a <- cowplot::get_legend(p_a + theme(legend.position = "bottom"))
legend_b <- cowplot::get_legend(p_b + theme(legend.position = "bottom",
                                            legend.direction = "horizontal"))
left_col <- cowplot::plot_grid(
  p_a + theme(legend.position = "none"), legend_a, ncol = 1, rel_heights = c(1, 0.18))
right_col <- cowplot::plot_grid(
  p_b + theme(legend.position = "none"), legend_b, ncol = 1, rel_heights = c(1, 0.08))
fig_body <- cowplot::plot_grid(left_col, right_col, ncol = 2, rel_widths = c(0.34, 0.66))

supp_fig <- cowplot::plot_grid(
  fig_body,
  cowplot::ggdraw() + cowplot::draw_label(
    paste0("Panel a: weighted vs. unweighted DISCO (points colored by unweighted quartile). ",
           "Panel b: linear regression adjusted for age, sex, clinic site, body weight, and height. ",
           "Error bars = 95% CI. \u03b2 = mean difference vs. Q1 in SD units of each outcome. ",
           "Quartiles derived from unweighted DISCO."),
    x = 0.01, y = 0.5, hjust = 0, vjust = 0.5, size = 6.5, color = "grey40"),
  ncol = 1, rel_heights = c(1, 0.07))

png(file.path(out_dir, "SuppFig_unweighted_DISCO.png"),
    width = 220/25.4, height = 150/25.4, units = "in", res = 300)
print(supp_fig); dev.off()

ggsave(file.path(out_dir, "SuppFig_unweighted_DISCO.pdf"), supp_fig,
       width = 220/25.4, height = 150/25.4, device = cairo_pdf)

readr::write_csv(std_df2, file.path(out_dir, "SuppFig_unweighted_DISCO_stdbeta.csv"))
cat("Saved: SuppFig_unweighted_DISCO.png / .pdf  +  _stdbeta.csv\n")
