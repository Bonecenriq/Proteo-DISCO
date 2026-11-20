# ==============================================================================
# STUDY       : SOMMA
# PROJECT     : SOMMA SomaLogic
# PROGRAM     : 09_survival_analysis.R
#
# PURPOSE     : Relate proteomic entropy to non-elective hospitalization, and
#               produce Figure 5 and Supplementary Figures 5 and 6.
#
# INPUTS      : - m (04_analysis_dataset.R)
#               - analysis_df (07_organ_entropy.R)
#               - Figure 2 helper objects (06_phenotype_associations.R)
#               - 09b_entropy_measures_organ_box.csv
#               - 05_entropy_measures.csv, 05z_entropy_measures_unweighted.csv
#               - 02_box_transformed_data_with_pheno.csv (medication count)
#
# OUTPUTS     : - Fig_hosp_KM_heatmap_BOX.png / .pdf     (Figure 5)
#               - Table_hosp_HR_BOX.csv
#               - Supp_Fig_hosp_HR_sensitivity_BOX.png / .pdf (Supplementary Figure 5)
#               - FigS_pheno_forest_meds.png
#               - FigS_hosp_HR_meds.png / .pdf
#               - FigS_meds_combined.png / .pdf          (Supplementary Figure 6)
#               - Table_hosp_HR_meds.csv
#
# DATA ACCESS : Individual-level SOMMA data are not included in this
#               repository. Access to the underlying study data is required
#               to reproduce the analyses.
#
# DEPENDENCIES: - 99_figure_setup.R
#               - 04_analysis_dataset.R, 06_phenotype_associations.R,
#                 07_organ_entropy.R
#               - survival, cowplot, scales, sandwich, lmtest
#
# VERSION     : 1.0
# DATE        : 2026-08-15
# AUTHOR      : Namki Hong
#
# REVISION HISTORY
# 08/15/26    Creation (split out of the combined analysis script)
#
# NOTES       : Cox proportional hazards models use the Efron method for ties with
#               robust standard errors. The base model is adjusted for age, sex,
#               clinic site, body weight and height; sensitivity models add GDF-15,
#               Cystatin C, or log medication count.
#
#               The final block (Supplementary Figure 6) reuses helper functions and
#               aesthetic objects defined in 06_phenotype_associations.R and earlier
#               in this script, and therefore requires those scripts to have been
#               run in the same R session.
#
#               The follow-up duration block reports median follow-up using the
#               reverse Kaplan-Meier estimator.
# ==============================================================================
source("code/99_figure_setup.R")

stopifnot(exists("m"), exists("analysis_df"))


# =============================================================================
# Figure 5 KM curve Forest plot & Supplementary Figure 5
# =============================================================================


dir.create(PATH_RESULTS, showWarnings=FALSE)

# ---- 0. Shared settings ----------------------------------------------------
organs <- c("blood","bone","brain","heart",
            "kidney","liver","lung",
            "lymphoid","skeletal")

base_cov <- "D1AGE2 + EEFEMALE + SITE_f + HWWGT + HWHGT"

q_colors_km <- c(
  "Q1 (lowest)"  = "#4575B4",  # blue
  "Q2"           = "#74C476",  # green
  "Q3"           = "#FD8D3C",  # orange
  "Q4 (highest)" = "#D73027"   # red
)

q_colors_km_dark <- c(
  "Q1 (lowest)"  = "#2C5AA0",  # deeper blue
  "Q2"           = "#41A85F",  # deeper green
  "Q3"           = "#E8730C",  # deeper orange
  "Q4 (highest)" = "#B81D13"   # deeper red
)


# =============================================================================
# 1. BOX file load analysis
# =============================================================================
organ_box_raw <- readr::read_csv(
  file.path(PATH_OUTPUT, "09b_entropy_measures_organ_box.csv"),
  show_col_types=FALSE
)

q_vars_box <- paste0("BOX_DISCO_JR_Q_", organs)

analysis_df <- meta_df |>
  dplyr::left_join(
    organ_box_raw |>
      dplyr::select(ID, all_of(q_vars_box)),
    by="ID"
  ) |>
  dplyr::mutate(
    SITE_f    = factor(SITE),
    GDF15_std = (SOGDF15 -
                   mean(SOGDF15, na.rm=TRUE)) /
      sd(SOGDF15, na.rm=TRUE),
    CysC_std  = (SOCYSC -
                   mean(SOCYSC, na.rm=TRUE)) /
      sd(SOCYSC, na.rm=TRUE)
  )

# quartile factor conversion
for (qv in q_vars_box) {
  analysis_df[[qv]] <- factor(
    analysis_df[[qv]], levels=1:4)
}

cat("Analysis n =", nrow(analysis_df), "\n")
cat("Events:", sum(analysis_df$HAEMERG,
                   na.rm=TRUE), "\n")


# =============================================================================
# 2. Helper function
# =============================================================================

# 2a. Full proteome HR  ( coxph fit)
extract_hr_full <- function(fit, model_name) {
  ci    <- summary(fit)$conf.int
  co    <- summary(fit)$coefficients
  q_idx <- grepl("^Q_fac[234]$", rownames(ci))
  data.frame(
    organ    = "full",
    model    = model_name,
    quartile = c("Q2","Q3","Q4"),
    HR       = ci[q_idx, "exp(coef)"],
    lo       = ci[q_idx, "lower .95"],
    hi       = ci[q_idx, "upper .95"],
    p        = co[q_idx, "Pr(>|z|)"]
  )
}

# 2b. Organ-specific HR (BOX quartile, efron+robust)
get_hr_all_q_box <- function(q_var, data,
                             extra="",
                             unadj=FALSE) {
  d <- data |>
    dplyr::filter(!is.na(HAFUTBLU),
                  !is.na(HAEMERG),
                  !is.na(.data[[q_var]]))
  if (nrow(d) < 50) return(NULL)
  
  fml <- if (unadj) {
    as.formula(paste0(
      "Surv(HAFUTBLU, HAEMERG) ~ ", q_var))
  } else {
    as.formula(paste0(
      "Surv(HAFUTBLU, HAEMERG) ~ ", q_var,
      " + ", base_cov,
      if (nchar(extra)>0)
        paste0(" + ", extra) else ""))
  }
  
  fit <- tryCatch(
    coxph(fml, data=d, ties="efron",
          robust=TRUE, id=ID),
    error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  ci  <- summary(fit)$conf.int
  co  <- summary(fit)$coefficients
  idx <- grepl(paste0(q_var,"[234]"),
               rownames(ci))
  if (sum(idx)==0) return(NULL)
  
  qs <- sub(paste0(".*",q_var), "",
            rownames(ci)[idx])
  data.frame(
    quartile = paste0("Q", qs),
    HR       = ci[idx, "exp(coef)"],
    lo       = ci[idx, "lower .95"],
    hi       = ci[idx, "upper .95"],
    p        = co[idx, "Pr(>|z|)"],
    n        = nrow(d),
    ev       = sum(d$HAEMERG, na.rm=TRUE)
  )
}

# 2c. HR data filtering function
tidy_hr <- function(df, organ_order_ref) {
  n_org <- length(organ_order_ref)
  df |>
    dplyr::mutate(
      organ_label = dplyr::recode(organ,
                                  "full"     = "Full Proteome",
                                  "blood"    = "Blood Vessel",
                                  "bone"     = "Bone",
                                  "brain"    = "Brain",
                                  "heart"    = "Heart",
                                  "kidney"   = "Kidney",
                                  "liver"    = "Liver",
                                  "lung"     = "Lung",
                                  "lymphoid" = "Lymphoid",
                                  "skeletal" = "Skeletal Muscle"
      ),
      sig       = p < 0.05,
      sig_label = dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ ""
      ),
      hr_main   = sprintf("%.2f%s",
                          HR, sig_label),
      hr_ci     = sprintf("(%.2f\u2013%.2f)",
                          lo, hi),
      fill_val  = dplyr::if_else(
        HR >= 1, HR, NA_real_),
      white_text = HR > 1.7
    ) |>
    dplyr::mutate(
      organ_label = factor(
        organ_label,
        levels=organ_order_ref),
      quartile    = factor(
        quartile, levels=c("Q2","Q3","Q4")),
      n_org       = n_org,
      y_screen    = n_org + 1 -
        as.numeric(organ_label),
      y_main      = y_screen + 0.14,
      y_ci        = y_screen - 0.22
    )
}

# =============================================================================
# 3. Organ-specific HR calculation (Unadjusted + Base)
# =============================================================================
cat("Computing organ-specific HR (BOX)...\n")

hr_organ_list <- list()

for (organ_nm in organs) {
  q_var <- paste0("BOX_DISCO_JR_Q_", organ_nm)
  for (model_nm in c("Unadjusted","Base")) {
    res <- get_hr_all_q_box(
      q_var, analysis_df,
      unadj=(model_nm=="Unadjusted"))
    if (!is.null(res)) {
      hr_organ_list[[
        length(hr_organ_list)+1]] <-
        data.frame(organ=organ_nm,
                   model=model_nm, res)
    }
  }
}

hr_organ_df <- dplyr::bind_rows(hr_organ_list)
cat("Organ HR rows:", nrow(hr_organ_df), "\n")

analysis_df$Q_fac <- factor(analysis_df$BOX_DISCO_JR_Q, levels = 1:4)

fit_unadj <- coxph(
  Surv(HAFUTBLU, HAEMERG) ~ Q_fac,
  data = analysis_df, ties = "efron", robust = TRUE, id = ID)

fit_base <- coxph(
  as.formula(paste0("Surv(HAFUTBLU, HAEMERG) ~ Q_fac + ", base_cov)),
  data = analysis_df, ties = "efron", robust = TRUE, id = ID)

fit_gdf <- coxph(
  as.formula(paste0("Surv(HAFUTBLU, HAEMERG) ~ Q_fac + ", base_cov, " + GDF15_std")),
  data = analysis_df, ties = "efron", robust = TRUE, id = ID)

fit_cysc <- coxph(
  as.formula(paste0("Surv(HAFUTBLU, HAEMERG) ~ Q_fac + ", base_cov, " + CysC_std")),
  data = analysis_df, ties = "efron", robust = TRUE, id = ID)

cat("Full-proteome Cox fits ready (unadj / base / gdf / cysc)\n")

# =============================================================================
# 4. Full proteome + Organ integration, sorting
# =============================================================================
hr_full <- dplyr::bind_rows(
  extract_hr_full(fit_unadj, "Unadjusted"),
  extract_hr_full(fit_base,  "Base")
)

hr_all <- dplyr::bind_rows(hr_full,
                           hr_organ_df)

# Organ order: Q4 Unadj HR descending
organ_order_heat <- hr_all |>
  dplyr::mutate(
    organ_label = dplyr::recode(organ,
                                "full"     = "Full Proteome",
                                "blood"    = "Blood Vessel",
                                "bone"     = "Bone",
                                "brain"    = "Brain",
                                "heart"    = "Heart",
                                "kidney"   = "Kidney",
                                "liver"    = "Liver",
                                "lung"     = "Lung",
                                "lymphoid" = "Lymphoid",
                                "skeletal" = "Skeletal Muscle")
  ) |>
  dplyr::filter(quartile=="Q4",
                model=="Unadjusted") |>
  dplyr::arrange(dplyr::desc(HR)) |>
  dplyr::pull(organ_label)

organ_order_heat <- c(
  "Full Proteome",
  organ_order_heat[
    organ_order_heat != "Full Proteome"]
)

cat("\nOrgan order (Q4 Unadj HR desc):\n")
print(organ_order_heat)

# data curation
heat_df <- hr_all |>
  tidy_hr(organ_order_heat) |>
  dplyr::mutate(
    model_label = factor(
      dplyr::if_else(model=="Unadjusted",
                     "Unadjusted","Adjusted"),
      levels=c("Unadjusted","Adjusted"))
  )

# yaxis face
y_faces_heat <- ifelse(
  rev(organ_order_heat)=="Full Proteome",
  "bold","plain")

# CSV save
heat_df |>
  dplyr::select(organ_label, model,
                quartile, HR, lo, hi,
                sig_label, sig) |>
  dplyr::arrange(model, organ_label,
                 quartile) |>
  readr::write_csv(
    file.path(PATH_RESULTS, "Table_hosp_HR_BOX.csv"))
cat("Saved: Table_hosp_HR_BOX.csv\n")

# =============================================================================
# 5. Heatmap common setting
# =============================================================================
heat_fill_scale <- ggplot2::scale_fill_gradientn(
  colours = c("#FEF0EE","#FBBCB4",
              "#F47E71","#D94E3A",
              "#A61C0E"),
  values  = scales::rescale(
    c(1.0,1.3,1.7,2.1,2.6)),
  name    = "Hazard ratio\n(vs Q1)",
  limits  = c(1.0,2.6),
  breaks  = c(1.0,1.5,2.0,2.5),
  na.value= "grey80"
)

heat_color_scale <- ggplot2::scale_color_manual(
  values=c("FALSE"="grey20","TRUE"="white"),
  guide ="none"
)

heat_theme <- ggplot2::theme_bw(base_size=9) +
  ggplot2::theme(
    axis.text.x      = ggplot2::element_text(
      size=8.5, color="black"),
    axis.title.x     = ggplot2::element_text(
      size=9),
    panel.grid       = ggplot2::element_blank(),
    panel.spacing    = ggplot2::unit(2,"mm"),
    plot.background  = ggplot2::element_rect(
      fill="white", color=NA),
    plot.margin      = ggplot2::margin(
      5,5,3,2,"mm")
  )

# common geom function
add_heat_geoms <- function(p) {
  p +
    ggplot2::geom_tile(
      ggplot2::aes(fill=fill_val),
      color="white", linewidth=0.6,
      width=0.95, height=0.95) +
    ggplot2::geom_tile(
      data=\(d) dplyr::filter(d, HR < 1),
      fill="grey80", color="white",
      linewidth=0.6, width=0.95,
      height=0.95) +
    # HR (main)
    ggplot2::geom_text(
      data=\(d) dplyr::filter(d, HR >= 1),
      ggplot2::aes(label=hr_main,
                   color=white_text,
                   y=y_main),
      size=3.0, fontface="bold",
      vjust=0.5) +
    # CI
    ggplot2::geom_text(
      data=\(d) dplyr::filter(d, HR >= 1),
      ggplot2::aes(label=hr_ci,
                   color=white_text,
                   y=y_ci),
      size=1.9, fontface="plain",
      vjust=0.5, alpha=0.85) +
    # HR < 1 (main)
    ggplot2::geom_text(
      data=\(d) dplyr::filter(d, HR < 1),
      ggplot2::aes(label=hr_main, y=y_main),
      color="grey40", size=3.0,
      fontface="bold", vjust=0.5) +
    # HR<1 CI
    ggplot2::geom_text(
      data=\(d) dplyr::filter(d, HR < 1),
      ggplot2::aes(label=hr_ci, y=y_ci),
      color="grey55", size=1.9,
      fontface="plain", vjust=0.5)
}

# =============================================================================
# 6. KM Curve 
# =============================================================================
surv_df2 <- analysis_df |>
  dplyr::filter(!is.na(HAFUTBLU), !is.na(HAEMERG), !is.na(BOX_DISCO_JR_Q)) |>
  dplyr::mutate(
    BOX_DISCO_JR_Q = factor(BOX_DISCO_JR_Q, levels = 1:4),
    Q_label = factor(
      dplyr::recode(as.character(BOX_DISCO_JR_Q),
                    "1" = "Q1 (lowest)", "2" = "Q2", "3" = "Q3", "4" = "Q4 (highest)"),
      levels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)")))


km_fit <- survfit(
  Surv(HAFUTBLU/365.25, HAEMERG) ~
    BOX_DISCO_JR_Q,
  data=surv_df2)

km_df <- data.frame(
  time   = km_fit$time,
  surv   = km_fit$surv,
  upper  = km_fit$upper,
  lower  = km_fit$lower,
  strata = rep(names(km_fit$strata),
               km_fit$strata)
) |>
  dplyr::mutate(
    Q = dplyr::case_when(
      grepl("=1",strata) ~ "Q1 (lowest)",
      grepl("=2",strata) ~ "Q2",
      grepl("=3",strata) ~ "Q3",
      grepl("=4",strata) ~ "Q4 (highest)"
    ),
    Q = factor(Q, levels=c(
      "Q1 (lowest)","Q2",
      "Q3","Q4 (highest)"))
  )

km_df0 <- data.frame(
  time=0, surv=1, upper=1, lower=1,
  strata=NA,
  Q=factor(
    c("Q1 (lowest)","Q2",
      "Q3","Q4 (highest)"),
    levels=c("Q1 (lowest)","Q2",
             "Q3","Q4 (highest)")))

km_plot_df <- dplyr::bind_rows(km_df0, km_df) |>
  dplyr::arrange(Q, time)

lr_test  <- survdiff(
  Surv(HAFUTBLU, HAEMERG) ~
    BOX_DISCO_JR_Q,
  data=surv_df2)
lr_p     <- 1 - pchisq(
  lr_test$chisq,
  df=length(lr_test$n)-1)
lr_label <- if (lr_p<0.001) "Log-rank p < 0.001" else
  sprintf("Log-rank p = %.3f", lr_p)

# Number at risk
risk_df <- lapply(0:4, function(t) {
  surv_df2 |>
    dplyr::group_by(Q_label) |>
    dplyr::summarise(
      n_risk=sum(HAFUTBLU/365.25 >= t,
                 na.rm=TRUE),
      .groups="drop") |>
    dplyr::mutate(time=t)
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    Q_label=factor(Q_label,
                   levels=rev(c(
                     "Q1 (lowest)","Q2",
                     "Q3","Q4 (highest)")))
  )

# KM plot
p_km <- ggplot2::ggplot(
  km_plot_df,
  ggplot2::aes(x=time, y=1-surv,
               color=Q, fill=Q)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin=1-upper,
                 ymax=1-lower),
    alpha=0.15, color=NA) +
  ggplot2::geom_step(linewidth=0.85) +
  ggplot2::annotate("text",
                    x=0.08, y=0.43,
                    label=lr_label,
                    size=2.6, hjust=0, color="grey30",
                    fontface="italic") +
  ggplot2::scale_color_manual(
    values=q_colors_km,
    name="Proteomic entropy\n(DISCO) quartile") +
  ggplot2::scale_fill_manual(
    values=q_colors_km,
    name="Proteomic entropy\n(DISCO) quartile") +
  ggplot2::scale_x_continuous(
    name="Follow-up (years)",
    breaks=0:4, limits=c(0,4.2)) +
  ggplot2::scale_y_continuous(
    name=paste0(
      "Cumulative incidence of\n",
      "non-elective hospitalization"),
    labels=scales::percent_format(accuracy=1),
    limits=c(0,0.48),
    breaks=seq(0,0.4,0.1)) +
  ggplot2::labs(subtitle="Unadjusted") +
  ggplot2::theme_bw(base_size=9) +
  ggplot2::theme(
    plot.subtitle    = ggplot2::element_text(
      size=7.5, color="grey40"),
    legend.position  = "none",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      color="grey92", linewidth=0.3),
    panel.background = ggplot2::element_rect(
      fill="white", color=NA),
    plot.background  = ggplot2::element_rect(
      fill="white", color=NA),
    plot.margin      = ggplot2::margin(
      5,5,3,3,"mm"))

# Number at risk plot
p_risk <- ggplot2::ggplot(
  risk_df,
  ggplot2::aes(x=time, y=Q_label,
               label=n_risk, color=Q_label)
) +
  ggplot2::geom_text(size=2.5,
                     fontface="bold") +
  ggplot2::scale_color_manual(
    values=rev(q_colors_km_dark),
    guide="none") +
  ggplot2::scale_x_continuous(
    breaks=0:4, limits=c(0,4.2)) +
  ggplot2::labs(x=NULL, y=NULL,
                title="Number at risk") +
  ggplot2::theme_void(base_size=7.5) +
  ggplot2::theme(
    plot.title   = ggplot2::element_text(
      size=6.5, hjust=0, color="grey40"),
    axis.text.y  = ggplot2::element_text(
      size=6.5, hjust=1,
      color=rev(unname(q_colors_km_dark)),
      margin=ggplot2::margin(r=3)),
    plot.background = ggplot2::element_rect(
      fill="white", color=NA),
    plot.margin  = ggplot2::margin(
      0,5,3,5,"mm"))

# KM legend
leg_km <- cowplot::get_legend(
  p_km +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.text      = ggplot2::element_text(
        size=7),
      legend.title     = ggplot2::element_text(
        size=7.5, face="bold"),
      legend.key.size  = ggplot2::unit(3.5,"mm")
    ) +
    ggplot2::scale_color_manual(
      values=q_colors_km,
      name="Proteomic entropy\n(DISCO) quartile") +
    ggplot2::scale_fill_manual(
      values=q_colors_km,
      name="Proteomic entropy\n(DISCO) quartile")
)

km_col <- cowplot::plot_grid(
  p_km, p_risk,
  ncol=1, rel_heights=c(4,1))

km_full <- cowplot::plot_grid(
  km_col, leg_km,
  ncol=1, rel_heights=c(1,0.13)) +
  ggplot2::theme(
    plot.background=ggplot2::element_rect(
      fill="white", color=NA))

# =============================================================================
# 7. Main Figure: KM + HR Heatmap (Unadjusted / Adjusted)
# =============================================================================
p_heat_main <- add_heat_geoms(
  ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x=quartile,
                 y=organ_label))
) +
  heat_fill_scale +
  heat_color_scale +
  ggplot2::facet_grid(model_label ~ .,
                      switch="y") +
  ggplot2::scale_x_discrete(
    name  ="Quartile vs Q1",
    labels=c("Q2","Q3","Q4 (highest)")
  ) +
  ggplot2::scale_y_discrete(
    name  =NULL,
    limits=rev(organ_order_heat)
  ) +
  ggplot2::labs(
    caption=paste0(
      "HR (95% CI). Bold = HR value. ",
      "* p<0.05, ** p<0.01, *** p<0.001. ",
      "Adjusted: age, sex, site, body weight, and height",
      "Grey = HR < 1.00 ",
      "Sorted by Q4 unadjusted HR.")
  ) +
  heat_theme +
  ggplot2::theme(
    strip.text.y     = ggplot2::element_text(
      size=9, face="bold", angle=0,
      margin=ggplot2::margin(l=3,r=3)),
    strip.background = ggplot2::element_rect(
      fill="grey92", color="grey70",
      linewidth=0.3),
    strip.placement  = "outside",
    axis.text.y      = ggplot2::element_text(
      size=8.5, color="black",
      face=y_faces_heat),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = ggplot2::element_text(
      size=7.5),
    legend.title     = ggplot2::element_text(
      size=8),
    legend.key.size  = ggplot2::unit(0.4,"cm"),
    legend.key.width = ggplot2::unit(1.8,"cm"),
    plot.caption     = ggplot2::element_text(
      size=6, color="grey40",
      margin=ggplot2::margin(t=4))
  )

fig_main <- cowplot::plot_grid(
  km_full, p_heat_main,
  ncol=2,
  rel_widths=c(1,1.6),
  labels=c("a","b"),
  label_size=11,
  label_fontface="bold"
) +
  ggplot2::theme(
    plot.background=ggplot2::element_rect(
      fill="white", color=NA))

ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_hosp_KM_heatmap_BOX.pdf"),
  plot=fig_main, width=260, height=170,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Fig_hosp_KM_heatmap_BOX.png"),
  plot=fig_main, width=260, height=170,
  units="mm", dpi=300)
cat("Saved: Fig_hosp_KM_heatmap_BOX\n")

# =============================================================================
# 8. Supplementary Figure 5. Sensitivity HR (Base+GDF15, Base+CysC)
# =============================================================================
cat("Computing sensitivity HR (BOX)...\n")

hr_organ_sens_list <- list()

for (organ_nm in organs) {
  q_var <- paste0("BOX_DISCO_JR_Q_", organ_nm)
  for (model_nm in c("Base + GDF-15",
                     "Base + CysC")) {
    extra <- dplyr::if_else(
      model_nm=="Base + GDF-15",
      "GDF15_std","CysC_std")
    res <- get_hr_all_q_box(
      q_var, analysis_df,
      extra=extra, unadj=FALSE)
    if (!is.null(res)) {
      hr_organ_sens_list[[
        length(hr_organ_sens_list)+1]] <-
        data.frame(organ=organ_nm,
                   model=model_nm, res)
    }
  }
}

hr_organ_sens_df <- dplyr::bind_rows(
  hr_organ_sens_list)

hr_full_sens <- dplyr::bind_rows(
  extract_hr_full(fit_gdf,  "Base + GDF-15"),
  extract_hr_full(fit_cysc, "Base + CysC")
)

heat_df_sens <- dplyr::bind_rows(
  hr_full_sens,
  hr_organ_sens_df
) |>
  tidy_hr(organ_order_heat) |>
  dplyr::mutate(
    model = factor(model,
                   levels=c("Base + GDF-15",
                            "Base + CysC"))
  )

cat("Sensitivity HR rows:",
    nrow(heat_df_sens), "\n")

# =============================================================================
# 9. Supp Figure: Sensitivity Heatmap (2panels)
# =============================================================================
make_heat_sens <- function(df, model_nm,
                           show_y=TRUE) {
  d <- df |> dplyr::filter(model==model_nm)
  
  p <- add_heat_geoms(
    ggplot2::ggplot(
      d,
      ggplot2::aes(x=quartile,
                   y=organ_label))
  ) +
    heat_fill_scale +
    heat_color_scale +
    ggplot2::scale_x_discrete(
      name  ="Quartile vs Q1",
      labels=c("Q2","Q3","Q4 (highest)")
    ) +
    ggplot2::scale_y_discrete(
      name  =NULL,
      limits=rev(organ_order_heat),
      labels=if(show_y)
        ggplot2::waiver()
      else
        setNames(
          rep("",length(organ_order_heat)),
          organ_order_heat)
    ) +
    ggplot2::labs(
      title = dplyr::case_when(
        model_nm=="Base + GDF-15" ~
          "a   Base + GDF-15 adjusted",
        model_nm=="Base + CysC" ~
          "b   Base + Cystatin-C adjusted"
      )
    ) +
    heat_theme +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(
        face="bold", size=9, hjust=0),
      axis.text.y     = ggplot2::element_text(
        size=8.5, color="black",
        face=if(show_y) y_faces_heat else
          rep("plain",
              length(organ_order_heat))),
      legend.position = "none"
    )
  p
}

p_gdf  <- make_heat_sens(
  heat_df_sens, "Base + GDF-15",
  show_y=TRUE)
p_cysc <- make_heat_sens(
  heat_df_sens, "Base + CysC",
  show_y=FALSE)

# Shared legend
leg_sens <- cowplot::get_legend(
  make_heat_sens(
    heat_df_sens, "Base + GDF-15",
    show_y=TRUE) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.text      = ggplot2::element_text(
        size=7.5),
      legend.title     = ggplot2::element_text(
        size=8),
      legend.key.size  = ggplot2::unit(0.4,"cm"),
      legend.key.width = ggplot2::unit(1.8,"cm"))
)

panels_sens <- cowplot::plot_grid(
  p_gdf, p_cysc,
  ncol=2, align="h", axis="tb",
  rel_widths=c(1.3,1))

fig_sens <- cowplot::plot_grid(
  panels_sens, leg_sens,
  ncol=1, rel_heights=c(1,0.07))

fig_sens <- cowplot::plot_grid(
  fig_sens,
  cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(
        "HR (95% CI). Bold = HR value. ",
        "* p<0.05, ** p<0.01, *** p<0.001. ",
        "Base model: adjusted for age, sex, site, body weight, and height",
        "GDF-15 and Cystatin-C: per SD ",
        "Grey = HR < 1.00 ",
        "Sorted by Q4 unadjusted HR."
      ),
      x=0.01, y=0.5, hjust=0, vjust=0.5,
      size=6.5, color="grey40"),
  ncol=1, rel_heights=c(1,0.05)
) +
  ggplot2::theme(
    plot.background=ggplot2::element_rect(
      fill="white", color=NA))

ggplot2::ggsave(
  file.path(PATH_RESULTS, "Supp_Fig_hosp_HR_sensitivity_BOX.pdf"),
  plot=fig_sens, width=260, height=170,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  file.path(PATH_RESULTS, "Supp_Fig_hosp_HR_sensitivity_BOX.png"),
  plot=fig_sens, width=260, height=170,
  units="mm", dpi=300)
cat("Saved: Supp_Fig_hosp_HR_sensitivity_BOX\n")

cat("\n=== Done ===\n")
cat("Output files:\n")
cat("  Fig_hosp_KM_heatmap_BOX.png\n")
cat("  Supp_Fig_hosp_HR_sensitivity_BOX.png\n")
cat("  Table_hosp_HR_BOX.csv\n")

# =============================================================
# Follow-up duration — Reverse KM
# =============================================================


surv_df <- m |>
  dplyr::filter(
    !is.na(HAFUTBLU),
    !is.na(HAEMERG)
  )

cat("n =", nrow(surv_df), "\n")
cat("Events:", sum(surv_df$HAEMERG), "\n")

# Reverse KM
fit_rev <- survival::survfit(
  survival::Surv(HAFUTBLU,
                 1 - HAEMERG) ~ 1,
  data=surv_df
)

fu <- summary(fit_rev)$table
cat("\n=== Median follow-up (Reverse KM) ===\n")
cat("Median:", round(fu["median"]/365.25, 2),
    "years\n")
cat("IQR:   ",
    round(fu["0.25"]/365.25, 2), "~",
    round(fu["0.75"]/365.25, 2),
    "years\n")
cat("Range: ",
    round(min(surv_df$HAFUTBLU)/365.25, 2),
    "~",
    round(max(surv_df$HAFUTBLU)/365.25, 2),
    "years\n")


# =============================================================
# DISCO (weighted) and DISCO (unweighted) comparison 
# =============================================================

# --- unweighted:  _UNW  ---
unw <- read_csv(file.path(PATH_OUTPUT, "05z_entropy_measures_unweighted.csv"), show_col_types = FALSE) |>
  distinct(ID, .keep_all = TRUE) |>
  rename_with(~ paste0(.x, "_UNW"), .cols = -ID)
# ID, BOX_DM_JRPCA_UNW, BOX_DM_JRSHRINK_UNW, BOX_DISCO_JR_UNW, ...

# --- weighted DISCO  ---
w <- read_csv(file.path(PATH_OUTPUT, "05_entropy_measures.csv"), show_col_types = FALSE) |>
  select(ID, BOX_DISCO_JR) |>          # weighted
  distinct(ID, .keep_all = TRUE)

# --- Join ---
df <- inner_join(w, unw, by = "ID")
cat("merged n =", nrow(df), "\n")

# --- correlation (weighted vs unweighted DISCO) ---
x <- df$BOX_DISCO_JR          # weighted
y <- df$BOX_DISCO_JR_UNW      # unweighted

cat("complete pairs n =", sum(complete.cases(x, y)), "\n")
print(cor.test(x, y, method = "pearson"))
print(cor.test(x, y, method = "spearman"))

# --- scatter plot ---
plot(x, y, pch = 16, col = "#4F2D2D",
     xlab = "BOX_DISCO_JR (weighted)",
     ylab = "BOX_DISCO_JR (unweighted)",
     main = sprintf("Pearson r = %.3f", cor(x, y, use = "complete.obs")))
abline(lm(y ~ x), col = "red", lwd = 2)


# =============================================================================
# Supplementary Figure 6
# SENSITIVITY: adding M1ANMEDS_LOG (log medication count) to the BASE model - Supplementary Figure 6
#   Figure a  — phenotype-association forest, SAME format as Figure 2 (panel b)
#               weighted-DISCO quartiles, base + M1ANMEDS_LOG
#   Figure b  — hospitalization Cox HR heatmap, SAME format as the Figure 5
#               sensitivity heatmap, single panel: base + M1ANMEDS_LOG
# -----------------------------------------------------------------------------
# PREREQUISITES (run 06_phenotype_associations.R and the Figure 5 blocks
#   From 06_phenotype_associations.R : get_std_beta(), make_forest_b(),
#                save_fig(), outcome_labels, outcome_order, domain_colors,
#                q_colors_main, q_shapes_v5, domain_label_df, `m`.
#   From this script (Figure 5 blocks) : get_hr_all_q_box(),
#                extract_hr_full(), tidy_hr(), add_heat_geoms(), heat_fill_scale,
#                heat_color_scale, heat_theme, organ_order_heat, y_faces_heat,
#                organs, base_cov, and `analysis_df` (07_organ_entropy.R).
#
#   [revision]
#     1) "VO2 peak" -> "VO2peak"
#     2) "Mitochondrial function" -> "Energetic capacity & efficiency"
# =============================================================================

dir.create(PATH_RESULTS, showWarnings = FALSE, recursive = TRUE)

# ---- 0. Load M1ANMEDS_LOG and merge into both analysis frames ---------------
meds <- readr::read_csv(file.path(PATH_DATA, "02_box_transformed_data_with_pheno.csv"), show_col_types = FALSE) |>
  dplyr::select(ID, M1ANMEDS_LOG) |>
  dplyr::distinct(ID, .keep_all = TRUE)

add_meds <- function(d) {
  if ("M1ANMEDS_LOG" %in% names(d)) d else dplyr::left_join(d, meds, by = "ID")
}

stopifnot(exists("m"), exists("analysis_df"))
m           <- add_meds(as.data.frame(m))
analysis_df <- add_meds(as.data.frame(analysis_df))
cat("M1ANMEDS_LOG merged. non-missing in m:",
    sum(!is.na(m$M1ANMEDS_LOG)), "/", nrow(m), "\n")

# ---- 0b. label/domain name
stopifnot(exists("outcome_labels"), exists("outcome_order"),
          exists("make_forest_b"), exists("save_fig"))

fix_vo2 <- function(x) sub("^VO(<sub>2</sub>|\u2082)[[:space:]]+peak", "VO\\1peak", x)
outcome_labels[] <- fix_vo2(outcome_labels)
outcome_order    <- fix_vo2(outcome_order)
cat("VO2 label in use:", outcome_labels[["VO2"]], "\n")

DOMAIN_ENERGY <- "Energetic capacity & efficiency"
domain_colors <- c("Physical performance" = "#0F6E56",
                   "Cognitive & frailty"  = "#534AB7")
domain_colors[DOMAIN_ENERGY] <- "#BA7517"
domain_label_df <- data.frame(
  y     = c(1.5, 3.5, 6.5),
  label = c("Energetic capacity\n& efficiency",
            "Cognitive &\nfrailty",
            "Physical\nperformance"),
  color = unname(domain_colors[c(DOMAIN_ENERGY, "Cognitive & frailty",
                                 "Physical performance")])
)
cat("3rd domain name:", DOMAIN_ENERGY, "\n")

# ---- 0c. get_std_beta 
get_std_beta_w <- function(outcome_var, data, extra_str = "",
                           q_var = "BOX_DISCO_JR_Q") {
  stopifnot(q_var %in% names(data))
  sd_out <- sd(data[[outcome_var]], na.rm = TRUE)
  data2  <- data |>
    dplyr::filter(!is.na(.data[[outcome_var]])) |>
    dplyr::mutate(
      Q_fac     = factor(.data[[q_var]], levels = 1:4),
      GDF15_std = (SOGDF15 - mean(SOGDF15, na.rm = TRUE)) / sd(SOGDF15, na.rm = TRUE),
      CysC_std  = (SOCYSC  - mean(SOCYSC,  na.rm = TRUE)) / sd(SOCYSC,  na.rm = TRUE)
    )
  base_cov_lm <- "D1AGE2 + EEFEMALE + SITE + HWWGT + HWHGT"
  fml <- as.formula(paste0(
    outcome_var, " ~ Q_fac + ", base_cov_lm,
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

# =============================================================================
# FIGURE a — phenotype forest (Fig 2 panel-b format), base + M1ANMEDS_LOG
# =============================================================================
outcomes_a <- list(
  VO2 = "TTPKVO2U", Speed = "NF400MPACE", Power = "LEPEAKPWR2", Steps = "ACSCFUAV",
  Digit = "DSCORR", VF = "FT0V2FN", OxPhos = "REMOXPHOS", CCR = "TTSSPKVO2"
)

# Build std with the SAME domain names used by domain_colors (correct colors).
std_med_a <- lapply(names(outcomes_a), function(nm) {
  d <- get_std_beta_w(outcomes_a[[nm]], m, extra_str = "M1ANMEDS_LOG")  # weighted DISCO + medication adjustment
  d$outcome <- nm; d
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    quartile_label = factor(paste0("Q", quartile), levels = c("Q4", "Q3", "Q2")),
    domain = dplyr::case_when(
      outcome %in% c("VO2", "Speed", "Power", "Steps") ~ "Physical performance",
      outcome %in% c("Digit", "VF")                    ~ "Cognitive & frailty",
      TRUE                                             ~ DOMAIN_ENERGY),  # <- CHANGE 2
    outcome_label = outcome_labels[outcome],       # <- VO2peak label applied automatically
    outcome_label = factor(outcome_label, levels = rev(outcome_order))
  )

cat("\n=== Figure a: std beta Q4 vs Q1 (base + M1ANMEDS_LOG) ===\n")
std_med_a |> dplyr::filter(quartile == 4) |>
  dplyr::select(outcome, beta, beta_lo, beta_hi, p_value) |>
  dplyr::mutate(dplyr::across(c(beta, beta_lo, beta_hi), \(x) round(x, 3)),
                p_value = round(p_value, 4)) |>
  as.data.frame() |> print()

fig_a <- make_forest_b(std_med_a, "a",
                       "Adjusted: age, sex, site, height, weight + medication count (log)")
save_fig(fig_a, file.path(PATH_RESULTS, "FigS_pheno_forest_meds.png"), w = 130, h = 140)

# =============================================================================
# FIGURE b — hospitalization Cox HR heatmap (Fig 5 sensitivity format)
#            single panel: base + M1ANMEDS_LOG
# =============================================================================
stopifnot(exists("organ_order_heat"), exists("add_heat_geoms"),
          exists("get_hr_all_q_box"), exists("extract_hr_full"), exists("tidy_hr"))

# full-proteome fit (Q_fac = whole-proteome DISCO quartile), base + M1ANMEDS_LOG
analysis_df$Q_fac <- factor(analysis_df$BOX_DISCO_JR_Q, levels = 1:4)
fit_med <- coxph(
  as.formula(paste0("Surv(HAFUTBLU, HAEMERG) ~ Q_fac + ", base_cov, " + M1ANMEDS_LOG")),
  data = analysis_df, ties = "efron", robust = TRUE, id = ID)

# organ-specific HRs, base + M1ANMEDS_LOG
hr_organ_med <- lapply(organs, function(o) {
  qv  <- paste0("BOX_DISCO_JR_Q_", o)
  res <- get_hr_all_q_box(qv, analysis_df, extra = "M1ANMEDS_LOG", unadj = FALSE)
  if (!is.null(res)) data.frame(organ = o, model = "Base + Medications", res)
}) |> dplyr::bind_rows()

heat_df_med <- dplyr::bind_rows(
  extract_hr_full(fit_med, "Base + Medications"),
  hr_organ_med
) |>
  tidy_hr(organ_order_heat) |>
  dplyr::mutate(model = factor(model, levels = "Base + Medications"))

cat("\n=== Figure b: Cox HR (Q4), base + M1ANMEDS_LOG ===\n")
heat_df_med |> dplyr::filter(quartile == "Q4") |>
  dplyr::select(organ_label, HR, lo, hi, sig_label) |>
  dplyr::mutate(dplyr::across(c(HR, lo, hi), \(x) round(x, 2))) |>
  as.data.frame() |> print()

# single heatmap panel (mirrors make_heat_sens with show_y = TRUE)
p_b_heat <- add_heat_geoms(
  ggplot2::ggplot(heat_df_med, ggplot2::aes(x = quartile, y = organ_label))
) +
  heat_fill_scale + heat_color_scale +
  ggplot2::scale_x_discrete(name = "Quartile vs Q1",
                            labels = c("Q2", "Q3", "Q4 (highest)")) +
  ggplot2::scale_y_discrete(name = NULL, limits = rev(organ_order_heat)) +
  ggplot2::labs(title = "b   Base + medication count (log) adjusted") +
  heat_theme +
  ggplot2::theme(
    plot.title  = ggplot2::element_text(face = "bold", size = 9, hjust = 0),
    axis.text.y = ggplot2::element_text(size = 8.5, color = "black", face = y_faces_heat),
    legend.position = "none")

leg_b <- cowplot::get_legend(
  p_b_heat + ggplot2::theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = ggplot2::element_text(size = 7.5),
    legend.title = ggplot2::element_text(size = 8),
    legend.key.size = ggplot2::unit(0.4, "cm"),
    legend.key.width = ggplot2::unit(1.8, "cm")))

fig_b <- cowplot::plot_grid(
  cowplot::plot_grid(p_b_heat, leg_b, ncol = 1, rel_heights = c(1, 0.08)),
  cowplot::ggdraw() + cowplot::draw_label(
    paste0("HR (95% CI). Bold = HR value. * p<0.05, ** p<0.01, *** p<0.001. ",
           "Base model: age, sex, site, body weight, height + medication count (log). ",
           "Grey = HR < 1.00. Sorted by Q4 unadjusted HR."),
    x = 0.01, y = 0.5, hjust = 0, vjust = 0.5, size = 6.5, color = "grey40"),
  ncol = 1, rel_heights = c(1, 0.06)) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))

ggplot2::ggsave(file.path(PATH_RESULTS, "FigS_hosp_HR_meds.pdf"), fig_b,
                width = 150, height = 170, units = "mm", device = cairo_pdf)
ggplot2::ggsave(file.path(PATH_RESULTS, "FigS_hosp_HR_meds.png"), fig_b,
                width = 150, height = 170, units = "mm", dpi = 300)

# CSV of the HR table
heat_df_med |>
  dplyr::select(organ_label, quartile, HR, lo, hi, p, sig_label) |>
  dplyr::arrange(organ_label, quartile) |>
  readr::write_csv(file.path(PATH_RESULTS, "Table_hosp_HR_meds.csv"))

# =============================================================================
# COMBINED — Figure a (forest) + Figure b (heatmap) in one panel
# =============================================================================
# left column: forest with its own quartile legend below
left_a <- cowplot::plot_grid(
  fig_a + ggplot2::theme(legend.position = "none"),
  cowplot::get_legend(fig_a + ggplot2::theme(
    legend.position = "bottom", legend.direction = "horizontal")),
  ncol = 1, rel_heights = c(1, 0.10))

# right column: heatmap with its own HR gradient legend below
right_b <- cowplot::plot_grid(
  p_b_heat, leg_b, ncol = 1, rel_heights = c(1, 0.08))

combined_body <- cowplot::plot_grid(
  left_a, right_b, ncol = 2, rel_widths = c(1, 1.15), align = "h", axis = "t")

fig_ab <- cowplot::plot_grid(
  combined_body,
  cowplot::ggdraw() + cowplot::draw_label(
    paste0("a: linear regression, standardized \u03b2 vs. Q1 (outcome SD units). ",
           "b: Cox HR (95% CI) for non-elective hospitalization; grey = HR < 1.00, ",
           "sorted by Q4 unadjusted HR. Both models adjusted for age, sex, site, ",
           "body weight, height, and medication count (log). * p<0.05, ** p<0.01, *** p<0.001."),
    x = 0.01, y = 0.5, hjust = 0, vjust = 0.5, size = 6.3, color = "grey40"),
  ncol = 1, rel_heights = c(1, 0.06)) +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))

ggplot2::ggsave(file.path(PATH_RESULTS, "FigS_meds_combined.pdf"), fig_ab,
                width = 300, height = 175, units = "mm", device = cairo_pdf)
ggplot2::ggsave(file.path(PATH_RESULTS, "FigS_meds_combined.png"), fig_ab,
                width = 300, height = 175, units = "mm", dpi = 300)

cat("\nSaved: FigS_pheno_forest_meds.png | FigS_hosp_HR_meds.png/.pdf | ",
    "FigS_meds_combined.png/.pdf | Table_hosp_HR_meds.csv\n")
