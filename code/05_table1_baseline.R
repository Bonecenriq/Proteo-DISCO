# ==============================================================================
# STUDY       : SOMMA
# PROJECT     : SOMMA SomaLogic
# PROGRAM     : 05_table1_baseline.R
#
# PURPOSE     : Produce Table 1 (baseline characteristics), overall and by
#               quartile of proteomic entropy (DISCO), with group comparisons
#               and Benjamini-Hochberg adjusted post-hoc tests.
#
# INPUTS      : - meta_df (04_analysis_dataset.R)
#
# OUTPUTS     : - Table1_baseline.csv
#               - Table1_baseline.png
#               - Table1_by_DISCO_quartile.csv
#               - Table1_by_DISCO_quartile_posthoc.csv
#               - Table1_by_DISCO_quartile.xlsx (if openxlsx or writexl is available)
#
# DATA ACCESS : Individual-level SOMMA data are not included in this
#               repository. Access to the underlying study data is required
#               to reproduce the analyses.
#
# DEPENDENCIES: - 99_figure_setup.R
#               - 04_analysis_dataset.R (object meta_df)
#               - Optional: openxlsx or writexl (Excel output), FSA (Dunn's test)
#
# VERSION     : 1.0
# DATE        : 2026-08-15
# AUTHOR      : Namki Hong
#
# REVISION HISTORY
# 08/15/26    Creation (split out of the combined analysis script)
#
# NOTES       : Continuous variables are summarised as mean +/- SD with one-way
#               ANOVA, except daily step count and Vigor-to-Frailty, which are
#               summarised as median [IQR] with the Kruskal-Wallis test.
#               Categorical variables are summarised as n (%) with the chi-square
#               test, or Fisher's exact test when any expected count is below 5.
#
#               The first block reproduces the original overall summary table; the
#               second block produces the final quartile-stratified Table 1.
# ==============================================================================
source("code/99_figure_setup.R")


# =============================================================================
# 1. OVERALL BASELINE SUMMARY TABLE
# =============================================================================

# helper functions
mn_sd <- function(x, d=1) {
  sprintf("%s \u00b1 %s",
          round(mean(x, na.rm=TRUE), d),
          round(sd(x, na.rm=TRUE), d))
}
n_pct <- function(x, val=1, d=1) {
  n   <- sum(x == val, na.rm=TRUE)
  tot <- sum(!is.na(x))
  sprintf("%d (%.1f%%)", n, 100*n/tot)
}

tbl1_data <- list(
  # Demographics
  "Clinical characteristics"                          =
    "Values",
  "Age, years"                 =
    mn_sd(meta_df$D1AGE2, 1),
  "Women, n (%)"               =
    n_pct(meta_df$EEFEMALE),
  "White race, n (%)"          =
    n_pct(meta_df$EERACE, 1),
  "Pittsburgh site, n (%)"     =
    n_pct(meta_df$SITE, 1),
  # Anthropometrics
  "Height, cm"                 =
    mn_sd(meta_df$HWHGTCM, 1),
  "Weight, kg"                 =
    mn_sd(meta_df$HWWGT, 1),
  "BMI, kg/m\u00b2"            =
    mn_sd(meta_df$HWBMI, 1),
  # Physical function
  "400m walking speed, m/s"    =
    mn_sd(meta_df$NF400MPACE, 3),
  "Leg peak power, W"          =
    mn_sd(meta_df$LEPEAKPWR2, 1),
  "Vigor-to-Frailty (0-12)"    =
    mn_sd(meta_df$FT0V2FN, 1),
  "Digit Symbol Substitution Test (0-133)"         =
    mn_sd(meta_df$DSCORR, 1),
  # Cardiorespiratory
  "VO\u2082 peak, mL/kg/min"   =
    mn_sd(meta_df$TTPKVO2, 1),
  "Cost-capacity ratio, %"        =
    mn_sd(meta_df$TTSSPKVO2, 1),
  "Max OxPhos, pmol/s/mg"      =
    sprintf("%s (n=%d)",
            mn_sd(meta_df$REMOXPHOS, 1),
            sum(!is.na(meta_df$REMOXPHOS))),
  # Physical activity
  "Daily step count"           =
    mn_sd(meta_df$ACSCFUAV, 0),
  # Biomarkers
  "GDF-15, log-transformed"    =
    mn_sd(meta_df$SOGDF15, 3),
  "Cystatin-C, log-transformed" =
    mn_sd(meta_df$SOCYSC, 3),
  # DISCO
  "DISCO score"                =
    mn_sd(meta_df$BOX_DISCO_JR, 3),
  "Quartile 1/2/3/4, n"        =
    paste(table(meta_df$BOX_DISCO_JR_Q),
          collapse=" / "),
  # Outcome
  "Non-elective hospitalization, n (%)" =
    sprintf("%d (%.1f%%)",
            sum(meta_df$HAEMERG, na.rm=TRUE),
            100*mean(meta_df$HAEMERG,
                     na.rm=TRUE))
)

tbl1_df <- data.frame(
  Characteristic = names(tbl1_data),
  Value          = unlist(tbl1_data),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Section headers
sections <- c(
  "Characteristics"      = "Demographics",
  "Height, cm"           = "Anthropometrics",
  "D3Cr muscle mass, kg" = "Body Composition",
  "Grip strength, kg"    = "Physical Function",
  "VO\u2082 peak, mL/kg/min" = "Cardiorespiratory & Metabolic",
  "Daily step count"     = "Physical Activity",
  "GDF-15, log-transformed" = "Circulating Biomarkers",
  "DISCO score"          = "Proteomic Entropy",
  "Non-elective hospitalization, n (%)" = "Clinical Outcome"
)

# Write CSV
readr::write_csv(tbl1_df,
                 file.path(PATH_RESULTS, "Table1_baseline.csv"))
cat("Saved: Table1_baseline.csv\n")

# Figure version (ggplot table)
p_tbl1 <- ggplot(
  tbl1_df,
  aes(x = 1, y = rev(seq_len(nrow(tbl1_df))))
) +
  geom_text(
    aes(label = Characteristic),
    x = 0.05, hjust = 0, size = 3.5
  ) +
  geom_text(
    aes(label = Value),
    x = 0.95, hjust = 1, size = 3.5
  ) +
  geom_hline(
    yintercept = nrow(tbl1_df) - 0.5,
    linewidth  = 0.8, color = "black"
  ) +
  geom_hline(
    yintercept = 0.5,
    linewidth  = 0.8, color = "black"
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = paste0(
    "Table 1. Baseline characteristics ",
    "(n = 838)")) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold", size = 11,
      hjust = 0,
      margin = margin(b=4))
  )

png(file.path(PATH_RESULTS, "Table1_baseline.png"),
    width  = 120 / 25.4,
    height = 200 / 25.4,
    units  = "in", res = 200)
print(p_tbl1)
dev.off()
cat("Saved: Table1_baseline.png\n")

# =============================================================================
# 2. TABLE 1 BY DISCO QUARTILE (final version)
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
.have <- function(p) requireNamespace(p, quietly = TRUE)
if (!.have("openxlsx") && !.have("writexl")) {
  try(install.packages("openxlsx"), silent = TRUE)
}
use_openxlsx <- .have("openxlsx")
use_writexl  <- !use_openxlsx && .have("writexl")
if (.have("FSA")) suppressMessages(require(FSA))   # optional: proper Dunn's test

out_dir <- PATH_RESULTS
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Data ----------------------------------------------------------------
# Expects `meta_df` (with BOX_DISCO_JR_Q = 1..4 and the baseline variables) to
# already exist in the session, i.e. after running 04_analysis_dataset.R.
# If starting fresh, uncomment the standalone loader below.
#
# library(readr); library(dplyr)
# pheno <- read_csv(file.path(PATH_DATA, "02_box_transformed_data_with_pheno.csv"), show_col_types = FALSE) |>
#   filter(COHORT == "Parent")
# ent   <- read_csv(file.path(PATH_OUTPUT, "05_entropy_measures.csv"), show_col_types = FALSE) |>
#   select(ID, BOX_DISCO_JR, BOX_DISCO_JR_Q) |> distinct(ID, .keep_all = TRUE)
# meta_df <- pheno |> select(-starts_with("seq.")) |> left_join(ent, by = "ID")

stopifnot(exists("meta_df"))
df      <- as.data.frame(meta_df)
grp_var <- "BOX_DISCO_JR_Q"
stopifnot(grp_var %in% names(df))

df  <- df[!is.na(df[[grp_var]]), ]
grp <- factor(df[[grp_var]], levels = sort(unique(df[[grp_var]])))
levels(grp) <- paste0("Q", levels(grp))
glv   <- levels(grp)
gsize <- as.integer(table(grp))
K     <- nlevels(grp)

# ---- 2. Variable specification ----------------------------------------------
# ty  : "cont" (continuous) | "cat" (categorical)
# test: "param"    -> mean (SD)     + ANOVA        + Tukey/pairwise-t (BH)   [DEFAULT]
#       "nonparam" -> median [IQR]  + Kruskal-Wallis + Dunn/pairwise-Wilcoxon (BH)
#       "auto"     -> decide by Shapiro-Wilk (kept available but not used here)
# lev : value counted as event for binary categorical n (%); omit for multi-level
# d   : decimals for continuous summaries

var_spec <- list(
  list(v = "D1AGE2",     lab = "Age, years",                          ty = "cont", test = "param",    d = 1),
  list(v = "EEFEMALE",   lab = "Women, n (%)",                        ty = "cat",  lev = 1),
  list(v = "EERACE",     lab = "White race, n (%)",                   ty = "cat",  lev = 1),
  list(v = "SITE",       lab = "Pittsburgh site, n (%)",              ty = "cat",  lev = 1),
  list(v = "HWHGTCM",    lab = "Height, cm",                          ty = "cont", test = "param",    d = 1),
  list(v = "HWWGT",      lab = "Weight, kg",                          ty = "cont", test = "param",    d = 1),
  list(v = "HWBMI",      lab = "BMI, kg/m2",                          ty = "cont", test = "param",    d = 1),
  list(v = "NF400MPACE", lab = "400m walking speed, m/s",             ty = "cont", test = "param",    d = 3),
  list(v = "LEPEAKPWR2", lab = "Leg peak power, W",                   ty = "cont", test = "param",    d = 1),
  list(v = "FT0V2FN",    lab = "Vigor-to-Frailty (0-12)",             ty = "cont", test = "nonparam", d = 1),
  list(v = "DSCORR",     lab = "DSST (0-133)",                        ty = "cont", test = "param",    d = 1),
  list(v = "TTPKVO2",    lab = "VO2 peak, mL/kg/min",                 ty = "cont", test = "param",    d = 1),
  list(v = "TTSSPKVO2",  lab = "Cost-capacity ratio, %",              ty = "cont", test = "param",    d = 1),
  list(v = "REMOXPHOS",  lab = "Max OxPhos, pmol/s/mg",               ty = "cont", test = "param",    d = 1),
  list(v = "ACSCFUAV",   lab = "Daily step count",                    ty = "cont", test = "nonparam", d = 0),
  list(v = "SOGDF15",    lab = "GDF-15, log",                         ty = "cont", test = "param",    d = 3),
  list(v = "SOCYSC",     lab = "Cystatin-C, log",                     ty = "cont", test = "param",    d = 3),
  list(v = "BOX_DISCO_JR", lab = "DISCO score",                       ty = "cont", test = "param",    d = 3),
  list(v = "HAEMERG",    lab = "Non-elective hospitalization, n (%)", ty = "cat",  lev = 1)
)
# keep only variables present in the data
var_spec <- Filter(function(s) s$v %in% names(df), var_spec)

# ---- 3. Helpers -------------------------------------------------------------
fmt_p <- function(p) ifelse(is.na(p), "-",
                            ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

is_normal <- function(x, g) {                    # TRUE only if normal in every group
  for (lv in levels(g)) {
    xi <- x[g == lv]; xi <- xi[!is.na(xi)]
    if (length(xi) >= 3 && length(xi) <= 5000 && stats::var(xi) > 0) {
      if (stats::shapiro.test(xi)$p.value < 0.05) return(FALSE)
    }
  }
  TRUE
}
mean_sd <- function(x, d) sprintf("%s \u00b1 %s",   # mean \u00b1 SD
                                  formatC(mean(x, na.rm = TRUE), format = "f", digits = d),
                                  formatC(sd(x,   na.rm = TRUE), format = "f", digits = d))
med_iqr <- function(x, d) { q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
sprintf("%s [%s, %s]", formatC(q[2], format = "f", digits = d),
        formatC(q[1], format = "f", digits = d),
        formatC(q[3], format = "f", digits = d)) }
n_pct   <- function(x, val) { n <- sum(x == val, na.rm = TRUE); tot <- sum(!is.na(x))
if (tot == 0) return("-"); sprintf("%d (%.1f%%)", n, 100 * n / tot) }

# continuous post-hoc -> long df of pairwise adjusted p-values (with direction)
posthoc_cont <- function(x, g, method, lab) {
  if (method == "param") {
    pm <- pairwise.t.test(x, g, p.adjust.method = "BH", pool.sd = TRUE)$p.value
    mname <- "ANOVA post-hoc: pairwise t (BH)"
  } else if (.have("FSA")) {
    dt <- FSA::dunnTest(x ~ g, method = "bh")$res     # Comparison, Z, P.unadj, P.adj
    res <- do.call(rbind, lapply(seq_len(nrow(dt)), function(i) {
      pr <- strsplit(dt$Comparison[i], " - ")[[1]]
      data.frame(Variable = lab, Comparison = paste(pr[1], "vs", pr[2]),
                 p_adj = dt$P.adj[i], stringsAsFactors = FALSE)
    }))
    res$Method <- "Kruskal-Wallis post-hoc: Dunn (BH)"
    return(res)
  } else {
    pm <- pairwise.wilcox.test(x, g, p.adjust.method = "BH", exact = FALSE)$p.value
    mname <- "Kruskal-Wallis post-hoc: pairwise Wilcoxon (BH)"
  }
  out <- list()
  for (i in rownames(pm)) for (j in colnames(pm)) {
    pij <- pm[i, j]
    if (!is.na(pij))
      out[[length(out) + 1]] <- data.frame(Variable = lab,
                                           Comparison = paste(i, "vs", j), p_adj = pij,
                                           Method = mname, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# direction label for a significant continuous pair, based on group centre
dir_lab <- function(x, g, a, b, method) {
  ca <- if (method == "param") mean(x[g == a], na.rm = TRUE) else median(x[g == a], na.rm = TRUE)
  cb <- if (method == "param") mean(x[g == b], na.rm = TRUE) else median(x[g == b], na.rm = TRUE)
  if (ca >= cb) paste0(a, ">", b) else paste0(a, "<", b)
}

# ---- 4. Build table ---------------------------------------------------------
rows <- list(); posthoc_all <- list()

for (s in var_spec) {
  x <- df[[s$v]]
  
  if (s$ty == "cont") {
    method <- if (is.null(s$test) || s$test == "auto")
      (if (is_normal(x, grp)) "param" else "nonparam") else s$test
    summ_fun <- if (method == "param") function(z) mean_sd(z, s$d) else function(z) med_iqr(z, s$d)
    
    overall <- summ_fun(x)
    cells   <- vapply(glv, function(lv) summ_fun(x[grp == lv]), character(1))
    
    if (method == "param") {
      p <- tryCatch(summary(aov(x ~ grp))[[1]][["Pr(>F)"]][1], error = function(e) NA)
      testname <- "ANOVA"
    } else {
      p <- tryCatch(kruskal.test(x ~ grp)$p.value, error = function(e) NA)
      testname <- "Kruskal-Wallis"
    }
    
    ph <- tryCatch(posthoc_cont(x, grp, method, s$lab), error = function(e) NULL)
    sig <- character(0)
    if (!is.null(ph)) {
      posthoc_all[[length(posthoc_all) + 1]] <- ph
      for (r in which(ph$p_adj < 0.05)) {
        pr <- strsplit(ph$Comparison[r], " vs ")[[1]]
        sig <- c(sig, dir_lab(x, grp, pr[1], pr[2], method))
      }
    }
    rows[[length(rows) + 1]] <- c(s$lab, overall, cells, fmt_p(p), testname, paste(sig, collapse = ", "))
    
  } else {  # categorical
    lev <- if (!is.null(s$lev)) s$lev else NA
    if (!is.na(lev)) {
      overall <- n_pct(x, lev)
      cells   <- vapply(glv, function(lv) n_pct(x[grp == lv], lev), character(1))
    } else {
      overall <- ""; cells <- rep("", K)   # multi-level: summarised in PostHoc only
    }
    tab <- table(grp, factor(x))
    ex  <- suppressWarnings(chisq.test(tab)$expected)
    if (any(ex < 5)) {
      p <- tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 1e4)$p.value, error = function(e) NA)
      testname <- "Fisher"
    } else {
      p <- suppressWarnings(chisq.test(tab)$p.value); testname <- "Chi-square"
    }
    # pairwise between quartiles
    combs <- combn(glv, 2); praw <- numeric(ncol(combs)); labs <- character(ncol(combs))
    for (k in seq_len(ncol(combs))) {
      a <- combs[1, k]; b <- combs[2, k]
      sub <- df[grp %in% c(a, b), , drop = FALSE]
      tt  <- table(droplevels(grp[grp %in% c(a, b)]), factor(sub[[s$v]]))
      exx <- suppressWarnings(chisq.test(tt)$expected)
      praw[k] <- tryCatch(
        if (any(exx < 5)) fisher.test(tt)$p.value else suppressWarnings(chisq.test(tt)$p.value),
        error = function(e) NA)
      labs[k] <- paste(a, "vs", b)
    }
    padj <- p.adjust(praw, "BH")
    posthoc_all[[length(posthoc_all) + 1]] <- data.frame(
      Variable = s$lab, Comparison = labs, p_adj = padj,
      Method = paste0(testname, " post-hoc: pairwise (BH)"), stringsAsFactors = FALSE)
    sig <- labs[!is.na(padj) & padj < 0.05]
    rows[[length(rows) + 1]] <- c(s$lab, overall, cells, fmt_p(p), testname, paste(sig, collapse = ", "))
  }
}

tab_df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
colnames(tab_df) <- c("Characteristic", sprintf("Overall (n=%d)", nrow(df)),
                      sprintf("%s (n=%d)", glv, gsize),
                      "P value", "Test", "Post hoc (p<0.05)")
rownames(tab_df) <- NULL

ph_df <- do.call(rbind, posthoc_all)
ph_df$p_adj_fmt <- fmt_p(ph_df$p_adj)
ph_df <- ph_df[, c("Variable", "Comparison", "Method", "p_adj_fmt")]
names(ph_df)[4] <- "p (BH-adjusted)"

# ---- 5. Write output --------------------------------------------------------
# Locale-proof UTF-8 CSV writer (keeps the "±" byte-correct on any OS/locale,
# incl. Excel opening it as UTF-8; add a BOM so Excel auto-detects encoding).
write_csv_utf8 <- function(dat, path, bom = TRUE) {
  q <- function(v) paste0("\"", gsub("\"", "\"\"", as.character(v)), "\"")
  lines <- c(paste(q(names(dat)), collapse = ","),
             apply(dat, 1, function(r) paste(q(r), collapse = ",")))
  con <- file(path, open = "wb")
  if (bom) writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)   # UTF-8 BOM for Excel
  writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
  close(con)
}
csv1 <- file.path(out_dir, "Table1_by_DISCO_quartile.csv")
csv2 <- file.path(out_dir, "Table1_by_DISCO_quartile_posthoc.csv")
write_csv_utf8(tab_df, csv1)
write_csv_utf8(ph_df,  csv2)

xlsx_path <- file.path(out_dir, "Table1_by_DISCO_quartile.xlsx")

if (use_openxlsx) {
  library(openxlsx)
  wb <- createWorkbook()
  hdr <- createStyle(textDecoration = "bold", fgFill = "#4F2D2D", fontColour = "white",
                     halign = "center", valign = "center", border = "TopBottom")
  chr <- createStyle(halign = "left"); ctr <- createStyle(halign = "center")
  note <- createStyle(fontSize = 9, textDecoration = "italic", fontColour = "#666666")
  
  addWorksheet(wb, "Table1")
  writeData(wb, "Table1", "Table 1. Baseline characteristics by DISCO quartile",
            startRow = 1, startCol = 1)
  addStyle(wb, "Table1", createStyle(textDecoration = "bold", fontSize = 12), 1, 1)
  writeData(wb, "Table1", tab_df, startRow = 3, headerStyle = hdr, borders = "rows",
            borderColour = "#DDDDDD")
  addStyle(wb, "Table1", ctr, rows = 4:(3 + nrow(tab_df)),
           cols = 2:ncol(tab_df), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, "Table1", chr, rows = 4:(3 + nrow(tab_df)), cols = 1, gridExpand = TRUE, stack = TRUE)
  nrow_note <- 3 + nrow(tab_df) + 2
  writeData(wb, "Table1", paste0(
    "Continuous variables are mean \u00b1 SD with one-way ANOVA, except daily step ",
    "count and Vigor-to-Frailty, which are median [IQR] with Kruskal-Wallis. ",
    "Categorical: n (%), Chi-square (Fisher if any expected count <5). ",
    "Post-hoc pairwise comparisons Benjamini-Hochberg adjusted; only pairs with ",
    "adjusted p<0.05 are listed (for continuous, direction of the higher group is shown)."),
    startRow = nrow_note, startCol = 1)
  addStyle(wb, "Table1", note, nrow_note, 1)
  setColWidths(wb, "Table1", cols = 1, widths = 34)
  setColWidths(wb, "Table1", cols = 2:(ncol(tab_df) - 1), widths = 18)
  setColWidths(wb, "Table1", cols = ncol(tab_df), widths = 40)
  freezePane(wb, "Table1", firstActiveRow = 4, firstActiveCol = 2)
  
  addWorksheet(wb, "PostHoc")
  writeData(wb, "PostHoc", ph_df, headerStyle = hdr, borders = "rows", borderColour = "#DDDDDD")
  setColWidths(wb, "PostHoc", cols = 1:4, widths = c(34, 16, 44, 16))
  freezePane(wb, "PostHoc", firstActiveRow = 2)
  
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat("Saved Excel:", xlsx_path, "\n")
  
} else if (use_writexl) {
  writexl::write_xlsx(list(Table1 = tab_df, PostHoc = ph_df), xlsx_path)
  cat("Saved Excel (writexl, unstyled):", xlsx_path, "\n")
} else {
  cat("Neither openxlsx nor writexl available; wrote CSVs instead:\n  ", csv1, "\n  ", csv2, "\n")
}

cat("Saved CSV backups:\n  ", csv1, "\n  ", csv2, "\n")
cat("\n--- Table 1 preview ---\n"); print(tab_df, row.names = FALSE)
