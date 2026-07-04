# =============================================================================
# SOMMA Proteomics — Complete Analysis Pipeline v14 Final
# =============================================================================

# =============================================================================
# 1. PACKAGES
# =============================================================================

library(tidyverse)
library(patchwork)
library(cowplot)
library(ggrepel)
library(RColorBrewer)
library(SomaScan.db)
library(limma)
library(clusterProfiler)
library(msigdbr)
library(cluster)
library(zoo)
library(survival)
library(sandwich)
library(lmtest)
library(GSVA)


# =============================================================================
# 2. DATA LOADING & PREPROCESSING
# =============================================================================

setwd("path/to/project")  # <- set this to your working directory (folder containing the input CSVs)
dir.create("results", showWarnings = FALSE)

# --- 2a. Protein expression (Box-Cox transformed, Parent n=838) ---
pheno_raw <- readr::read_csv(
  "02_box_transformed_data_with_pheno.csv",
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

# --- 2b. DISCO entropy merge ---
entropy_raw <- readr::read_csv(
  "05_entropy_measures.csv",
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

# --- 2c. Sort by BOX_DISCO_JR ---
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
     file = "results/checkpoint_section2.RData")
cat("Checkpoint saved.\n")


# =============================================================================
# 3. PROTEIN NAME MAPPING + DEA (continuous DISCO)
# =============================================================================

# --- 3a. SomaScan.db mapping ---
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

# --- 3b. limma DEA (continuous BOX_DISCO_JR) ---
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
  "results/DEA_continuous_BOX_DISCO.csv"
)
cat("Saved: DEA_continuous_BOX_DISCO.csv\n")


# =============================================================================
# 4. Figure 1b. TABLE (Baseline Characteristics) -> Changes as table 1
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

# Save CSV
readr::write_csv(tbl1_df,
                 "results/Table1_baseline.csv")
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

png("results/Table1_baseline.png",
    width  = 120 / 25.4,
    height = 200 / 25.4,
    units  = "in", res = 200)
print(p_tbl1)
dev.off()
cat("Saved: Table1_baseline.png\n")

# =============================================================================
# TABLE 1 (Baseline Characteristics) -> final version
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
.have <- function(p) requireNamespace(p, quietly = TRUE)
if (!.have("openxlsx") && !.have("writexl")) {
  try(install.packages("openxlsx"), silent = TRUE)
}
use_openxlsx <- .have("openxlsx")
use_writexl  <- !use_openxlsx && .have("writexl")
if (.have("FSA")) suppressMessages(require(FSA))   # optional: proper Dunn's test

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Data ----------------------------------------------------------------
# Expects `meta_df` (with BOX_DISCO_JR_Q = 1..4 and the baseline variables) to
# already exist in the session, i.e. after running Sections 1-2 of the main
# pipeline. If starting fresh, uncomment the standalone loader below.
#
# library(readr); library(dplyr)
# pheno <- read_csv("02_box_transformed_data_with_pheno.csv", show_col_types = FALSE) |>
#   filter(COHORT == "Parent")
# ent   <- read_csv("05_entropy_measures.csv", show_col_types = FALSE) |>
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


# =============================================================================
# FIGURE 2: PHENOTYPE ASSOCIATION & SUPPLEMENTARY FIGURE 1
# =============================================================================

# --- 5a. get_std_beta function ---
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

outcome_labels <- c(
  VO2    = "VO\u2082 peak (mL/min)",
  Speed  = "Walking speed (m/s)",
  Power  = "Leg peak power (W)",
  Steps  = "Daily step count",
  Digit  = "DSST score (0-133)\n(Lower = poorer cognition)",
  VF     = "Vigor-to-Frailty (0-12)\n(Higher = more frail)",
  OxPhos = "Max OxPhos (pmol/s/mg)",
  CCR    = "Cost-capacity ratio (%)",
  GDF15  = "GDF-15 (log)",
  CysC   = "Cystatin-C (log)"
)

outcome_order <- c(
  "VO\u2082 peak (mL/min)",
  "Walking speed (m/s)",
  "Leg peak power (W)",
  "Daily step count",
  "DSST score (0-133)\n(Lower = poorer cognition)",
  "Vigor-to-Frailty (0-12)\n(Higher = more frail)",
  "Max OxPhos (pmol/s/mg)",
  "Cost-capacity ratio (%)"
)

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
        "Mitochondrial function",
      TRUE ~ "Biomarker"),
    # ── Fix: character first, then factor conversion ──────────────
    outcome_label = outcome_labels[outcome],  # character
    outcome_label = factor(
      outcome_label,
      levels = c(
        rev(outcome_order),     # existing 8 (for panel b)
        "GDF-15 (log)",         # added
        "Cystatin-C (log)"      # added
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
  "Physical performance"       = "#0F6E56",
  "Cognitive & frailty"   = "#534AB7",
  "Mitochondrial function" = "#BA7517"
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

domain_label_df <- data.frame(
  y     = c(1.5, 3.5, 6.5),
  label = c("Mitochondrial\nfunction",
            "Cognitive &\nfrailty",
            "Physical\nperformance"),
  color = unname(domain_colors[c(
    "Mitochondrial function",
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

# ── edit: std_bio as.character() add ────────────────────────
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

# --- Panel b: Forest (physical/cog/metabolic) ---
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
    aes(x=0.82, y=y, label=label,
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
  labs(title="b",
       x="Adjusted difference vs. Q1 (outcome SD units)",
       y=NULL) +
  coord_cartesian(clip="off") +
  theme_bw(base_size=9) +
  theme(
    plot.title         = element_text(
      face="bold", size=10, hjust=0),
    axis.text.y        = element_text(size=8.5),
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

# --- Panel c: Biomarker (revised) ---
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
    axis.text.y        = element_text(size=8.5),
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

fig2_v8 <- cowplot::plot_grid(
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

png("results/Fig2_v8.png",
    width=220/25.4, height=150/25.4,
    units="in", res=300)
print(fig2_v8)
dev.off()
cat("Saved: results/Fig2_v8.png\n")

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
        TRUE ~ "Mitochondrial function"),
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
      aes(x=0.82, y=y, label=label,
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
      axis.text.y        = element_text(
        size=8.5),
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
         "results/Fig2_panel_b_base.png")
save_fig(fig_b_v2,
         "results/Fig2_panel_b_GDF15.png")
save_fig(fig_b_v3,
         "results/Fig2_panel_b_CysC.png")

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

png("results/Fig2_sensitivity_3panel.png",
    width=330/25.4, height=130/25.4,
    units="in", res=300)
print(fig2_sensitivity)
dev.off()
cat("Saved: results/Fig2_sensitivity_3panel.png\n")

std_v1 <- make_std_df("",          "Base")
std_v2 <- make_std_df("GDF15_std", "Base + GDF-15")
std_v3 <- make_std_df("CysC_std",  "Base + Cystatin-C")

# v1
fig_b_v1 <- make_forest_b(std_v1, NULL,
                          "Adjusted: age, sex, site, height, weight")
fig_b_v2 <- make_forest_b(std_v2, "a",
                          "Adjusted: base + GDF-15 (continuous)")
fig_b_v3 <- make_forest_b(std_v3, "b",
                          "Adjusted: base + Cystatin-C (continuous)")

# 
save_fig(fig_b_v2,
         "results/Fig2_panel_b_GDF15.png")
save_fig(fig_b_v3,
         "results/Fig2_panel_b_CysC.png")

# 
leg_shared <- cowplot::get_legend(
  fig_b_v1 + theme(
    legend.position="bottom",
    legend.direction="horizontal"))

# GDF-15 / Cystatin-C adujsted panels
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

png("results/FigS_sensitivity_2panel.png",
    width=220/25.4, height=130/25.4,
    units="in", res=300)
print(fig2_sensitivity)
dev.off()
cat("Saved: results/FigS_sensitivity_2panel.png\n")

# =============================================================================
# Supplementary Figure 2: Age-unweighted DISCO 
# =============================================================================

# =============================================================================
# SUPPLEMENTARY FIGURE — Unweighted DISCO
#   Panel a: correlation between weighted and unweighted DISCO
#   Panel b: association of UNWEIGHTED-DISCO quartiles with cross-sectional
#            outcomes — SAME forest format as main-text Figure 2b
# Only difference vs. Fig 2b: quartiles come from unweighted DISCO.
# =============================================================================

library(dplyr); library(ggplot2); library(cowplot); library(readr)

out_dir <- "results"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 0. Merge unweighted DISCO into the analysis data `m` -------------------
# `m` is the same analysis data frame used for Figure 2 (weighted BOX_DISCO_JR
# + all outcomes). We add the unweighted columns with a _UNW tag.
unw <- readr::read_csv("05z_entropy_measures_unweighted.csv", show_col_types = FALSE) |>
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
outcome_labels <- c(
  VO2 = "VO\u2082 peak (mL/min)", Speed = "Walking speed (m/s)",
  Power = "Leg peak power (W)", Steps = "Daily step count",
  Digit = "DSST score (0-133)\n(Lower = poorer cognition)",
  VF = "Vigor-to-Frailty (0-12)\n(Higher = more frail)",
  OxPhos = "Max OxPhos (pmol/s/mg)", CCR = "Cost-capacity ratio (%)",
  GDF15 = "GDF-15 (log)", CysC = "Cystatin-C (log)"
)
outcome_order <- c(
  "VO\u2082 peak (mL/min)", "Walking speed (m/s)", "Leg peak power (W)",
  "Daily step count", "DSST score (0-133)\n(Lower = poorer cognition)",
  "Vigor-to-Frailty (0-12)\n(Higher = more frail)",
  "Max OxPhos (pmol/s/mg)", "Cost-capacity ratio (%)"
)

std_df2 <- lapply(names(outcomes_all), function(nm) {
  d <- get_std_beta(outcomes_all[[nm]], m_unw); d$outcome <- nm; d
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    quartile_label = factor(paste0("Q", quartile), levels = c("Q4", "Q3", "Q2")),
    domain = dplyr::case_when(
      outcome %in% c("VO2", "Speed", "Power", "Steps") ~ "Physical performance",
      outcome %in% c("Digit", "VF")                    ~ "Cognitive & frailty",
      outcome %in% c("OxPhos", "CCR")                  ~ "Mitochondrial function",
      TRUE                                             ~ "Biomarker"),
    outcome_label = outcome_labels[outcome],
    outcome_label = factor(outcome_label,
                           levels = c(rev(outcome_order), "GDF-15 (log)", "Cystatin-C (log)"))
  )

cat("=== Unweighted-DISCO Std beta, Q4 vs Q1 ===\n")
std_df2 |> dplyr::filter(quartile == 4) |>
  dplyr::select(outcome, beta, beta_lo, beta_hi, p_value) |>
  dplyr::mutate(dplyr::across(c(beta, beta_lo, beta_hi), round, 3),
                p_value = round(p_value, 4)) |>
  as.data.frame() |> print()

# ---- 3. Aesthetics (identical to Figure 2) ----------------------------------
domain_colors <- c("Physical performance" = "#0F6E56",
                   "Cognitive & frailty"  = "#534AB7",
                   "Mitochondrial function" = "#BA7517")
q_colors_main <- c("Q2" = "#74C476", "Q3" = "#FD8D3C", "Q4" = "#D73027")
q_shapes_v5   <- c("Q2" = 21, "Q3" = 24, "Q4" = 23)
q_colors_scatter <- c("1" = "#4575B4", "2" = "#74C476", "3" = "#FD8D3C", "4" = "#D73027")

domain_label_df <- data.frame(
  y = c(1.5, 3.5, 6.5),
  label = c("Mitochondrial\nfunction", "Cognitive &\nfrailty", "Physical\nperformance"),
  color = unname(domain_colors[c("Mitochondrial function", "Cognitive & frailty",
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
  geom_text(data = domain_label_df, aes(x = 0.82, y = y, label = label, color = I(color)),
            hjust = 1, size = 2.2, lineheight = 0.85, fontface = "bold", inherit.aes = FALSE) +
  scale_color_manual(values = domain_colors, guide = "none") +
  shared_shape_b + shared_fill_b +
  guides(shape = guide_legend(nrow = 1, override.aes = list(size = 2.5)),
         fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), limits = c(-0.9, 0.9),
                     labels = function(x) sprintf("%.1f", x)) +
  labs(title = "b", x = "Adjusted difference vs. Q1 (outcome SD units)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
        axis.text.y = element_text(size = 8.5),
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


# =============================================================================
# FIGURE 3: ORGAN SPECIFIC ENTROPY AND OUTCOMES & SUPPLEMENTARY FIGURE 3
# =============================================================================
library(tidyverse)
library(survival)

# ── 0. Common settings ──────────────────────────────────────────────
organs <- c("blood","bone","brain","heart",
            "kidney","liver","lung",
            "lymphoid","skeletal")

organ_order <- c(
  "Full\nProteome",
  "Blood\nVessel","Bone","Brain","Heart",
  "Kidney","Liver","Lung","Lymphoid",
  "Skeletal\nMuscle"
)

outcome_order_fig <- c(
  "VO\u2082 peak (mL/min)",
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
  VO2    = "VO\u2082 peak (mL/min)",
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

dir.create("results", showWarnings=FALSE)

# ── 1. BOX  ─────────────────────
organ_box_raw <- readr::read_csv(
  "09b_entropy_measures_organ_box.csv",
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
            "VO2 peak (mL/min)" =
              "VO\u2082 peak (mL/min)"),
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
        "VO2 peak (mL/min)" =
          "VO\u2082 peak (mL/min)")
    ),
  results_cont_df |>
    dplyr::filter(organ != "full") |>
    dplyr::mutate(
      organ_label   = recode_organ(organ),
      outcome_label = dplyr::recode(
        as.character(outcome_label),
        "VO2 peak (mL/min)" =
          "VO\u2082 peak (mL/min)")
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
        "VO2 peak (mL/min)" =
          "VO\u2082 peak (mL/min)")
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
    axis.text.y      = ggplot2::element_text(
      size=8.5, color="black"),
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
  "results/Fig_organ_dot_Q2Q3Q4_BOX.pdf",
  plot=p_dot_q234, width=220, height=140,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  "results/Fig_organ_dot_Q2Q3Q4_BOX.png",
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
    axis.text.y      = ggplot2::element_text(
      size=8.5, color="black"),
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
  "results/Fig_organ_dot_continuous_BOX.pdf",
  plot=p_dot_cont, width=200, height=130,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  "results/Fig_organ_dot_continuous_BOX.png",
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
          "results/Table_Fig3b_organ_outcome_continuous.csv",
          row.names = FALSE)
cat("Saved: Table_Fig3b_organ_outcome_continuous.csv\n")


# =============================================================================
# 3b. Adjusted continuous-DISCO DEA (limma)  ->  Supplementary Table S2
#   Model : aptamer_z ~ DISCO_z + age + sex + site   (age/sex/site adjusted)
#   beta  : per 1-SD DISCO, in aptamer SD units (standardized coefficient)
#   Output: results/DEA_lm_BOX_DISCO_adj.csv
#   Prereqs (Sections 2-3): prot_df, meta_df, mapping_tbl, anno_tbl
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

# --- annotate (reuse Section 3 mapping/annotation) ---
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
  "results/DEA_lm_BOX_DISCO_adj.csv"
)
cat("Saved: DEA_lm_BOX_DISCO_adj.csv\n")

# --- sanity check against the manuscript thresholds ---
n_up <- sum(res_adj_anno$FDR < 0.001 & res_adj_anno$beta >  0.2, na.rm = TRUE)
n_dn <- sum(res_adj_anno$FDR < 0.001 & res_adj_anno$beta < -0.2, na.rm = TRUE)
cat(sprintf("|beta|>0.2 & FDR<0.001  ->  UP: %d | DOWN: %d\n", n_up, n_dn))

# =============================================================
# Figure 4. Volcano plot and bubble plot
# =============================================================

library(tidyverse)
library(ggrepel)
library(patchwork)

col_up      <- "#C0392B"
col_dn      <- "#2471A3"
col_up_dark <- "#7B241C"
col_dn_dark <- "#1A5276"
col_ns      <- "grey80"


# =============================================================
# 1. DEA data load
# =============================================================
dea_df <- readr::read_csv(
  "results/DEA_lm_BOX_DISCO_adj.csv",
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
# 3c. Hallmark GSEA on continuous DISCO  ->  Supplementary Table S3
#   Ranking : limma moderated t-statistic for the DISCO term (from 3b / S2),
#             collapsed to gene level (strongest aptamer = max |t| per SYMBOL).
#   Sets    : 50 MSigDB Hallmark gene sets (msigdbr, Homo sapiens, "H").
#   Method  : clusterProfiler::GSEA (pre-ranked), BH FDR across Hallmark sets.
#   Output  : results/GSEA_Hallmark_DISCO.csv   (all 50 pathways retained)
#   Prereqs : res_adj_anno from the S2 block (needs columns SYMBOL and t).
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

readr::write_csv(gsea_out, "results/GSEA_Hallmark_DISCO.csv")
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
  "results/GSEA_Hallmark_DISCO.csv",
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
  "results/Fig_GSEA_Hallmark_DISCO.pdf",
  plot=p_gsea, width=200, height=180,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  "results/Fig_GSEA_Hallmark_DISCO.png",
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
  "results/Fig_volcano_GSEA_combined.pdf",
  plot  = fig_combined,
  width = 380, height = 175,
  units = "mm", device = cairo_pdf)
ggplot2::ggsave(
  "results/Fig_volcano_GSEA_combined.png",
  plot  = fig_combined,
  width = 380, height = 175,
  units = "mm", dpi   = 300)
cat("Saved: Fig_volcano_GSEA_combined\n")



# =============================================================
# Supplementary Figure 4. Organ-specific Entropy and hallmark pathways
# =============================================================

library(tidyverse)
library(scales)

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
  "results/DEA_lm_BOX_DISCO_adj.csv",
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
  "results/Table_organ_protein_beta.csv",
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
  "results/Supp_Fig_organ_protein_heatmap.pdf",
  plot=p_heat_prot, width=190, height=250,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  "results/Supp_Fig_organ_protein_heatmap.png",
  plot=p_heat_prot, width=190, height=250,
  units="mm", dpi=300)
cat("Saved: Supp_Fig_organ_protein_heatmap\n")


# =============================================================
# SUPP FIG B: Organ × Pathway GSEA Heatmap
# =============================================================

# ── B-1. data load ─────────────────────────────────────────
gsea_full <- readr::read_csv(
  "results/GSEA_Hallmark_DISCO.csv",
  show_col_types=FALSE
)

gsea_organ <- readr::read_csv(
  "results/GSEA_organ_Hallmark.csv",
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
  "results/Supp_Fig_GSEA_organ_heatmap.pdf",
  plot=p_heat_gsea, width=220, height=180,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  "results/Supp_Fig_GSEA_organ_heatmap.png",
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
  "results/Supp_Fig_organ_combined.pdf",
  plot  = supp_combined,
  width = 420, height = 260,
  units = "mm", device = cairo_pdf)
ggplot2::ggsave(
  "results/Supp_Fig_organ_combined.png",
  plot  = supp_combined,
  width = 420, height = 260,
  units = "mm", dpi   = 300)
cat("Saved: Supp_Fig_organ_combined\n")



# =============================================================================
# Figure 5 KM curve Forest plot & Supplementary Figure 5
# =============================================================================

library(tidyverse)
library(survival)
library(cowplot)
library(scales)

dir.create("results", showWarnings=FALSE)

# ── 0. Common settings ──────────────────────────────────────────────
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
  "09b_entropy_measures_organ_box.csv",
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
    "results/Table_hosp_HR_BOX.csv")
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
    # HR main
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
    # HR<1 main
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
  "results/Fig_hosp_KM_heatmap_BOX.pdf",
  plot=fig_main, width=260, height=170,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  "results/Fig_hosp_KM_heatmap_BOX.png",
  plot=fig_main, width=260, height=170,
  units="mm", dpi=300)
cat("Saved: Fig_hosp_KM_heatmap_BOX\n")

# =============================================================================
# 8. Sensitivity HR (Base+GDF15, Base+CysC)
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

# shared legend
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
  "results/Supp_Fig_hosp_HR_sensitivity_BOX.pdf",
  plot=fig_sens, width=260, height=170,
  units="mm", device=cairo_pdf)
ggplot2::ggsave(
  "results/Supp_Fig_hosp_HR_sensitivity_BOX.png",
  plot=fig_sens, width=260, height=170,
  units="mm", dpi=300)
cat("Saved: Supp_Fig_hosp_HR_sensitivity_BOX\n")

cat("\n=== done ===\n")
cat("output files:\n")
cat("  Fig_hosp_KM_heatmap_BOX.png\n")
cat("  Supp_Fig_hosp_HR_sensitivity_BOX.png\n")
cat("  Table_hosp_HR_BOX.csv\n")

# =============================================================
# Follow-up duration — Reverse KM
# =============================================================

library(survival)

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
# DISCO (weighted) and DISCO (unweighted) comparison - Supplementary Figure 2
# =============================================================
library(readr); library(dplyr)

# --- unweighted:  _UNW  ---
unw <- read_csv("05z_entropy_measures_unweighted.csv", show_col_types = FALSE) |>
  distinct(ID, .keep_all = TRUE) |>
  rename_with(~ paste0(.x, "_UNW"), .cols = -ID)
# ID, BOX_DM_JRPCA_UNW, BOX_DM_JRSHRINK_UNW, BOX_DISCO_JR_UNW, ...

# --- weighted DISCO  ---
w <- read_csv("05_entropy_measures.csv", show_col_types = FALSE) |>
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
# SENSITIVITY: adding M1ANMEDS_LOG (log medication count) to the BASE model - Supplementary Figure 6
#   Figure a  — phenotype-association forest, SAME format as Figure 2 (panel b)
#               weighted-DISCO quartiles, base + M1ANMEDS_LOG
#   Figure b  — hospitalization Cox HR heatmap, SAME format as the Figure 5
#               sensitivity heatmap, single panel: base + M1ANMEDS_LOG
# -----------------------------------------------------------------------------
# PREREQUISITES (run the Figure 2 and Figure 5 scripts first, in this session):
#   From Fig 2 : get_std_beta(), make_forest_b(), save_fig(), outcome_labels,
#                outcome_order, domain_colors, q_colors_main, q_shapes_v5,
#                domain_label_df, and the data frame `m`.
#   From Fig 5 : get_hr_all_q_box(), extract_hr_full(), tidy_hr(),
#                add_heat_geoms(), heat_fill_scale, heat_color_scale, heat_theme,
#                organ_order_heat, y_faces_heat, organs, base_cov, `analysis_df`.
# =============================================================================

library(dplyr); library(ggplot2); library(cowplot); library(readr); library(survival)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ---- 0. Load M1ANMEDS_LOG and merge into both analysis frames ---------------
meds <- readr::read_csv("02_box_transformed_data_with_pheno.csv", show_col_types = FALSE) |>
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

# =============================================================================
# FIGURE a — phenotype forest (Fig 2 panel-b format), base + M1ANMEDS_LOG
# =============================================================================
outcomes_a <- list(
  VO2 = "TTPKVO2U", Speed = "NF400MPACE", Power = "LEPEAKPWR2", Steps = "ACSCFUAV",
  Digit = "DSCORR", VF = "FT0V2FN", OxPhos = "REMOXPHOS", CCR = "TTSSPKVO2"
)

# Build std with the SAME domain names used by domain_colors (correct colors).
std_med_a <- lapply(names(outcomes_a), function(nm) {
  d <- get_std_beta(outcomes_a[[nm]], m, extra_str = "M1ANMEDS_LOG")   # <- extra covariate
  d$outcome <- nm; d
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    quartile_label = factor(paste0("Q", quartile), levels = c("Q4", "Q3", "Q2")),
    domain = dplyr::case_when(
      outcome %in% c("VO2", "Speed", "Power", "Steps") ~ "Physical performance",
      outcome %in% c("Digit", "VF")                    ~ "Cognitive & frailty",
      TRUE                                             ~ "Mitochondrial function"),
    outcome_label = outcome_labels[outcome],
    outcome_label = factor(outcome_label, levels = rev(outcome_order))
  )

cat("\n=== Figure a: std beta Q4 vs Q1 (base + M1ANMEDS_LOG) ===\n")
std_med_a |> dplyr::filter(quartile == 4) |>
  dplyr::select(outcome, beta, beta_lo, beta_hi, p_value) |>
  dplyr::mutate(dplyr::across(c(beta, beta_lo, beta_hi), round, 3),
                p_value = round(p_value, 4)) |>
  as.data.frame() |> print()

fig_a <- make_forest_b(std_med_a, "a",
                       "Adjusted: age, sex, site, height, weight + medication count (log)")

save_fig(fig_a, "results/FigS_pheno_forest_meds.png", w = 130, h = 140)

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
  dplyr::mutate(dplyr::across(c(HR, lo, hi), round, 2)) |>
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

ggplot2::ggsave("results/FigS_hosp_HR_meds.pdf", fig_b,
                width = 150, height = 170, units = "mm", device = cairo_pdf)
ggplot2::ggsave("results/FigS_hosp_HR_meds.png", fig_b,
                width = 150, height = 170, units = "mm", dpi = 300)

# CSV of the HR table
heat_df_med |>
  dplyr::select(organ_label, quartile, HR, lo, hi, p, sig_label) |>
  dplyr::arrange(organ_label, quartile) |>
  readr::write_csv("results/Table_hosp_HR_meds.csv")

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

ggplot2::ggsave("results/FigS_meds_combined.pdf", fig_ab,
                width = 300, height = 175, units = "mm", device = cairo_pdf)
ggplot2::ggsave("results/FigS_meds_combined.png", fig_ab,
                width = 300, height = 175, units = "mm", dpi = 300)

cat("\nSaved: FigS_pheno_forest_meds.png | FigS_hosp_HR_meds.png/.pdf | ",
    "FigS_meds_combined.png/.pdf | Table_hosp_HR_meds.csv\n")


