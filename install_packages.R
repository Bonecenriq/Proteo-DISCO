# ==============================================================================
# install_packages.R
# Installs all R packages required by the Proteo-DISCO pipeline.
# Usage (from the project root):  Rscript install_packages.R
# ==============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- CRAN packages -----------------------------------------------------------
cran_pkgs <- c(
  # data handling
  "tidyverse", "data.table", "assertthat", "haven", "plyr", "readxl",
  "openxlsx", "writexl", "DescTools", "moments", "future", "future.apply",
  # SomaScan I/O
  "SomaDataIO",
  # statistics
  "MASS", "glmnet", "MuMIn", "emmeans", "survival", "survminer", "timeROC",
  "sandwich", "lmtest", "FSA", "msigdbr",
  # visualization
  "ggpubr", "ggrepel", "ggtext", "cowplot", "patchwork", "scales",
  "corrplot", "ggcorrplot", "pheatmap", "RColorBrewer"
)

# ---- Bioconductor packages ---------------------------------------------------
bioc_pkgs <- c(
  "limma", "fgsea", "GSVA", "clusterProfiler", "ComplexHeatmap",
  "AnnotationDbi", "SomaScan.db"
)

install_if_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) install.packages(missing)
}

install_if_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) BiocManager::install(missing, update = FALSE, ask = FALSE)
}

install_if_missing_cran(cran_pkgs)
install_if_missing_bioc(bioc_pkgs)

# ---- Check -------------------------------------------------------------------
all_pkgs  <- c(cran_pkgs, bioc_pkgs)
still_missing <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(still_missing) == 0) {
  message("All ", length(all_pkgs), " packages are installed.")
} else {
  warning("These packages could not be installed: ",
          paste(still_missing, collapse = ", "))
}
