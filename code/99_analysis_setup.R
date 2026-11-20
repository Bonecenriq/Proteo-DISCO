# ==============================================================================
# STUDY       : SOMMA
# PROJECT     : SOMMA SomaLogic
# PROGRAM     : 99_analysis_setup.R
#
# PURPOSE     : Define project paths, package dependencies, global analysis
#               parameters, metadata, and reusable functions for the SOMMA
#               SomaLogic proteomic entropy analyses.
#
# DATA ACCESS : Individual-level SOMMA data are not included in this repository.
#               Reproduction of the analyses requires authorized access to the
#               underlying study data.
#
#               Users should configure PROJECT_DIR below to point to the local
#               root directory containing the data/, output/, figures/, and
#               code/ directories.
#
# DEPENDENCIES:
#               - R >= 4.0
#               - R packages listed in required_packages below
#
# VERSION     : 1.0
# DATE        : 2025-11-24
# AUTHOR      : Tyler Mansfield
#
# REVISION HISTORY
# 09/30/25    Creation
# 01/16/26    Removed inverse-normal transformation from primary analyses
# 01/30/26    Set Box-Cox transformation as the primary transformation
#
# NOTES       : This file is sourced by the analysis scripts and is not intended
#               to be run independently.
#
#               Local paths and access to individual-level data are intentionally
#               excluded from the repository.
# ==============================================================================

# ==============================================================================
# PATH CONFIGURATION
# ==============================================================================
# Set this path to the local directory containing the SOMMA analysis files.
#
# The directory should contain:
#   data/
#   output/
#   figures/
#   code/
#
# IMPORTANT:
# Individual-level study data are not included in the GitHub repository.
# Users must have authorized access to the underlying SOMMA data and configure
# this path to point to their local data directory.

PROJECT_DIR <- "/path/to/SOMMA_SomaLogic"

PATH_DATA   <- file.path(PROJECT_DIR, "data")
PATH_OUTPUT <- file.path(PROJECT_DIR, "output")
PATH_FIGURE <- file.path(PROJECT_DIR, "figures")
PATH_CODE   <- file.path(PROJECT_DIR, "code")

# Load custom analysis functions
source(file.path(PATH_CODE, "99_entropy_methods.R"))

# Check that required directories exist
required_dirs <- c(
  PATH_DATA,
  PATH_OUTPUT,
  PATH_FIGURE,
  PATH_CODE
)

missing_dirs <- required_dirs[!dir.exists(required_dirs)]

if (length(missing_dirs) > 0) {
  stop(
    "The following directories do not exist:\n",
    paste(missing_dirs, collapse = "\n"),
    "\n\nPlease update PROJECT_DIR in 99_analysis_setup.R."
  )
}

# ==============================================================================
# PACKAGES
# ==============================================================================

required_packages <- c(
  "tidyverse",
  "SomaDataIO",
  "openxlsx",
  "readxl",
  "corrplot",
  "ggcorrplot",
  "data.table",
  "assertthat",
  "ggpubr",
  "moments",
  "future",
  "future.apply",
  "ggrepel",
  "DescTools",
  "glmnet",
  "pheatmap",
  "MuMIn",
  "MASS",
  "ComplexHeatmap",
  "survival",
  "survminer",
  "timeROC",
  "fgsea",
  "emmeans",
  "limma",
  "GSVA"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following R packages are required but not installed:\n",
    paste(missing_packages, collapse = "\n"),
    "\n\nPlease install the missing packages before running the analysis."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))


### GLOBAL VARS  -------------------------------------------------
data_transform <- "box"
data_transforms <- c("log", "box", "invnorm")

entropy_methods <- c("DM_JRPCA", "DM_JRSHRINK", "DISCO_JR")
entropy_methods_labs <- c("M-dist (PCA)", "M-dist (Shrink)", "DISCO")

entropy_vars <- toupper(outer(data_transforms, entropy_methods, paste, sep = "_")) %>% t() %>% as.vector()
entropy_vars_labs_short <- outer(entropy_methods_labs, data_transforms, paste, sep = ", ") %>% as.vector()
entropy_vars_labs_short <- entropy_methods_labs
rm(entropy_methods_labs)

# ==============================================================================
# OUTCOME DEFINITIONS
# ==============================================================================

# Primary outcomes.
primary_outcomes <- c(
  "D1AGE2",
  "TTPKVO2U",
  "NF400MPACE",
  "LEPEAKPWR2"
)

primary_outcomes_labs_short <- c(
  "Age",
  "Absolute VO2 Peak",
  "Walking Speed",
  "Leg Power"
)

primary_outcomes_labs_short_with_units <- c(
  "Age, years",
  "Absolute VO2 Peak, mL/min",
  "Walking Speed, m/s",
  "Leg Power, watts"
)

primary_outcomes_labs_long <- c(
  "Age, years",
  "Absolute VO2 Peak, mL/min",
  "Walking Speed, m/s",
  "Leg Power, watts"
)


# Secondary outcomes.
sec_outcomes <- c(
  "ACSCFUAVK",
  "DSCORR",
  "REMOXPHOS",
  "TTSSPKVO2",
  "FT0V2FN",
  "SOGDF15SD",
  "SOCYSCSD"
)

sec_outcomes_labs_short <- c(
  "Daily Step Count",
  "Digit Symbol Substitution",
  "Max Oxphos",
  "Cost Capacity Ratio",
  "Vigor to Frailty",
  "GDF-15 (SomaLogic)",
  "Cystatin-C (SomaLogic)"
)

sec_outcomes_labs_short_with_units <- c(
  "Daily Step Count, thousand",
  "Digit Symbol Substitution, number correct",
  "Max Oxphos, pmol/(s*mg)",
  "Cost Capacity Ratio, percent",
  "Vigor to Frailty, 0-12",
  "GDF-15 (SomaLogic), SD",
  "Cystatin-C (SomaLogic), SD"
)

sec_outcomes_labs_long <- c(
  "Daily Step Count, thousand steps per day",
  "Digit Symbol Substitution, number correct",
  "Maximal Complex I&II Supported Oxphos, pmol/(s*mg)",
  "Cost Capacity Ratio of VO2 from Steady Speed and Peak, percent",
  "Vigor to Frailty, 0-12",
  "GDF-15 from SomaLogic, SD of log(RFU)",
  "Cystatin-C from SomaLogic, SD of log(RFU)"
)


# All outcomes.
all_outcomes <- c(
  primary_outcomes,
  sec_outcomes
)

all_outcomes_labs_short <- c(
  primary_outcomes_labs_short,
  sec_outcomes_labs_short
)

all_outcomes_labs_short_with_units <- c(
  primary_outcomes_labs_short_with_units,
  sec_outcomes_labs_short_with_units
)

all_outcomes_labs_long <- c(
  primary_outcomes_labs_long,
  sec_outcomes_labs_long
)


### Aptamer metadata --------------------
aptamer_meta <- read_excel(file.path(PATH_DATA, "SO000_Annotations.xlsx"), 
                           sheet = "Analyte Annotations") %>%
  as.data.table()


# ==============================================================================
# TRANSFORMATION FUNCTIONS
# ==============================================================================

# Regularized Box-Cox transformation
#
# The Box-Cox transformation parameter (lambda) is estimated by maximizing the
# profile log-likelihood over a grid from -10 to 10 in increments of 0.5.
#
# Before estimating lambda, values are divided by their median to improve
# numerical stability. The estimated lambda is then multiplied by the
# regularization parameter.
#
# regularize = 1:
#   No regularization; use the maximum-likelihood estimate directly.
#
# regularize < 1:
#   Shrinks the estimated lambda toward zero (the log transformation).
#
# The transformation requires strictly positive input values.

boxcox_transform <- function(X, regularize = 1)
{
  assertthat::assert_that(
    is.numeric(X),
    all(X > 0, na.rm = TRUE)
  )
  
  assertthat::assert_that(
    is.numeric(regularize),
    regularize >= 0,
    regularize <= 1
  )
  
  # Scale values for numerical stability.
  X <- X / median(X, na.rm = TRUE)
  
  n <- length(X)
  
  # Candidate Box-Cox parameters.
  lambda_grid <- seq(-10, 10, by = 0.5)
  
  # Profile log-likelihood for each candidate lambda.
  loglike <- sapply(
    lambda_grid,
    function(lambda) {
      
      Z <- if (lambda == 0) {
        log(X)
      } else {
        (X^lambda - 1) / lambda
      }
      
      - (n / 2) * log(mean((Z - mean(Z))^2)) +
        (lambda - 1) * sum(log(X))
    }
  )
  
  # Estimate lambda and apply regularization toward zero.
  lambda_hat <- lambda_grid[which.max(loglike)] * regularize
  
  # Apply the regularized Box-Cox transformation.
  if (lambda_hat == 0) {
    log(X)
  } else {
    (X^lambda_hat - 1) / lambda_hat
  }
}


# Inverse-normal transformation
#
# Applies a rank-based inverse-normal transformation using the Blom
# approximation. Ties are assigned their average rank.

invnorm_transform <- function(x)
{
  r <- rank(
    x,
    ties.method = "average",
    na.last = "keep"
  )
  
  n <- sum(!is.na(x))
  
  p <- (r - 3 / 8) / (n + 1 / 4)
  
  qnorm(p)
}