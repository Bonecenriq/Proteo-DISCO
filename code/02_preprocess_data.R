# ==============================================================================
# STUDY       : SOMMA
# PROGRAM     : 02_preprocess_data.R
#
# PURPOSE     : Pre-process phenotypic and SomaLogic proteomic data for the
#               SOMMA SomaLogic entropy analyses.
#
# INPUTS      : - SOMMA Parent phenotypic data
#               - SOMMA Jr phenotypic data
#               - SomaLogic ADAT file
#               - SomaLogic sample annotation file
#               - Hospitalization endpoint data
#               - Death endpoint data
#               - SomaLogic aptamer QC metadata
#
# OUTPUTS     : - Log2-transformed and parent-standardized proteomic data
#               - Regularized Box-Cox-transformed and parent-standardized
#                 proteomic data
#               - Inverse-normal-transformed proteomic data
#               - Unscaled log2-transformed proteomic data
#               - Untransformed proteomic data
#               - Phenotypic and outcome variables merged with each
#                 proteomic dataset
#
# DATA ACCESS : Individual-level SOMMA data are not included in this
#               repository. Access to the underlying study data is required
#               to reproduce the analyses.
#
# DEPENDENCIES:
#               - R
#               - haven
#               - data.table
#               - dplyr
#               - plyr
#               - stringr
#               - readxl
#               - SomaDataIO (or package providing read_adat)
#               - Additional analysis functions and objects defined in
#                 99_analysis_setup.R
#
# VERSION     : 1.0
# DATE        : 2025-11-24
# AUTHOR      : Tyler Mansfield
#
# REVISION HISTORY
# 11/24/25    Creation (TM)
# 01/30/26    Removed GDF-15 and Cystatin C from entropy dataset and instead
#             treated them as outcomes.
# 02/16/26    Added non-elective hospitalization event status and follow-up time.
#
# NOTES       : Local file paths are specified separately from the analysis
#               code and should be configured for the user's environment.
# ==============================================================================

source("code/99_analysis_setup.R")




### Pre-process data ----------------------------------------------------------
# Read in phenotypic traits
pheno_data_parent <- haven::read_sas(
  file.path(PATH_DATA, "prot_subset.sas7bdat")
) %>%
  mutate(
    SEX = factor(EEFEMALE, 0:1, c("Men", "Women")),
    SITE = factor(SITE),
    COHORT = "Parent"
  ) %>%
  as.data.table()

pheno_data_jr <- haven::read_sas(
  file.path(PATH_DATA, "protjr_subset.sas7bdat")
) %>%
  mutate(
    SEX = factor(EEFEMALE, 0:1, c("Men", "Women")),
    SITE = factor(SITE),
    COHORT = "Jr"
  ) %>%
  as.data.table()

pheno_data <- plyr::rbind.fill(pheno_data_parent, pheno_data_jr) %>%
  as.data.table()

# ==============================================================================
# SOMALOGIC DATA
# ==============================================================================

my_adat <- read_adat(
  file.path(PATH_DATA, "SO000_2025_05.adat")
)
# Subset to those in the proteomics dataset
somma_subset <- my_adat %>% 
  dplyr::filter(str_detect(SampleId, "^SOMMA")) %>% # Only keep SOMMA Participants
  pull(SampleId) %>%
  str_remove_all("SOMMA")


pheno_data <- pheno_data[ID %in% somma_subset]

# Merge with Sample meta data
sample_meta <- read_excel(file.path(PATH_DATA, "SO000_Annotations.xlsx"),
                          sheet = "Sample Annotations") %>%
  as.data.table() %>%
  .[, ID := as.numeric(str_remove(SampleId, "SOMMA"))]

pheno_data <- merge(pheno_data, sample_meta, by = "ID")
rm(sample_meta)

# Make some standardized covariates and outcomes
pheno_data[COHORT == "Parent", D1AGE2_STD := scale(D1AGE2)]
pheno_data[COHORT == "Parent", D1AGE2_5YR := (D1AGE2 - 70) / 5]
pheno_data[COHORT == "Parent", HWWGT_STD := scale(HWWGT)]
pheno_data[COHORT == "Parent", HWHGT_STD := scale(HWHGT)]
pheno_data[COHORT == "Parent", WGTHGT_IXN := HWWGT_STD*HWHGT_STD]
pheno_data[COHORT == "Parent", MRMAFATTOT_STD := scale(MRMAFATTOT)]
pheno_data[COHORT == "Parent", FATHGT_IXN := MRMAFATTOT_STD*HWHGT_STD]
pheno_data[COHORT == "Parent", M1ANMEDS_LOG := log10(1 + M1ANMEDS)]

# Standardized outcomes just in parent
for (outcome in all_outcomes)
{
  if (!(outcome %in% names(pheno_data))) {next}
  pheno_data[COHORT == "Parent", paste0(outcome, "_PARENT_STD") := scale(get(outcome))]
}

### ADD IN HOSPITALIZATIONS ---------------------------------------------------

hosp_data <- haven::read_sas(file.path(PATH_DATA, "ha999_2026_01_W1.sas7bdat")) %>%
  as.data.table()

hosp_data_long <- melt(
  hosp_data,
  id.vars = c("ID", "SITE", "HAHOSP", "HANUM"),
  measure.vars = patterns(
    H0TTEH   = "^H0TTEH[1-6]$",
    H0TTEUH = "^H0TTEUH[1-6]$",
    H0LOSTAY = "^H0LOSTAY[1-6]$",
    H0ELECT   = "^H0ELECT[1-6]$",
    H0EMERG = "^H0EMERG[1-6]$",
    H0ICU = "^H0ICU[1-6]$",
    H0ICUV = "^H0ICUV[1-6]$",
    H0CHF = "^H0CHF[1-6]$",
    H0DPCD= "^H0DPCD[1-6]$",
    H0PCD = "^H0PCD[1-6]$"
  ),
  variable.name = "event_num"
) %>%
  .[!is.na(H0TTEH)]

# Subset to non-elective hospitalizations
hosp_data_long <- hosp_data_long[H0EMERG == 1]

# Merge in event status and followup time
death_data <- haven::read_sas(file.path(PATH_DATA, "ef999_2026_01_w1.sas7bdat")) %>%
  as.data.table()
pheno_data <- pheno_data[, HAEMERG := ifelse(ID %in% hosp_data_long$ID, 1, 0)] %>%
  merge(hosp_data[, .(ID, HAFUTBLU)], by = "ID", all.x = T) %>%
  merge(death_data[, .(ID, EF48STAT, DADEAD)], by = "ID", all.x = T)


### Proteomic Data Transformation --------------------------------------------------
# Remove GDF-15 (seq.4374.45) and Cystatin C (seq.2609.59) from entropy dataset and instead put them as outcomes
pheno_data_additions <- my_adat %>% 
  dplyr::filter(SampleType == "Sample" & SampleId %in% pheno_data$SampleId) %>%
  dplyr::select(c("SampleId", "seq.4374.45", "seq.2609.59")) %>%
  dplyr::mutate(ID = as.numeric(str_remove_all(SampleId, "SOMMA")),
                SOGDF15 = log10(seq.4374.45),
                # Cystatin-C has some outliers: Winsorize
                SOCYSC =  Winsorize(log10(seq.2609.59), val = quantile(log10(seq.2609.59), probs = c(0.01, 0.99)))) %>%
  dplyr::mutate(SOGDF15SD = scale(SOGDF15), SOCYSCSD = scale(SOCYSC)) %>% 
  dplyr::select(ID, SOGDF15, SOGDF15SD, SOCYSC, SOCYSCSD)

pheno_data <- merge(pheno_data, 
                    pheno_data_additions, 
                    by = "ID",
                    all.x = T) %>%
  as.data.table() %>%
  .[COHORT == "Parent", SOGDF15SD_PARENT_STD := scale(SOGDF15)] %>%
  .[COHORT == "Parent", SOCYSCSD_PARENT_STD := scale(SOCYSC)]
rm(pheno_data_additions)
  
# We will use the log2 transformation on all aptamers
sample_data <- my_adat %>% 
  dplyr::filter(SampleType == "Sample" & SampleId %in% pheno_data$SampleId) %>%
  dplyr::select(starts_with("seq")) %>%
  dplyr::select(-c("seq.4374.45", "seq.2609.59")) # Remove GDF-15 (seq.4374.45) and Cystatin C (seq.2609.59)

sample_data_id <- my_adat %>% 
  dplyr::filter(SampleType == "Sample" & SampleId %in% pheno_data$SampleId) %>%
  dplyr::pull("SampleId")


# ============================================================================
# LOG2 TRANSFORM
# ============================================================================
log_transformed_data <- copy(sample_data) %>%
  log2() %>%
  as.data.table()

# Add IDs
log_transformed_data[, ID := str_remove(sample_data_id, "SOMMA")]

# Merge cohort information
log_transformed_data <- merge(
  log_transformed_data,
  pheno_data[, .(ID = as.character(ID), COHORT)],
  by = "ID"
)

# STANDARDIZE USING PARENT COHORT MEAN/SD
protein_vars <- setdiff(names(log_transformed_data), c("ID", "COHORT"))

# Matrix form is MUCH faster
prot_mat <- as.matrix(log_transformed_data[, ..protein_vars])

# Parent rows
parent_ix <- log_transformed_data$COHORT == "Parent"

# Parent reference statistics
parent_means <- colMeans(prot_mat[parent_ix, ], na.rm = TRUE)

parent_sds <- apply(prot_mat[parent_ix, ],
                    2,
                    sd,
                    na.rm = TRUE)

# Vectorized standardization
prot_mat <- sweep(prot_mat, 2, parent_means, "-")
prot_mat <- sweep(prot_mat, 2, parent_sds, "/")

# Put back
log_transformed_data[, (protein_vars) := as.data.table(prot_mat)]

# ============================================================================
# REGULARIZED BOX-COX TRANSFORMATION
# ============================================================================

box_transformed_data <- as.data.table(sample_data)

# Apply Box-Cox transform
num_cols <- names(box_transformed_data)

for (col in num_cols)
{
  box_transformed_data[[col]] <-
    boxcox_transform(box_transformed_data[[col]],
                     regularize = 0.5)
}

# Add IDs
box_transformed_data[, ID := str_remove(sample_data_id, "SOMMA")]

# Merge cohort information
box_transformed_data <- merge(
  box_transformed_data,
  pheno_data[, .(ID = as.character(ID), COHORT)],
  by = "ID"
)

# STANDARDIZE USING PARENT COHORT MEAN/SD
prot_mat <- as.matrix(box_transformed_data[, ..protein_vars])

# Parent rows
parent_ix <- box_transformed_data$COHORT == "Parent"

# Parent reference statistics
parent_means <- colMeans(prot_mat[parent_ix, ], na.rm = TRUE)

parent_sds <- apply(prot_mat[parent_ix, ],
                    2,
                    sd,
                    na.rm = TRUE)

# Vectorized standardization
prot_mat <- sweep(prot_mat, 2, parent_means, "-")
prot_mat <- sweep(prot_mat, 2, parent_sds, "/")

# Put back
box_transformed_data[, (protein_vars) := as.data.table(prot_mat)]

# ============================================================================
# INVERSE NORMAL TRANSFORMATION
# ============================================================================

invnorm_transformed_data <- as.data.table(sample_data)
num_cols <- names(invnorm_transformed_data)
for (col in num_cols) {
  invnorm_transformed_data[[col]] <- invnorm_transform(invnorm_transformed_data[[col]])
}
invnorm_transformed_data$ID <- str_remove(sample_data_id, "SOMMA")


#### Three other scales
# Log2 transformed: no SCALE
log_transformed_data_noscale <- copy(sample_data) %>%
  log2() %>%
  as.data.table()
log_transformed_data_noscale$ID <- str_remove(sample_data_id, "SOMMA")

no_transform_data <- copy(sample_data) %>%
  as.data.table()
no_transform_data$ID <- str_remove(sample_data_id, "SOMMA")

### Proteomic Data Aptamer Exclusions --------------------------------------------------
aptamer_include <- aptamer_meta %>%
  .[`Total CV Plasma` < 0.20] %>%
  .[`F-Statistic Plasma` > 2.43] %>%
  .[`Fraction below LoDB` < 0.05] %>%
  .[`Signal to Noise Ratio` > 2.5] %>%
  .[Shaky == 0] %>%
  .[!(AptName %in% c("seq.4374.45", "seq.2609.59"))]

log_transformed_data <- log_transformed_data[, mget(c("ID", aptamer_include$AptName))]
box_transformed_data <- box_transformed_data[, mget(c("ID", aptamer_include$AptName))]
invnorm_transformed_data <- invnorm_transformed_data[, mget(c("ID", aptamer_include$AptName))]
log_transformed_data_noscale <- log_transformed_data_noscale[, mget(c("ID", aptamer_include$AptName))]
no_transform_data <- no_transform_data[, mget(c("ID", aptamer_include$AptName))]
pheno_data$ID <- as.character(pheno_data$ID)



# Combine all into one large data.table
fwrite(merge(pheno_data, log_transformed_data, by = "ID"), 
       file.path(PATH_DATA, "02_log_transformed_data_with_pheno.csv"))
fwrite(merge(pheno_data, box_transformed_data, by = "ID"), 
       file.path(PATH_DATA, "02_box_transformed_data_with_pheno.csv"))
fwrite(merge(pheno_data, invnorm_transformed_data, by = "ID"), 
       file.path(PATH_DATA, "02_invnorm_transformed_data_with_pheno.csv"))
fwrite(merge(pheno_data, log_transformed_data_noscale, by = "ID"), 
       file.path(PATH_DATA, "02_log2noscale_transformed_data_with_pheno.csv"))
fwrite(merge(pheno_data, no_transform_data, by = "ID"), 
       file.path(PATH_DATA, "02_no_transform_data_with_pheno.csv"))