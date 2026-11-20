# ==============================================================================
# STUDY       : SOMMA
# PROJECT     : SOMMA SomaLogic
# PROGRAM     : 99_entropy_methods.R
#
# PURPOSE     : Reusable functions for calculating proteomic entropy measures
#               used in the SOMMA SomaLogic analyses.
#
# METHODS     : - DISCO (Distance of Covariance)
#               - Optimized DISCO implementation using rank-one covariance
#                 updates
#
# DATA ACCESS : These functions operate on individual-level study data that are
#               not included in this repository. Authorized access to the
#               underlying SOMMA data is required to reproduce the analyses.
#
# DEPENDENCIES:
#               - data.table
#               - assertthat
#
# VERSION     : 1.0
# DATE        : 2025-12-01
# AUTHOR      : Tyler Mansfield
#
# REVISION HISTORY
# 12/01/25    Creation
#
# NOTES       : This file contains reusable functions and is sourced by the
#               analysis scripts. Package dependencies are loaded centrally
#               in 99_analysis_setup.R.
#
#               The DISCO_fast() implementation is mathematically equivalent
#               to DISCO() for the supported Pearson-correlation implementation
#               and is intended for substantially faster computation on large
#               proteomic datasets.
# ==============================================================================

# Package dependencies are loaded centrally in 99_analysis_setup.R.

# ==============================================================================
# DISCO
# ==============================================================================

#' Calculate DISCO (Distance of Covariance)
#'
#' Calculates the Distance of Covariance (DISCO) for each target individual
#' relative to a reference population. DISCO quantifies the extent to which
#' adding an individual's biomarker profile to the reference population
#' perturbs the correlation structure of the biomarker network.
#'
#' @param ref_dt A data.frame or data.table containing biomarker measurements
#'   for the reference population.
#' @param target_dt A data.frame or data.table containing biomarker measurements
#'   for target individuals.
#' @param cols Character vector specifying the biomarker columns used to
#'   construct the correlation network.
#' @param age_assoc Named numeric vector containing the age association for
#'   each biomarker. These associations are used to calculate pairwise
#'   biomarker weights.
#' @param method Correlation method passed to stats::cor(). Default is
#'   "pearson".
#' @param log10_output Logical indicating whether DISCO values should be
#'   log10-transformed. Default is TRUE.
#'
#' @return A numeric vector containing one DISCO value for each row of
#'   target_dt.

DISCO <- function(ref_dt, target_dt, cols, age_assoc, 
                     method = "pearson", log10_output = TRUE) {
  assertthat::assert_that(all(cols %in% names(ref_dt)))
  assertthat::assert_that(all(cols %in% names(target_dt)))
  assertthat::assert_that(all(cols %in% names(age_assoc)))
  
  if (!data.table::is.data.table(ref_dt)) ref_dt <- data.table::as.data.table(ref_dt)
  if (!data.table::is.data.table(target_dt)) target_dt <- data.table::as.data.table(target_dt)
  
  X_ref <- as.matrix(ref_dt[, ..cols])
  X_tgt <- as.matrix(target_dt[, ..cols])
  
  n_ref <- nrow(X_ref)
  m <- ncol(X_ref)
  
  # Reference correlation matrix (fixed baseline)
  R_ref <- stats::cor(X_ref, method = method)
  
  # Weight matrix (paper uses a product-based weighting tied to age associations)
  w <- abs(outer(age_assoc[cols], age_assoc[cols]))
  diag(w) <- 0
  
  # We'll sum only unique pairs (i < j), so normalize weights over i<j
  w_lower <- w[lower.tri(w)]
  w_lower <- w_lower / sum(w_lower)
  
  disco_vals <- numeric(nrow(X_tgt))
  
  for (i in seq_len(nrow(X_tgt))) {
    X_pert <- rbind(X_ref, X_tgt[i, , drop = FALSE])
    R_pert <- stats::cor(X_pert, method = method)
    
    diff_sq <- (R_pert - R_ref)^2
    disco_raw <- (n_ref^2) * sum(diff_sq[lower.tri(diff_sq)] * w_lower)
    
    disco_vals[i] <- disco_raw
  }
  
  if (log10_output) disco_vals <- log10(disco_vals)
  disco_vals
}


#' Calculate DISCO using an optimized Pearson-correlation implementation
#'
#' This function calculates the Distance of Covariance (DISCO) for each target
#' individual relative to a reference population. DISCO quantifies the extent
#' to which adding an individual's biomarker profile to the reference
#' population perturbs the correlation structure of the biomarker network.
#'
#' This implementation is optimized for large datasets. Rather than
#' recalculating the full correlation matrix for each target individual, it
#' uses a rank-one update of the centered cross-product matrix to calculate
#' the perturbed covariance and correlation matrices.
#'
#' The optimized implementation currently supports Pearson correlation only.
#' It assumes that the biomarker matrix contains finite, non-missing values
#' and that all reference biomarkers have non-zero variance.
#'
#' @param ref_dt A data.frame or data.table containing biomarker measurements
#'   for the reference population used to define the baseline correlation
#'   structure.
#' @param target_dt A data.frame or data.table containing biomarker measurements
#'   for the individuals whose DISCO values are to be calculated.
#' @param cols Character vector specifying the biomarker columns used to
#'   construct the correlation network.
#' @param age_assoc Named numeric vector containing the age association for
#'   each biomarker. Pairwise weights are calculated as the absolute product
#'   of the age associations.
#' @param method Character string specifying the correlation method. The
#'   optimized implementation currently supports "pearson" only.
#' @param log10_output Logical indicating whether DISCO values should be
#'   log10-transformed. Default is TRUE.
#' @param show_progress Logical indicating whether a progress bar should be
#'   displayed. Default is TRUE.
#'
#' @return A numeric vector with one DISCO value for each row of target_dt.
#'   Larger values indicate greater perturbation of the reference correlation
#'   structure.

DISCO_fast <- function(ref_dt,
                       target_dt,
                       cols,
                       age_assoc,
                       method = "pearson",
                       log10_output = TRUE,
                       show_progress = TRUE) {
  
  assertthat::assert_that(all(cols %in% names(ref_dt)))
  assertthat::assert_that(all(cols %in% names(target_dt)))
  assertthat::assert_that(all(cols %in% names(age_assoc)))
  
  if (!data.table::is.data.table(ref_dt)) ref_dt <- data.table::as.data.table(ref_dt)
  if (!data.table::is.data.table(target_dt)) target_dt <- data.table::as.data.table(target_dt)
  
  # Extract matrices
  X_ref <- as.matrix(ref_dt[, ..cols])
  X_tgt <- as.matrix(target_dt[, ..cols])
  
  n_ref <- nrow(X_ref)
  p <- ncol(X_ref)
  if (n_ref < 3) stop("Reference group must have at least 3 samples.")
  if (method != "pearson") stop("This optimized version currently supports method = 'pearson' only.")
  
  # --- Reference mean and scatter (M2) ---
  # mu_ref: p-vector
  mu_ref <- colMeans(X_ref)
  
  # Centered reference data
  Xc_ref <- sweep(X_ref, 2, mu_ref, FUN = "-")
  
  # M2_ref = sum_{i=1..n} (x_i - mu)(x_i - mu)^T  (p x p)
  # We only keep lower triangle + diagonal to avoid big allocations later
  # Crossprod gives the full matrix once here (acceptable once)
  M2_ref <- crossprod(Xc_ref)  # p x p
  
  # Reference covariance and correlation
  cov_ref <- M2_ref / (n_ref - 1)
  var_ref <- diag(cov_ref)
  sd_ref <- sqrt(var_ref)
  
  # Protect against zero-variance features
  if (any(sd_ref == 0)) {
    zero_vars <- cols[sd_ref == 0]
    stop("Some reference variables have zero variance: ", paste(zero_vars, collapse = ", "),
         ". Remove or filter these features.")
  }
  
  # Lower-triangle index vectors (fixed)
  lower_idx <- which(lower.tri(M2_ref), arr.ind = TRUE)
  li <- lower_idx[, 1]
  lj <- lower_idx[, 2]
  
  # Reference correlation lower triangle as vector
  cov_ref_lower <- cov_ref[lower.tri(cov_ref)]
  denom_ref_lower <- sd_ref[li] * sd_ref[lj]
  R_ref_lower <- cov_ref_lower / denom_ref_lower
  
  # --- Weights (lower triangle), normalized consistently with the sum we take ---
  a <- as.numeric(age_assoc[cols])
  w_lower <- abs(a[li] * a[lj])
  if (all(w_lower == 0)) stop("All weights are zero. Check age_assoc values for selected cols.")
  w_lower <- w_lower / sum(w_lower)
  
  # Store diagonal of M2_ref (needed for updated variances)
  M2_ref_diag <- diag(M2_ref)
  
  # Pre-store lower triangle of M2_ref (so we can update it without forming outer products)
  M2_ref_lower <- M2_ref[lower.tri(M2_ref)]
  
  # Output
  n_tgt <- nrow(X_tgt)
  disco_vals <- numeric(n_tgt)
  
  # Progress bar
  pb <- NULL
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = n_tgt, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  # --- Loop over target samples ---
  # We are adding one sample to the reference set; the new sample size is n_ref + 1.
  # Welford multivariate update:
  # delta  = x - mu
  # mu_new = mu + delta/(n+1)
  # delta2 = x - mu_new
  # M2_new = M2 + outer(delta, delta2)
  # Then cov_new = M2_new / ( (n+1) - 1 ) = M2_new / n_ref
  denom_cov_new <- n_ref  # since (n_ref + 1 - 1) = n_ref
  
  for (i in seq_len(n_tgt)) {
    x <- X_tgt[i, ]
    
    delta <- x - mu_ref
    mu_new <- mu_ref + delta / (n_ref + 1)
    delta2 <- x - mu_new
    
    # Update diagonal variances cheaply: diag(M2_new) = diag(M2_ref) + delta * delta2
    M2_new_diag <- M2_ref_diag + (delta * delta2)
    var_new <- M2_new_diag / denom_cov_new
    sd_new <- sqrt(var_new)
    
    # Again protect against zeros (shouldn’t happen if ref was OK, but can if denom pushes tiny)
    if (any(sd_new == 0)) {
      disco_vals[i] <- NA_real_
      if (show_progress) utils::setTxtProgressBar(pb, i)
      next
    }
    
    # Update only lower triangle of M2 using rank-one formula on indices:
    # M2_new_lower[k] = M2_ref_lower[k] + delta[li[k]] * delta2[lj[k]]
    M2_new_lower <- M2_ref_lower + (delta[li] * delta2[lj])
    
    # Convert to covariance lower triangle
    cov_new_lower <- M2_new_lower / denom_cov_new
    
    # Convert to correlation lower triangle
    denom_new_lower <- sd_new[li] * sd_new[lj]
    R_new_lower <- cov_new_lower / denom_new_lower
    
    # DISCO: n^2 * sum_{i<j} w_ij (R_new - R_ref)^2
    diff2 <- (R_new_lower - R_ref_lower)^2
    disco_vals[i] <- (n_ref^2) * sum(diff2 * w_lower)
    
    if (show_progress) utils::setTxtProgressBar(pb, i)
  }
  
  if (log10_output) disco_vals <- log10(disco_vals)
  disco_vals
}

