###############################################################################
# calculate empirical p values ================================================
###############################################################################

#' Calculate empirical p-values using absolute values
#'
#' Computes empirical p-values by comparing the absolute values of
#' observed statistics against a null distribution. Uses continuity correction
#' to avoid p-values of exactly zero.
#'
#' @param .stat Numeric vector of observed statistics.
#' @param .null Numeric vector representing the null distribution.
#'
#' @return Numeric vector of empirical p-values, same length as `.stat`.
#'   Values represent the probability of observing statistics as extreme or
#'   more extreme than the observed values under the null hypothesis.
#'
#' @details
#' The function applies a continuity correction by adding 1 to both the
#' numerator and denominator: `(sum(abs(.null) >= abs(.stat)) + 1) /
#' (length(.null) + 1)`. This prevents p-values of exactly zero and provides
#' more conservative estimates.
#'
#' @examples
#' \dontrun{
#' observed <- c(0.5, -0.8, 0.3)
#' null_dist <- rnorm(1000, mean = 0, sd = 0.2)
#' p_values <- empirical_p(observed, null_dist)
#' }
empirical_p <- function(.stat, .null) {
  # we want to get both negative and positive loading genes
  abs_stat <- abs(.stat)
  abs_null <- abs(.null)

  vapply(
    abs_stat,
    function(s) {
      (sum(abs_null >= s) + 1) / (length(abs_null) + 1) # +1 for continuity correction
    },
    numeric(1)
  )
}

###############################################################################
# generate null loadings ======================================================
###############################################################################

#' Generate null loadings for jackstraw procedure
#'
#' Randomly permutes a subset of features (rows) and performs PCA to generate
#' null loadings for statistical testing. This is a core component of the
#' jackstraw procedure for identifying statistically significant features in
#' PCA.
#'
#' @param x A numeric matrix with features as rows and samples as columns.
#' @param s Integer. Number of features to randomly sample and permute.
#' @param ndim Integer. Number of principal components to compute. If `NULL`,
#'   all components are computed.
#' @param center Logical. Whether to center the data (subtract column means)
#'   before PCA. Default is `TRUE`.
#' @param scale. Logical. Whether to scale the data (divide by column standard
#'   deviations) before PCA. Default is `TRUE`.
#'
#' @return A data frame containing the null loadings for the permuted features.
#'   Rows correspond to the `s` randomly selected features, columns correspond
#'   to principal components.
#'
#' @examples
#' \dontrun{
#' # Generate null loadings for 10 randomly selected features
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' null_loadings <- js_generate_null_loadings(
#'   data,
#'   s = 10,
#'   ndim = 5,
#'   center = TRUE,
#'   scale. = TRUE
#' )
#' }
js_generate_null_loadings <- function(
  x,
  s,
  ndim = NULL,
  center = TRUE,
  scale. = TRUE
) {
  idx <- sample(nrow(x), size = s)
  js_dat <- x
  for (i in idx) {
    js_dat[i, ] <- sample(js_dat[i, ])
  }
  js_pca <- prcomp(t(js_dat), center = center, scale. = scale., rank. = ndim)
  null_loadings <- js_pca$rotation[idx, ] |>
    as.data.frame()
  rownames(null_loadings) <- NULL
  return(null_loadings)
}

###############################################################################
# run PCA with jackstraw procedure ============================================
###############################################################################

#' Run PCA with jackstraw significance testing
#'
#' Performs principal component analysis and uses the jackstraw procedure to
#' identify features with statistically significant loadings. The method
#' generates a null distribution by randomly permuting subsets of features
#' across multiple iterations.
#'
#' @param x A numeric matrix with features as rows and samples as columns.
#' @param s Integer. Number of features to randomly sample and permute in each
#'   jackstraw iteration.
#' @param B Integer. Number of jackstraw iterations to perform for building
#'   the null distribution.
#' @param ndim Integer. Number of principal components to compute. If `NULL`,
#'   all components are computed.
#' @param center Logical. Whether to center the data (subtract column means)
#'   before PCA. Default is `TRUE`.
#' @param scale. Logical. Whether to scale the data (divide by column standard
#'   deviations) before PCA. Default is `TRUE`.
#'
#' @return A list containing:
#'   \item{pca_out}{The result of `prcomp()` on the original data.}
#'   \item{pval}{A matrix of empirical p-values with the same dimensions as
#'     the PCA loadings matrix.}
#'   \item{nulldist}{A data frame containing all null loadings from the
#'     jackstraw iterations, with an `iter` column indicating the iteration
#'     number.}
#'
#' @details
#' The jackstraw procedure works by:
#' 1. Computing PCA on the original data
#' 2. For each of `B` iterations, randomly selecting `s` features, permuting
#'    their values across samples, and computing PCA loadings
#' 3. Comparing observed loadings to the null distribution using empirical
#'    p-values based on absolute values
#'
#' Both the original PCA and jackstraw iterations use the same centering and
#' scaling parameters to ensure comparable results.
#'
#' @examples
#' \dontrun{
#' # Run jackstraw on a gene expression matrix
#' data <- matrix(rnorm(5000), nrow = 500, ncol = 10)
#' rownames(data) <- paste0("Gene", 1:500)
#'
#' results <- js_run_pca_with_jackstraw(
#'   data,
#'   s = 50,
#'   B = 100,
#'   ndim = 5,
#'   center = TRUE,
#'   scale. = TRUE
#' )
js_run_pca_with_jackstraw <- function(
  x,
  s,
  B,
  ndim = NULL,
  center = TRUE,
  scale. = TRUE
) {
  message("Running PCA on original data")
  pca_out <- prcomp(t(x), center = center, scale. = scale., rank. = ndim)
  loadings <- pca_out$rotation

  message("Generating null distribution of loadings")
  nulldist <- purrr::map(
    seq_len(B),
    ~ js_generate_null_loadings(x, s, ndim, center = center, scale. = scale.),
    .progress = TRUE
  ) |>
    dplyr::bind_rows(.id = "iter")

  message("Calculating empirical p values")
  empp <- loadings
  for (pc in colnames(empp)) {
    empp[, pc] <- empirical_p(empp[, pc], nulldist[[pc]])
  }

  return(list(pca_out = pca_out, pval = empp, nulldist = nulldist))
}

# end =========================================================================
