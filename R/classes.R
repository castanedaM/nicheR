#' nicheR_ellipsoid Class Constructor
#'
#' Internal helper to create a nicheR_ellipsoid object.
#'
#' @param dimensions The number of variables/dimensions.
#' @param var_names Names of the original variables.
#' @param centroid The numeric center vector.
#' @param cov_matrix The covariance matrix.
#' @param Sigma_inv The precision matrix.
#' @param chol_Sigma The Cholesky decomposition.
#' @param eigen List of eigenvectors and eigenvalues.
#' @param cl Confidence level.
#' @param chi2_cutoff Chi-squared quantile.
#' @param semi_axes_lengths Radii of the ellipsoid.
#' @param axes_coordinates List of vertex matrices.
#' @param volume Calculated hyper-volume.
#' @param cov_limits Axis-aligned limits.
#'
#' @return
#' An object of class \code{nicheR_ellipsoid} with the fields described above.
#' 
#' @keywords internal

new_nicheR_ellipsoid <- function(dimensions, var_names, centroid, cov_matrix,
                                 Sigma_inv, chol_Sigma, eigen, cl,
                                 chi2_cutoff, semi_axes_lengths,
                                 axes_coordinates, volume, cov_limits) {
  
  structure(
    list(
      dimensions = dimensions,
      var_names = var_names,
      centroid = centroid,
      cov_matrix = cov_matrix,
      Sigma_inv = Sigma_inv,
      chol_Sigma = chol_Sigma,
      eigen = eigen,
      cl = cl,
      chi2_cutoff = chi2_cutoff,
      semi_axes_lengths = semi_axes_lengths,
      axes_coordinates = axes_coordinates,
      volume = volume,
      cov_limits = cov_limits
    ),
    class = "nicheR_ellipsoid"
  )
}



#' nicheR_community Class Constructor
#'
#' @description
#' Helper function to construct a \code{nicheR_community} object.
#'
#' @param ellipse_community A list of \code{nicheR_ellipsoid} objects.
#' @param reference The \code{nicheR_ellipsoid} object used as a template.
#' @param pattern Character. The generation pattern: "random", "nested", or
#'   "conserved".
#' @param n Integer. Number of ellipses in the community.
#' @param smallest_proportion Numeric. Minimum scaling factor used.
#' @param largest_proportion Numeric. Maximum scaling factor used
#'   (NA if not applicable).
#' @param bias Numeric. Bias exponent used (NA if not applicable).
#' @param seed Integer. Seed used for reproducibility (NA if not applicable).
#'
#' @return
#' An object of class \code{nicheR_community}. The object
#' aggregates a collection of ellipsoids, the reference niche they were
#' derived from, and the metadata of the generation process.
#'
#' @keywords internal
#'
#' @export

new_nicheR_community <- function(ellipse_community, reference, pattern, n,
                                 smallest_proportion, largest_proportion = NA,
                                 bias = NA, seed = NA) {

  details <- data.frame(
    pattern = pattern,
    n = n,
    cl = reference$cl,
    smallest_proportion = smallest_proportion,
    largest_proportion = largest_proportion,
    bias = bias,
    seed = seed
  )
  
  structure(
    list(
      details = details,
      reference = reference,
      ellipse_community = ellipse_community
    ),
    class = "nicheR_community"
  )
}