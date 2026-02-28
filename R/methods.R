#' Print a nicheR Ellipsoid Object
#'
#' Provides a concise summary of a \code{nicheR_ellipsoid} object created by
#' \code{\link{build_ellipsoid}}. The printed output includes dimensionality,
#' chi-square cutoff, centroid, covariance matrix, principal semi-axis lengths,
#' axis endpoints, and ellipsoid volume.
#'
#' @param x A \code{nicheR_ellipsoid} object.
#' @param digits Integer. Number of decimal places used when printing numeric
#'   values. Default is 3.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This is an S3 method for objects of class \code{"nicheR_ellipsoid"}.
#' The function formats and rounds key quantities for readability but does
#' not modify the underlying object.
#'
#' @return
#' The input object \code{x}, returned invisibly.
#'
#' @seealso \code{\link{build_ellipsoid}}
#'
#' @method print nicheR_ellipsoid
#' @export

print.nicheR_ellipsoid <- function(x, digits = 3, ...) {

  cat("nicheR Ellipsoid Object\n")
  cat("----------------------\n")

  cat("Dimensions:        ", x$dimensions, "D\n", sep = "")
  cat("Chi-square cutoff: ", round(x$chi2_cutoff, digits), "\n", sep = "")

  cat("Centroid (mu):     ",
      paste(round(x$centroid, digits), collapse = ", "),
      "\n", sep = "")

  cat("\nCovariance matrix:\n")
  print(round(x$cov_matrix, digits))

  cat("\nEllipsoid semi-axis lengths:\n  ",
      paste(round(x$semi_axes_lengths, digits), collapse = ", "),
      "\n", sep = "")

  cat("\nEllipsoid axis endpoints:\n")

  for(i in seq_len(x$dimensions)){
    cat(ifelse(i == 1, "", "\n"), " Axis ", i, ":\n", sep = "")
    print(round(x$axes_coordinates[[i]], digits))
  }

  cat("\nEllipsoid volume:\n")
  print(round(x$volume, digits))

  invisible(x)
}