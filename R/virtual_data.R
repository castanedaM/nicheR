#' Generate data within based on a ellipsoidal niche
#'
#' @description
#' Simulates \code{n} random points from a multivariate normal distribution
#' defined by the centroid and covariance matrix of a \code{nicheR_ellipsoid}
#' object.
#'
#' @param object A \code{nicheR_ellipsoid} object containing at least
#'   \code{centroid} and \code{cov_matrix}.
#' @param n Integer. The number of virtual points to generate. Default = 100.
#'
#' @details
#' The function uses eigen-decomposition to transform standard normal
#' variables into the coordinate system defined by the ellipsoid's
#' covariance structure.
#'
#' @return
#' A matrix with \code{n} rows and columns corresponding to the
#' dimensions of the input ellipsoid.
#'
#' @export

virtual_data <- function(object, n = 100) {
  # Detecting potential errors
  if (missing(object)) stop("Argument 'object' is required.")
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Argument 'object' must be of class 'nicheR_ellipsoid'.")
  }

  # Ellipsoid features
  centroid <- object$centroid
  cov_matrix <- object$cov_matrix

  # Number of dimensions
  p <- length(centroid)

  # Eigen-decomposition of the covariance matrix
  eS <- eigen(cov_matrix, symmetric = TRUE)
  ev <- eS$values

  # Generate standard normal samples
  vdata <- matrix(rnorm(p * n), nrow = n)

  # Transform using the Square Root of Sigma (V * L^0.5)
  vdata <- drop(centroid) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(vdata)

  # Handle dimension names
  rownames(vdata) <- colnames(cov_matrix)

  # Return the generated data
  return(t(vdata))
}
