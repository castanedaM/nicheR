#' Calculate n-dimensional ellipsoid metrics
#'
#' This function computes the geometric and statistical characteristics of an
#' n-dimensional ellipsoid based on a covariance matrix and a centroid.
#' It calculates the semi-axes lengths, the coordinates of the vertices,
#' the volume, and checks for Spectral Positive Definiteness (SPD).
#'
#' @param cov_matrix A square, numeric covariance matrix representing the
#'   shape and orientation of the niche.
#' @param centroid A numeric vector representing the center of the ellipsoid
#'   (e.g., the mean of the variables).
#' @param cl A numeric value between 0 and 1 representing the confidence
#'   level (e.g., 0.95 for a 95\% confidence ellipsoid).
#' @param verbose Logical; if \code{TRUE}, the function prints progress messages
#'   to the console. Default is \code{TRUE}.
#'
#' @return
#' An object of class \code{\link{nicheR_ellipsoid}}. With all elements
#' computed for the ellipsoid produced by the input parameters.
#'
#' @export
#' @importFrom stats qchisq
#'
#' @examples
#' cm <- matrix(c(11.11, 0, 0, 17777.78), nrow = 2)
#' ctr <- c(20, 600)
#' ellie <- ellipsoid_calculator(cm, ctr, 0.95)

ellipsoid_calculator <- function(cov_matrix,
                                 centroid,
                                 cl,
                                 verbose = TRUE) {
  if (missing(cov_matrix) || missing(centroid) || missing(cl)) {
    stop("All arguments ('cov_matrix', 'centroid', 'cl') must be provided.")
  }
  if (!is.matrix(cov_matrix) || !is.numeric(cov_matrix)) {
    stop("'cov_matrix' must be a numeric matrix.")
  }
  if (nrow(cov_matrix) != ncol(cov_matrix)) {
    stop("'cov_matrix' must be a square matrix.")
  }
  if (!is.numeric(centroid) || length(centroid) != ncol(cov_matrix)) {
    stop("'centroid' must be a numeric vector, 
    with length equal to the number of columns in 'cov_matrix'.")
  }
  if (!is.numeric(cl) || cl <= 0 || cl >= 1) {
    stop("'cl' must be a numeric value between 0 and 1.")
  }

  verbose_message <- function(...) if(isTRUE(verbose)) message(...)

  verbose_message("Step: computing ellipsoid metrics...\n")

  # SPD + inverse
  chol_Sigma <- tryCatch(chol(cov_matrix), error = function(e) NULL)
  if (is.null(chol_Sigma)) stop("Updated covariance matrix is not SPD.")
  Sigma_inv <- chol2inv(chol_Sigma)

  # Cutoff
  chi2_cutoff <- stats::qchisq(cl, df = ncol(cov_matrix))

  # Eigen + axes
  eig <- eigen(cov_matrix, symmetric = TRUE)

  semi_axes_lengths <- sqrt(pmax(eig$values, 0) * chi2_cutoff)
  names(semi_axes_lengths) <- paste0("axis_", seq_along(semi_axes_lengths))

  # Axis coordinates
  axes_coordinates <- lapply(seq_len(ncol(cov_matrix)), function(i) {
    current_length <- semi_axes_lengths[i]
    
    ## Create the matrix of vertex coordinates for the current axis
    rbind(vertex_a = centroid - current_length * eig$vectors[, i],
          vertex_b = centroid + current_length * eig$vectors[, i])
  })
  names(axes_coordinates) <- paste0("axis_", seq_along(axes_coordinates))

  volume <- ellipsoid_volume(n_dimensions = ncol(cov_matrix),
                             semi_axes_lengths = semi_axes_lengths)

  cov_limits <- covariance_limits(cov_matrix)

  verbose_message("Done: updated ellipsoidal niche boundary. For remianing limits see out$cov_limits_remaining\n")

  new_nicheR_ellipsoid(
    dimensions = ncol(cov_matrix),
    var_names = colnames(cov_matrix),
    centroid = centroid,
    cov_matrix = cov_matrix,
    Sigma_inv = Sigma_inv,
    chol_Sigma = chol_Sigma,
    eigen = list(vectors = eig$vectors, values = eig$values),
    cl = cl,
    chi2_cutoff = chi2_cutoff,
    semi_axes_lengths = semi_axes_lengths,
    axes_coordinates = axes_coordinates,
    volume = volume,
    cov_limits = cov_limits
  )
}
