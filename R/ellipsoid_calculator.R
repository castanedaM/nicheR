#' Calculate n-dimensional ellipsoid metrics
#'
#' Computes geometric and probabilistic metrics for an n-dimensional ellipsoid
#' defined by a centroid and covariance matrix, including semi-axis lengths,
#' axis vertices, and hypervolume for a chi-square confidence contour.
#'
#' @details
#' The ellipsoid boundary is defined by the constant Mahalanobis distance contour:
#' \deqn{(x - \mu)^\top \Sigma^{-1} (x - \mu) = c^2,}
#' where \eqn{\mu} is the centroid, \eqn{\Sigma} is the covariance matrix, and
#' \eqn{c^2 = \chi^2_{n}(\mathrm{cl})} is the chi-square cutoff with \eqn{n}
#' degrees of freedom.
#'
#' The covariance matrix must be symmetric positive definite (SPD). The inverse
#' covariance is computed via the Cholesky factorization. Semi-axis lengths are
#' computed from covariance eigenvalues \eqn{\lambda_i} as:
#' \deqn{a_i = \sqrt{\lambda_i c^2}.}
#'
#' Axis vertices are computed along each eigenvector direction as
#' \eqn{\mu \pm a_i v_i}. Hypervolume is computed with
#' \code{\link{ellipsoid_volume}}, and covariance-derived limits with
#' \code{\link{covariance_limits}}.
#'
#' @param cov_matrix A square, numeric covariance matrix \eqn{\Sigma}. Must be SPD.
#'   Row/column names (if provided) are used as variable names in the output.
#' @param centroid Numeric vector giving the centroid \eqn{\mu}. Must have length
#'   equal to \code{ncol(cov_matrix)}.
#' @param cl Numeric confidence level in (0, 1). Used to compute the chi-square
#'   cutoff defining the ellipsoid contour.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'
#' @return
#' An object of class \code{"nicheR_ellipsoid"} created by
#' \code{\link{new_nicheR_ellipsoid}}, containing ellipsoid geometry and
#' associated quantities (e.g., centroid, covariance matrix, chi-square cutoff,
#' semi-axis lengths, axis vertex coordinates, volume, and covariance limits).
#'
#' @examples
#' cm <- matrix(c(11.11, 0,
#'                0, 17777.78), nrow = 2, byrow = TRUE)
#' colnames(cm) <- rownames(cm) <- c("var1", "var2")
#' ctr <- c(20, 600)
#' ell <- ellipsoid_calculator(cov_matrix = cm, centroid = ctr, cl = 0.95, verbose = FALSE)
#'
#' @importFrom stats qchisq
#' @export
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

   verbose_message(verbose, "Step: computing ellipsoid metrics...\n")

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
    vec <- eig$vectors[, i] * semi_axes_lengths[i]

    ## Create the matrix of vertex coordinates for the current axis
    rbind(vertex_a = centroid - vec, vertex_b = centroid + vec)
  })
  names(axes_coordinates) <- paste0("axis_", seq_along(axes_coordinates))

  volume <- ellipsoid_volume(n_dimensions = ncol(cov_matrix),
                             semi_axes_lengths = semi_axes_lengths)

  cov_limits <- covariance_limits(cov_matrix)

  verbose_message(verbose, "Done: updated ellipsoidal niche metrics")

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
