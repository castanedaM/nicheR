#' Build a Probabilistic Ellipsoidal Niche
#'
#' Constructs a multivariate normal (MVN) ellipsoid in environmental space,
#' defined by a centroid, marginal tolerances, and optional correlation-based
#' tilt among dimensions. The ellipsoid represents a probability contour of
#' constant Mahalanobis distance.
#'
#' The niche is defined as:
#' \deqn{X \sim \mathcal{N}(\mu, \Sigma)}
#'
#' where the covariance matrix is decomposed as:
#' \deqn{\Sigma = D R D}
#'
#' with \eqn{D} representing marginal standard deviations and \eqn{R} the
#' correlation (tilt) matrix.
#'
#' Two probability levels are used:
#'
#' \describe{
#'   \item{range_level}{Probability mass represented by the user-supplied
#'   minimum and maximum values. Used only to infer marginal standard deviations
#'   when defining the multivariate normal niche.}
#'
#'   \item{level}{Probability mass enclosed by the ellipsoid contour, defined by
#'   a chi-square quantile of the Mahalanobis distance.}
#' }
#'
#' @param range Input defining the niche centroid and marginal tolerances.
#'   Supported formats include:
#'   \itemize{
#'     \item data.frame with columns \code{min} and \code{max}
#'     \item list with \code{mu} and \code{sd}, or \code{mean} and \code{sd}
#'     \item list with \code{mu} and covariance matrix \code{Sigma}
#'     \item data.frame or matrix of observations (columns = dimensions)
#'   }
#'
#' @param range_level Numeric between 0 and 1. Probability mass represented by
#'   the supplied minimum and maximum values when \code{range} is defined by
#'   bounds. Default is 0.95.
#'
#' @param level Numeric between 0 and 1. Probability mass enclosed by the
#'   ellipsoid contour. Default is 0.95.
#'
#' @param cor_tilt Optional correlation structure controlling ellipsoid
#'   orientation. Can be:
#'   \itemize{
#'     \item \code{NULL} for no tilt (independent dimensions)
#'     \item a \eqn{p \times p} correlation matrix
#'     \item a matrix or data.frame of observations used to estimate correlation
#'   }
#'
#' @param n_points Integer. Number of surface points used to approximate the
#'   ellipsoid boundary for visualization. Default is 100.
#'
#' @return An object of class \code{"nicheR_ellipsoid"} containing:
#' \itemize{
#'   \item \code{mu}: niche centroid
#'   \item \code{Sigma}: covariance matrix
#'   \item \code{Sigma_inv}: inverse covariance matrix
#'   \item \code{axes_sd}: marginal standard deviations
#'   \item \code{axes_boundary}: semi-axis lengths at the specified \code{level}
#'   \item \code{boundary_points}: surface coordinates for plotting
#'   \item \code{c2}: chi-square cutoff
#' }
#'
#' @details
#' The ellipsoid boundary corresponds to:
#' \deqn{(x - \mu)^T \Sigma^{-1} (x - \mu) \le \chi^2_p(\text{level})}
#'
#' This formulation allows separation of:
#' \itemize{
#'   \item marginal niche breadth (\code{range})
#'   \item correlation structure (\code{cor_tilt})
#'   \item probability contour (\code{level})
#' }
#'
#' @seealso \code{\link{parse_range}}, \code{\link{parse_cor_tilt}}
#'
#' @export
build_ellipsoid <- function(range,
                            range_level = 0.95,  # distribution of individual dimensions
                            level = 0.95,        # distribution of the niche area
                            cor_tilt = NULL,
                            n_points = 100) {

  # Validate inputs
  stopifnot(is.numeric(range_level), length(range_level) == 1, range_level > 0, range_level < 1)
  stopifnot(is.numeric(level), length(level) == 1, level > 0, level < 1)
  stopifnot(is.numeric(n_points), length(n_points) == 1, n_points >= 10)

  # 1) Parse range into mu and marginal sd
  parsed <- parse_range(range, range_level = range_level)
  mu <- parsed$mu
  sd_vec <- parsed$sd
  p <- length(mu)

  # 2) Convert cor_tilt into correlation matrix R
  R <- parse_cor_tilt(cor_tilt = cor_tilt, p = p)

  # 3) Build covariance matrix: Sigma = D R D
  D <- diag(sd_vec, nrow = p, ncol = p)
  Sigma <- D %*% R %*% D

  # 4) Validate Sigma is symmetric positive definite
  Sigma <- (Sigma + t(Sigma)) / 2

  chol_Sigma <- tryCatch(chol(Sigma), error = function(e) NULL)

  if (is.null(chol_Sigma)) {
    stop(
      "Sigma is not positive definite.\n",
      "Check cor_tilt and implied standard deviations from range."
    )
  }

  Sigma_inv <- chol2inv(chol_Sigma)

  # 5) Chi square cutoff for the ellipsoid boundary at level
  c2 <- stats::qchisq(level, df = p)

  # 6) Eigen for axes and rotation
  eig <- eigen(Sigma, symmetric = TRUE)
  axes_boundary <- sqrt(eig$values * c2)

  volume <- ellipsoid_volume(n_dimensions = p,
                       semi_axes_length = axes_boundary)

  # 7) Surface points for plotting
  boundary_points <- ellipsoid_surface_points(mu, Sigma, c2, n_points)

  out <- list(
    p = p,
    mu = as.numeric(mu),
    Sigma = Sigma,
    Sigma_inv = Sigma_inv,
    chol_Sigma = chol_Sigma,
    eigen = list(vectors = eig$vectors, values = eig$values),
    range_level = range_level,
    level = level,
    c2 = c2,
    axes_sd = as.numeric(sd_vec),
    axes_boundary = as.numeric(axes_boundary),
    volume = volume,
    boundary_points = boundary_points
  )

  class(out) <- c("nicheR_ellipsoid", "list")

  out
}


#' Parse Niche Range Input
#'
#' Converts user-defined niche range inputs into a centroid vector and marginal
#' standard deviations suitable for constructing a multivariate normal niche.
#'
#' @param range Input defining the niche center and tolerances. Supported formats:
#' \itemize{
#'   \item data.frame with columns \code{min} and \code{max}
#'   \item list with \code{mu} and \code{sd}
#'   \item list with \code{mean} and \code{sd}
#'   \item list with \code{mu} and covariance matrix \code{Sigma}
#'   \item matrix or data.frame of observations
#' }
#'
#' @param range_level Numeric between 0 and 1. Probability mass represented by
#'   the supplied minimum and maximum values when bounds are provided.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{mu}: centroid vector
#'   \item \code{sd}: marginal standard deviations
#'   \item \code{p}: number of dimensions
#'   \item \code{input_type}: detected input format
#' }
#'
#' @details
#' When \code{min} and \code{max} values are provided, they are interpreted as a
#' central probability interval of a normal distribution with probability
#' \code{range_level}. Standard deviations are derived using the corresponding
#' normal quantile.
#'
#' @keywords internal
parse_range <- function(range, range_level = 0.95) {

  if (!is.numeric(range_level) || length(range_level) != 1 || range_level <= 0 || range_level >= 1) {
    stop("range_level must be a single number between 0 and 1 (exclusive).")
  }

  # Option A: data frame with min and max
  if (is.data.frame(range) && all(c("min", "max") %in% names(range))) {

    mins <- as.numeric(range$min)
    maxs <- as.numeric(range$max)

    if (length(mins) == 0) stop("range has no rows.")
    if (length(mins) != length(maxs)) stop("min and max must have the same length.")
    if (any(!is.finite(mins)) || any(!is.finite(maxs))) stop("range contains non-finite values in min/max.")
    if (any(maxs <= mins)) stop("Each max must be greater than min.")
    mu_vec <- (mins + maxs) / 2

    z <- stats::qnorm((1 + range_level) / 2)
    if (!is.finite(z) || z <= 0) stop("Invalid range_level produced non-finite z score.")

    sd_vec <- (maxs - mins) / (2 * z)
    if (any(!is.finite(sd_vec)) || any(sd_vec <= 0)) stop("Derived sd contains non-positive or non-finite values.")

    return(list(mu = mu_vec, sd = sd_vec, p = length(mu_vec), input_type = "min_max"))
  }

  # Option B: list input
  if (is.list(range)) {

    if (!is.null(range$mu) && !is.null(range$sd)) {
      mu_vec <- as.numeric(range$mu)
      sd_vec <- as.numeric(range$sd)

      if (length(mu_vec) == 0) stop("range$mu is empty.")
      if (length(sd_vec) != length(mu_vec)) stop("range$sd must be the same length as range$mu.")
      if (any(!is.finite(mu_vec)) || any(!is.finite(sd_vec))) stop("range$mu or range$sd contains non-finite values.")
      if (any(sd_vec <= 0)) stop("range$sd must be positive.")

      return(list(mu = mu_vec, sd = sd_vec, p = length(mu_vec), input_type = "mu_sd"))
    }

    if (!is.null(range$mean) && !is.null(range$sd)) {
      mu_vec <- as.numeric(range$mean)
      sd_vec <- as.numeric(range$sd)

      if (length(mu_vec) == 0) stop("range$mean is empty.")
      if (length(sd_vec) != length(mu_vec)) stop("range$sd must be the same length as range$mean.")
      if (any(!is.finite(mu_vec)) || any(!is.finite(sd_vec))) stop("range$mean or range$sd contains non-finite values.")
      if (any(sd_vec <= 0)) stop("range$sd must be positive.")

      return(list(mu = mu_vec, sd = sd_vec, p = length(mu_vec), input_type = "mean_sd"))
    }

    if (!is.null(range$mu) && !is.null(range$Sigma)) {
      mu_vec <- as.numeric(range$mu)
      Sigma <- as.matrix(range$Sigma)

      if (length(mu_vec) == 0) stop("range$mu is empty.")
      if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) stop("range$Sigma must be a square matrix.")
      if (nrow(Sigma) != length(mu_vec)) stop("range$Sigma dimensions must match length of range$mu.")
      if (any(!is.finite(mu_vec)) || any(!is.finite(Sigma))) stop("range$mu or range$Sigma contains non-finite values.")

      sd_vec <- sqrt(diag(Sigma))
      if (any(!is.finite(sd_vec)) || any(sd_vec <= 0)) stop("Sigma implies non-positive or non-finite marginal SDs.")

      return(list(mu = mu_vec, sd = sd_vec, p = length(mu_vec), input_type = "mu_Sigma"))
    }
  }

  # Option C: observations
  if (is.data.frame(range) || is.matrix(range)) {

    X <- as.matrix(range)

    if (nrow(X) < 2) stop("range as observations must have at least 2 rows.")
    if (ncol(X) < 1) stop("range as observations must have at least 1 column.")
    if (any(!is.finite(X))) stop("range observations contain non-finite values.")

    mu_vec <- colMeans(X)
    sd_vec <- apply(X, 2, stats::sd)

    if (any(!is.finite(mu_vec)) || any(!is.finite(sd_vec))) stop("Computed mu/sd from observations contains non-finite values.")
    if (any(sd_vec <= 0)) stop("At least one column has zero (or negative) SD; cannot build ellipsoid.")

    return(list(mu = as.numeric(mu_vec), sd = as.numeric(sd_vec), p = length(mu_vec), input_type = "observations"))
  }

  stop(
    "Unsupported range format. Provide one of:\n",
    "  (A) data.frame with columns min and max;\n",
    "  (B) list(mu, sd) or list(mean, sd) or list(mu, Sigma);\n",
    "  (C) data.frame/matrix of observations (cols = dimensions)."
  )
}


# Helpers for correlation tilt (orientation)
make_identity_cor <- function(p) {
  diag(1, p)
}

suggest_cor_tilt <- function(p, sd_vec) {
  # TO DO: later implement logic to generate candidate cor_tilt matrices that
  # keep Sigma positive definite for the given sd_vec. Also send message to user
  # instead of having the user  to restaar for him to send option
  NULL
}


#' Parse Correlation Tilt Matrix
#'
#' Converts user-supplied correlation information into a valid correlation
#' matrix controlling ellipsoid orientation in environmental space.
#'
#' @param cor_tilt Correlation structure describing variable co-variation.
#' Can be:
#' \itemize{
#' \item \code{NULL} for no tilt (identity matrix)
#' \item a \eqn{p \times p} correlation matrix
#' \item a matrix or data.frame of observations used to estimate correlation
#' }
#'
#' @param p Integer. Number of environmental dimensions.
#'
#' @return A \eqn{p \times p} correlation matrix.
#'
#' @details
#' The correlation matrix defines ellipsoid orientation but does not affect
#' marginal niche breadth, which is controlled separately by standard deviations.
#'
#' @keywords internal
parse_cor_tilt <- function(cor_tilt, p) {

  # Default: no tilt (independent dimensions)
  if (missing(cor_tilt) || is.null(cor_tilt)) {
    suggest_cor_tilt(p = p, sd_vec = NULL)
    return(make_identity_cor(p))
  }

  # Case 1: user provides correlation matrix directly
  if (is.matrix(cor_tilt)) {

    if (!all(dim(cor_tilt) == c(p, p)))
      stop("cor_tilt must be a p x p matrix.")

    R <- (cor_tilt + t(cor_tilt)) / 2

    if (any(!is.finite(R)))
      stop("cor_tilt contains non-finite values.")

    if (any(abs(diag(R) - 1) > 1e-6))
      stop("cor_tilt must have 1s on the diagonal (correlation matrix).")

    return(R)
  }

  # Case 2: observations used to estimate correlation
  # The user gives a cloud of points that reflects how variables tend to vary together.
  if (is.data.frame(cor_tilt)) cor_tilt <- as.matrix(cor_tilt)

  if (is.matrix(cor_tilt)) {

    if (ncol(cor_tilt) != p)
      stop("cor_tilt as observations must have p columns.")

    if (nrow(cor_tilt) < 3)
      stop("cor_tilt must have at least 3 rows to estimate correlation.")

    if (any(!is.finite(cor_tilt)))
      stop("cor_tilt observations contain non-finite values.")

    R <- stats::cor(cor_tilt, use = "pairwise.complete.obs")
    R <- (R + t(R)) / 2
    diag(R) <- 1

    return(R)
  }

  stop(
    "Unsupported cor_tilt format. Provide:\n",
    "  • NULL (no tilt),\n",
    "  • a p x p correlation matrix, or\n",
    "  • a matrix/data.frame of observations with p columns."
  )
}



#' Compute Ellipsoidal Niche Volume
#'
#' Calculates the geometric volume (or area in 2D) of a p-dimensional ellipsoid
#' defined by its semi-axis lengths.
#'
#' For an ellipsoid with semi-axes \eqn{a_1, a_2, \ldots, a_p}, the volume is:
#'
#' \deqn{
#' V_p = V_{\text{n},p} \prod_{i=1}^{p} a_i
#' }
#'
#' where the volume of the p-dimensional n is:
#'
#' \deqn{
#' V_{\text{n},p} =
#' \frac{\pi^{p/2}}{\Gamma\left(\frac{p}{2} + 1\right)}
#' }
#'
#' In the context of probabilistic niche modeling, the semi-axis lengths are
#' typically defined as:
#'
#' \deqn{
#' a_i = \sqrt{\lambda_i \, c^2}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\lambda_i} are the eigenvalues of the covariance matrix
#'   \eqn{\Sigma}
#'   \item \eqn{c^2 = \chi^2_p(\text{level})} is the chi-square cutoff defining the
#'   probability contour
#' }
#'
#' This formulation allows ellipsoid volume to be computed directly from either
#' covariance matrices or axis lengths derived from them.
#'
#' @param n_dimensions Integer. Number of dimensions (\eqn{p}).
#' @param semi_axes_length Numeric vector of length \code{n_dimensions}
#'   containing the ellipsoid semi-axis lengths.
#'
#' @return Numeric. Geometric volume of the ellipsoid.
#'
#' @details
#' \strong{Interpretation:}
#'
#' \itemize{
#'   \item In 2D, the returned value is area.
#'   \item In 3D, the returned value is volume.
#'   \item In higher dimensions, the result is a hypervolume measured in the
#'     product of environmental units.
#' }
#'
#' The volume depends on the selected probability contour (\code{level}) and
#' therefore increases monotonically with increasing \code{level}.
#'
#' @seealso
#' \code{\link{build_ellipsoid}},
#' \code{\link{ellipsoid_surface_points}}
#'
#' @export
ellipsoid_volume <- function(n_dimensions, semi_axes_length) {

  # ---- validation ----
  if (missing(n_dimensions))
    stop("Argument 'n_dimensions' must be defined.")

  if (missing(semi_axes_length))
    stop("Argument 'semi_axes_length' must be defined.")

  if (!is.numeric(n_dimensions) || length(n_dimensions) != 1)
    stop("'n_dimensions' must be a single numeric value.")

  if (!is.numeric(semi_axes_length))
    stop("'semi_axes_length' must be numeric.")

  if (length(semi_axes_length) != n_dimensions)
    stop("'semi_axes_length' must have length equal to 'n_dimensions'.")

  if (any(semi_axes_length <= 0))
    stop("All semi-axis lengths must be positive.")

  # ---- p-dimensional n volume ----
  unit_ball_volume <- pi^(n_dimensions / 2) /
    gamma(n_dimensions / 2 + 1)


  # ---- ellipsoid volume ----
  volume <- unit_ball_volume * prod(semi_axes_length)

  return(volume)
}



#' Generate Ellipsoid Surface Coordinates
#'
#' Computes boundary coordinates of a 2D or 3D ellipsoid defined by a covariance
#' matrix and chi-square probability contour.
#'
#' @param mu Numeric vector. Ellipsoid centroid.
#' @param Sigma Covariance matrix.
#' @param c2 Chi-square cutoff defining the ellipsoid contour.
#' @param n_points Number of points used to approximate the surface.
#'
#' @return A data.frame containing surface coordinates:
#' \itemize{
#'   \item \code{x1, x2} for 2D
#'   \item \code{x1, x2, x3} for 3D
#' }
#'
#' @details
#' Surface points are generated by rotating and scaling unit circles or spheres
#' using the eigen decomposition of the covariance matrix.
#'
#' @keywords internal
ellipsoid_surface_points <- function(mu,
                                     Sigma,
                                     c2,
                                     n_points) {

  p <- length(mu)

  # eigen decomposition gives rotation + axis lengths
  eig <- eigen(Sigma, symmetric = TRUE)

  # semi-axis lengths at chi-square level
  axes_boundary <- sqrt(eig$values * c2)

  # rotation matrix
  R <- eig$vectors

  if (p == 2) {

    t <- seq(0, 2 * pi, length.out = n_points)

    unit_ellipse <- rbind(
      axes_boundary[1] * cos(t),
      axes_boundary[2] * sin(t)
    )

    rotated_points <- t(R %*% unit_ellipse)

    surface <- sweep(rotated_points, 2, mu, "+")

    surface <- as.data.frame(surface)
    colnames(surface) <- c("x1", "x2")

    return(surface)
  }

  if (p == 3) {

    grid <- expand.grid(
      u = seq(0, 2 * pi, length.out = n_points),
      v = seq(0, pi, length.out = n_points)
    )

    unit_ellipsoid <- rbind(
      axes_boundary[1] * sin(grid$v) * cos(grid$u),
      axes_boundary[2] * sin(grid$v) * sin(grid$u),
      axes_boundary[3] * cos(grid$v)
    )

    rotated_points <- t(R %*% unit_ellipsoid)

    surface <- sweep(rotated_points, 2, mu, "+")

    surface <- as.data.frame(surface)
    colnames(surface) <- c("x1", "x2", "x3")

    return(surface)
  }

  stop("Surface points only implemented for 2D or 3D.")
}


#' @export
print.nicheR_ellipsoid <- function(x, digits = 3, ...) {

  cat("nicheR Ellipsoid Object\n")
  cat("-----------------------\n")

  cat(" Dimensions:      ", x$p, "D\n", sep = "")
  cat(" Centroid (mu):   ",
      paste(round(x$mu, digits), collapse = ", "), "\n", sep = "")
  cat(" Marginal SDs:    ",
      paste(round(x$axes_sd, digits), collapse = ", "), "\n", sep = "")

  cat(" Range level:     ", round(x$range_level, digits), "\n", sep = "")
  cat(" Ellipsoid level: ", round(x$level, digits), "\n", sep = "")
  cat(" Chi-square c²:   ", round(x$c2, digits), "\n", sep = "")

  cat("\n Ellipsoid semi-axes at level:\n  ",
      paste(round(x$axes_boundary, digits), collapse = ", "), "\n", sep = "")

  cat("\n Correlation tilt (R):\n")
  print(round(cov2cor(x$Sigma), digits))

  cat("\n Covariance matrix (Sigma):\n")
  print(round(x$Sigma, digits))

  cat("\n Inverse covariance matrix (Sigma⁻¹):\n")
  print(round(x$Sigma_inv, digits))

  cat("\n Volume:\n")
  print(round(x$volume, digits))

  cat("\n Boundary points: ",
      nrow(x$boundary_points),
      ifelse(x$p == 2, " (x, y)\n", " (x, y, z)\n"),
      sep = "")

  invisible(x)
}
