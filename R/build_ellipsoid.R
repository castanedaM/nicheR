#' Build a 2D or 3D Ellipsoid Object
#'
#' Creates a structured object representing a rotated ellipsoid. This object
#' contains all necessary matrices for geometric and statistical calculations,
#' as well as surface points for plotting. This function is a core component for
#' defining virtual niches in environmental space.
#'
#' @param center A numeric vector of length 2 or 3 defining the ellipsoid's center.
#'   This vector can be named (e.g., `c(Temp = 15, Precip = 800)`) to link
#'   dimensions to environmental variables.
#' @param axes A numeric vector with the same length as `center`, defining the
#'   semi-axis lengths *before* rotation. These values represent the standard
#'   deviations of the ellipsoid along its principal axes. Must be positive.
#' @param angles A numeric vector of rotation angles in radians. For 2D, a single
#'   angle is sufficient. For 3D, a vector of three angles `c(ax, ay, az)`
#'   is required, corresponding to rotations around the x, y, and z axes,
#'   respectively.
#' @param n_points An integer specifying the resolution for the plotted
#'   surface. A higher value results in a smoother surface but increases
#'   computation time and file size.
#'
#' @return An S3 object of class `ellipsoid`, which is a list containing:
#'   \item{center}{The numeric vector defining the ellipsoid's center.}
#'   \item{axes}{The numeric vector of semi-axis lengths.}
#'   \item{angles}{The numeric vector of rotation angles.}
#'   \item{dimen}{The dimensionality of the ellipsoid (2 or 3).}
#'   \item{R}{The rotation matrix used to orient the ellipsoid.}
#'   \item{Sigma}{The covariance matrix of the ellipsoid, based on `axes` and
#'     `R`.}
#'   \item{Sigma_inv}{The inverse of the covariance matrix, useful for
#'     Mahalanobis distance calculations.}
#'   \item{surface}{A data frame of points representing the ellipsoid's surface,
#'     suitable for plotting.}
#'
#' @details The function internally calculates the covariance matrix (`Sigma`)
#'   and its inverse, which are fundamental for computing Mahalanobis distance
#'   to determine if a point falls within the ellipsoid. The `surface` points
#'   are generated as a set of points on a unit sphere (in 3D) or circle (in 2D)
#'   that are then scaled by `axes` and rotated by `R` before being shifted to
#'   the `center`.
#'
#' @family ellipsoid functions
#' @export
build_ellps <- function(center = c(x = 0, y = 0),
                            axes = c(1, 1),
                            angles = 0,
                            n_points = 100) {

  dimen <- length(center)

  # --- 1. Input Validation ---
  if (!(dimen %in% c(2, 3))) {
    stop("'center' must have exactly 2 or 3 elements.")
  }

  if (length(axes) != dimen) {
    stop("'center' and 'axes' must have the same number of elements.")
  }

  if (!(length(angles) %in% c(1, dimen))) {
    stop(sprintf(
      "For a %sD shape, 'angles' must be a single value or a vector of length %s.",
      dimen, dimen
    ))
  }

  # Warning when using a single vector
  if (length(angles) == 1) {
    warning(sprintf("All axes will be applied the same angle specified: %s",
                    angles))
    angles <- rep(angles, dimen)  # repeat to match number of dimensions
  }

  # --- 2. Rotation Matrix (R) ---
  if (dimen == 2) {
    theta <- angles[1]
    R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  } else {
    ax <- angles[1]; ay <- angles[2]; az <- angles[3]
    Rx <- matrix(c(1, 0, 0, 0, cos(ax), -sin(ax), 0, sin(ax), cos(ax)), 3, 3)
    Ry <- matrix(c(cos(ay), 0, sin(ay), 0, 1, 0, -sin(ay), 0, cos(ay)), 3, 3)
    Rz <- matrix(c(cos(az), -sin(az), 0, sin(az), cos(az), 0, 0, 0, 1), 3, 3)
    R <- Rz %*% Ry %*% Rx
  }

  # --- 3. Shape (Covariance) Matrix (Sigma) ---
  D_sq      <- diag(axes^2)
  Sigma     <- R %*% D_sq %*% t(R)
  Sigma_inv <- solve(Sigma)

  # Assign variable names to match users named variables, but they wont have
  # names if the user does not provide names
  var_names <- names(center)
  if (!is.null(var_names)) {
    colnames(Sigma) <- rownames(Sigma) <- var_names
    colnames(Sigma_inv) <- rownames(Sigma_inv) <- var_names
  }

  # --- 4. Surface Point Generation ---
  if (dimen == 2) {
    t <- seq(0, 2 * pi, length.out = n_points)
    unit_ellipse <- rbind(axes[1] * cos(t), axes[2] * sin(t))
    rotated_points <- t(R %*% unit_ellipse)
    surface <- sweep(rotated_points, 2, center, "+")
    surface <- as.data.frame(surface)
    colnames(surface) <- c("x", "y")
  } else { # 3D
    grid <- expand.grid(u = seq(0, 2 * pi, length.out = n_points),
                        v = seq(0, pi, length.out = n_points))

    unit_ellipsoid <- rbind(axes[1] * sin(grid$v) * cos(grid$u),
                            axes[2] * sin(grid$v) * sin(grid$u),
                            axes[3] * cos(grid$v))

    rotated_points <- t(R %*% unit_ellipsoid)

    surface <- sweep(rotated_points, 2, center, "+")

    surface <- as.data.frame(surface)
    colnames(surface) <- c("x", "y", "z")
  }

  structure(
    list(center = center,
         axes = axes,
         angles = angles,
         dimen = dimen,
         R = R,
         Sigma = Sigma,
         Sigma_inv = Sigma_inv,
         surface = surface),
    class = "ellipsoid"
  )
}
