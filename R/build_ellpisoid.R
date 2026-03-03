#' Build a probabilistic ellipsoidal niche from ranges
#'
#' Builds an ellipsoidal niche in multivariate environmental space using a
#' multivariate normal (MVN) contour defined by a constant Mahalanobis distance.
#' The ellipsoid is parameterized from user-provided variable ranges by deriving
#' a centroid and marginal standard deviations, and assuming a diagonal covariance
#' matrix.
#'
#' @details
#' \code{range} must be a 2-row \code{matrix} or \code{data.frame} with variables
#' in columns. Rows represent lower and upper bounds for each variable (row order
#' may be min/max or max/min). The centroid is computed as:
#' \deqn{\mu_i = (m_i + M_i)/2.}
#'
#' Marginal standard deviations are derived assuming bounds represent
#' approximately \eqn{\pm 3} standard deviations:
#' \deqn{\sigma_i = (M_i - m_i)/6,}
#' and a diagonal covariance matrix is assumed:
#' \deqn{\Sigma = \mathrm{diag}(\sigma_1^2,\dots,\sigma_n^2).}
#'
#' The ellipsoid contour is defined using a chi-square cutoff
#' \eqn{c^2 = \chi^2_{n}(\mathrm{cl})}, where \eqn{n} is the number of variables.
#'
#' @param range A 2-row \code{matrix} or \code{data.frame} of bounds, with variables
#'   as columns. Rows may be ordered as min/max or max/min. Column names are required
#'   and used as variable names.
#' @param cl Numeric confidence level in (0, 1). Used to compute the chi-square
#'   cutoff defining the ellipsoid contour.
#' @param verbose Logical; if \code{TRUE}, prints brief progress messages.
#'
#' @return
#' An object of class \code{"nicheR_ellipsoid"} produced by
#' \code{\link{ellipsoid_calculator}} (via \code{\link{new_nicheR_ellipsoid}}),
#' containing ellipsoid geometry and associated quantities (e.g., centroid,
#' covariance matrix, chi-square cutoff, semi-axis lengths, axis vertex coordinates,
#' volume, and covariance limits).
#'
#' @examples
#' rng <- data.frame(bio1 = c(10, 20),
#'                   bio2 = c(20, 30))
#' ell <- build_ellipsoid(range = rng, cl = 0.95, verbose = FALSE)
#' print(ell)
#'
#' @export
build_ellipsoid <- function(range,
                            cl = 0.99,
                            verbose = TRUE){

  verbose_message(verbose, "Starting: building ellipsoidal niche from ranges...\n")

  # Input checks -------------------------------------------------------------

  if(!(is.data.frame(range) || is.matrix(range))){
    stop("range must be a data.frame or matrix.")
  }

  range <- as.matrix(range)

  if(is.null(colnames(range))){
    stop("range must have column names (variable names).")
  }

  var_names <- colnames(range)

  if(nrow(range) != 2L){
    stop("range must have exactly 2 rows (min/max) and variables as columns.")
  }

  if(!is.numeric(cl) || length(cl) != 1L || !is.finite(cl) || cl <= 0 || cl >= 1){
    stop("cl must be a single finite number strictly between 0 and 1.")
  }

  # Range parse --------------------------------------------------------------

  r1 <- as.numeric(range[1, ])
  r2 <- as.numeric(range[2, ])

  if(any(!is.finite(r1)) || any(!is.finite(r2))){
    stop("range contains non-finite values.")
  }

  if(all(r1 < r2)){
    mins <- r1
    maxs <- r2
  }else if(all(r2 < r1)){
    mins <- r2
    maxs <- r1
  }else{
    stop("Each variable must have max > min. Please check range formatting.")
  }

  # Center and marginal SDs --------------------------------------------------

  mu_vec <- (mins + maxs)/2
  sd_vec <- (maxs - mins)/6

  if(any(!is.finite(mu_vec)) || any(!is.finite(sd_vec))){
    stop("Derived mu_vec/sd_vec contain non-finite values.")
  }

  if(any(sd_vec <= 0)){
    stop("Derived sd_vec must be positive for all variables (check mins/maxs).")
  }

  # Covariance handling ------------------------------------------------------


  verbose_message(verbose, "Step: computing covariance matrix...\n")

  cov_matrix <- diag(sd_vec^2, nrow = length(mu_vec), ncol = length(mu_vec))
  rownames(cov_matrix) <- var_names
  colnames(cov_matrix) <- var_names
  names(mu_vec) <- var_names

  verbose_message(verbose, "Done: created ellipsoidal niche.\n")

  out <- ellipsoid_calculator(cov_matrix = cov_matrix,
                              centroid = mu_vec, cl = cl,
                              verbose = verbose)
  out
}
