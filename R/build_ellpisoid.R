#' Build a probabilistic ellipsoidal niche from ranges
#'
#' Constructs an ellipsoidal niche in multivariate environmental space using a
#' multivariate normal (MVN) contour of constant Mahalanobis distance. The niche
#' is defined by a centroid (midpoints of user provided ranges), marginal spreads
#' (derived from ranges), and an optional user supplied variance covariance matrix.
#'
#' @details
#' This function expects \code{range} as a 2 row matrix or data.frame with variables in
#' columns. The two rows represent lower and upper bounds (order can be min max
#' or max min). The centroid is computed as the midpoint of each variable range:
#' \deqn{\mu_i = (m_i + M_i)/2.}
#'
#' If \code{cov_matrix} is not provided, marginal standard deviations are derived from
#' ranges using the rule of thumb that the provided bounds represent roughly
#' \eqn{\pm 3} standard deviations:
#' \deqn{\sigma_i = (M_i - m_i)/6,}
#' and a diagonal covariance matrix is assumed:
#' \deqn{\Sigma = \mathrm{diag}(\sigma_1^2,\dots,\sigma_n^2).}
#'
#' The ellipsoid boundary is defined by:
#' \deqn{(x - \mu)^\top \Sigma^{-1} (x - \mu) = c^2,}
#' where \eqn{c^2} is the chi square cutoff \eqn{\chi^2_{n}(\mathrm{cl})} with
#' \eqn{n} degrees of freedom.
#'
#' Principal semi axis lengths are computed from the eigenvalues of \eqn{\Sigma}:
#' \deqn{a_i = \sqrt{\lambda_i c^2}.}
#'
#' @param range A 2 row \code{matrix} or \code{data.frame} of bounds, with variables as
#'   columns. Rows may be ordered as min max or max min. Column names must be
#'   provided and are stored as \code{var_names}.
#' @param cl Numeric confidence level in (0, 1). Used to set the chi square cutoff
#'   that determines the ellipsoid contour.
#' @param cov_matrix Optional variance covariance matrix \eqn{\Sigma} for the
#'   ellipsoid (must be \code{dimensions x dimensions}, symmetric positive definite).
#'   If provided, it is used directly (no rescaling from \code{range}).
#' @param verbose Logical; if \code{TRUE}, print brief progress messages.
#'
#' @return
#' An object of class \code{"nicheR_ellipsoid"} with components:
#' \itemize{
#'   \item \code{dimensions}: number of variables (dimensions).
#'   \item \code{var_names}: variable names (from \code{colnames(range)}).
#'   \item \code{centroid}: centroid vector \eqn{\mu}.
#'   \item \code{cov_matrix}: covariance matrix \eqn{\Sigma}.
#'   \item \code{Sigma_inv}: inverse covariance \eqn{\Sigma^{-1}}.
#'   \item \code{chol_Sigma}: Cholesky factor of \eqn{\Sigma}.
#'   \item \code{eigen}: list with \code{vectors} and \code{values} from eigen decomposition.
#'   \item \code{cl}: confidence level used.
#'   \item \code{chi2_cutoff}: chi square cutoff \eqn{c^2}.
#'   \item \code{truncation_level}: truncation level stored (set equal to \code{cl}).
#'   \item \code{axes_sd}: derived marginal standard deviations (from \code{range}).
#'   \item \code{semi_axes_lengths}: principal semi axis lengths \eqn{a_i}.
#'   \item \code{axis_points}: list of principal axis endpoint coordinates (neg and pos).
#'   \item \code{volume}: ellipsoid hypervolume at the specified contour.
#'   \item \code{cov_limits}: output from \code{covariance_limits()}.
#' }
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

  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)

  verbose_message("Starting: building ellipsoidal niche from ranges...\n")

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


  verbose_message("Step: computing covariance matrix...\n")

  cov_matrix <- diag(sd_vec^2, nrow = length(mu_vec), ncol = length(mu_vec))
  rownames(cov_matrix) <- var_names
  colnames(cov_matrix) <- var_names
  names(mu_vec) <- var_names

  verbose_message("Done: created ellipsoidal niche.\n")

  out <- ellipsoid_calculator(cov_matrix = cov_matrix,
                              centroid = mu_vec, cl = cl,
                              verbose = verbose)
  out
}