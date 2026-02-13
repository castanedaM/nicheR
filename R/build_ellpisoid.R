#' Build a Probabilistic Ellipsoidal Niche from Ranges
#'
#' Constructs an ellipsoidal niche in multivariate environmental space using a
#' multivariate normal (MVN) contour of constant Mahalanobis distance. The niche
#' is defined by a centroid (midpoints of user-provided ranges), per-variable
#' marginal spreads (derived from ranges), and an optional user-supplied
#' variance–covariance matrix.
#'
#' @details
#' This function expects `range` as a 2-row matrix/data.frame with variables in
#' columns. The two rows represent lower and upper bounds (order can be min/max
#' or max/min). The centroid is computed as the midpoint of each variable's
#' range:
#' \deqn{\mu_i = (m_i + M_i)/2.}
#'
#' If `cov_matrix` is not provided, marginal standard deviations are derived from
#' ranges using the rule-of-thumb that the provided bounds represent roughly
#' \eqn{\pm 3} standard deviations:
#' \deqn{\sigma_i = (M_i - m_i)/6,}
#' and a diagonal covariance matrix is assumed:
#' \deqn{\Sigma = \mathrm{diag}(\sigma_1^2,\dots,\sigma_n^2).}
#'
#' The ellipsoid boundary is defined by:
#' \deqn{(x - \mu)^\top \Sigma^{-1} (x - \mu) = c^2,}
#' where \eqn{c^2} is the chi-square cutoff \eqn{\chi^2_{n}(\mathrm{cl})} with
#' \eqn{n} degrees of freedom.
#'
#' Principal semi-axis lengths are computed from the eigenvalues of \eqn{\Sigma}:
#' \deqn{a_i = \sqrt{\lambda_i c^2}.}
#'
#' @param range A 2-row `matrix` or `data.frame` of bounds, with variables as
#'   columns. Rows may be ordered as min/max or max/min. Column names must be
#'   provided (variable names).
#' @param cl Numeric confidence level in (0, 1). Used to set the chi-square
#'   cutoff that determines the ellipsoid contour (not per-variable limits).
#' @param cov_matrix Optional variance–covariance matrix \eqn{\Sigma} for the
#'   ellipsoid (must be `dimensions x dimensions`, symmetric positive definite).
#'   If provided, it is used directly (no rescaling from `range`).
#' @param verbose Logical. If `TRUE`, prints brief progress messages.
#'
#' @return An object of class `"nicheR_ellipsoid"` with components:
#' \itemize{
#'   \item `dimensions`: number of variables (dimensions).
#'   \item `mu_vec`: centroid vector \eqn{\mu}.
#'   \item `cov_matrix`: covariance matrix \eqn{\Sigma}.
#'   \item `Sigma_inv`: inverse covariance \eqn{\Sigma^{-1}}.
#'   \item `chol_Sigma`: Cholesky factor of \eqn{\Sigma}.
#'   \item `eigen`: list with `vectors` and `values` from eigen-decomposition.
#'   \item `cl`: confidence level used.
#'   \item `chi2_cutoff`: chi-square cutoff \eqn{c^2}.
#'   \item `axes_sd`: derived marginal standard deviations (from `range`).
#'   \item `semi_axes_lengths`: principal semi-axis lengths \eqn{a_i}.
#'   \item `axis_points`: list of principal-axis endpoint coordinates (LB/UB).
#'   \item `volume`: ellipsoid hypervolume at the specified contour.
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
                            cov_matrix = NULL,
                            verbose = TRUE){


# TO DO: add the limit/threshold for trucations, NA mahalanobis, 0 in
# suitability. Do not binaries, everything outside the ellipsoid NA. Add the
# mahalanobis value and suitability correspond to the limit.

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

  if(!is.null(cov_matrix)){

    if(is.data.frame(cov_matrix)){
      cov_matrix <- as.matrix(cov_matrix)
    }

    if(!is.matrix(cov_matrix)){
      stop("cov_matrix must be a matrix (variance–covariance matrix).")
    }

    if(is.null(colnames(cov_matrix)) || is.null(rownames(cov_matrix))){
      stop("cov_matrix must have row and column names matching range column names.")
    }

    if(!all(colnames(cov_matrix) == var_names) ||
       !all(rownames(cov_matrix) == var_names)){
      stop("cov_matrix row/column names must match range column names (same order).")
    }
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

  dimensions <- length(mu_vec)

  # Covariance handling ------------------------------------------------------

  if(is.null(cov_matrix)){

    verbose_message("Step: computing covariance matrix...\n")

    cov_matrix <- diag(sd_vec^2, nrow = dimensions, ncol = dimensions)
    rownames(cov_matrix) <- var_names
    colnames(cov_matrix) <- var_names

  }else{

    verbose_message("Step: Checking covariance matrix...\n")

    if(any(dim(cov_matrix) != c(dimensions, dimensions))){
      stop("cov_matrix must be a square matrix with dimensions x dimensions.")
    }

    if(any(!is.finite(cov_matrix))){
      stop("cov_matrix contains non-finite values.")
    }

    if(!isTRUE(all.equal(cov_matrix, t(cov_matrix), tolerance = 1e-6))){
      stop("cov_matrix must be symmetric. You can use (cov_matrix + t(cov_matrix))/2 to make symmetric.")
    }
  }

  chol_Sigma <- tryCatch(chol(cov_matrix), error = function(e) NULL)

  if(is.null(chol_Sigma)){
    stop("cov_matrix must be symmetric positive definite (SPD).")
  }

  verbose_message("Step: computing safe covariance limits... see out$cov_limits for options\n")

  cov_limits <- covariance_limits(cov_matrix)
  Sigma_inv <- chol2inv(chol_Sigma)

  # Ellipsoid metrics --------------------------------------------------------

  verbose_message("Step: computing ellipsoid metrics...\n")

  chi2_cutoff <- stats::qchisq(cl, df = dimensions)

  eig <- eigen(cov_matrix, symmetric = TRUE)

  vals <- pmax(eig$values, 0)
  semi_axes_lengths <- sqrt(vals * chi2_cutoff)

  axis_points <- lapply(seq_len(dimensions), function(i){
    list(
      neg = mu_vec - semi_axes_lengths[i] * eig$vectors[, i],
      pos = mu_vec + semi_axes_lengths[i] * eig$vectors[, i]
    )
  })

  volume <- ellipsoid_volume(
    n_dimensions = dimensions,
    semi_axes_lengths = semi_axes_lengths
  )

  verbose_message("Done: created ellipsoidal niche.\n")

  out <- list(
    dimensions = dimensions,
    var_names = var_names,
    centroid = mu_vec,
    cov_matrix = cov_matrix,
    Sigma_inv = Sigma_inv,
    chol_Sigma = chol_Sigma,
    eigen = list(vectors = eig$vectors, values = eig$values),
    cl = cl,
    chi2_cutoff = chi2_cutoff,
    axes_sd = sd_vec,
    semi_axes_lengths = as.numeric(semi_axes_lengths),
    axis_points = axis_points,
    volume = volume,
    cov_limits = cov_limits
  )

  class(out) <- "nicheR_ellipsoid"

  out
}

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
print.nicheR_ellipsoid <- function(x, digits = 3, ...){

  cat("nicheR Ellipsoid Object\n")
  cat("----------------------\n")

  cat("Dimensions:        ", x$dimensions, "D\n", sep = "")
  cat("Chi-square cutoff: ", round(x$chi2_cutoff, digits), "\n", sep = "")

  cat("Centroid (mu):     ",
      paste(round(x$centroid, digits), collapse = ", "),
      "\n", sep = "")

  cat("\nCovariance matrix:\n")
  print(round(x$cov_matrix, digits))

  cat("\nSemi-axis lengths:\n  ",
      paste(round(x$semi_axes_lengths, digits), collapse = ", "),
      "\n", sep = "")

  cat("\nPrincipal axis endpoints (neg → pos):\n")

  for(i in seq_len(x$dimensions)){
    cat(" Axis", i, ":\n", sep = "")
    cat("   neg: ",
        paste(round(x$axis_points[[i]]$neg, digits), collapse = ", "),
        "\n", sep = "")
    cat("   pos: ",
        paste(round(x$axis_points[[i]]$pos, digits), collapse = ", "),
        "\n", sep = "")
  }

  cat("\nEllipsoid volume:\n")
  print(round(x$volume, digits))

  invisible(x)
}

