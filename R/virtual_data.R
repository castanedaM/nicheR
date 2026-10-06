#' Generate data based on a ellipsoidal niche
#'
#' @description
#' Simulates \code{n} random points from a multivariate normal distribution
#' defined by the centroid and covariance matrix of a \code{nicheR_ellipsoid}
#' object.
#'
#' @usage
#' virtual_data(object, n = 100, truncate = FALSE, effect = "direct", seed = 1)
#'
#' @param object A \code{nicheR_ellipsoid} object containing at least
#'   \code{centroid} and \code{cov_matrix}.
#' @param n Integer. The number of virtual points to generate. Default = 100.
#' @param truncate Logical. If \code{TRUE}, points are constrained
#'   within the confidence limit (\code{cl}) defined in the object.
#' @param effect Character. The distribution pattern of points.
#'   \code{"direct"} (default) creates a concentration near the centroid.
#'   \code{"inverse"} creates higher density towards the edges.
#'   \code{"uniform"} distributes points evenly throughout the ellipsoid volume.
#'   Note: \code{"inverse"} and \code{"uniform"} require \code{truncate = TRUE}.
#' @param seed Integer. Random seed for reproducibility. Default = 1.
#'   Set to \code{NULL} for no seeding.
#'
#' @details
#' When \code{truncate = FALSE}, the function generates points from a standard
#' multivariate normal distribution defined by the ellipsoid's centroid and
#' covariance matrix, without any constraints on their location. The function
#' uses eigen-decomposition to transform standard normal variables into the
#' coordinate system defined by the ellipsoid's covariance structure.
#'
#' When \code{truncate = TRUE}, every point falls inside the ellipsoid, that
#' is, where the squared Mahalanobis distance \eqn{Md \le} \code{chi2_cutoff}.
#' How the points are distributed inside it depends on the \code{effect}
#' argument:
#' \itemize{
#'   \item \code{"direct"}: Points follow the multivariate normal density
#'   (\eqn{\exp(-0.5 \times Md)}) truncated at the ellipsoid boundary, which
#'   clusters them near the centroid. They are drawn from the normal
#'   distribution itself, and the draws that fall outside are discarded.
#'   \item \code{"inverse"}: Points follow the complement of the normal density
#'   (\eqn{1 - \exp(-0.5 \times Md)}), which pushes them toward the edges.
#'   Candidates are drawn uniformly inside the ellipsoid and each one is kept
#'   with a probability proportional to that weight.
#'   \item \code{"uniform"}: All locations within the ellipsoid are equally
#'   likely, resulting in a uniform distribution through its volume.
#' }
#'
#' For \code{"inverse"} and \code{"uniform"}, candidates are drawn uniformly
#' within a bounding box around the ellipsoid, the centroid plus and minus
#' \eqn{\sqrt{diag(\Sigma) \times}} \code{chi2_cutoff}, and those outside
#' the ellipsoid are removed.
#'
#' @return
#' A matrix with \code{n} rows and columns corresponding to the
#' environmental variables (dimensions) of the input \code{object}.
#'
#' @importFrom stats runif rnorm pchisq
#'
#' @examples
#' # Loading data
#' ## Reference niche
#' data("ref_ellipse", package = "nicheR")
#'
#' # Generate virtual data from the reference niche
#' vdata_direct <- virtual_data(ref_ellipse, n = 100, effect = "direct")
#' vdata_inverse <- virtual_data(ref_ellipse, n = 100,
#'                               effect = "inverse", truncate = TRUE)
#'
#' # Check a sample of the generated data
#' head(vdata_direct)
#' head(vdata_inverse)
#' @export
virtual_data <- function(object,
                         n = 100,
                         truncate = FALSE,
                         effect = "direct",
                         seed = 1) {
  # Detecting potential errors
  if (missing(object)) {
    stop("Argument 'object' is required.")
  }
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Argument 'object' must be of class 'nicheR_ellipsoid'.")
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("Argument 'n' must be an integer > 0.")
  }
  if (!is.logical(truncate)) {
    stop("Argument 'truncate' must be a logical.")
  }
  if (!is.character(effect)) {
    stop("Argument 'effect' must be a character.")
  }
  if (!effect %in% c("direct", "inverse", "uniform")) {
    stop("Argument 'effect' must be 'direct', 'inverse', or 'uniform'.")
  }
  if (effect %in% c("inverse", "uniform") && !truncate) {
    stop("Effect 'inverse' and 'uniform' only possible when 'truncate = TRUE'.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Ellipsoid features
  centroid <- object$centroid
  cov_matrix <- object$cov_matrix
  p <- length(centroid)

  # Eigen-decomposition of the covariance matrix
  es <- object$eigen
  ev <- object$eigen$values

  if (truncate) {
    # Chi-squared cutoff that defines the ellipsoid boundary
    conf_cutoff <- object$chi2_cutoff

    final_points <- matrix(nrow = 0, ncol = p)

    if (effect == "direct") {
      # Truncated multivariate normal. Points are drawn from the normal
      # itself and those outside the ellipsoid are discarded, so what is kept
      # follows the normal density inside the boundary in any number of
      # dimensions.

      ## Square root of Sigma (V * L^0.5), as in the untruncated case
      to_cov <- es$vectors %*% diag(sqrt(pmax(ev, 0)), p)

      ## Share of normal draws expected to fall inside the ellipsoid
      p_inside <- pchisq(conf_cutoff, df = p)

      while (nrow(final_points) < n) {
        ## Batch size, with a margin so one pass is usually enough
        batch_size <- ceiling((n - nrow(final_points)) /
                                max(p_inside, 0.05) * 1.1) + 10L

        ## Standard normal draws
        z <- matrix(rnorm(p * batch_size), nrow = batch_size)

        ## For a standard normal draw, the squared Mahalanobis distance of
        ## the transformed point is the squared length of the draw
        inside <- rowSums(z^2) <= conf_cutoff

        ### Safety check: if no points are inside, skip to next iteration
        if (!any(inside)) next

        ## Transform to the ellipsoid and add to our collection
        pts <- sweep(z[inside, , drop = FALSE] %*% t(to_cov), 2L,
                     centroid, "+")
        final_points <- rbind(final_points, pts)
      }

    } else {
      # Uniform candidates inside the ellipsoid, for "uniform" and "inverse"

      ## Get the range for each variable across all axes
      half <- sqrt(diag(cov_matrix) * conf_cutoff)
      v_min <- centroid - half
      v_max <- centroid + half

      inv_cov <- object$Sigma_inv

      ## Largest weight "inverse" can take, reached on the boundary
      w_max <- 1 - exp(-0.5 * conf_cutoff)

      while (nrow(final_points) < n) {
        ## Batch size
        batch_size <- ceiling((n - nrow(final_points)) / 0.1)

        ## Generate uniform random points in the space defined by the axes
        v_raw_cube <- mapply(runif, n = rep(batch_size, p),
                             min = v_min, max = v_max)

        ## Mahalanobis distance to points
        diffs <- sweep(v_raw_cube, 2L, centroid, "-")
        d2 <- rowSums((diffs %*% inv_cov) * diffs)

        ## Get points and distances within the confidence cutoff
        inside <- d2 <= conf_cutoff

        ### Safety check: if no points are inside, skip to next iteration
        if (!any(inside)) next

        v_raw_cube <- v_raw_cube[inside, , drop = FALSE]
        d2 <- d2[inside]

        if (effect == "uniform") {
          ## Every candidate inside has the same weight. The draw is kept as
          ## it was, so a given seed still returns the same points.
          weights <- rep(1, nrow(v_raw_cube))

          keep <- sample(seq_len(nrow(v_raw_cube)),
                         size = min(nrow(v_raw_cube), n - nrow(final_points)),
                         prob = weights, replace = FALSE)

        } else {
          ## Inverse. Each candidate is kept with probability equal to its
          ## weight over the largest weight, so the kept points follow the
          ## weight itself. Picking a fixed number from a small pool does
          ## not, since it keeps nearly every candidate whatever its weight.
          weights <- 1 - exp(-0.5 * d2)

          keep <- which(runif(nrow(v_raw_cube)) < weights / w_max)
        }

        ## Add the sampled points to our collection
        final_points <- rbind(final_points, v_raw_cube[keep, , drop = FALSE])
      }
    }

    # Trim to exactly n and ensure names
    final_points <- final_points[1:n, , drop = FALSE]
    colnames(final_points) <- colnames(cov_matrix)
    return(final_points)

  } else {
    # Generate standard normal samples
    vdata <- matrix(rnorm(p * n), nrow = n)

    # Transform using the Square Root of Sigma (V * L^0.5)
    vdata <- drop(centroid) + es$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
      t(vdata)

    # Handle dimension names
    rownames(vdata) <- colnames(cov_matrix)

    # Return the generated data
    return(t(vdata))
  }
}
