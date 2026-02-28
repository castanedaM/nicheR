#' Generate data based on a ellipsoidal niche
#'
#' @description
#' Simulates \code{n} random points from a multivariate normal distribution
#' defined by the centroid and covariance matrix of a \code{nicheR_ellipsoid}
#' object.
#'
#' @param object A \code{nicheR_ellipsoid} object containing at least
#'   \code{centroid} and \code{cov_matrix}.
#' @param n Integer. The number of virtual points to generate. Default = 100.
#' @param truncate Logical. If \code{TRUE} (default), points are constrained
#'   within the confidence limit (\code{cl}) defined in the object.
#' @param effect Character. The distribution pattern of points.
#'   \code{"direct"} (default) creates a concentration near the centroid.
#'   \code{"inverse"} creates higher density towards the edges.
#'   \code{"uniform"} distributes points evenly throughout the ellipsoid volume.
#'   Note: \code{"inverse"} and \code{"uniform"} require \code{truncate = TRUE}.
#' @param seed Integer. Random seed for reproducibility. Default = 1.
#'   Set to NULL for no seeding.
#'
#' @details
#' When \code{truncate = FALSE}, the function generates points from a standard
#' multivariate normal distribution defined by the ellipsoid's centroid and
#' covariance matrix, without any constraints on their location. The function
#' uses eigen-decomposition to transform standard normal variables into the
#' coordinate system defined by the ellipsoid's covariance structure.
#'
#' When \code{truncate = TRUE}, the function generates candidate points
#' uniformly distributed within a bounding box (hyper-cube) defined by the
#' ellipsoid's \code{axes_coordinates}. Points falling outside the ellipsoid
#' (where Mahalanobis distance $Md >$ \code{chi2_cutoff}) are removed.
#'
#' From this filtered pool, \code{n} points are selected using weighted random
#' sampling without replacement. The weights are determined by the \code{effect}
#' argument:
#' \itemize{
#'   \item \code{"direct"}: Weights are proportional to the multivariate normal
#'   density ($\exp(-0.5 \cdot Md)$), clustering points near the centroid.
#'   \item \code{"inverse"}: Weights are proportional to the complement of the
#'   normal density ($1 - \exp(-0.5 \cdot Md)$), pushing points toward the edges.
#'   \item \code{"uniform"}: All points within the ellipsoid have equal weight,
#'   resulting in a uniform spatial distribution.
#' }
#'
#' @return
#' A matrix with \code{n} rows and columns corresponding to the
#' environmental variables (dimensions) of the input \code{object}.
#'
#' @export

virtual_data <- function(object,
                         n = 100,
                         truncate = TRUE,
                         effect = "direct",
                         seed = 1) {
  # Detecting potential errors
  if (missing(object)) {
    stop("Argument 'object' is required.")
  }
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Argument 'object' must be of class 'nicheR_ellipsoid'.")
  }
  if (!is.numeric(n) && n >= 1) {
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
    # Get the range for each variable across all axes
    all_coords <- do.call(rbind, object$axes_coordinates)
    v_min <- apply(all_coords, 2, min)
    v_max <- apply(all_coords, 2, max)

    # Truncating with Mahalanobis distance and chi-squared cutoff
    conf_cutoff <- object$chi2_cutoff
    inv_cov <- object$Sigma_inv
    
    final_points <- matrix(nrow = 0, ncol = p)
    
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

      ## Multivariate normal from the Mahalanobis distance
      mvnd <- exp(-0.5 * d2)
      
      if (effect == "direct") {
        weights <- mvnd
      } else if (effect == "inverse") {
        weights <- 1 - mvnd
      } else {
        weights <- rep(1, nrow(v_raw_cube)) # Uniform
      }
      
      ## Sample points based on weights, ensuring we don't exceed n
      keep <- sample(seq_len(nrow(v_raw_cube)),
                     size = min(nrow(v_raw_cube), n - nrow(final_points)),
                     prob = weights, replace = FALSE)
      
      ## Add the sampled points to our collection
      final_points <- rbind(final_points, v_raw_cube[keep, , drop = FALSE])
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