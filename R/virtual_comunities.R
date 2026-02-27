#' Generate random ellipses constrained by a point cloud and a reference ellipse
#'
#' @description
#' Creates n random ellipses with centroids sampled from an irregular point
#' cloud. Covariance matrices are built using random rotations and scaled
#' eigenvalues restricted by user-defined limits.
#'
#' @param object A nicheR_ellipsoid object used as a reference ellipse, and
#'   containing at least \code{covariance_matrix} and \code{cl}.
#' @param background Matrix or Dataframe. The 2D point cloud (coordinates)
#'   used to select random centroids.
#' @param n Integer. Number of ellipses to generate.
#' @param smaller_proportion Numeric. Minimum scaling factor for the covariance.
#'   Must be between 0 and 1. Default = 0.1. This controls how much smaller
#'   the new ellipses can be compared to the reference.
#' @param uniform_centroids Logical. If TRUE, centroids are sampled more
#'   uniformly across the background using a grid-based thinning approach.
#'   Default = TRUE.
#' @param resolution Integer. Number of cells per side in the grid to deal with
#'   point density variation across background.
#' @param seed Integer. Random seed for reproducibility. Default = 1.
#' Set to NULL for no seeding.
#'
#' @return A list of length \code{n} with ellipse features.
#' 
#' @export

random_ellipses <- function(object,
                            background,
                            n = 10,
                            smaller_proportion = 0.1,
                            uniform_centroids = TRUE,
                            resolution = 50,
                            seed = 1) {
  
  # Input validation
  if (missing(object)) stop("Argument 'object' is required.")
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Argument 'object' must be of class 'nicheR_ellipsoid'.")
  }
  if (missing(background)) stop("Argument 'background' is required.")
  if (!is.null(seed)) set.seed(seed)
  
  # Extract reference covariance and level
  background <- as.matrix(background)
  ref_cov <- object$cov_matrix
  ref_vars <- diag(ref_cov)
  ref_level <- object$level
  
  # Sample centroids
  if (uniform_centroids) {
    ## Create a grid and assign points to grid cells
    x_range <- seq(min(background[, 1]), max(background[, 1]),
                   length.out = resolution)
    y_range <- seq(min(background[, 2]), max(background[, 2]),
                   length.out = resolution)
    grid_id <- paste(findInterval(background[, 1], x_range),
                     findInterval(background[, 2], y_range), sep = "_")
    
    thin_idx <- tapply(seq_len(nrow(background)), grid_id, function(x) {
      if(length(x) == 1) return(x) else return(sample(x, 1))
    })
    
    idx <- sample(thin_idx, n, replace = TRUE)
  } else {
    idx <- sample(seq_len(nrow(background)), n, replace = TRUE)
  }
  
  ## Extract centroids for the selected indices
  centroids <- background[idx, , drop = FALSE]
  
  # Generate Ellipses
  rand_ellipses <- lapply(1:n, function(i) {
    ## Randomly scale variances between smaller_proportion and 1.0 of original
    ## This keeps the new variances within the "largest possible" limit
    v_scales <- runif(2, min = smaller_proportion, max = 1)
    new_vars <- ref_vars * v_scales
    sds <- sqrt(new_vars)
    
    ## Calculate covariance limits based on new variances
    max_cov <- sds[1] * sds[2]
    
    ## Pick a random covariance
    ## To avoid being too close to the edge use a multiplier of 0.9
    new_cov_val <- runif(1, min = -max_cov * 0.9, max = max_cov * 0.9)
    
    ## Reconstruct the variance-covariance matrix
    new_varcov <- matrix(c(new_vars[1], new_cov_val,
                           new_cov_val, new_vars[2]), nrow = 2)
    
    ## Centroid
    cent <- centroids[i, ]

    ## Name covariance matrix and centroid for downstream use
    colnames(new_varcov) <- rownames(new_varcov) <- colnames(ref_cov)
    names(cent) <- colnames(ref_cov)
    
    ellipsoid_calculator(cov_matrix = new_varcov, centroid = cent,
                         cl = ref_level)
  })
  
  return(rand_ellipses)
}



#' Generate nested ellipses based on a reference ellipse
#'
#' @description
#' Creates a sequence of nested ellipses by scaling the covariance matrix of a
#' reference ellipse. The distribution of the nested ellipses can be controlled
#' using a bias exponent to cluster them toward the border or the centroid.
#'
#' @param object An object of class "nicheR_ellipsoid"
#'    describing an initial ellipse. Must contain \code{centroid},
#'    \code{cov_matrix}, and \code{cl}.
#' @param n Integer. Number of nested ellipses to generate. Default is 10.
#' @param smaller_proportion Numeric scalar in \code{(0, 1)}. The scale of the
#'    smallest ellipse relative to the original. Default is \code{0.1}.
#' @param bias Numeric. An exponent controlling the spacing of the nested
#'    ellipses.
#'    \itemize{
#'      \item \code{bias = 1}: Linear spacing (default).
#'      \item \code{0 < bias < 1}: Clusters ellipses toward the **border**
#'      (outer original ellipse).
#'      \item \code{bias > 1}: Clusters ellipses toward the **centroid**
#'      (inner smallest ellipse).
#'    }
#'
#' @return
#' A list of \code{n} "nicheR_ellipsoid" objects. The largest ellipse
#' corresponds to the original reference, and the smallest is that scaled by
#' \code{smaller_proportion}.
#'
#' @details
#' The function generates a sequence of scale factors $k$ using the formula:
#' $k_i = \text{smaller\_proportion} + (1 - \text{smaller\_proportion}) \ times t_i^{\text{bias}}$,
#' where $t_i$ is a linear sequence from 1 down to 0.
#' 
#' @export

nested_ellipses <- function(object,
                            n = 10,
                            smaller_proportion = 0.1,
                            bias = 1) {
  # Basic checks
  if (missing(object)) stop("Argument 'object' is required.")
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Argument 'object' must be of class 'nicheR_ellipsoid'.")
  }
  if (is.null(object$centroid) ||
      is.null(object$cov_matrix) ||
      is.null(object$cl)) {
    stop("'object' must contain 'centroid', 'cov_matrix', and 'cl'.")
  }
  if (!is.numeric(bias) || length(bias) != 1L || bias <= 0) {
    stop("'bias' must be a single numeric value > 0.")
  }
  
  # Generate biased scale factors
  # t goes from 1 to 0. Raising it to 'bias' shifts the distribution
  t <- seq(1, 0, length.out = n)
  scales <- smaller_proportion + (1 - smaller_proportion) * (t^bias)
  
  # Generate nested ellipses
  ell_list <- lapply(scales, function(k) {
    scaled_cov <- k * object$cov_matrix
    ellipsoid_calculator(cov_matrix = scaled_cov,
                         centroid = object$centroid,
                         cl = object$cl)
  })
  
  return(ell_list)
}