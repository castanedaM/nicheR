#' Generate random ellipses constrained by a point cloud and a reference ellipse
#'
#' @description
#' Creates n random ellipses with centroids sampled from an irregular point
#' cloud. Covariance matrices are built using random rotations and scaled
#' eigenvalues restricted by user-defined limits.
#'
#' @param object A nicheR_ellipsoid object used as a reference ellipse
#'   (the biggest to be generated), and containing at least
#'   \code{covariance_matrix} and \code{cl}.
#' @param background Matrix or Dataframe. The 2D point cloud (coordinates)
#'   used to select random centroids.
#' @param n Integer. Number of ellipses to generate.
#' @param smallest_proportion Numeric. Minimum scaling factor for the variance.
#'   Must be between 0 and 1. Default = 0.1. This controls how much smaller
#'   the new ellipses can be compared to the reference.
#' @param largest_proportion Numeric. Maximum scaling factor for the variance
#'   relative to the reference. Default = 1.0. This controls how much larger
#'   the new ellipses can be compared to the reference. Values larger than 1
#'   will result in ellipses that exceed the reference size.
#' @param thin_background Logical. If TRUE, centroids are sampled more
#'   uniformly across the background using a grid-based thinning approach.
#'   Default = FALSE.
#' @param resolution Integer. Number of cells per side in the grid to deal with
#'   point density variation across background. Default = 50.
#' @param seed Integer. Random seed for reproducibility. Default = 1.
#' Set to NULL for no seeding.
#'
#' @return
#' An object of class \code{nicheR_community} containing the generated
#' ellipses, the reference object, and generation metadata.
#'
#' @export

random_ellipses <- function(object,
                            background,
                            n = 10,
                            smallest_proportion = 0.1,
                            largest_proportion = 1.0,
                            thin_background = FALSE,
                            resolution = 50,
                            seed = 1) {
  
  # Input validation
  if (missing(object)) {
    stop("Argument 'object' is required.")
  }
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Argument 'object' must be of class 'nicheR_ellipsoid'.")
  }
  if (missing(background)) {
    stop("Argument 'background' is required.")
  }
  if (!identical(colnames(background), colnames(object$cov_matrix))) {
    stop("'background' column names must match those of 'object$cov_matrix'.")
  }
  if (smallest_proportion <= 0 || smallest_proportion >= 1) {
    stop("'smallest_proportion' must be a single numeric value, >0, <1.")
  }
  if (largest_proportion < smallest_proportion) {
    stop("'largest_proportion' must be a single value > 'smallest_proportion'.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Extract reference covariance and level
  background <- as.matrix(background)
  ref_cov <- object$cov_matrix
  ref_vars <- diag(ref_cov)
  ref_level <- object$cl
  
  # Sample centroids
  if (thin_background) {
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
    ## Randomly scale variances between smallest_proportion and 1.0 of original
    ## This keeps the new variances within the "largest possible" limit
    v_scales <- runif(2, min = smallest_proportion, max = largest_proportion)
    new_vars <- ref_vars * v_scales
    sds <- sqrt(new_vars)
    
    ## Pick a random covariance
    ## To avoid being too close to the edge use a multiplier of 0.9
    rho <- runif(1, min = -0.9, max = 0.9)
    new_cov_val <- rho * sds[1] * sds[2]
    
    ## Reconstruct the variance-covariance matrix
    new_varcov <- matrix(c(new_vars[1], new_cov_val,
                           new_cov_val, new_vars[2]), nrow = 2)
    
    ## Centroid
    cent <- centroids[i, ]

    ## Name covariance matrix and centroid for downstream use
    colnames(new_varcov) <- rownames(new_varcov) <- colnames(ref_cov)
    names(cent) <- colnames(ref_cov)
    
    ellipsoid_calculator(cov_matrix = new_varcov, centroid = cent,
                         cl = ref_level, verbose = FALSE)
  })

  return(new_nicheR_community(
    ellipse_community = rand_ellipses,
    reference = object,
    pattern = "random",
    n = n,
    smallest_proportion = smallest_proportion,
    largest_proportion = largest_proportion,
    thin_background = thin_background,
    resolution = resolution,
    seed = seed
  ))
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
#' @param smallest_proportion Numeric scalar in \code{(0, 1)}. The scale of the
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
#' An object of class \code{nicheR_community} containing the generated
#' ellipses, the reference object, and generation metadata.
#'
#' @details
#' The largest ellipse corresponds to the original reference, and the
#' smallest is that scaled by \code{smallest_proportion}.
#'
#' The function generates a sequence of scale factors $k$ using the formula:
#' $k_i = \text{smaller\_proportion} + (1 - \text{smaller\_proportion}) \ times t_i^{\text{bias}}$,
#' where $t_i$ is a linear sequence from 1 down to 0.#' A list of \code{n} "nicheR_ellipsoid" objects.
#' 
#' @export

nested_ellipses <- function(object,
                            n = 10,
                            smallest_proportion = 0.1,
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
  scales <- smallest_proportion + (1 - smallest_proportion) * (t^bias)
  
  # Generate nested ellipses
  ell_list <- lapply(scales, function(k) {
    scaled_cov <- (k^2) * object$cov_matrix
    ellipsoid_calculator(cov_matrix = scaled_cov,
                         centroid = object$centroid,
                         cl = object$cl, verbose = FALSE)
  })
  
  return(new_nicheR_community(
    ellipse_community = ell_list,
    reference = object,
    pattern = "nested",
    n = n,
    smallest_proportion = smallest_proportion,
    bias = bias
  ))
}




#' Generate ellipses via multivariate normal biased sampling
#'
#' @description
#' Creates a set of ellipses with centroids sampled from a background, biased
#' by their proximity to the centroid to a reference niche. Includes an option
#' to thin the background to reduce centroid sampling bias due to point-density.
#'
#' @param object A nicheR_ellipsoid object used as the reference. This is will
#'    be considered the "largest" ellipse to be generated.
#' @param background Matrix or Dataframe. The 2D point cloud (coordinates)
#'    used to select centroids for the ellipses.
#' @param n Integer. Number of ellipses to generate. Default = 10.
#' @param smallest_proportion Numeric scalar in \code{(0, 1)}. The scale of the
#'    smallest ellipse relative to the original. Default is \code{0.1}.
#' @param largest_proportion Numeric. Maximum scaling factor for the variance
#'   relative to the reference. Default = 1.0. This controls how much larger
#'   the new ellipses can be compared to the reference. Values larger than 1
#'   will result in ellipses that exceed the reference size.
#' @param thin_background Logical. If TRUE, centroids are sampled more
#'   uniformly across the background using a grid-based thinning approach.
#'   Default = FALSE.
#' @param resolution Integer. Number of cells per side in the grid to deal with
#'   point density variation across background. Default = 100.
#' @param seed Integer. Random seed for reproducibility. Default = 1.
#' Set to NULL for no seeding.
#'
#' @return
#' An object of class \code{nicheR_community} containing the generated
#' ellipses, the reference object, and generation metadata.
#'
#' @details
#' Ellipses are generated to simulate a community of niches with varying
#' degrees of similarity to the reference. The distribution of the generated
#' ellipses is influenced by the proximity to the reference and the density
#' of the background points.
#'
#' @export

conserved_ellipses <- function(object,
                               background,
                               n = 10,
                               smallest_proportion = 0.1,
                               largest_proportion = 1.0,
                               thin_background = FALSE,
                               resolution = 100,
                               seed = 1) {

  # Input validation
  if (missing(object)) {
    stop("Argument 'object' is required.")
  }
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Argument 'object' must be of class 'nicheR_ellipsoid'.")
  }
  if (missing(background)) {
    stop("Argument 'background' is required.")
  }
  if (!identical(colnames(background), colnames(object$cov_matrix))) {
    stop("'background' column names must match those of 'object$cov_matrix'.")
  }
  if (smallest_proportion <= 0 || smallest_proportion >= 1) {
    stop("'smallest_proportion' must be a single numeric value, >0, <1.")
  }
  if (largest_proportion < smallest_proportion) {
    stop("'largest_proportion' must be a single value > 'smallest_proportion'.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Extract reference metrics
  background <- as.matrix(background)
  ref_cov  <- object$cov_matrix
  ref_vars <- diag(ref_cov)
  ref_level <- object$cl
  
  # Background boundaries for reflection logic
  bg_mvnd <- predict.nicheR_ellipsoid(object, background, verbose = FALSE)
  bg_prob <- bg_mvnd$suitability  # this is MVND at each background point

  # Handle Background Thinning if asked
  if (thin_background) {
    ## Create grid indices
    x_range <- seq(min(background[, 1]), max(background[, 1]),
                   length.out = resolution)
    y_range <- seq(min(background[, 2]), max(background[, 2]),
                   length.out = resolution)
    
    grid_id <- paste(findInterval(background[, 1], x_range),
                     findInterval(background[, 2], y_range), sep = "_")
    
    ## Select one representative index per grid cell
    thin_idx <- tapply(seq_len(nrow(background)), grid_id, function(x) {
      if(length(x) == 1) return(x) else return(sample(x, 1))
    })
    
    ## Subset both background and the associated probabilities
    background <- background[thin_idx, , drop = FALSE]
    bg_prob <- bg_prob[thin_idx]
  }
  
  # Sample centroids with bias toward the reference centroid
  if (sum(bg_prob) > 0) {
    cent_samp <- sample(seq_len(nrow(background)), size = n,
                        prob = bg_prob, replace = TRUE)
    centroids <- background[cent_samp, , drop = FALSE]
  } else {
    stop("'background' has points too far from the reference in 'object'.")
  }
  
  # Generate Ellipses
  results <- lapply(1:n, function(i) {
    ## Randomly scale variances between smallest_proportion and 1.0 of original
    ## This keeps the new variances within the "largest possible" limit
    v_scales <- runif(2, min = smallest_proportion, max = largest_proportion)
    new_vars <- ref_vars * v_scales
    sds <- sqrt(new_vars)
    
    ## Pick a random covariance
    ## To avoid being too close to the edge use a multiplier of 0.9
    rho <- runif(1, min = -0.9, max = 0.9)
    new_cov_val <- rho * sds[1] * sds[2]
    
    ## Reconstruct the variance-covariance matrix
    new_varcov <- matrix(c(new_vars[1], new_cov_val,
                           new_cov_val, new_vars[2]), nrow = 2)

    ## Centroid
    cent <- centroids[i, ]

    ## Name covariance matrix and centroid for downstream use
    colnames(new_varcov) <- rownames(new_varcov) <- colnames(ref_cov)
    names(cent) <- colnames(ref_cov)

    ellipsoid_calculator(cov_matrix = new_varcov, centroid = cent,
                         cl = ref_level, verbose = FALSE)
  })

  return(new_nicheR_community(
    ellipse_community = results,
    reference = object,
    pattern = "conserved",
    n = n,
    smallest_proportion = smallest_proportion,
    largest_proportion = largest_proportion,
    thin_background = thin_background,
    resolution = resolution,
    seed = seed
  ))
}