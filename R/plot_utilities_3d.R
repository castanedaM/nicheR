#' Plot a nicheR ellipsoid in 3D environmental space
#'
#' @description
#' Creates an interactive 3D plot of a \code{nicheR_ellipsoid} object with
#' support for background points or suitability surfaces.
#'
#' @usage
#' plot_ellipsoid_3d(object, dim = c(1, 2, 3), wire = FALSE, aspect = TRUE,
#'                   background = NULL, prediction = NULL, col_layer = NULL,
#'                   pal = hcl.colors(100, palette = "Viridis"),
#'                   rev_pal = FALSE, bg_sample = NULL, col_ell = "#8b0000",
#'                   alpha_ell = 1, alpha_bg = 1, col_bg = "#8A8A8A",
#'                   fixed_lims = NULL, xlab = NULL, ylab = NULL,
#'                   zlab = NULL, ...)
#'
#' @param object A \code{nicheR_ellipsoid} object constructed with at least
#'    3 dimensions.
#' @param dim Integer vector of length 3. Indices of dimensions to plot.
#' @param wire Logical. If \code{TRUE}, plots wireframe. The default,
#'    \code{FALSE}, plots a shaded volume.
#' @param aspect Logical. If \code{TRUE}, maintains aspect ratio (1:1:1).
#' @param background Optional data frame/matrix of background points.
#'    This argument has priority over \code{prediction} if both are provided.
#' @param prediction Optional data frame/matrix for prediction surfaces.
#' @param col_layer Character or \code{NULL}. Column in \code{prediction}
#'    to use for coloring.
#' @param pal Color palette function or character vector.
#' @param rev_pal Logical. If \code{TRUE}, reverses the palette.
#' @param bg_sample Integer or \code{NULL}. Subsample size for large data.
#' @param col_ell Color of the ellipsoid. Default is \code{"#000000"}.
#' @param alpha_ell Transparency of the ellipsoid. Default is \code{1}.
#'    Not applied to wireframe mode.
#' @param alpha_bg Transparency of background points. Default is \code{1}.
#'    Also applied to prediction points.
#' @param col_bg Color for background points.
#' @param fixed_lims Named list with \code{xlim}, \code{ylim}, and \code{zlim}.
#' @param xlab x-axis label. The default, \code{NULL}, uses ellipsoid object
#'    variable names, if any found.
#' @param ylab y-axis label. The default, \code{NULL}, uses ellipsoid object
#'    variable names, if any found.
#' @param zlab z-axis label. The default, \code{NULL}, uses ellipsoid object
#'    variable names, if any found.
#' @param ... Additional graphical parameters.
#'
#' @return A 3d plot of the ellipsoid in E-space, shown in viewer.
#'
#' @importFrom grDevices hcl.colors adjustcolor
#'
#' @export
#'
#' @examples
#' # Building an ellipsoid
#' ## Define ranges for three variables
#' range <- data.frame(bio_1  = c(22, 32),
#'                     bio_12 = c(800, 4200),
#'                     bio_15 = c(45, 115))
#'
#' ## Build the ellipsoid
#' ell5 <- build_ellipsoid(ranges = range)
#' ell5$cov_limits
#'
#' ell5 <- update_ellipsoid_covariance(ell5, c("bio_1-bio_12" = 200,
#'                                             "bio_1-bio_15" = 0,
#'                                             "bio_12-bio_15" = -3000))
#'
#' # Plot the ellipsoid in 3D
#' \donttest{
#' if(requireNamespace("rgl", quietly = TRUE)){
#'   plot_ellipsoid_3d(ell5)
#'  }
#' }
plot_ellipsoid_3d <- function(object,
                              dim = c(1, 2, 3),
                              wire = FALSE,
                              aspect = TRUE,
                              background = NULL,
                              prediction = NULL,
                              col_layer = NULL,
                              pal = hcl.colors(100, palette = "Viridis"),
                              rev_pal = FALSE,
                              bg_sample = NULL,
                              col_ell = "#8b0000",
                              alpha_ell = 1,
                              alpha_bg = 1,
                              col_bg = "#8A8A8A",
                              fixed_lims = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              zlab = NULL,
                              ...) {

  # Argument Checks
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' is required, please install it before trying again.")
  }
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Provide a valid 'nicheR_ellipsoid' object.")
  }
  if (length(dim) != 3) {
    stop("'dim' must be a numeric vector of length 3.")
  }

  # Extract variable names for the specified dimensions
  vars <- object$var_names[dim]

  # Axes labels
  xlab <- if (!is.null(xlab)) xlab else vars[1]
  ylab <- if (!is.null(ylab)) ylab else vars[2]
  zlab <- if (!is.null(zlab)) zlab else vars[3]

  # Handle Data (Background or Prediction)
  ## Prioritize background if both are provided
  plot_data <- if (!is.null(background)) background else prediction

  if (!is.null(plot_data)) {
    ## Subsampling logic
    if (!is.null(bg_sample) && nrow(plot_data) > bg_sample) {
      idx <- sample(seq_len(nrow(plot_data)), bg_sample)
      plot_data <- plot_data[idx, , drop = FALSE]
    }

    ## Point colors
    if (!is.null(col_layer) && is.null(background)) {
      if (is.function(pal)) pal <- pal(100)
      if (isTRUE(rev_pal)) pal <- rev(pal)

      col_vals <- plot_data[[col_layer]]

      # Map internal points to palette; zeros/NAs get col_bg
      is_outside <- is.na(col_vals) | col_vals == 0

      pt_colors <- rep(col_bg, nrow(plot_data))

      # Internal mapping utility
      if (any(!is_outside)) {
        pt_colors[!is_outside] <- pal[map_to_pal(col_vals[!is_outside],
                                                 pal_len = length(pal))]
      }
    } else {
      pt_colors <- col_bg
    }

    # Setup 3D Scene
    rgl::plot3d(plot_data[, vars], col = pt_colors, alpha = alpha_bg,
                xlab = xlab, ylab = ylab, zlab = zlab,
                xlim = fixed_lims$xlim, ylim = fixed_lims$ylim,
                zlim = fixed_lims$zlim, aspect = aspect, box = FALSE, ...)


  } else {
    # Initialize 3D scene using ellipsoid centroid to ensure axes are drawn
    rgl::plot3d(t(object$centroid[dim]), type = "n",
                xlab = xlab, ylab = ylab, zlab = zlab,
                xlim = fixed_lims$xlim, ylim = fixed_lims$ylim,
                zlim = fixed_lims$zlim, aspect = aspect, box = FALSE, ...)
  }

  # Add Ellipsoid
  ell_mesh <- rgl::ellipse3d(x = object$cov_matrix[dim, dim],
                             centre = object$centroid[dim],
                             t = sqrt(object$chi2_cutoff))

  if (isTRUE(wire)) {
    rgl::wire3d(ell_mesh, col = col_ell)
  } else {
    rgl::shade3d(ell_mesh, col = col_ell, alpha = alpha_ell)
  }
}

#' Add data to an existing 3D E-space plot
#'
#' @usage
#' add_data_3d(data, dim = c(1, 2, 3), col_layer = NULL, alpha = 1, ...)
#'
#' @param data A data frame or matrix containing the points.
#' @param dim Integer vector of length 3. Indices of dimensions to plot.
#' @param col_layer Character or \code{NULL}. Column for coloring.
#' @param alpha Transparency for the points. Default is \code{1}.
#' @param ... Additional arguments passed to \code{rgl::points3d}.
#'
#' @return adds data to an existing 3D e-space plot
#'
#' @export
#'
#' @examples
#' # Building an ellipsoid
#' ## Define ranges for three variables
#' range <- data.frame(bio_1  = c(22, 32),
#'                     bio_12 = c(800, 4200),
#'                     bio_15 = c(45, 115))
#'
#' ## Build the ellipsoid
#' ell5 <- build_ellipsoid(ranges = range)
#' ell5$cov_limits
#'
#' ell5u <- update_ellipsoid_covariance(ell5, c("bio_1-bio_12" = 200,
#'                                              "bio_1-bio_15" = 0,
#'                                              "bio_12-bio_15" = -3000))
#' \donttest{
#' if(requireNamespace("rgl", quietly = TRUE)){
#'
#' # Plot the ellipsoid in 3D
#' plot_ellipsoid_3d(ell5u)
#'
#' # Add background points
#' add_data_3d(back_data[, c(3, 7, 10)], col = "#8A8A8A")
#'
#'  }
#' }
add_data_3d <- function(data,
                        dim = c(1, 2, 3),
                        col_layer = NULL,
                        alpha = 1,
                        ...) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' is required, please install it before trying again.")
  }
  if (missing(data) || !is.data.frame(data) && !is.matrix(data)) {
    stop("Provide a valid data frame or matrix for 'data'.")
  }

  # variable selection
  vars <- colnames(data)[dim]

  # Lean point addition
  rgl::points3d(data[, vars], alpha = alpha, ...)
}



#' Add an ellipsoid to an existing 3D E-space plot
#'
#' @usage
#' add_ellipsoid_3d(object, dim = c(1, 2, 3), wire = FALSE, col_ell = "#800000",
#'                   alpha_ell = 1, ...)
#'
#' @param object A \code{nicheR_ellipsoid} object.
#' @param dim Integer vector of length 3.
#' @param wire Logical. If \code{TRUE}, plots wireframe, otherwise plots a
#'    shaded volume. Default is \code{FALSE}.
#' @param col_ell Color of the ellipsoid. Default is \code{"#000000"}.
#' @param alpha_ell Transparency of the ellipsoid. Default is
#'    \code{1}. Not applied if \code{wire = TRUE}.
#' @param ... Additional arguments passed to \code{rgl::wire3d} or
#'    \code{rgl::shade3d}.
#'
#' @return adds a new ellipsoid object to the 3d rendering of the 3D plot
#'
#' @export
#'
#' @examples
#' # Building an ellipsoid
#' ## Define ranges for three variables
#' range <- data.frame(bio_1  = c(22, 32),
#'                     bio_12 = c(800, 4200),
#'                     bio_15 = c(45, 115))
#'
#' ## Build the ellipsoid
#' ell5 <- build_ellipsoid(ranges = range)
#' ell5$cov_limits
#'
#' ell5u <- update_ellipsoid_covariance(ell5, c("bio_1-bio_12" = 200,
#'                                              "bio_1-bio_15" = 0,
#'                                              "bio_12-bio_15" = -3000))
#'\donttest{
#' if(requireNamespace("rgl", quietly = TRUE)){
#'
#' # Plot the ellipsoid in 3D
#' plot_ellipsoid_3d(ell5u)
#'
#' # Add the original ellipsoid as a wireframe
#' add_ellipsoid_3d(ell5, wire = TRUE, col = "#0e008b")
#'
#'   }
#' }
add_ellipsoid_3d <- function(object,
                             dim = c(1, 2, 3),
                             wire = FALSE,
                             col_ell = "#800000",
                             alpha_ell = 1,
                             ...) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' is required, please install it before trying again.")
  }
  if (!inherits(object, "nicheR_ellipsoid")) {
    stop("Provide a valid 'nicheR_ellipsoid' object.")
  }
  if (length(dim) != 3) {
    stop("'dim' must be a numeric vector of length 3.")
  }

  ell_mesh <- rgl::ellipse3d(x = object$cov_matrix[dim, dim],
                             centre = object$centroid[dim],
                             t = sqrt(object$chi2_cutoff))

  if (isTRUE(wire)) {
    rgl::wire3d(ell_mesh, col = col_ell, ...)
  } else {
    rgl::shade3d(ell_mesh, col = col_ell, alpha = alpha_ell, ...)
  }
}
