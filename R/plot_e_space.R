#' Plot environmental space with optional ellipsoid overlays
#'
#' Produces pairwise views of a 3D environmental space with optional overlays of a
#' virtual niche ellipsoid and occurrence points. In 2D mode it returns a 3-panel
#' grid of pairwise scatterplots with projected ellipse boundaries. In 3D mode it
#' returns an interactive scatter3d of background environments with an optional
#' ellipsoid wireframe and centroids.
#'
#' @param env_bg A data.frame of background environments with at least three numeric
#'   predictor columns. These columns must contain the variables referenced by `x`,
#'   `y`, and `z`.
#' @param x,y,z Column specifications for the three predictors to display. Each may
#'   be a single column name or a single 1-based integer index into `env_bg`.
#' @param labels Character vector of length 3 giving axis labels for the x, y, and z
#'   variables in display order. Defaults to `c("ENV 1", "ENV 2", "ENV 3")`.
#' @param n_bg Positive number giving the maximum number of background rows to plot.
#'   If `nrow(env_bg) > n_bg`, a random subset of size `n_bg` is drawn using `rand_seed`.
#' @param niche Optional ellipsoid-like object describing the niche in E space.
#'   The validator requires `niche$center` and `niche$axes` to be numeric length 3.
#'   In 2D mode, projected ellipses are drawn using `build_ellipsoid`. In 3D mode,
#'   if present, an ellipsoid wireframe is read from `niche$surface` and the centroid
#'   is plotted from `niche$center`.
#' @param show.pts.in Logical. If `TRUE` and `niche` is provided, points inside the
#'   ellipsoid are computed with `get_suitable_env` and overplotted.
#' @param occ_pts Optional data.frame of occurrence points that includes the same
#'   predictor columns used for `x`, `y`, and `z`. These are overplotted if supplied.
#' @param rand_seed Integer used to seed the background downsampling.
#' @param show.occ.density Logical. If `TRUE` and `occ_pts` is provided, adds marginal
#'   density panels for each variable in 2D mode. Ignored in 3D mode.
#' @param plot.3d Logical. If `TRUE`, returns an interactive `plotly` 3D scatter with
#'   optional ellipsoid wireframe. If `FALSE`, returns a `ggpubr` grid of 2D panels.
#'
#' @details
#' - Column selectors `x`, `y`, and `z` are resolved to names during validation and are
#'   then used consistently across plotting layers.
#' - When `show.pts.in = TRUE`, a message clarifies that some points may appear outside
#'   a 2D ellipse boundary due to dimensionality. This does not indicate an error in the
#'   function or plotting.
#' - Large `n_bg` values trigger warnings since plotting may be slow.
#'
#' @return
#' - If `plot.3d = TRUE`: a `plotly` object.
#' - If `plot.3d = FALSE`: a `ggpubr` object containing arranged `ggplot2` panels.
#'
#' @section Validation:
#' Inputs are validated by `validate_plot_e_space_args()` before plotting. See that
#' helper for exact checks and error messages.
#'
#' @seealso [validate_plot_e_space_args()], [build_ellipsoid()], [get_suitable_env()]
#'
#' @examples
#' \dontrun{
#' # Minimal 2D example with three numeric columns
#' set.seed(1)
#' env_bg <- data.frame(A = rnorm(5e4), B = rnorm(5e4), C = rnorm(5e4))
#'
#' # 2D pairwise view
#' p2d <- plot_e_space(env_bg, x = "A", y = "B", z = "C",
#'                     labels = c("Temp", "Salinity", "SST"),
#'                     n_bg = 10000, plot.3d = FALSE)
#'
#' # 3D interactive background only
#' p3d <- plot_e_space(env_bg, 1, 2, 3, n_bg = 20000, plot.3d = TRUE)
#' }
#' @export


library(plotly)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

plot_e_space <- function(env_bg,
                         x, y, z,
                         labels = c("ENV 1", "ENV 2", "ENV 3"),
                         n_bg = 10000,
                         niche = NULL,
                         show.pts.in = FALSE,
                         occ_pts = NULL,
                         rand_seed = 1234,
                         show.occ.density = FALSE, # only for 2D plots
                         plot.3d = FALSE){

  # -- 0. Validate (no rlang in helper) --
  v <- validate_plot_e_space_args(env_bg, x, y, z,
                                  labels, n_bg, niche, show.pts.in,
                                  occ_pts, show.occ.density)

  # Downsample info + action
  if (nrow(env_bg) > n_bg) {
    message(sprintf("Sampling %d of %d rows from 'env_bg' for plotting.", n_bg, nrow(env_bg)))
    set.seed(rand_seed)
    env_bg <- env_bg[sample.int(nrow(env_bg), size = n_bg, replace = FALSE), ]
  }

  # --- 1. Base scatter layers ---
  if (isTRUE(plot.3d)) {
    return_plot <- plotly::plot_ly(data = env_bg,
                                   x = env_bg[[if (is.numeric(x)) names(env_bg)[x] else x]],
                                   y = env_bg[[if (is.numeric(y)) names(env_bg)[y] else y]],
                                   z = env_bg[[if (is.numeric(z)) names(env_bg)[z] else z]],
                                   type = "scatter3d",
                                   mode = "markers",
                                   marker = list(color = "lightgrey", size = 2),
                                   name = "Background Environments") %>%
      plotly::layout(
        title = list(text = "Background Environments (E-space)"),
        scene = list(
          xaxis = list(title = list(text = labels[1])),
          yaxis = list(title = list(text = labels[2])),
          zaxis = list(title = list(text = labels[3]))
        ),
        legend = list(x = 0.05, y = 0.95)
      )

  } else {

    # Use .data pronoun inside aes as you had; validator already resolved names
    p_main_y_x <- ggplot2::ggplot(env_bg, ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]],
                                                   y = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]])) +
      ggplot2::geom_point(alpha = 0.5, color = "grey", pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    p_main_z_x <- ggplot2::ggplot(env_bg, ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]],
                                                   y = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]])) +
      ggplot2::geom_point(alpha = 0.5, color = "grey", pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    p_main_z_y <- ggplot2::ggplot(env_bg, ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]],
                                                   y = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]])) +
      ggplot2::geom_point(alpha = 0.5, color = "grey", pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    x_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[1])) + ggplot2::xlab(NULL)
    y_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[2])) + ggplot2::xlab(NULL)
    z_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[3])) + ggplot2::xlab(NULL)

    return_plot <- ggpubr::ggarrange(
      x_name, p_main_y_x, p_main_z_x,
      NULL,   y_name,    p_main_z_y,
      NULL,   NULL,      z_name,
      ncol = 3, nrow = 3,
      widths = c(0.15, 0.425, 0.425),
      heights = c(0.45, 0.45, 0.15)
    )
  }

  # --- 2. Ellipsoid overlays and extras ---
  if (!is.null(niche)) {
    if(isTRUE(plot.3d)){

      return_plot <- return_plot %>%
        add_trace(data = niche$surface,
                  x = niche$surface[[names(niche$surface)[1]]],
                  y = niche$surface[[names(niche$surface)[2]]],
                  z = niche$surface[[names(niche$surface)[3]]],
                  type ="scatter3d", mode="lines",
                  line =list(color="blue"),
                  name ="Niche Boundary", inherit = FALSE) %>%
        add_markers(x = niche$center[1], y = niche$center[2], z = niche$center[3],
                    marker = list(color = 'red', size = 5),
                    name = "Niche Centroid") %>%
        plotly::layout(
          title = "Virtual Niche Boundary in E-space",
          legend = list(x = 0.05, y = 0.95)
        )

      if (isTRUE(show.pts.in)) {
        # Use base subsetting with resolved names from validator
        pts_in <- get_suitable_env(niche = FN_1,
                                   env_bg = env_bg[, v$col_names, drop = FALSE],
                                   out = "data.frame")

        return_plot <- return_plot %>%
          add_markers(data = pts_in,
                      x = pts_in[[names(pts_in)[1]]],
                      y = pts_in[[names(pts_in)[2]]],
                      z = pts_in[[names(pts_in)[3]]],
                      marker = list(color="darkgreen", size = 3),
                      name = "Suitable Environments", inherit = FALSE) %>%
          plotly::layout(
            title = "Virtual Niche Suitable Environment in E-space",
            legend = list(x = 0.05, y = 0.95)
          )
      }

      if (!is.null(occ_pts)) {

        return_plot <- return_plot %>%
          add_markers(data = occ_pts,
                      x = occ_pts[[if (is.numeric(x)) names(occ_pts)[x] else x]],
                      y = occ_pts[[if (is.numeric(y)) names(occ_pts)[y] else y]],
                      z = occ_pts[[if (is.numeric(z)) names(occ_pts)[z] else z]],
                      marker = list(color="orange", size = 3),
                      name = "Sampled Occurrences", inherit = FALSE) %>%
          plotly::layout(
            title = "Virtual Niche and Sampled Occurrences in E-space",
            legend = list(x = 0.05, y = 0.95)
          )
      }


    }else{

      center_y_x <- c(niche$center[2], niche$center[1])
      axes_y_x   <- c(niche$axes[2],   niche$axes[1])
      angles_y_x <- c(niche$angles[2], niche$angles[1])

      center_z_x <- c(niche$center[3], niche$center[1])
      axes_z_x   <- c(niche$axes[3],   niche$axes[1])
      angles_z_x <- c(niche$angles[3], niche$angles[1])

      center_z_y <- c(niche$center[3], niche$center[2])
      axes_z_y   <- c(niche$axes[3],   niche$axes[2])
      angles_z_y <- c(niche$angles[3], niche$angles[2])


      # 2D ellipse projections in each plane
      ell2d_y_x <- build_ellipsoid(center = center_y_x,
                                   axes = axes_y_x,
                                   angles = angles_y_x)
      ell2d_z_x <- build_ellipsoid(center = center_z_x,
                                   axes = axes_z_x,
                                   angles = angles_z_x)
      ell2d_z_y <- build_ellipsoid(center = center_z_y,
                                   axes = axes_z_y,
                                   angles = angles_z_y)

      ell_y_x <- p_main_y_x +
        geom_path(data = ell2d_y_x$surface, mapping = aes(x, y),
                  color = "royalblue", size = 0.5) +
        geom_segment(aes(x = ell2d_y_x$center[1] - ell2d_y_x$axes[1],
                         xend = ell2d_y_x$center[1] + ell2d_y_x$axes[1],
                         y = ell2d_y_x$center[2], yend = ell2d_y_x$center[2]),
                     color = "paleturquoise", linetype = "dashed") +
        geom_segment(aes(y = ell2d_y_x$center[2] - ell2d_y_x$axes[2],
                         yend = ell2d_y_x$center[2] + ell2d_y_x$axes[2],
                         x = ell2d_y_x$center[1], xend = ell2d_y_x$center[1]),
                     color = "paleturquoise", linetype = "dashed") +
        geom_point(aes(x = ell2d_y_x$center[1], y = ell2d_y_x$center[2]),
                   color = "tomato", size = 2)

      ell_z_x <- p_main_z_x +
        geom_path(data = ell2d_z_x$surface, mapping = aes(x, y),
                  color = "royalblue", size = 0.5) +
        geom_segment(aes(x = ell2d_z_x$center[1] - ell2d_z_x$axes[1],
                         xend = ell2d_z_x$center[1] + ell2d_z_x$axes[1],
                         y = ell2d_z_x$center[2], yend = ell2d_z_x$center[2]),
                     color = "paleturquoise", linetype = "dashed") +
        geom_segment(aes(y = ell2d_z_x$center[2] - ell2d_z_x$axes[2],
                         yend = ell2d_z_x$center[2] + ell2d_z_x$axes[2],
                         x = ell2d_z_x$center[1], xend = ell2d_z_x$center[1]),
                     color = "paleturquoise", linetype = "dashed") +
        geom_point(aes(x = ell2d_z_x$center[1], y = ell2d_z_x$center[2]),
                   color = "tomato", size = 2)

      ell_z_y <- p_main_z_y +
        geom_path(data = ell2d_z_y$surface, mapping = aes(x, y),
                  color = "royalblue", size = 0.5) +
        geom_segment(aes(x = ell2d_z_y$center[1] - ell2d_z_y$axes[1],
                         xend = ell2d_z_y$center[1] + ell2d_z_y$axes[1],
                         y = ell2d_z_y$center[2], yend = ell2d_z_y$center[2]),
                     color = "paleturquoise", linetype = "dashed") +
        geom_segment(aes(y = ell2d_z_y$center[2] - ell2d_z_y$axes[2],
                         yend = ell2d_z_y$center[2] + ell2d_z_y$axes[2],
                         x = ell2d_z_y$center[1], xend = ell2d_z_y$center[1]),
                     color = "paleturquoise", linetype = "dashed") +
        geom_point(aes(x = ell2d_z_y$center[1], y = ell2d_z_y$center[2]),
                   color = "tomato", size = 2)

      return_plot <- ggpubr::ggarrange(
        x_name, ell_y_x, ell_z_x,
        NULL,   y_name,  ell_z_y,
        NULL,   NULL,    z_name,
        ncol = 3, nrow = 3,
        widths = c(0.15, 0.425, 0.425),
        heights = c(0.45, 0.45, 0.15)
      )

      if (isTRUE(show.pts.in)) {
        message("Note: Some points may fall outside the ellipse boundary due to their dimensionality. This does not indicate an error in the function or plotting.")

        # Use base subsetting with resolved names from validator
        pts_in <- get_suitable_env(niche = niche,
                                   env_bg = env_bg[, v$col_names, drop = FALSE],
                                   out = "data.frame")

        ell_y_x <- ell_y_x +
          ggplot2::geom_point(data = pts_in,
                              ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]],
                                           y = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]]),
                              color = "darkolivegreen3", size = 0.5)
        ell_z_x <- ell_z_x +
          ggplot2::geom_point(data = pts_in,
                              ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]],
                                           y = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]]),
                              color = "darkolivegreen3", size = 0.5)
        ell_z_y <- ell_z_y +
          ggplot2::geom_point(data = pts_in,
                              ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]],
                                           y = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]]),
                              color = "darkolivegreen3", size = 0.5)

        return_plot <- ggpubr::ggarrange(
          x_name, ell_y_x, ell_z_x,
          NULL,   y_name,  ell_z_y,
          NULL,   NULL,    z_name,
          ncol = 3, nrow = 3,
          widths = c(0.15, 0.425, 0.425),
          heights = c(0.45, 0.45, 0.15)
        )
      }


      if (!is.null(occ_pts)) {
        ell_y_x <- ell_y_x +
          geom_point(data = occ_pts, aes(x = .data[[y]], y = .data[[x]]),
                     color = "darkorange", size = 0.5)
        ell_z_x <- ell_z_x +
          geom_point(data = occ_pts, aes(x = .data[[z]], y = .data[[x]]),
                     color = "darkorange", size = 0.5)
        ell_z_y <- ell_z_y +
          geom_point(data = occ_pts, aes(x = .data[[z]], y = .data[[y]]),
                     color = "darkorange", size = 0.5)

        return_plot <- ggpubr::ggarrange(
          x_name, ell_y_x, ell_z_x,
          NULL,   y_name,  ell_z_y,
          NULL,   NULL,    z_name,
          ncol = 3, nrow = 3,
          widths = c(0.15, 0.425, 0.425),
          heights = c(0.45, 0.45, 0.15)
        )

        if (isTRUE(show.occ.density)) {
          # Use base range() to avoid tidyselect here
          rng_z <- range(env_bg[[ if (is.numeric(z)) names(env_bg)[z] else z ]], na.rm = TRUE)
          rng_y <- range(env_bg[[ if (is.numeric(y)) names(env_bg)[y] else y ]], na.rm = TRUE)
          rng_x <- range(env_bg[[ if (is.numeric(x)) names(env_bg)[x] else x ]], na.rm = TRUE)

          env_z_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]])) +
            ggplot2::geom_density(fill = "darkorange", alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_z) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank())

          env_y_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]])) +
            ggplot2::geom_density(fill = "darkorange", alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_y) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank())

          env_x_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]])) +
            ggplot2::geom_density(fill = "darkorange", alpha = 0.6) +
            ggplot2::coord_flip() +
            ggplot2::scale_x_continuous(limits = rng_x) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank())

          env_y_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]])) +
            ggplot2::geom_density(fill = "darkorange", alpha = 0.6) +
            ggplot2::coord_flip() +
            ggplot2::scale_x_continuous(limits = rng_y) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank())

          return_plot <- ggpubr::ggarrange(
            NULL,      env_y_top, env_z_top, NULL,
            x_name,    ell_y_x,   ell_z_x,   env_x_right,
            NULL,      y_name,    ell_z_y,   env_y_right,
            NULL,      NULL,      z_name,    NULL,
            ncol = 4, nrow = 4,
            widths  = c(0.1, 0.4, 0.4, 0.1),
            heights = c(0.1, 0.4, 0.4, 0.1)
          )
        }
      }
    }
  }

  return(return_plot)
}

#' Validate arguments for plot_e_space
#'
#' Performs input checks for \code{plot_e_space} and emits informative errors
#' or warnings. Ensures:
#' \itemize{
#'   \item \code{env_bg} is a data frame with at least three numeric predictor columns.
#'   \item \code{x}, \code{y}, \code{z} resolve to distinct, numeric columns in \code{env_bg}.
#'   \item \code{labels} is a character vector of length 3.
#'   \item \code{n_bg} is a single, positive number (large values trigger a warning).
#'   \item Optional objects (\code{niche}, \code{occ_pts}) have the required structure.
#' }
#'
#' @param env_bg A data.frame with at least three numeric predictor columns.
#' @param x,y,z Column names or single integer indices in \code{env_bg} for the
#'   three predictors, in the same order used to define the niche or ellipsoid.
#' @param labels Character vector of length 3 used as axis labels.
#' @param n_bg Positive number indicating the maximum number of background rows
#'   to plot. Values much larger than typical point clouds may slow plotting.
#' @param niche Optional list with required elements \code{center} and \code{axes}
#'   (both numeric vectors of length 3), and \emph{recommended} element \code{angles}
#'   (numeric vector of length 3). If \code{angles} is missing, 2D projections in the
#'   main function will fail.
#' @param show.pts.in Logical. If TRUE, \code{niche} must be provided so that
#'   points inside the ellipsoid can be computed.
#' @param occ_pts Optional data.frame containing the same three predictor columns
#'   referenced by \code{x}, \code{y}, and \code{z}; columns must be numeric.
#' @param show.occ.density Logical. If TRUE, \code{occ_pts} must be provided.
#'
#' @return An invisible list with element \code{col_names}, the resolved predictor
#'   column names in the order \code{c(x, y, z)}.
#'
#' @details This validator does not evaluate plotting mode. In the main function,
#'   density panels are available only for 2D plots. Provide \code{angles} in
#'   \code{niche} to enable 2D ellipsoid projections.
#'
#' @keywords internal
#' @seealso \code{\link{plot_e_space}}
#' @export
validate_plot_e_space_args <- function(env_bg, x, y, z,
                                       labels, n_bg,
                                       niche, show.pts.in,
                                       occ_pts, show.occ.density) {
  # --- env_bg ---
  if (!is.data.frame(env_bg)) stop("'env_bg' must be a data.frame.")
  if (ncol(env_bg) < 3) stop("'env_bg' must have at least 3 columns.")
  if (nrow(env_bg) < 1) stop("'env_bg' has no rows to plot.")

  # Helper: resolve a single spec (name or index) to a valid column name
  resolve_one <- function(spec) {
    if (is.character(spec) && length(spec) == 1) {
      if (!spec %in% names(env_bg)) {
        stop(sprintf("Column '%s' not found in 'env_bg'.", spec))
      }
      spec
    } else if (is.numeric(spec) && length(spec) == 1) {
      idx <- as.integer(spec)
      if (is.na(idx) || idx < 1 || idx > ncol(env_bg)) {
        stop("'x', 'y', and 'z' indices must be within the number of columns in 'env_bg'.")
      }
      names(env_bg)[idx]
    } else {
      stop("'x', 'y', and 'z' must be column names or single integer indices.")
    }
  }

  col_x <- resolve_one(x)
  col_y <- resolve_one(y)
  col_z <- resolve_one(z)
  col_names <- c(col_x, col_y, col_z)

  # Distinctness
  if (any(duplicated(col_names))) {
    stop("'x', 'y', and 'z' must refer to three distinct columns in 'env_bg'.")
  }

  # Numeric columns
  non_numeric <- col_names[!vapply(env_bg[, col_names, drop = FALSE], is.numeric, logical(1))]
  if (length(non_numeric)) {
    stop(sprintf("These columns must be numeric in 'env_bg': %s",
                 paste(non_numeric, collapse = ", ")))
  }

  # labels
  if (!(is.character(labels) && length(labels) == 3 && all(!is.na(labels)))) {
    stop("'labels' must be a non-NA character vector of length 3.")
  }

  # n_bg
  if (!(length(n_bg) == 1 && is.finite(n_bg) && n_bg > 0)) {
    stop("'n_bg' must be a single positive number.")
  }
  if (n_bg > 100000) {
    warning("'n_bg' is very large. Plotting may be slow.", call. = FALSE, immediate. = TRUE)
  } else if (n_bg > 10000) {
    warning("Selected number of background points is large and may slow plotting.",
            call. = FALSE, immediate. = TRUE)
  }
  if (is.numeric(n_bg) && n_bg %% 1 != 0) {
    warning("'n_bg' is not an integer; it will be truncated when sampling.", call. = FALSE, immediate. = TRUE)
  }
  if (nrow(env_bg) < n_bg) {
    warning(sprintf("'n_bg' (%d) exceeds available rows in 'env_bg' (%d). All rows will be used.",
                    as.integer(n_bg), nrow(env_bg)), call. = FALSE, immediate. = TRUE)
  }

  # niche object if provided
  if (!is.null(niche)) {
    need <- c("center", "axes")
    miss <- setdiff(need, names(niche))
    if (length(miss)) {
      stop(sprintf("'niche' is missing required fields: %s", paste(miss, collapse = ", ")))
    }
    if (!(is.numeric(niche$center) && length(niche$center) == 3)) {
      stop("'niche$center' must be numeric of length 3.")
    }
    if (!(is.numeric(niche$axes) && length(niche$axes) == 3)) {
      stop("'niche$axes' must be numeric of length 3.")
    }
    if (any(!is.finite(niche$axes)) || any(niche$axes <= 0)) {
      stop("'niche$axes' must contain finite, positive values.")
    }
    # Angles are required by the main function for 2D projections
    if (!("angles" %in% names(niche))) {
      stop("'niche$angles' is required for 2D projections in 'plot_e_space'. Provide a numeric vector of length 3.")
    }
    if (!(is.numeric(niche$angles) && length(niche$angles) == 3)) {
      stop("'niche$angles' must be numeric of length 3.")
    }
    if (any(!is.finite(niche$angles))) {
      stop("'niche$angles' must contain finite values.")
    }
  }

  # show.pts.in requires niche
  if (isTRUE(show.pts.in) && is.null(niche)) {
    warning("'show.pts.in' is TRUE but 'niche' is NULL. Points inside cannot be computed.",
            call. = FALSE, immediate. = TRUE)
  }

  # occ_pts and density flags
  if (!is.null(occ_pts)) {
    if (!is.data.frame(occ_pts)) stop("'occ_pts' must be a data.frame.")
    missing_occ <- setdiff(col_names, names(occ_pts))
    if (length(missing_occ)) {
      stop(sprintf("'occ_pts' is missing required columns: %s", paste(missing_occ, collapse = ", ")))
    }
    non_num_occ <- col_names[!vapply(occ_pts[, col_names, drop = FALSE], is.numeric, logical(1))]
    if (length(non_num_occ)) {
      stop(sprintf("These 'occ_pts' columns must be numeric: %s", paste(non_num_occ, collapse = ", ")))
    }
    if (nrow(occ_pts) < 1) warning("'occ_pts' has zero rows; nothing to plot.", call. = FALSE, immediate. = TRUE)
  } else if (isTRUE(show.occ.density)) {
    warning("'show.occ.density' is TRUE but 'occ_pts' is NULL. Density panels will be skipped.",
            call. = FALSE, immediate. = TRUE)
  }

  invisible(list(col_names = col_names))
}



