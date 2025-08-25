#' Plot Environmental Space with Optional Ellipsoid Overlays
#'
#' Produces pairwise views of a 3D environmental space with optional overlays of a
#' virtual niche ellipsoid and occurrence points. In 2D mode, it returns a grid
#' of pairwise scatterplots with projected ellipse boundaries. In 3D mode, it
#' returns an interactive `plotly` scatterplot.
#'
#' @param env_bg A `data.frame` of background environments with at least three
#'   numeric predictor columns. These columns must contain the variables
#'   referenced by `x`, `y`, and `z`.
#' @param x,y,z Column specifications for the three predictors to display. Each may
#'   be a single column name (character string) or a single 1-based integer index
#'   into `env_bg`.
#' @param labels Character vector of length 3 giving axis labels for the x, y,
#'   and z variables in display order. Defaults to `c("ENV 1", "ENV 2", "ENV 3")`.
#' @param n_bg Positive integer giving the maximum number of background rows to plot.
#'   If `nrow(env_bg)` is greater than `n_bg`, a random subset of size `n_bg` is
#'   drawn. Using a large `n_bg` may slow plotting.
#' @param niche Optional object of class `ellipsoid` describing the niche. If
#'   provided, its boundary and center will be plotted. For 2D plots, the
#'   object must contain `niche$angles`.
#' @param show.pts.in Logical. If `TRUE` and `niche` is provided, points from
#'   `env_bg` that fall inside the ellipsoid are highlighted. A warning is
#'   issued if `niche` is not provided.
#' @param occ_pts Optional `data.frame` of occurrence points that includes the
#'   same predictor columns used for `x`, `y`, and `z`. These are overplotted
#'   if supplied.
#' @param rand_seed Integer used to set the random number generator seed for
#'   reproducible background downsampling.
#' @param show.occ.density Logical. If `TRUE` and `occ_pts` is provided, this
#'   adds marginal density panels for each variable. This is only supported in
#'   2D plots (`plot.3d = FALSE`).
#' @param plot.3d Logical. If `TRUE`, returns an interactive `plotly` 3D scatter
#'   plot. If `FALSE` (the default), returns a static `ggpubr` grid of 2D panels.
#'
#' @return
#' - If `plot.3d = TRUE`: A `plotly` object.
#' - If `plot.3d = FALSE`: A `ggpubr` object containing arranged `ggplot2` panels.
#'
#' @details This function is a powerful visualization tool for understanding the
#'   relationship between a species' niche, its occurrences, and the available
#'   environmental space. The 2D views are especially useful for publication-quality
#'   figures, while the 3D interactive plot is great for data exploration. Note that
#'   in 2D plots, points from a 3D space may appear to fall outside the projected
#'   2D ellipse boundary. This is expected and does not indicate an error.
#'
#' @family plotting functions
#' @seealso [validate_plot_e_space_args()], [build_ellps()], [get_suitable_env()]
#' @importFrom RColorBrewer brewer.pal
#' @export
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
        pts_in <- get_suitable_environment(niche = FN_1,
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


    } else{ # 2D plots

      center_y_x <- c(niche$center[2], niche$center[1])
      axes_y_x   <- c(niche$axes[2],   niche$axes[1])

      center_z_x <- c(niche$center[3], niche$center[1])
      axes_z_x   <- c(niche$axes[3],   niche$axes[1])

      center_z_y <- c(niche$center[3], niche$center[2])
      axes_z_y   <- c(niche$axes[3],   niche$axes[2])

      ell2d_y_x <- build_ellps(center = center_y_x, axes = axes_y_x, angles = c(0, 0))
      ell2d_z_x <- build_ellps(center = center_z_x, axes = axes_z_x, angles = c(0, 0))
      ell2d_z_y <- build_ellps(center = center_z_y, axes = axes_z_y, angles = c(0, 0))

      ell_y_x <- p_main_y_x +
        ggplot2::geom_path(data = ell2d_y_x$surface, aes(x, y), color = "royalblue", linewidth = 0.5) +
        ggplot2::annotate("segment",
                          x = ell2d_y_x$center[1] - ell2d_y_x$axes[1],
                          xend = ell2d_y_x$center[1] + ell2d_y_x$axes[1],
                          y = ell2d_y_x$center[2], yend = ell2d_y_x$center[2],
                          color = "paleturquoise", linetype = "dashed") +
        ggplot2::annotate("segment",
                          y = ell2d_y_x$center[2] - ell2d_y_x$axes[2],
                          yend = ell2d_y_x$center[2] + ell2d_y_x$axes[2],
                          x = ell2d_y_x$center[1], xend = ell2d_y_x$center[1],
                          color = "paleturquoise", linetype = "dashed") +
        ggplot2::annotate("point", x = ell2d_y_x$center[1], y = ell2d_y_x$center[2],
                          color = "tomato", size = 2)

      ell_z_x <- p_main_z_x +
        ggplot2::geom_path(data = ell2d_z_x$surface, aes(x, y), color = "royalblue", linewidth = 0.5) +
        ggplot2::annotate("segment",
                          x = ell2d_z_x$center[1] - ell2d_z_x$axes[1],
                          xend = ell2d_z_x$center[1] + ell2d_z_x$axes[1],
                          y = ell2d_z_x$center[2], yend = ell2d_z_x$center[2],
                          color = "paleturquoise", linetype = "dashed") +
        ggplot2::annotate("segment",
                          y = ell2d_z_x$center[2] - ell2d_z_x$axes[2],
                          yend = ell2d_z_x$center[2] + ell2d_z_x$axes[2],
                          x = ell2d_z_x$center[1], xend = ell2d_z_x$center[1],
                          color = "paleturquoise", linetype = "dashed") +
        ggplot2::annotate("point", x = ell2d_z_x$center[1], y = ell2d_z_x$center[2],
                          color = "tomato", size = 2)

      ell_z_y <- p_main_z_y +
        ggplot2::geom_path(data = ell2d_z_y$surface, aes(x, y), color = "royalblue", linewidth = 0.5) +
        ggplot2::annotate("segment",
                          x = ell2d_z_y$center[1] - ell2d_z_y$axes[1],
                          xend = ell2d_z_y$center[1] + ell2d_z_y$axes[1],
                          y = ell2d_z_y$center[2], yend = ell2d_z_y$center[2],
                          color = "paleturquoise", linetype = "dashed") +
        ggplot2::annotate("segment",
                          y = ell2d_z_y$center[2] - ell2d_z_y$axes[2],
                          yend = ell2d_z_y$center[2] + ell2d_z_y$axes[2],
                          x = ell2d_z_y$center[1], xend = ell2d_z_y$center[1],
                          color = "paleturquoise", linetype = "dashed") +
        ggplot2::annotate("point", x = ell2d_z_y$center[1], y = ell2d_z_y$center[2],
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
        # Use base subsetting with resolved names from validator
        pts_in <- get_suitable_env(niche = niche, # Corrected FN_1 to niche
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
          ggplot2::geom_point(data = occ_pts, ggplot2::aes(x = .data[[y]], y = .data[[x]]),
                              color = "darkorange", size = 0.5)
        ell_z_x <- ell_z_x +
          ggplot2::geom_point(data = occ_pts, ggplot2::aes(x = .data[[z]], y = .data[[x]]),
                              color = "darkorange", size = 0.5)
        ell_z_y <- ell_z_y +
          ggplot2::geom_point(data = occ_pts, ggplot2::aes(x = .data[[z]], y = .data[[y]]),
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
          # Use base range() to avoid tidyselect here
          rng_z <- range(env_bg[[ if (is.numeric(z)) names(env_bg)[z] else z ]], na.rm = TRUE)
          rng_y <- range(env_bg[[ if (is.numeric(y)) names(env_bg)[y] else y ]], na.rm = TRUE)
          rng_x <- range(env_bg[[ if (is.numeric(x)) names(env_bg)[x] else x ]], na.rm = TRUE)

          env_z_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]])) +
            ggplot2::geom_density(fill = "darkorange", alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_z) +
            ggplot2::scale_y_continuous(n.breaks = 3) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank())

          env_y_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]])) +
            ggplot2::geom_density(fill = "darkorange", alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_y) +
            ggplot2::scale_y_continuous(n.breaks = 3) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank())

          env_x_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]])) +
            ggplot2::geom_density(fill = "darkorange", alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_x) +
            ggplot2::coord_flip() +
            ggplot2::scale_y_continuous(n.breaks = 3) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank())

          env_y_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]])) +
            ggplot2::geom_density(fill = "darkorange", alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_y) +
            ggplot2::coord_flip() +
            ggplot2::scale_y_continuous(n.breaks = 3) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank())

          return_plot <- ggpubr::ggarrange(
            NULL, env_y_top, env_z_top, NULL,
            x_name, ell_y_x, ell_z_x, env_x_right,
            NULL, y_name, ell_z_y, env_y_right,
            NULL, NULL, z_name, NULL,
            ncol = 4, nrow = 4,
            widths = c(0.1, 0.4, 0.4, 0.1),
            heights = c(0.1, 0.4, 0.4, 0.1)
          )
        }
      }
    }
  }

  return(return_plot)
}
