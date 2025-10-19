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
                         plot.3d = FALSE,
                         colors = NULL,
                         palette = "default"){

  # -- 0. Validate (no rlang in helper) --
  v <- validate_plot_e_space_args(env_bg, x, y, z,
                                  labels, n_bg, niche, show.pts.in,
                                  occ_pts, show.occ.density)

  # Set colors for legend and plot
  palettes <- list(

    default = list(
      bg           = "#9093A2FF",
      ellipsoid    = "#2A363BFF",
      centroid     = "#D72000FF",
      tolerance    = "#EE6100FF",
      suitable_env = "#FED789FF",
      occ          = "#B4BF3AFF"),

    palette2 = list(
      bg           = "#9CA9BAFF",
      ellipsoid    = "#3D619DFF",
      centroid     = "#345084FF",
      tolerance    = "#693829FF",
      suitable_env = "#CFB267FF",
      occ          = "#A56A3EFF"),


    palette3 = list(
      bg           = "#C8CCC6FF",
      ellipsoid    = "#023743FF",
      centroid     = "#72874EFF",
      tolerance    = "#476F84FF",
      suitable_env = "#FED789FF",
      occ          = "#A4BED5FF"),

    palette4 = list(
      bg           = "#C0D1CEFF",
      ellipsoid    = "#859B6CFF",
      centroid     = "#B74954FF",
      tolerance    = "#A99364FF",
      suitable_env = "#C2DDB2FF",
      occ          = "#EBA49EFF"),

    palette5 = list(
      bg           = "#A89F8EFF",
      ellipsoid    = "#7887A4FF",
      centroid     = "#A8CDECFF",
      tolerance    = "#682C37FF",
      suitable_env = "#F6955EFF",
      occ          = "#9B6981FF"),


    palette6 = list(
      bg           = "#D3D4D8FF",
      ellipsoid    = "#731A12FF",
      centroid     = "#F2D43DFF",
      tolerance    = "#3F858CFF",
      suitable_env = "#D9814EFF",
      occ          = "#707322FF")

  )


  # pick palette
  base_colors <- palettes[[palette]]


  if (is.null(colors)) {
    colors <- base_colors
  } else {
    if (is.null(names(colors))) {
      # user gave unnamed vector -> assign in order
      names(colors) <- names(base_colors)[seq_along(colors)]
      if(!is.list(colors)){ colors <- as.list(colors)}

      message("No names detected in 'colors'. ",
              "Using the order provided: ",
              paste(names(colors), collapse = ", "),
              ". Remaining colors filled with defaults.\n",
              "If you want to change a specific object, use a named list",
              "Available options are: ",
              paste(names(base_colors), collapse = ", "), ".\n",
              "Example: colors = list(ellipsoid = 'black', bg = '#72874EFF')")
    }

    if(!is.list(colors)){ colors <- as.list(colors)}

    # merge with defaults (keeps user overrides, fills missing)
    colors <- modifyList(base_colors, colors)
  }




  # Check options for  legend
  ## Toggle these on/off
  opts <- list(
    background_point = !is.null(n_bg),
    trace_line       = TRUE,
    centroid_point   = TRUE,
    tolerance_range_line  = TRUE,
    suitable_point   = show.pts.in,
    occurrence_point = !is.null(occ_pts)
  )

  ## 1) Legend spec (base R data.frame)
  legend_items <- data.frame(
    id       = c("background_point","trace_line","centroid_point",
                 "tolerance_range_line","suitable_point","occurrence_point"),
    type     = c("point","line","point","line","point","point"),
    label    = c("Background environments","Niche boundary","Niche centroid",
                 "Tolerance ranges","Suitable evironments","Occurrence"),
    color    = c(colors[["bg"]], colors[["ellipsoid"]], colors[["centroid"]],
                 colors[["tolerance"]], colors[["suitable_env"]], colors[["occ"]]),
    linetype = c(NA, 1, NA, 2, NA, NA),
    stringsAsFactors = FALSE
  )

  ## Mark active rows using the opts list
  active <- logical(nrow(legend_items))
  for (i in seq_len(nrow(legend_items))) {
    active[i] <- isTRUE(opts[[ legend_items$id[i] ]])
  }

  legend_items <- legend_items[active, , drop = FALSE]

  legend_plot <- NULL

  # Early exit if nothing to show
  if (nrow(legend_items) == 0) {
    legend_plot <- ggplot() + theme_void()
  } else {
    # 2) Auto layout
    top_y   <- 2        # top anchor
    spacing <- 0.25      # vertical gap between rows
    x_point <- 0.00
    x_text  <- 0.1
    x0_line <- -0.2
    x1_line <-  0.00

    legend_items <- legend_items %>%
      mutate(
        row = dplyr::row_number(),
        y   = top_y - (row - 1) * spacing,
        x_point = x_point,
        x_text  = x_text,
        x0_line = x0_line,
        x1_line = x1_line
      )

    # 3) Build plot with per-row aesthetics (use *identity* scales)
    legend_base <- ggplot() +
      coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, 2.5), clip = "off") +
      theme_void() +
      theme(legend.position = "none")

    legend_plot <- legend_base +
      # points
      geom_point(
        data = filter(legend_items, type == "point"),
        aes(x = x_point, y = y, colour = color),
        size = 2
      ) +
      # lines
      geom_segment(
        data = filter(legend_items, type == "line"),
        aes(x = x0_line, xend = x1_line, y = y, yend = y,
            colour = color, linetype = linetype),
        linewidth = 0.5
      ) +
      # labels
      geom_text(
        data = legend_items,
        aes(x = x_text, y = y, label = label),
        hjust = 0
      ) +
      scale_colour_identity() +
      scale_linetype_identity()
  }



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
      ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    p_main_z_x <- ggplot2::ggplot(env_bg, ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]],
                                                       y = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]])) +
      ggplot2::geom_point(alpha = 0.5, color =  colors[["bg"]], pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    p_main_z_y <- ggplot2::ggplot(env_bg, ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]],
                                                       y = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]])) +
      ggplot2::geom_point(alpha = 0.5, color =  colors[["bg"]], pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    x_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[1]),
                         fontface = "bold")
    y_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[2]),
                         fontface = "bold")
    z_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[3]),
                         fontface = "bold")

    return_plot <- ggpubr::ggarrange(
      x_name, p_main_y_x, p_main_z_x,
      legend_plot,   y_name,    p_main_z_y,
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
                  line =list(color= colors[["ellipsoid"]]),
                  name ="Niche Boundary", inherit = FALSE) %>%
        add_markers(x = niche$center[1], y = niche$center[2], z = niche$center[3],
                    marker = list(color =  colors[["centroid"]], size = 5),
                    name = "Niche Centroid") %>%
        plotly::layout(
          title = "Virtual Niche Boundary in E-space",
          legend = list(x = 0.05, y = 0.95)
        )

      if (isTRUE(show.pts.in)) {
        # Use base subsetting with resolved names from validator
        pts_in <- get_suitable_env(niche = niche,
                                   env_bg = env_bg[, v$col_names, drop = FALSE],
                                   out = "data.frame")

        return_plot <- return_plot %>%
          add_markers(data = pts_in,
                      x = pts_in[[names(pts_in)[1]]],
                      y = pts_in[[names(pts_in)[2]]],
                      z = pts_in[[names(pts_in)[3]]],
                      marker = list(color= colors[["suitable_env"]], size = 3),
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
                      marker = list(color= colors[["occ"]], size = 3),
                      name = "Sampled Occurrences", inherit = FALSE) %>%
          plotly::layout(
            title = "Virtual Niche and Sampled Occurrences in E-space",
            legend = list(x = 0.05, y = 0.95)
          )
      }


    } else{ # 2D plots

      center_y_x <- c(niche$center[2], niche$center[1])
      axes_y_x   <- c(niche$axes[2],   niche$axes[1])
      angle_y_x <- c(niche$angles[2], niche$angles[1])

      center_z_x <- c(niche$center[3], niche$center[1])
      axes_z_x   <- c(niche$axes[3],   niche$axes[1])
      angle_z_x <- c(niche$angles[3], niche$angles[1])

      center_z_y <- c(niche$center[3], niche$center[2])
      axes_z_y   <- c(niche$axes[3],   niche$axes[2])
      angle_z_y <- c(niche$angles[3], niche$angles[2])

      ell2d_y_x <- build_ellps(center = center_y_x, axes = axes_y_x, angles = angle_y_x)
      ell2d_z_x <- build_ellps(center = center_z_x, axes = axes_z_x, angles = angle_z_x)
      ell2d_z_y <- build_ellps(center = center_z_y, axes = axes_z_y, angles = angle_z_y)

      ell_y_x <- p_main_y_x
      ell_z_x <- p_main_z_x
      ell_z_y <- p_main_z_y


      if (any(niche$angles != 0)) {
        angle_warn <- ggplot2::ggplot() + ggplot2::theme_void() +
          ggplot2::geom_text(ggplot2::aes(0, 0,
                             label = "Note: The ellipsoid is angled, its shape may appear distorted and some points may fall outside due to dimensionality."),
                             size = 2)
      } else {
        angle_warn <- NULL
      }

      if (isTRUE(show.pts.in)) {
        # Use base subsetting with resolved names from validator
        pts_in <- get_suitable_env(niche = niche,
                                   env_bg = env_bg[, v$col_names, drop = FALSE],
                                   out = "data.frame")

        ell_y_x <- ell_y_x +
          ggplot2::geom_point(data = pts_in,
                              ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]],
                                           y = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]]),
                              color = colors[["suitable_env"]], size = 0.5)
        ell_z_x <- ell_z_x +
          ggplot2::geom_point(data = pts_in,
                              ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]],
                                           y = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]]),
                              color = colors[["suitable_env"]], size = 0.5)
        ell_z_y <- ell_z_y +
          ggplot2::geom_point(data = pts_in,
                              ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]],
                                           y = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]]),
                              color = colors[["suitable_env"]], size = 0.5)

      }


      if (!is.null(occ_pts)) {
        ell_y_x <- ell_y_x +
          ggplot2::geom_point(data = occ_pts, ggplot2::aes(x = .data[[y]], y = .data[[x]]),
                              color =  colors[["occ"]], size = 0.5)
        ell_z_x <- ell_z_x +
          ggplot2::geom_point(data = occ_pts, ggplot2::aes(x = .data[[z]], y = .data[[x]]),
                              color =  colors[["occ"]], size = 0.5)
        ell_z_y <- ell_z_y +
          ggplot2::geom_point(data = occ_pts, ggplot2::aes(x = .data[[z]], y = .data[[y]]),
                              color =  colors[["occ"]], size = 0.5)
      }


      ell_y_x <- ell_y_x +
        ggplot2::geom_path(data = ell2d_y_x$surface, aes(x, y), color =  colors[["ellipsoid"]], linewidth = 0.5) +
        ggplot2::annotate("segment",
                          x = ell2d_y_x$center[1] - ell2d_y_x$axes[1],
                          xend = ell2d_y_x$center[1] + ell2d_y_x$axes[1],
                          y = ell2d_y_x$center[2], yend = ell2d_y_x$center[2],
                          color =  colors[["tolerance"]], linetype = "dashed") +
        ggplot2::annotate("segment",
                          y = ell2d_y_x$center[2] - ell2d_y_x$axes[2],
                          yend = ell2d_y_x$center[2] + ell2d_y_x$axes[2],
                          x = ell2d_y_x$center[1], xend = ell2d_y_x$center[1],
                          color =  colors[["tolerance"]], linetype = "dashed") +
        ggplot2::annotate("point", x = ell2d_y_x$center[1], y = ell2d_y_x$center[2],
                          color =  colors[["centroid"]], size = 2)

      ell_z_x <- ell_z_x +
        ggplot2::geom_path(data = ell2d_z_x$surface, aes(x, y), color =  colors[["ellipsoid"]], linewidth = 0.5) +
        ggplot2::annotate("segment",
                          x = ell2d_z_x$center[1] - ell2d_z_x$axes[1],
                          xend = ell2d_z_x$center[1] + ell2d_z_x$axes[1],
                          y = ell2d_z_x$center[2], yend = ell2d_z_x$center[2],
                          color =  colors[["tolerance"]], linetype = "dashed") +
        ggplot2::annotate("segment",
                          y = ell2d_z_x$center[2] - ell2d_z_x$axes[2],
                          yend = ell2d_z_x$center[2] + ell2d_z_x$axes[2],
                          x = ell2d_z_x$center[1], xend = ell2d_z_x$center[1],
                          color =  colors[["tolerance"]], linetype = "dashed") +
        ggplot2::annotate("point", x = ell2d_z_x$center[1], y = ell2d_z_x$center[2],
                          color =  colors[["centroid"]], size = 2)

      ell_z_y <- ell_z_y +
        ggplot2::geom_path(data = ell2d_z_y$surface, aes(x, y), color =  colors[["ellipsoid"]], linewidth = 0.5) +
        ggplot2::annotate("segment",
                          x = ell2d_z_y$center[1] - ell2d_z_y$axes[1],
                          xend = ell2d_z_y$center[1] + ell2d_z_y$axes[1],
                          y = ell2d_z_y$center[2], yend = ell2d_z_y$center[2],
                          color =  colors[["tolerance"]], linetype = "dashed") +
        ggplot2::annotate("segment",
                          y = ell2d_z_y$center[2] - ell2d_z_y$axes[2],
                          yend = ell2d_z_y$center[2] + ell2d_z_y$axes[2],
                          x = ell2d_z_y$center[1], xend = ell2d_z_y$center[1],
                          color =  colors[["tolerance"]], linetype = "dashed") +
        ggplot2::annotate("point", x = ell2d_z_y$center[1], y = ell2d_z_y$center[2],
                          color =  colors[["centroid"]], size = 2)


      return_plot <- ggpubr::ggarrange(x_name,        ell_y_x,   ell_z_x,
                                       legend_plot,   y_name,    ell_z_y,
                                       NULL,          angle_warn,      z_name,
                                       ncol = 3, nrow = 3,
                                       widths = c(0.15, 0.425, 0.425),
                                       heights = c(0.45, 0.45, 0.15))

      if (isTRUE(show.occ.density)) {
        if (!is.null(occ_pts)) {
          # Use base range() to avoid tidyselect here
          # Use base range() to avoid tidyselect here
          rng_z <- range(env_bg[[ if (is.numeric(z)) names(env_bg)[z] else z ]], na.rm = TRUE)
          rng_y <- range(env_bg[[ if (is.numeric(y)) names(env_bg)[y] else y ]], na.rm = TRUE)
          rng_x <- range(env_bg[[ if (is.numeric(x)) names(env_bg)[x] else x ]], na.rm = TRUE)

          env_z_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(z)) names(env_bg)[z] else z ]])) +
            ggplot2::geom_density(fill =  colors[["occ"]], alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_z) +
            ggplot2::scale_y_continuous(n.breaks = 3) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_text(size = 5))

          env_y_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]])) +
            ggplot2::geom_density(fill =  colors[["occ"]], alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_y) +
            ggplot2::scale_y_continuous(n.breaks = 3) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_text(size = 5))

          env_x_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(x)) names(env_bg)[x] else x ]])) +
            ggplot2::geom_density(fill =  colors[["occ"]], alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_x) +
            ggplot2::coord_flip() +
            ggplot2::scale_y_continuous(n.breaks = 3) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(size = 5))

          env_y_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[ if (is.numeric(y)) names(env_bg)[y] else y ]])) +
            ggplot2::geom_density(fill =  colors[["occ"]], alpha = 0.6) +
            ggplot2::scale_x_continuous(limits = rng_y) +
            ggplot2::coord_flip() +
            ggplot2::scale_y_continuous(n.breaks = 3) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(size = 5))

          return_plot <- ggpubr::ggarrange(
            NULL, env_y_top, env_z_top, NULL,
            x_name, ell_y_x, ell_z_x, env_x_right,
            legend_plot, y_name, ell_z_y, env_y_right,
            NULL, angle_warn, z_name, NULL,
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

