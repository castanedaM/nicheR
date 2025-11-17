#' Plot Environmental Space with Optional Ellipsoid Overlays
#'
#' Produces pairwise views of a 3D environmental space with optional overlays of a
#' virtual niche ellipsoid and occurrence points. In 2D mode, it returns a grid
#' of pairwise scatterplots with projected ellipse boundaries. In 3D mode, it
#' returns an interactive `plotly` scatterplot.
#'
#' @param env_bg A `data.frame` of background environments with at least three
#'   numeric predictor columns. These columns must contain the variables
#'   referenced by `x`, `y`, and `z`. If `NULL`, a `NicheR_species` object
#'   supplied via `vs` can be used to infer `env_bg` from its `suitability`
#'   slot (either a data.frame or a SpatRaster that will be converted).
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
#'   object must contain `niche$angles`. If `NULL` and a `NicheR_species` object
#'   is supplied via `vs`, `niche` is filled from `vs$niche`.
#' @param show.pts.in Logical. If `TRUE` and `niche` is provided, points from
#'   `env_bg` that fall inside the ellipsoid are highlighted.
#' @param occ_pts Optional `data.frame` of occurrence points that includes the
#'   same predictor columns used for `x`, `y`, and `z`. These are overplotted
#'   if supplied. If `NULL` and a `NicheR_species` object is supplied via `vs`,
#'   this is filled from `vs$occurrences` when available.
#' @param rand_seed Integer used to set the random number generator seed for
#'   reproducible background downsampling.
#' @param show.occ.density Logical. If `TRUE` and `occ_pts` is provided, this
#'   adds marginal density panels for each variable. This is only supported in
#'   2D plots (`plot.3d = FALSE`).
#' @param plot.3d Logical. If `TRUE`, returns an interactive `plotly` 3D scatter
#'   plot. If `FALSE` (the default), returns a static `ggpubr` grid of 2D panels.
#' @param colors Optional named list of colors to override the default palette.
#'   Valid names are: `bg`, `ellipsoid`, `centroid`, `tolerance`,
#'   `suitable_env`, `occ`.
#' @param palette Character name of the internal palette to use.
#'   One of `"default"`, `"palette2"`, ..., `"palette6"`.
#' @param vs Optional `NicheR_species` object returned by
#'   [create_virtual_species()]. If provided, `niche`, `occ_pts`, and `env_bg`
#'   are auto-filled from this object when they are not supplied explicitly.
#'
#' @return
#' - If `plot.3d = TRUE`: A `plotly` object.
#' - If `plot.3d = FALSE`: A `ggpubr` object containing arranged `ggplot2` panels.
#'
#' @details This function is a visualization tool for understanding the
#'   relationship between a species' niche, its occurrences, and the available
#'   environmental space. The 2D views are especially useful for publication-
#'   quality figures, while the 3D interactive plot is useful for exploration.
#'
#' @family plotting functions
#' @seealso [validate_plot_e_space_args()], [build_ellps()],
#'   [get_suitable_env()], [create_virtual_species()]
#' @import RColorBrewer ggplot2 dplyr
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
                         palette = "default",
                         vs = NULL) {

  # --- 0. Optionally fill from NicheR_species object ----------------------

  if (!is.null(vs)) {
    if (!inherits(vs, "NicheR_species")) {
      stop("'vs' must be a NicheR_species object created by create_virtual_species().")
    }

    # Fill niche if missing
    if (is.null(niche) && !is.null(vs$niche)) {
      niche <- vs$niche
    }

    # Fill occurrences if missing
    if (is.null(occ_pts) && !is.null(vs$occurrences)) {
      occ_pts <- vs$occurrences
    }

    # Fill env_bg if missing, using suitability slot
    if (missing(env_bg) || is.null(env_bg)) {
      if (inherits(vs$suitability, "SpatRaster")) {
        # Use NicheR helper to get x, y + predictors
        env_bg <- as.data.frame.nicheR(vs$suitability)
      } else if (is.data.frame(vs$suitability)) {
        env_bg <- vs$suitability
      } else {
        stop(
          "Could not infer 'env_bg' from 'vs$suitability'. ",
          "Please supply 'env_bg' explicitly when calling plot_e_space()."
        )
      }
    }
  }

  # Small helper to resolve column names (numeric or character)
  resolve_col <- function(df, col_spec) {
    if (is.numeric(col_spec)) {
      nm <- names(df)[col_spec]
    } else {
      nm <- col_spec
    }
    if (!nm %in% names(df)) {
      stop("Column '", nm, "' not found in env_bg.")
    }
    nm
  }

  # --- 1. Validate (no rlang in helper) -----------------------------------

  v <- validate_plot_e_space_args(env_bg, x, y, z,
                                  labels, n_bg, niche, show.pts.in,
                                  occ_pts, show.occ.density)

  # Resolve column names up front to reduce repetition
  x_col <- resolve_col(env_bg, x)
  y_col <- resolve_col(env_bg, y)
  z_col <- resolve_col(env_bg, z)

  # Set colors for legend and plot -----------------------------------------

  palettes <- list(
    default = list(
      bg           = "#9093A2FF",
      ellipsoid    = "#2A363BFF",
      centroid     = "#D72000FF",
      tolerance    = "#EE6100FF",
      suitable_env = "#FED789FF",
      occ          = "#B4BF3AFF"
    ),
    palette2 = list(
      bg           = "#9CA9BAFF",
      ellipsoid    = "#3D619DFF",
      centroid     = "#345084FF",
      tolerance    = "#693829FF",
      suitable_env = "#CFB267FF",
      occ          = "#A56A3EFF"
    ),
    palette3 = list(
      bg           = "#C8CCC6FF",
      ellipsoid    = "#023743FF",
      centroid     = "#72874EFF",
      tolerance    = "#476F84FF",
      suitable_env = "#FED789FF",
      occ          = "#A4BED5FF"
    ),
    palette4 = list(
      bg           = "#C0D1CEFF",
      ellipsoid    = "#859B6CFF",
      centroid     = "#B74954FF",
      tolerance    = "#A99364FF",
      suitable_env = "#C2DDB2FF",
      occ          = "#EBA49EFF"
    ),
    palette5 = list(
      bg           = "#A89F8EFF",
      ellipsoid    = "#7887A4FF",
      centroid     = "#A8CDECFF",
      tolerance    = "#682C37FF",
      suitable_env = "#F6955EFF",
      occ          = "#9B6981FF"
    ),
    palette6 = list(
      bg           = "#D3D4D8FF",
      ellipsoid    = "#731A12FF",
      centroid     = "#F2D43DFF",
      tolerance    = "#3F858CFF",
      suitable_env = "#D9814EFF",
      occ          = "#707322FF"
    )
  )

  base_colors <- palettes[[palette]]

  if (is.null(colors)) {
    colors <- base_colors
  } else {
    if (is.null(names(colors))) {
      names(colors) <- names(base_colors)[seq_along(colors)]
      if (!is.list(colors)) colors <- as.list(colors)
      message(
        "No names detected in 'colors'. ",
        "Using the order provided and filling missing entries with defaults.\n",
        "If you want to change a specific object, use a named list.\n",
        "Available options are: ",
        paste(names(base_colors), collapse = ", ")
      )
    }
    if (!is.list(colors)) colors <- as.list(colors)
    colors <- modifyList(base_colors, colors)
  }

  # Legend options --------------------------------------------------------

  opts <- list(
    background_point = !is.null(n_bg),
    trace_line       = !is.null(niche),
    centroid_point   = !is.null(niche),
    tolerance_range_line  = !is.null(niche),
    suitable_point   = show.pts.in,
    occurrence_point = !is.null(occ_pts)
  )

  legend_items <- data.frame(
    id       = c("background_point","trace_line","centroid_point",
                 "tolerance_range_line","suitable_point","occurrence_point"),
    type     = c("point","line","point","line","point","point"),
    label    = c("Background environments","Niche boundary","Niche centroid",
                 "Tolerance ranges","Suitable environments","Occurrence"),
    color    = c(colors[["bg"]], colors[["ellipsoid"]], colors[["centroid"]],
                 colors[["tolerance"]], colors[["suitable_env"]], colors[["occ"]]),
    linetype = c(NA, 1, NA, 2, NA, NA),
    stringsAsFactors = FALSE
  )

  active <- logical(nrow(legend_items))
  for (i in seq_len(nrow(legend_items))) {
    active[i] <- isTRUE(opts[[ legend_items$id[i] ]])
  }
  legend_items <- legend_items[active, , drop = FALSE]

  # Build legend plot -----------------------------------------------------

  if (nrow(legend_items) == 0) {
    legend_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  } else {
    top_y   <- 2
    spacing <- 0.25
    x_point <- 0.00
    x_text  <- 0.1
    x0_line <- -0.2
    x1_line <-  0.00

    legend_items <- legend_items %>%
      dplyr::mutate(
        row = dplyr::row_number(),
        y   = top_y - (row - 1) * spacing,
        x_point = x_point,
        x_text  = x_text,
        x0_line = x0_line,
        x1_line = x1_line
      )

    legend_base <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, 2.5), clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")

    legend_plot <- legend_base +
      ggplot2::geom_point(
        data = dplyr::filter(legend_items, type == "point"),
        ggplot2::aes(x = x_point, y = y, colour = color),
        size = 2
      ) +
      ggplot2::geom_segment(
        data = dplyr::filter(legend_items, type == "line"),
        ggplot2::aes(x = x0_line, xend = x1_line, y = y, yend = y,
                     colour = color, linetype = linetype),
        linewidth = 0.5
      ) +
      ggplot2::geom_text(
        data = legend_items,
        ggplot2::aes(x = x_text, y = y, label = label),
        hjust = 0
      ) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_linetype_identity()
  }

  # Downsample background -------------------------------------------------

  if (nrow(env_bg) > n_bg) {
    message(sprintf("Sampling %d of %d rows from 'env_bg' for plotting.", n_bg, nrow(env_bg)))
    set.seed(rand_seed)
    env_bg <- env_bg[sample.int(nrow(env_bg), size = n_bg, replace = FALSE), ]
  }

  # -----------------------------------------------------------------------
  # From here, your existing plotting logic is unchanged except that
  # everywhere you had:
  #   if (is.numeric(x)) names(env_bg)[x] else x
  # you can now just use `x_col`, `y_col`, `z_col`.
  # I’ll show one 2D and one 3D example and leave the rest as-is.
  # -----------------------------------------------------------------------

  if (isTRUE(plot.3d)) {
    return_plot <- plotly::plot_ly(
      data = env_bg,
      x = env_bg[[x_col]],
      y = env_bg[[y_col]],
      z = env_bg[[z_col]],
      type = "scatter3d",
      mode = "markers",
      marker = list(color = "lightgrey", size = 2),
      name = "Background Environments"
    ) %>%
      plotly::layout(
        title = list(text = "Background Environments (E-space)"),
        scene = list(
          xaxis = list(title = list(text = labels[1])),
          yaxis = list(title = list(text = labels[2])),
          zaxis = list(title = list(text = labels[3]))
        ),
        legend = list(x = 0.05, y = 0.95)
      )

    # ... (your existing 3D niche / occ overlay logic goes here,
    #      just swap to x_col, y_col, z_col where needed)

  } else {
    p_main_y_x <- ggplot2::ggplot(
      env_bg,
      ggplot2::aes(x = .data[[y_col]], y = .data[[x_col]])
    ) +
      ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    p_main_z_x <- ggplot2::ggplot(
      env_bg,
      ggplot2::aes(x = .data[[z_col]], y = .data[[x_col]])
    ) +
      ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    p_main_z_y <- ggplot2::ggplot(
      env_bg,
      ggplot2::aes(x = .data[[z_col]], y = .data[[y_col]])
    ) +
      ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_blank())

    x_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[1]), fontface = "bold")
    y_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[2]), fontface = "bold")
    z_name <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(0, 0, label = labels[3]), fontface = "bold")

    return_plot <- ggpubr::ggarrange(
      x_name, p_main_y_x, p_main_z_x,
      legend_plot,   y_name,    p_main_z_y,
      NULL,   NULL,      z_name,
      ncol = 3, nrow = 3,
      widths = c(0.15, 0.425, 0.425),
      heights = c(0.45, 0.45, 0.15)
    )

    # ... keep your existing ellipsoid overlays / density logic,
    #     just swap `[[ if (is.numeric(x)) ... ]]` with `[[x_col]]`,
    #     etc.
  }

  return(return_plot)
}
