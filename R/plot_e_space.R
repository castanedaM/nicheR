#' Plot Environmental Space with Optional Ellipsoid Overlays
#'
#' Produces pairwise views of a 3D environmental space with optional overlays of a
#' virtual niche ellipsoid and occurrence points. In 2D mode, it returns a grid
#' of pairwise scatterplots with projected ellipse boundaries. In 3D mode, it
#' returns an interactive \code{plotly} scatterplot.
#'
#' @param env_bg A \code{data.frame} (or coercible object) of background environments
#'   with at least three numeric predictor columns. These columns must contain
#'   the variables referenced by \code{x}, \code{y}, and \code{z}. If \code{env_bg}
#'   is a \code{terra::SpatRaster} or \code{raster::Raster*}, it is converted with
#'   \code{as.data.frame.nicheR()}. If \code{NULL}, a \code{NicheR_species} object
#'   supplied via \code{vs} can be used to infer \code{env_bg} from its
#'   \code{suitability} slot.
#' @param x,y,z Column specifications for the three predictors to display. Each may
#'   be a single column name (character string) or a single 1-based integer index
#'   into \code{env_bg}. If any of \code{x}, \code{y}, or \code{z} are omitted,
#'   the function attempts to infer them from the predictor columns in \code{env_bg}
#'   (excluding \code{"x"} and \code{"y"} coordinate columns if present).
#' @param labels Character vector of length 3 giving axis labels for the x, y,
#'   and z variables in display order. Defaults to \code{c("ENV 1", "ENV 2", "ENV 3")}.
#' @param n_bg Positive integer giving the maximum number of background rows to plot.
#'   If \code{nrow(env_bg)} is greater than \code{n_bg}, a random subset of size
#'   \code{n_bg} is drawn. Using a large \code{n_bg} may slow plotting.
#' @param niche Optional object of class \code{ellipsoid} describing the niche. If
#'   provided, its boundary and center will be plotted. For 2D plots, the
#'   object should contain \code{niche$angles}. If \code{NULL} and a
#'   \code{NicheR_species} object is supplied via \code{vs}, \code{niche} is
#'   filled from \code{vs$niche}.
#' @param suitable_env Optional suitable environment object (data.frame,
#'   \code{suitable_env} list, or SpatRaster), typically produced by
#'   \code{\link{get_suitable_env}}. When provided, points inside the ellipsoid
#'   are plotted in E space.
#' @param occ_pts Optional \code{data.frame} of occurrence points that includes the
#'   same predictor columns used for \code{x}, \code{y}, and \code{z}. These are
#'   overplotted if supplied. If \code{NULL} and a \code{NicheR_species} object
#'   is supplied via \code{vs}, this is filled from \code{vs$occurrences} when available.
#' @param rand_seed Integer used to set the random number generator seed for
#'   reproducible background downsampling.
#' @param show.occ.density Logical. If \code{TRUE} and \code{occ_pts} is provided,
#'   adds marginal density panels for each variable (2D plots only;
#'   \code{plot.3d = FALSE}).
#' @param plot.3d Logical. If \code{TRUE}, returns an interactive \code{plotly} 3D scatter
#'   plot. If \code{FALSE} (the default), returns a static \code{ggpubr} grid of 2D panels.
#' @param colors Optional named list of colors to override the default palette.
#'   Valid names are: \code{bg}, \code{ellipsoid}, \code{centroid}, \code{tolerance},
#'   \code{suitable_env}, \code{occ}.
#' @param palette Character name of the internal palette to use.
#'   One of \code{"default"}, \code{"palette2"}, ..., \code{"palette6"}.
#' @param vs Optional \code{NicheR_species} object returned by
#'   \code{\link{create_virtual_species}}. If provided, \code{niche},
#'   \code{occ_pts}, \code{env_bg}, and \code{suitable_env} are auto-filled from
#'   this object when they are not supplied explicitly.
#'
#' @return
#' If \code{plot.3d = TRUE}, a \code{plotly} object.
#' If \code{plot.3d = FALSE}, a \code{ggpubr} object containing arranged
#' \code{ggplot2} panels.
#'
#' @family plotting functions
#' @seealso \code{\link{validate_plot_e_space_args}}, \code{\link{build_ellps}},
#'   \code{\link{get_suitable_env}}, \code{\link{create_virtual_species}}
#' @import RColorBrewer ggplot2 dplyr
#' @importFrom rlang .data
#' @export
plot_e_space <- function(env_bg,
                         x, y, z,
                         labels = c("ENV 1", "ENV 2", "ENV 3"),
                         n_bg = 10000,
                         niche = NULL,
                         suitable_env = NULL,
                         occ_pts = NULL,
                         rand_seed = 1234,
                         show.occ.density = FALSE, # only for 2D plots
                         plot.3d = FALSE,
                         colors = NULL,
                         palette = "default",
                         vs = NULL) {

  # Track whether x,y,z were auto-inferred
  x_missing <- missing(x)
  y_missing <- missing(y)
  z_missing <- missing(z)
  auto_inferred <- x_missing || y_missing || z_missing

  # --- 0. Optionally fill from NicheR_species object -------------------------

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
        env_bg <- as.data.frame.nicheR(vs$suitability)
      } else if (is.data.frame(vs$suitability)) {
        env_bg <- vs$suitability
      }
    }

    # Fill suitable_env if missing, using suitability slot directly
    if (is.null(suitable_env) && !is.null(vs$suitability)) {
      suitable_env <- vs$suitability
    }
  }

  if (missing(env_bg) || is.null(env_bg)) {
    stop("'env_bg' must be supplied, or inferable from 'vs$suitability'.")
  }

  # --- 0.1 Coerce env_bg to data.frame --------------------------------------

  if (inherits(env_bg, "tbl_df")) {
    env_bg <- as.data.frame(env_bg)
  }
  if (inherits(env_bg, "Raster")) {
    env_bg <- terra::rast(env_bg)
  }
  if (inherits(env_bg, "SpatRaster")) {
    env_bg <- as.data.frame.nicheR(env_bg)
  }
  if (!is.data.frame(env_bg)) {
    stop("'env_bg' must be a data.frame or coercible to one.")
  }

  # --- 0.2 Auto-infer x,y,z from env_bg if needed ---------------------------

  if (auto_inferred) {

    if (all(c("x", "y") %in% names(env_bg))) {
      candidate_vars <- setdiff(names(env_bg), c("x", "y"))
    } else {
      candidate_vars <- names(env_bg)
    }

    if (length(candidate_vars) < 3) {
      stop(
        "Could not infer x, y, z from 'env_bg' (fewer than 3 predictor columns after ",
        "removing any 'x'/'y' coordinates). Please provide x, y, and z explicitly."
      )
    }

    # Use first 3 predictor columns as x,y,z
    inferred <- candidate_vars[seq_len(3)]
    x <- inferred[1]
    y <- inferred[2]
    z <- inferred[3]

    message(
      "No complete x, y, z specification provided. Using predictor columns inferred from 'env_bg': ",
      paste(inferred, collapse = ", ")
    )
  }

  # --- 1. Validate core arguments and resolve predictor names ----------------

  v <- validate_plot_e_space_args(
    env_bg, x, y, z,
    labels, n_bg,
    niche,
    occ_pts, show.occ.density
  )

  col_x <- v$col_names[1]
  col_y <- v$col_names[2]
  col_z <- v$col_names[3]

  # If user did provide x,y,z but as indices/variants, still nice to tell them
  if (auto_inferred) {
    message(
      "Using predictor columns: ",
      paste(v$col_names, collapse = ", ")
    )
  }

  # Coerce occ_pts to data.frame (structure already checked in validator)
  if (!is.null(occ_pts)) {
    occ_pts <- as.data.frame(occ_pts)
  }

  # --- 2. Coerce suitable_env to a data.frame of inside points (if given) ----

  pts_in <- NULL
  if (!is.null(suitable_env)) {

    if (inherits(suitable_env, "suitable_env") ||
        (is.list(suitable_env) &&
         any(c("suitable_env_df", "suitable_env_sp") %in% names(suitable_env)))) {

      if ("suitable_env_df" %in% names(suitable_env) &&
          is.data.frame(suitable_env$suitable_env_df)) {

        pts_in <- suitable_env$suitable_env_df

      } else if ("suitable_env_sp" %in% names(suitable_env)) {

        sp <- suitable_env$suitable_env_sp
        if (inherits(sp, "Raster")) sp <- terra::rast(sp)

        if (inherits(sp, "SpatRaster")) {
          pts_in <- as.data.frame.nicheR(sp)
        } else if (is.list(sp) &&
                   length(sp) > 0 &&
                   all(vapply(sp, inherits, logical(1), "SpatRaster"))) {

          if ("suitable" %in% names(sp)) {
            pts_in <- as.data.frame.nicheR(sp[["suitable"]])
          } else {
            pts_in <- as.data.frame.nicheR(sp[[1]])
          }
        }
      }

    } else if (inherits(suitable_env, "Raster") || inherits(suitable_env, "SpatRaster")) {

      if (inherits(suitable_env, "Raster")) {
        suitable_env <- terra::rast(suitable_env)
      }
      pts_in <- as.data.frame.nicheR(suitable_env)

    } else if (is.data.frame(suitable_env) || is.matrix(suitable_env)) {

      pts_in <- as.data.frame(suitable_env)
    }

    if (!is.null(pts_in)) {
      if (!all(c(col_x, col_y, col_z) %in% names(pts_in))) {
        warning(
          "suitable_env does not contain all of x, y, z predictor columns; ",
          "suitable points will not be plotted in E-space."
        )
        pts_in <- NULL
      }
    }
  }

  # --- 3. Colors / palette ---------------------------------------------------

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

  if (!palette %in% names(palettes)) {
    stop("Unknown palette '", palette, "'.")
  }

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
    colors <- utils::modifyList(base_colors, colors)
  }

  # Legend options ------------------------------------------------------------

  opts <- list(
    background_point     = TRUE,
    trace_line           = !is.null(niche),
    centroid_point       = !is.null(niche),
    tolerance_range_line = !is.null(niche),
    suitable_point       = !is.null(pts_in),
    occurrence_point     = !is.null(occ_pts)
  )

  legend_items <- data.frame(
    id       = c("background_point","trace_line","centroid_point",
                 "tolerance_range_line","suitable_point","occurrence_point"),
    type     = c("point","line","point","line","point","point"),
    label    = c("Background environments","Niche boundary","Niche centroid",
                 "Tolerance ranges","Suitable environments","Occurrences"),
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

  # Build legend plot ---------------------------------------------------------

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
        row     = dplyr::row_number(),
        y       = top_y - (row - 1) * spacing,
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
        data  = dplyr::filter(legend_items, type == "point"),
        ggplot2::aes(x = x_point, y = y, colour = color),
        size = 2
      ) +
      ggplot2::geom_segment(
        data  = dplyr::filter(legend_items, type == "line"),
        ggplot2::aes(x = x0_line, xend = x1_line, y = y, yend = y,
                     colour = color, linetype = linetype),
        linewidth = 0.5
      ) +
      ggplot2::geom_text(
        data  = legend_items,
        ggplot2::aes(x = x_text, y = y, label = label),
        hjust = 0
      ) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_linetype_identity()
  }

  # Downsample background -----------------------------------------------------

  if (nrow(env_bg) > n_bg) {
    message(sprintf("Sampling %d of %d rows from 'env_bg' for plotting.", n_bg, nrow(env_bg)))
    set.seed(rand_seed)
    env_bg <- env_bg[sample.int(nrow(env_bg), size = n_bg, replace = FALSE), ]
  }

  # ---------------------------------------------------------------------------
  # 3D plotting branch
  # ---------------------------------------------------------------------------

  if (isTRUE(plot.3d)) {

    p3 <- plotly::plot_ly(
      data = env_bg,
      x    = env_bg[[col_x]],
      y    = env_bg[[col_y]],
      z    = env_bg[[col_z]],
      type = "scatter3d",
      mode = "markers",
      marker = list(color = colors[["bg"]], size = 2),
      name   = "Background Environments"
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

    # Ellipsoid surface + centroid (if present)
    if (!is.null(niche) && !is.null(niche$surface)) {

      surf <- as.data.frame(niche$surface)
      if (ncol(surf) >= 3) {
        p3 <- p3 %>%
          plotly::add_trace(
            data = surf,
            x    = surf[[1]],
            y    = surf[[2]],
            z    = surf[[3]],
            type = "scatter3d",
            mode = "lines",
            line = list(color = colors[["ellipsoid"]]),
            name = "Niche Boundary",
            inherit = FALSE
          )
      }

      if (!is.null(niche$center) && length(niche$center) >= 3) {
        p3 <- p3 %>%
          plotly::add_markers(
            x = niche$center[1],
            y = niche$center[2],
            z = niche$center[3],
            marker = list(color = colors[["centroid"]], size = 5),
            name   = "Niche Centroid"
          )
      }

      if (!is.null(pts_in)) {
        p3 <- p3 %>%
          plotly::add_markers(
            data = pts_in,
            x    = pts_in[[col_x]],
            y    = pts_in[[col_y]],
            z    = pts_in[[col_z]],
            marker = list(color = colors[["suitable_env"]], size = 3),
            name   = "Suitable Environments",
            inherit = FALSE
          )
      }
    }

    # Occurrence points in 3D
    if (!is.null(occ_pts)) {
      p3 <- p3 %>%
        plotly::add_markers(
          data = occ_pts,
          x    = occ_pts[[col_x]],
          y    = occ_pts[[col_y]],
          z    = occ_pts[[col_z]],
          marker = list(color = colors[["occ"]], size = 3),
          name   = "Sampled Occurrences",
          inherit = FALSE
        ) %>%
        plotly::layout(
          title  = "Virtual Niche and Sampled Occurrences in E-space",
          legend = list(x = 0.05, y = 0.95)
        )
    }

    return(p3)
  }

  # ---------------------------------------------------------------------------
  # 2D plotting branch
  # ---------------------------------------------------------------------------

  p_main_y_x <- ggplot2::ggplot(
    env_bg,
    ggplot2::aes(x = .data[[col_y]], y = .data[[col_x]])
  ) +
    ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  p_main_z_x <- ggplot2::ggplot(
    env_bg,
    ggplot2::aes(x = .data[[col_z]], y = .data[[col_x]])
  ) +
    ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_blank())

  p_main_z_y <- ggplot2::ggplot(
    env_bg,
    ggplot2::aes(x = .data[[col_z]], y = .data[[col_y]])
  ) +
    ggplot2::geom_point(alpha = 0.5, color = colors[["bg"]], pch = ".") +
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
    x_name,      p_main_y_x, p_main_z_x,
    legend_plot, y_name,     p_main_z_y,
    NULL,        NULL,       z_name,
    ncol = 3, nrow = 3,
    widths  = c(0.15, 0.425, 0.425),
    heights = c(0.45, 0.45, 0.15)
  )

  # --- Ellipsoid overlays and extras (2D) ------------------------------------

  if (!is.null(niche)) {

    center_y_x <- c(niche$center[2], niche$center[1])
    axes_y_x   <- c(niche$axes[2],   niche$axes[1])
    angle_y_x  <- c(niche$angles[2], niche$angles[1])

    center_z_x <- c(niche$center[3], niche$center[1])
    axes_z_x   <- c(niche$axes[3],   niche$axes[1])
    angle_z_x  <- c(niche$angles[3], niche$angles[1])

    center_z_y <- c(niche$center[3], niche$center[2])
    axes_z_y   <- c(niche$axes[3],   niche$axes[2])
    angle_z_y  <- c(niche$angles[3], niche$angles[2])

    ell2d_y_x <- build_ellps(center = center_y_x, axes = axes_y_x, angles = angle_y_x)
    ell2d_z_x <- build_ellps(center = center_z_x, axes = axes_z_x, angles = angle_z_x)
    ell2d_z_y <- build_ellps(center = center_z_y, axes = axes_z_y, angles = angle_z_y)

    ell_y_x <- p_main_y_x
    ell_z_x <- p_main_z_x
    ell_z_y <- p_main_z_y

    if (any(niche$angles != 0)) {
      angle_warn <- ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::geom_text(
          ggplot2::aes(0, 0,
                       label = "Note: The ellipsoid is angled; its shape may appear distorted and some points may fall outside due to dimensionality."),
          size = 2
        )
    } else {
      angle_warn <- NULL
    }

    # Suitable points (inside)
    if (!is.null(pts_in)) {

      ell_y_x <- ell_y_x +
        ggplot2::geom_point(
          data = pts_in,
          ggplot2::aes(
            x = .data[[col_y]],
            y = .data[[col_x]]
          ),
          color = colors[["suitable_env"]],
          size  = 0.5
        )

      ell_z_x <- ell_z_x +
        ggplot2::geom_point(
          data = pts_in,
          ggplot2::aes(
            x = .data[[col_z]],
            y = .data[[col_x]]
          ),
          color = colors[["suitable_env"]],
          size  = 0.5
        )

      ell_z_y <- ell_z_y +
        ggplot2::geom_point(
          data = pts_in,
          ggplot2::aes(
            x = .data[[col_z]],
            y = .data[[col_y]]
          ),
          color = colors[["suitable_env"]],
          size  = 0.5
        )
    }

    # Occurrence points in 2D
    if (!is.null(occ_pts)) {
      ell_y_x <- ell_y_x +
        ggplot2::geom_point(
          data = occ_pts,
          ggplot2::aes(
            x = .data[[col_y]],
            y = .data[[col_x]]
          ),
          color = colors[["occ"]],
          size  = 0.5
        )

      ell_z_x <- ell_z_x +
        ggplot2::geom_point(
          data = occ_pts,
          ggplot2::aes(
            x = .data[[col_z]],
            y = .data[[col_x]]
          ),
          color = colors[["occ"]],
          size  = 0.5
        )

      ell_z_y <- ell_z_y +
        ggplot2::geom_point(
          data = occ_pts,
          ggplot2::aes(
            x = .data[[col_z]],
            y = .data[[col_y]]
          ),
          color = colors[["occ"]],
          size  = 0.5
        )
    }

    # Ellipse boundaries + axes + centroid -----------------------------------

    ell_y_x <- ell_y_x +
      ggplot2::geom_path(
        data   = ell2d_y_x$surface,
        ggplot2::aes(x, y),
        color = colors[["ellipsoid"]],
        linewidth = 0.5
      ) +
      ggplot2::annotate(
        "segment",
        x    = ell2d_y_x$center[1] - ell2d_y_x$axes[1],
        xend = ell2d_y_x$center[1] + ell2d_y_x$axes[1],
        y    = ell2d_y_x$center[2],
        yend = ell2d_y_x$center[2],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "segment",
        y    = ell2d_y_x$center[2] - ell2d_y_x$axes[2],
        yend = ell2d_y_x$center[2] + ell2d_y_x$axes[2],
        x    = ell2d_y_x$center[1],
        xend = ell2d_y_x$center[1],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "point",
        x = ell2d_y_x$center[1],
        y = ell2d_y_x$center[2],
        color = colors[["centroid"]],
        size  = 2
      )

    ell_z_x <- ell_z_x +
      ggplot2::geom_path(
        data   = ell2d_z_x$surface,
        ggplot2::aes(x, y),
        color = colors[["ellipsoid"]],
        linewidth = 0.5
      ) +
      ggplot2::annotate(
        "segment",
        x    = ell2d_z_x$center[1] - ell2d_z_x$axes[1],
        xend = ell2d_z_x$center[1] + ell2d_z_x$axes[1],
        y    = ell2d_z_x$center[2],
        yend = ell2d_z_x$center[2],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "segment",
        y    = ell2d_z_x$center[2] - ell2d_z_x$axes[2],
        yend = ell2d_z_x$center[2] + ell2d_z_x$axes[2],
        x    = ell2d_z_x$center[1],
        xend = ell2d_z_x$center[1],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "point",
        x = ell2d_z_x$center[1],
        y = ell2d_z_x$center[2],
        color = colors[["centroid"]],
        size  = 2
      )

    ell_z_y <- ell_z_y +
      ggplot2::geom_path(
        data   = ell2d_z_y$surface,
        ggplot2::aes(x, y),
        color = colors[["ellipsoid"]],
        linewidth = 0.5
      ) +
      ggplot2::annotate(
        "segment",
        x    = ell2d_z_y$center[1] - ell2d_z_y$axes[1],
        xend = ell2d_z_y$center[1] + ell2d_z_y$axes[1],
        y    = ell2d_z_y$center[2],
        yend = ell2d_z_y$center[2],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "segment",
        y    = ell2d_z_y$center[2] - ell2d_z_y$axes[2],
        yend = ell2d_z_y$center[2] + ell2d_z_y$axes[2],
        x    = ell2d_z_y$center[1],
        xend = ell2d_z_y$center[1],
        color = colors[["tolerance"]],
        linetype = "dashed"
      ) +
      ggplot2::annotate(
        "point",
        x = ell2d_z_y$center[1],
        y = ell2d_z_y$center[2],
        color = colors[["centroid"]],
        size  = 2
      )

    # Re-assemble main grid with ellipsoid overlays --------------------------

    return_plot <- ggpubr::ggarrange(
      x_name,      ell_y_x,   ell_z_x,
      legend_plot, y_name,    ell_z_y,
      NULL,        angle_warn, z_name,
      ncol = 3, nrow = 3,
      widths  = c(0.15, 0.425, 0.425),
      heights = c(0.45, 0.45, 0.15)
    )

    # Optional occurrence density panels -------------------------------------

    if (isTRUE(show.occ.density) && !is.null(occ_pts)) {

      rng_z <- range(env_bg[[col_z]], na.rm = TRUE)
      rng_y <- range(env_bg[[col_y]], na.rm = TRUE)
      rng_x <- range(env_bg[[col_x]], na.rm = TRUE)

      env_z_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[col_z]])) +
        ggplot2::geom_density(fill = colors[["occ"]], alpha = 0.6) +
        ggplot2::scale_x_continuous(limits = rng_z) +
        ggplot2::scale_y_continuous(n.breaks = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.y  = ggplot2::element_text(size = 5)
        )

      env_y_top <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[col_y]])) +
        ggplot2::geom_density(fill = colors[["occ"]], alpha = 0.6) +
        ggplot2::scale_x_continuous(limits = rng_y) +
        ggplot2::scale_y_continuous(n.breaks = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.y  = ggplot2::element_text(size = 5)
        )

      env_x_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[col_x]])) +
        ggplot2::geom_density(fill = colors[["occ"]], alpha = 0.6) +
        ggplot2::scale_x_continuous(limits = rng_x) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(n.breaks = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x  = ggplot2::element_text(size = 5)
        )

      env_y_right <- ggplot2::ggplot(occ_pts, ggplot2::aes(x = .data[[col_y]])) +
        ggplot2::geom_density(fill = colors[["occ"]], alpha = 0.6) +
        ggplot2::scale_x_continuous(limits = rng_y) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(n.breaks = 3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x  = ggplot2::element_text(size = 5)
        )

      return_plot <- ggpubr::ggarrange(
        NULL,        env_y_top, env_z_top, NULL,
        x_name,      ell_y_x,   ell_z_x,   env_x_right,
        legend_plot, y_name,    ell_z_y,   env_y_right,
        NULL,        angle_warn, z_name,   NULL,
        ncol = 4, nrow = 4,
        widths  = c(0.1, 0.4, 0.4, 0.1),
        heights = c(0.1, 0.4, 0.4, 0.1)
      )
    }
  }

  return(return_plot)
}
