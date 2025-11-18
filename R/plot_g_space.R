#' Plot Geographic Space with Optional Suitability/Distance and Occurrences
#'
#' Draws a basemap and overlays (a) suitable environments as tiles,
#' (b) an optional continuous distance-to-centroid surface, and/or (c)
#' occurrence points. The map extent is primarily inferred from \code{env_bg}
#' (if available), otherwise from the coordinates in \code{suitable_env} and
#' finally from \code{occ_pts}, always with a small buffer where needed.
#'
#' @param env_bg Optional background environment object used when
#'   \code{suitable_env} is not supplied and suitability must be computed.
#'   Typically a \code{terra::SpatRaster} (or \code{raster::Raster*}).
#' @param suitable_env Optional precomputed suitable environment object. Can be:
#'   \itemize{
#'     \item a \code{data.frame} with columns \code{x}, \code{y}, and optionally \code{dist_sq},
#'     \item a \code{"suitable_env"} object returned by \code{\link{get_suitable_env}}
#'           with a \code{suitable_env_df} or \code{suitable_env_sp} component,
#'     \item a \code{terra::SpatRaster} or a list of \code{terra::SpatRaster} objects.
#'   }
#'   When a data.frame is used, it is assumed that each row corresponds to a
#'   suitable cell (i.e. presence-only). When \code{dist_sq} is present, a
#'   distance panel can also be drawn.
#' @param occ_pts Optional \code{data.frame} of occurrences with columns
#'   \code{x}, \code{y} (assumed longitude/latitude, WGS84). Plotted as points
#'   if supplied. If \code{vs} is supplied and \code{occ_pts} is \code{NULL},
#'   the function will try to use \code{vs$occurrences}.
#' @param show.occ.density Logical (currently unused; reserved for future
#'   density panels).
#' @param colors Optional named list to override aesthetics. Recognized names:
#'   \code{bg} (basemap fill), \code{suitable_env} (suitable tiles),
#'   \code{occ_fill}, \code{occ_stroke}, and \code{dist} (RColorBrewer palette
#'   name for distance tiles).
#' @param palette Character palette key. One of \code{"default"},
#'   \code{"palette2"}, or \code{"palette3"}.
#' @param surface Character; which surface(s) to plot. One of:
#'   \itemize{
#'     \item \code{"both"}: show suitability tiles and (if available) distance tiles,
#'     \item \code{"suit"}: show only suitability tiles,
#'     \item \code{"dist"}: show only the distance-to-centroid surface
#'           (requires a \code{dist_sq} column).
#'   }
#' @param vs Optional object of class \code{"NicheR_species"} created by
#'   \code{\link{create_virtual_species}}. When provided, any of
#'   \code{env_bg}, \code{niche}, \code{occ_pts}, or \code{suitable_env} that
#'   are \code{NULL} will be filled from this object if possible.
#' @param niche Optional object of class \code{"ellipsoid"} (from
#'   \code{\link{build_ellps}}). Required when suitability is computed internally
#'   (i.e., when \code{suitable_env} is \code{NULL}).
#'
#' @return A \code{ggpubr} object produced by \code{\link[ggpubr]{ggarrange}},
#'   containing one or two map panels plus a matching legend panel.
#'
#' @family plotting functions
#' @seealso \code{\link{build_ellps}}, \code{\link{get_suitable_env}},
#'   \code{\link{create_virtual_species}}
#'
#' @export
plot_g_space <- function(env_bg = NULL,
                         suitable_env = NULL,
                         occ_pts = NULL,
                         show.occ.density = FALSE,
                         colors = NULL,
                         palette = "default",
                         surface = c("both", "suit", "dist"),
                         vs = NULL,
                         niche = NULL) {

  surface <- match.arg(surface)

  # ---- 0. Pull components from NicheR_species if provided -----------------
  if (!is.null(vs)) {
    if (!inherits(vs, "NicheR_species")) {
      stop("'vs' must be an object of class 'NicheR_species' created by create_virtual_species().")
    }

    if (is.null(env_bg) && !is.null(vs$call_args) && "env_bg" %in% names(vs$call_args)) {
      env_bg <- vs$call_args$env_bg
    }
    if (is.null(niche) && !is.null(vs$niche)) {
      niche <- vs$niche
    }
    if (is.null(occ_pts) && !is.null(vs$occurrences)) {
      occ_pts <- vs$occurrences
    }
    if (is.null(suitable_env) && !is.null(vs$suitability)) {
      suitable_env <- vs$suitability
    }
  }

  # ---- 1. Palette / color handling ----------------------------------------
  palettes <- list(
    default = list(
      bg           = "#F0F0F0FF",
      suitable_env = "#FED789FF",
      occ_fill     = "#B4BF3AFF",
      occ_stroke   = "#B4BF3AFF",
      dist         = "YlOrRd"   # brewer palette name
    ),
    palette2 = list(
      bg           = "#E0ECF4FF",
      suitable_env = "#9ECAE1FF",
      occ_fill     = "#08519CFF",
      occ_stroke   = "#08306BFF",
      dist         = "YlGnBu"
    ),
    palette3 = list(
      bg           = "#F0F0F0FF",
      suitable_env = "#BDBDBDFF",
      occ_fill     = "#252525FF",
      occ_stroke   = "#000000FF",
      dist         = "OrRd"
    )
  )

  if (!palette %in% names(palettes)) {
    stop("Unknown palette '", palette, "'. Use one of: ",
         paste(names(palettes), collapse = ", "), ".")
  }

  base_colors <- palettes[[palette]]

  if (is.null(colors)) {
    default_colors <- base_colors
  } else {
    if (is.null(names(colors))) {
      names(colors) <- names(base_colors)[seq_along(colors)]
      if (!is.list(colors)) colors <- as.list(colors)
      message(
        "No names detected in 'colors'. Using the provided order and filling ",
        "missing entries with defaults.\n",
        "Named options are: ",
        paste(names(base_colors), collapse = ", "), "."
      )
    }
    if (!is.list(colors)) colors <- as.list(colors)
    default_colors <- utils::modifyList(base_colors, colors)
  }

  # ---- 2. Occurrence points validation ------------------------------------

  if (!is.null(occ_pts) && !all(c("x", "y") %in% names(occ_pts))) {
    stop("`occ_pts` must have columns named 'x' and 'y'.")
  }

  # ---- 3. Standardize suitable_env into a data.frame ----------------------

  suitable_df <- NULL

  if (!is.null(suitable_env)) {

    # Case 1: plain data.frame / matrix
    if (is.data.frame(suitable_env) || is.matrix(suitable_env)) {

      suitable_df <- as.data.frame(suitable_env)

      # Case 2: "suitable_env" style object
    } else if (is.list(suitable_env) &&
               any(c("suitable_env_df", "suitable_env_sp") %in% names(suitable_env))) {

      if ("suitable_env_df" %in% names(suitable_env) &&
          is.data.frame(suitable_env$suitable_env_df)) {

        suitable_df <- suitable_env$suitable_env_df

      } else if ("suitable_env_sp" %in% names(suitable_env)) {

        sp <- suitable_env$suitable_env_sp
        if (inherits(sp, "Raster")) {
          sp <- terra::rast(sp)
        }

        if (inherits(sp, "SpatRaster")) {
          rast_stack <- sp
        } else if (is.list(sp) &&
                   length(sp) > 0 &&
                   all(vapply(sp, inherits, logical(1), "SpatRaster"))) {

          rast_stack <- do.call(c, sp)
        } else {
          stop(
            "'suitable_env$suitable_env_sp' must be a SpatRaster or a list of SpatRasters ",
            "when no 'suitable_env_df' is present."
          )
        }

        suitable_df <- as.data.frame.nicheR(rast_stack)
      } else {
        stop(
          "`suitable_env` list must contain either 'suitable_env_df' or 'suitable_env_sp'."
        )
      }

      # Case 3: Raster / SpatRaster / list of SpatRasters directly
    } else {

      if (inherits(suitable_env, "Raster")) {
        suitable_env <- terra::rast(suitable_env)
      }

      if (inherits(suitable_env, "SpatRaster")) {

        rast_stack <- suitable_env

      } else if (is.list(suitable_env) &&
                 length(suitable_env) > 0 &&
                 all(vapply(suitable_env, inherits, logical(1), "SpatRaster"))) {

        rast_stack <- do.call(c, suitable_env)

      } else {
        stop(
          "`suitable_env` must be either:\n",
          "  * a data.frame with columns 'x', 'y' (and optionally 'dist_sq'), or\n",
          "  * a 'suitable_env' object with 'suitable_env_df' or 'suitable_env_sp', or\n",
          "  * a SpatRaster / list of SpatRasters.\n"
        )
      }

      suitable_df <- as.data.frame.nicheR(rast_stack)
    }

    # sanity check for coordinates
    if (!all(c("x", "y") %in% names(suitable_df))) {
      stop("`suitable_env` (after coercion) must contain 'x' and 'y' columns.")
    }

  } else {
    # no suitable_env provided; fall back to internal computation if possible
    if (!is.null(niche) && !is.null(env_bg)) {
      suitable_df <- get_suitable_env(
        niche     = niche,
        env_bg    = env_bg,
        out.suit  = "data.frame",
        distances = TRUE,
        verbose   = FALSE
      )
      message(
        "Suitability was not supplied; computed internally via get_suitable_env().\n",
        "For repeated plots, consider calling get_suitable_env() once and passing ",
        "the result to `suitable_env`."
      )
    } else {
      suitable_df <- NULL
    }
  }

  has_suitable <- !is.null(suitable_df)
  has_distance <- has_suitable && ("dist_sq" %in% names(suitable_df))

  if (surface == "dist" && !has_distance) {
    stop(
      "surface = 'dist' requested, but 'suitable_env' does not contain a 'dist_sq' column.\n",
      "Run get_suitable_env(..., distances = TRUE) or provide a data.frame with 'dist_sq'."
    )
  }

  # ---- 4. Legend configuration --------------------------------------------

  opts <- list(
    background_point = TRUE,
    suitable_point   = has_suitable,
    occurrence_point = !is.null(occ_pts)
  )

  legend_items <- data.frame(
    id    = c("background_point","suitable_point","occurrence_point"),
    type  = c("point","point","point"),
    label = c("Background environments","Suitable environments","Occurrence"),
    color = c(default_colors[["bg"]],
      default_colors[["suitable_env"]],
      default_colors[["occ_fill"]]
    ),
    stringsAsFactors = FALSE
  )

  active <- logical(nrow(legend_items))
  for (i in seq_len(nrow(legend_items))) {
    active[i] <- isTRUE(opts[[ legend_items$id[i] ]])
  }
  legend_items <- legend_items[active, , drop = FALSE]

  if (nrow(legend_items) == 0) {
    legend_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  } else {
    top_y   <- 2
    spacing <- 0.25
    x_point <- 0.00

    legend_items <- legend_items %>%
      dplyr::mutate(
        row     = dplyr::row_number(),
        y       = top_y - (row - 1) * spacing,
        x_point = x_point,
        x_text  = 0.01
      )

    legend_base <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 2.5), clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none")

    legend_plot <- legend_base +
      ggplot2::geom_point(
        data  = dplyr::filter(legend_items, type == "point"),
        ggplot2::aes(x = x_point, y = y, colour = color),
        size  = 2
      ) +
      ggplot2::geom_text(
        data  = legend_items,
        ggplot2::aes(x = x_text, y = y, label = label),
        hjust = 0
      ) +
      ggplot2::scale_colour_identity()
  }

  # ---- 5. Compute map extent: env_bg -> suitable_env -> occ_pts -----------

  x_lim <- y_lim <- NULL

  # 5A. Try env_bg first
  if (!is.null(env_bg)) {
    if (inherits(env_bg, "Raster")) {
      env_bg_r <- terra::rast(env_bg)
    } else if (inherits(env_bg, "SpatRaster")) {
      env_bg_r <- env_bg
    } else {
      env_bg_r <- NULL
    }

    if (!is.null(env_bg_r)) {
      ex <- terra::ext(env_bg_r)

      if (all(is.finite(as.vector(ex)))) {
        x_lim <- c(ex[1], ex[2])
        y_lim <- c(ex[3], ex[4])
      }
    } else if (is.data.frame(env_bg) || is.matrix(env_bg)) {
      env_bg_df <- as.data.frame(env_bg)
      if (all(c("x", "y") %in% names(env_bg_df))) {
        coords <- env_bg_df[, c("x", "y")]
        coords <- coords[stats::complete.cases(coords), , drop = FALSE]
        if (nrow(coords) > 0) {
          x_rng <- range(coords$x, na.rm = TRUE)
          y_rng <- range(coords$y, na.rm = TRUE)
          dx    <- x_rng[2] - x_rng[1]
          dy    <- y_rng[2] - y_rng[1]
          if (!is.finite(dx) || dx <= 0) dx <- 1
          if (!is.finite(dy) || dy <= 0) dy <- 1
          pad_x <- 0.05 * dx
          pad_y <- 0.05 * dy
          x_lim <- c(x_rng[1] - pad_x, x_rng[2] + pad_x)
          y_lim <- c(y_rng[1] - pad_y, y_rng[2] + pad_y)
        }
      }
    }
  }

  # 5B. If env_bg did not give an extent, use suitable_env coords
  if (is.null(x_lim) || is.null(y_lim)) {
    coord_sources <- list()
    if (has_suitable) {
      coord_sources[[length(coord_sources) + 1]] <- suitable_df[, c("x", "y")]
    }
    if (!is.null(occ_pts) && all(c("x", "y") %in% names(occ_pts))) {
      coord_sources[[length(coord_sources) + 1]] <- occ_pts[, c("x", "y")]
    }

    if (length(coord_sources) > 0) {
      coords <- do.call(rbind, coord_sources)
      coords <- coords[stats::complete.cases(coords), , drop = FALSE]

      if (nrow(coords) > 0) {
        x_rng <- range(coords$x, na.rm = TRUE)
        y_rng <- range(coords$y, na.rm = TRUE)
        dx    <- x_rng[2] - x_rng[1]
        dy    <- y_rng[2] - y_rng[1]
        if (!is.finite(dx) || dx <= 0) dx <- 1
        if (!is.finite(dy) || dy <= 0) dy <- 1
        pad_x <- 0.1 * dx
        pad_y <- 0.1 * dy
        x_lim <- c(x_rng[1] - pad_x, x_rng[2] + pad_x)
        y_lim <- c(y_rng[1] - pad_y, y_rng[2] + pad_y)
      }
    }
  }

  # ---- 6. Basemap ----------------------------------------------------------

  world <- ggplot2::map_data("world")

  base_map <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = default_colors[["bg"]]
    ) +
    ggplot2::xlab("Longitude") +
    ggplot2::ylab("Latitude") +
    ggplot2::theme_bw()

  if (!is.null(x_lim) && !is.null(y_lim)) {
    base_map <- base_map +
      ggplot2::coord_quickmap(xlim = x_lim, ylim = y_lim, expand = FALSE)
  } else {
    base_map <- base_map +
      ggplot2::coord_quickmap()
  }

  # ---- 7. Occurrence points (convert to sf) --------------------------------

  occ_pts_sp <- NULL
  if (!is.null(occ_pts)) {
    occ_pts_sp <- sf::st_as_sf(
      occ_pts[, c("x", "y")],
      coords = c("x", "y"),
      crs    = 4326
    )
  }

  add_occ <- function(p) {
    if (!is.null(occ_pts_sp)) {
      p <- p +
        ggplot2::geom_sf(
          data  = occ_pts_sp,
          ggplot2::aes(geometry = geometry),
          color = default_colors[["occ_stroke"]],
          fill  = default_colors[["occ_fill"]],
          pch   = 21,
          size  = 0.75
        )
    }
    p
  }

  # ---- 8. Build one or two map panels -------------------------------------

  map_list <- list()

  # Panel 1: binary suitable (presence-only tiles)
  if (has_suitable && surface %in% c("both", "suit")) {
    p_suit <- base_map +
      ggplot2::geom_tile(
        data = suitable_df,
        ggplot2::aes(x = x, y = y),
        fill = default_colors[["suitable_env"]]
      ) +
      ggplot2::ggtitle("Suitable Environments")

    p_suit <- add_occ(p_suit)
    map_list[["suitable"]] <- p_suit
  }

  # Panel 2: distance surface
  if (has_distance && surface %in% c("both", "dist")) {
    p_dist <- base_map +
      ggplot2::geom_tile(
        data = suitable_df,
        ggplot2::aes(x = x, y = y, fill = dist_sq)
      ) +
      ggplot2::scale_fill_distiller(
        name    = "Distance to centroid",
        palette = default_colors[["dist"]]
      ) +
      ggplot2::ggtitle("Distance surface") +
      ggplot2::theme(legend.position = "bottom")

    p_dist <- add_occ(p_dist)
    map_list[["dist_sq"]] <- p_dist
  }

  # If no suitability at all, just show basemap + occurrences
  if (length(map_list) == 0L) {
    p_base <- add_occ(base_map)
    return(ggpubr::ggarrange(p_base, legend_plot, widths = c(0.7, 0.3)))
  }

  # ---- 9. Arrange panels + legend -----------------------------------------

  if (length(map_list) == 1L) {
    main_panel <- map_list[[1]]

    return(ggpubr::ggarrange(main_panel,
                             legend_plot,
                             ncol    = 1,
                             heights = c(0.7, 0.3)))

  } else {

    main_panel <- ggpubr::ggarrange(
      plotlist = map_list,
      nrow     = length(map_list),
      labels   = NULL,
      heights  = c(0.45, 0.55)
    )

    return(ggpubr::ggarrange(main_panel,
                             legend_plot,
                             ncol    = 1,
                             heights = c(0.8, 0.2)))
  }

}
