#' Plot Geographic Space with Optional Suitability/Distance and Occurrences
#'
#' Draws a world basemap and overlays (a) suitable environments as tiles
#' derived from a niche ellipsoid, (b) a continuous distance-to-centroid
#' surface, and/or (c) occurrence points. A compact custom legend is built
#' to match the layers shown. If both `show.suitable` and `show.distance`
#' are `TRUE`, the function disables `show.distance` and proceeds with
#' `show.suitable = TRUE`.
#'
#' @param env_bg A background environment object used by
#'   [get_suitable_env()]. Typically a `terra::SpatRaster` (or `raster::Raster*`)
#'   providing the environmental layers and extent. If `vs` is supplied and
#'   `env_bg` is `NULL`, the function will try to use `vs$call_args$env_bg`.
#' @param n_bg Integer (currently unused; reserved for future downsampling of
#'   the background grid).
#' @param niche Optional object of class `"ellipsoid"` (from [build_ellps()]).
#'   Required when `show.suitable = TRUE` or `show.distance = TRUE`. If `vs`
#'   is supplied and `niche` is `NULL`, the function will try to use `vs$niche`.
#' @param show.suitable Logical. If `TRUE`, fills suitable environments
#'   (inside the ellipsoid) as tiles.
#' @param show.distance Logical. If `TRUE`, fills tiles by squared distance to
#'   the ellipsoid centroid (`dist_sq`). If both flags are `TRUE`, distance is
#'   turned off and suitability is shown (with a message).
#' @param occ_pts Optional `data.frame` of occurrences with columns `x`, `y`
#'   (assumed longitude/latitude, WGS84). Plotted as points if supplied. If
#'   `vs` is supplied and `occ_pts` is `NULL`, the function will try to use
#'   `vs$occurrences`.
#' @param show.occ.density Logical (currently unused; reserved for future
#'   density panels).
#' @param colors Optional named list to override aesthetics. Recognized names:
#'   `bg` (basemap fill), `suitable_env` (suitable tiles), `occ_fill`,
#'   `occ_stroke`, and `dist` (RColorBrewer palette name for distance tiles).
#'   Example: `list(suitable_env = "#66CC99", dist = "YlOrRd")`.
#' @param palette Character palette key. One of `"default"`, `"palette2"`,
#'   or `"palette3"`.
#' @param vs Optional object of class `"NicheR_species"` created by
#'   [create_virtual_species()]. When provided, any of `env_bg`, `niche`, and
#'   `occ_pts` that are `NULL` will be filled from this object if possible.
#'
#' @return A `ggpubr` object produced by [ggpubr::ggarrange()], containing the
#'   map panel and a matching legend panel.
#'
#' @details
#' The suitability/distance surface is obtained via
#' `get_suitable_env(niche, env_bg, out = "data.frame", distances = TRUE)`,
#' which returns columns `x`, `y`, and `dist_sq`. When `show.distance = TRUE`,
#' a continuous legend is provided by `scale_fill_distiller()`. The custom
#' legend on the right reflects only the discrete layers actually drawn
#' (background, suitable, occurrences).
#'
#' If a `NicheR_species` object is provided via `vs`, the function:
#' \itemize{
#'   \item uses `vs$call_args$env_bg` as `env_bg` when `env_bg` is `NULL`,
#'   \item uses `vs$niche` as `niche` when `niche` is `NULL`,
#'   \item uses `vs$occurrences` as `occ_pts` when `occ_pts` is `NULL`.
#' }
#' Explicit arguments always take precedence over values inferred from `vs`.
#'
#' @section Notes:
#' - `niche` must be provided (directly or via `vs`) if either `show.suitable`
#'   or `show.distance` is `TRUE`.
#' - If both `show.suitable` and `show.distance` are `TRUE`, the function sets
#'   `show.distance <- FALSE` and informs the user.
#' - `occ_pts` must include columns named `x` and `y`.
#'
#' @family plotting functions
#' @seealso [build_ellps()], [get_suitable_env()], [create_virtual_species()]
#'
#' @export
plot_g_space <- function(env_bg,
                         n_bg = 10000,
                         niche = NULL,
                         show.suitable = FALSE,
                         show.distance = FALSE,
                         occ_pts = NULL,
                         show.occ.density = FALSE,
                         colors = NULL,
                         palette = "default",
                         vs = NULL) {

  # ---- 0. Pull components from NicheR_species if provided -----------------
  if (!is.null(vs)) {
    if (!inherits(vs, "NicheR_species")) {
      stop("'vs' must be an object of class 'NicheR_species' created by create_virtual_species().")
    }

    # env_bg: try call_args$env_bg if user didn't pass env_bg explicitly
    if (missing(env_bg) || is.null(env_bg)) {
      if (!is.null(vs$call_args) && "env_bg" %in% names(vs$call_args)) {
        env_bg <- vs$call_args$env_bg
      }
    }

    # niche: use vs$niche if not supplied explicitly
    if (is.null(niche) && !is.null(vs$niche)) {
      niche <- vs$niche
    }

    # occ_pts: use vs$occurrences if not supplied explicitly
    if (is.null(occ_pts) && !is.null(vs$occurrences)) {
      occ_pts <- vs$occurrences
    }
  }

  # ---- 1. Palette / color handling ----------------------------------------
  palettes <- list(
    default = list(
      bg           = "#FED789FF",
      suitable_env = "#B4BF3AFF",
      occ_fill     = "black",
      occ_stroke   = "black",
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
    # Allow unnamed list/vector: assign in base_colors order
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
    default_colors <- modifyList(base_colors, colors)
  }

  # ---- 2. Input validation & auto-fixes -----------------------------------

  if (missing(env_bg) || is.null(env_bg)) {
    stop("`env_bg` must be provided (directly or via `vs`) to set extent / background.")
  }

  # If both TRUE, turn distance off
  if (isTRUE(show.suitable) && isTRUE(show.distance)) {
    show.distance <- FALSE
    message("Both `show.suitable` and `show.distance` were TRUE. ",
            "Setting `show.distance = FALSE` to avoid conflicts.")
  }

  # If either suitable or distance is requested, need a niche object
  if ((isTRUE(show.suitable) || isTRUE(show.distance)) && is.null(niche)) {
    stop("`niche` must not be NULL when `show.suitable` or `show.distance` is TRUE ",
         "(provide it directly or via `vs`).")
  }

  # Occurrence points validation
  if (!is.null(occ_pts)) {
    if (!all(c("x", "y") %in% names(occ_pts))) {
      stop("`occ_pts` must have columns named 'x' and 'y'.")
    }
  }

  # ---- 3. Legend configuration --------------------------------------------

  opts <- list(
    background_point = TRUE,
    suitable_point   = show.suitable,
    occurrence_point = !is.null(occ_pts)
  )

  legend_items <- data.frame(
    id       = c("background_point","suitable_point","occurrence_point"),
    type     = c("point","point","point"),
    label    = c("Background environments","Suitable environments","Occurrence"),
    color    = c(
      default_colors[["bg"]],
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

  # ---- 4. Basemap ----------------------------------------------------------

  world <- ggplot2::map_data("world")

  return_plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = default_colors[["bg"]]
    ) +
    ggplot2::xlab("Longitude") +
    ggplot2::ylab("Latitude") +
    ggplot2::theme_bw()

  # ---- 5. Occurrence points (convert to sf) --------------------------------

  if (!is.null(occ_pts)) {
    occ_pts_sp <- sf::st_as_sf(
      occ_pts[, c("x", "y")],
      coords = c("x", "y"),
      crs    = 4326
    )
  }

  # ---- 6. Suitable / distance surfaces ------------------------------------

  if (isTRUE(show.suitable) || isTRUE(show.distance)) {
    # env_bg can be Raster* or SpatRaster; get_suitable_env handles both
    suitable_g_space <- get_suitable_env(
      niche     = niche,
      env_bg    = env_bg,
      out       = "data.frame",
      distances = TRUE
    )
  }

  if (isTRUE(show.distance)) {
    return_plot <- return_plot +
      ggplot2::geom_tile(
        data = suitable_g_space,
        ggplot2::aes(x = x, y = y, fill = dist_sq)
      ) +
      ggplot2::scale_fill_distiller(
        name    = "Distance to centroid",
        palette = default_colors[["dist"]]
      ) +
      ggplot2::theme(
        legend.position = "bottom"
      )
  }

  if (isTRUE(show.suitable)) {
    return_plot <- return_plot +
      ggplot2::geom_tile(
        data = suitable_g_space,
        ggplot2::aes(x = x, y = y),
        fill = default_colors[["suitable_env"]]
      )
  }

  if (!is.null(occ_pts)) {
    return_plot <- return_plot +
      ggplot2::geom_sf(
        data  = occ_pts_sp,
        ggplot2::aes(geometry = geometry),
        color = default_colors[["occ_stroke"]],
        fill  = default_colors[["occ_fill"]],
        pch   = 21,
        size  = 0.75
      )
  }

  ggpubr::ggarrange(return_plot, legend_plot, widths = c(0.7, 0.3))
}
