#' Plot geographic space with suitability, distance, and occurrences
#'
#' Draws a basemap and overlays any combination of:
#' \itemize{
#'   \item binary suitable environments (tiles),
#'   \item an optional continuous distance-to-centroid surface,
#'   \item occurrence points.
#' }
#'
#' The map extent is primarily inferred from \code{env_bg} (if available),
#' otherwise from the coordinates in \code{suitable_env} and finally from
#' \code{occ_pts}, always with a small buffer where needed.
#'
#' When a \code{NicheR_species} object (\code{vs}) is supplied, any of
#' \code{env_bg}, \code{suitable_env}, \code{occ_pts}, or \code{niche} that are
#' \code{NULL} are filled using the generic getters:
#' \code{nr_get_env()}, \code{nr_get_suitable_all()}, \code{nr_get_occ()}, and
#' \code{nr_get_niche()}.
#'
#' @param env_bg Optional background environment used mainly to define the map
#'   extent. Can be:
#'   \itemize{
#'     \item a \code{terra::SpatRaster} or \code{raster::Raster*} object,
#'     \item a \code{data.frame} or \code{matrix} with columns \code{x}, \code{y},
#'       and optionally additional predictors.
#'   }
#'
#' @param suitable_env Optional suitability object, used when plotting
#'   binary suitability and/or distance-to-centroid surfaces. This can be:
#'   \itemize{
#'     \item a \code{data.frame} or \code{matrix} with columns \code{x}, \code{y}
#'       (and optionally \code{dist_sq}),
#'     \item a \code{"suitable_env"} object returned by \code{\link{get_suitable_env}},
#'     \item a \code{NicheR_species} suitability object inside \code{vs},
#'     \item a \code{terra::SpatRaster} or a list of \code{SpatRaster} objects
#'       containing layers named \code{"suitable"} and/or \code{"dist_sq"}.
#'   }
#'   For memory efficiency, the most efficient input is a data frame of suitable
#'   cells (with \code{x}, \code{y}, and optionally \code{dist_sq}). If only
#'   rasters are available, they are converted to a data frame via
#'   \code{as.data.frame.nicheR()}, which may be slow for large rasters.
#'
#' @param occ_pts Optional data frame of occurrence points with columns
#'   \code{x}, \code{y} (assumed longitude/latitude in WGS84). If omitted and
#'   \code{vs} is provided, occurrences are retrieved via \code{nr_get_occ(vs)}.
#'
#' @param show.occ.density Reserved for future use (currently ignored).
#'
#' @param colors Optional named list to override aesthetics. Recognized names:
#'   \code{bg} (basemap fill), \code{suitable_env} (suitable tiles),
#'   \code{occ_fill}, \code{occ_stroke}, and \code{dist} (RColorBrewer palette
#'   name for the distance tiles).
#'
#' @param palette Character palette key. One of \code{"default"},
#'   \code{"palette2"}, or \code{"palette3"}.
#'
#' @param show.in.plot Character vector specifying what to show in the plot.
#'   This controls which surfaces are drawn and whether occurrences are
#'   overlaid. It accepts both abbreviations and full names:
#'   \itemize{
#'     \item \code{"suit"} or \code{"suitability"}: plot the binary suitable environments.
#'     \item \code{"dist"} or \code{"distance"}: plot the distance-to-centroid surface.
#'     \item \code{"occ"} or \code{"occurrences"}: overlay occurrence points
#'           on all panels (if available).
#'     \item \code{"bg"} or \code{"background"}: only affects messaging; the
#'           basemap is always drawn.
#'     \item \code{"none"}: no suitability/distance surfaces (only basemap and,
#'           if requested via \code{"occ"}, occurrences).
#'   }
#'
#'   \strong{Combinations} are allowed, for example:
#'   \itemize{
#'     \item \code{show.in.plot = c("suit", "occ")}:
#'           one panel with suitable environments and occurrences.
#'     \item \code{show.in.plot = c("suit", "dist", "occ")}:
#'           two panels, one for suitability and one for distance, both with occurrences.
#'     \item \code{show.in.plot = "dist"}:
#'           one panel with the distance surface only.
#'   }
#'
#'   If \code{show.in.plot} is \code{NULL}, the function chooses a default:
#'   \itemize{
#'     \item If a suitability object is available, defaults to
#'           \code{c("suitability", "occurrences")} when occurrences exist,
#'           otherwise just \code{"suitability"}.
#'     \item If no suitability is available but occurrences exist, defaults to
#'           \code{"occurrences"}.
#'     \item Otherwise defaults to \code{"none"}.
#'   }
#'
#'   A brief message is printed indicating the resolved mode and reminding the
#'   user to set \code{show.in.plot} explicitly if a different combination is desired.
#'
#' @param vs Optional object of class \code{"NicheR_species"} created by
#'   \code{\link{create_virtual_species}}. When provided, any of \code{env_bg},
#'   \code{suitable_env}, \code{occ_pts}, or \code{niche} that are \code{NULL}
#'   will be filled from this object via the \code{nr_get_*()} helpers.
#'
#' @param niche Optional object of class \code{"ellipsoid"} (from
#'   \code{\link{build_ellps}}). Currently not required for plotting but kept
#'   for possible future extensions.
#'
#' @return A \code{ggpubr} object produced by \code{\link[ggpubr]{ggarrange}},
#'   containing one or two map panels plus a matching legend panel.
#'
#' @family plotting functions
#' @seealso \code{\link{build_ellps}}, \code{\link{get_suitable_env}},
#'   \code{\link{create_virtual_species}}, \code{\link{nr_get}}
#'
#' @export
plot_g_space <- function(env_bg = NULL,
                         suitable_env = NULL,
                         occ_pts = NULL,
                         show.occ.density = FALSE,
                         colors = NULL,
                         palette = "default",
                         show.in.plot = NULL,
                         vs = NULL,
                         niche = NULL) {


  ## 0. Pull components from NicheR_species (using nr_get_*)

  if (!is.null(vs)) {
    if (!inherits(vs, "NicheR_species")) {
      stop("'vs' must be an object of class 'NicheR_species' created by create_virtual_species().")
    }

    if (is.null(env_bg)) {
      env_from_vs <- nr_get_env(vs)
      if (!is.null(env_from_vs)) {
        env_bg <- env_from_vs
        message("Using env_bg retrieved from 'vs' via nr_get_env().")
      }
    }

    if (is.null(niche)) {
      niche_from_vs <- nr_get_niche(vs)
      if (!is.null(niche_from_vs)) {
        niche <- niche_from_vs
      }
    }

    if (is.null(occ_pts)) {
      occ_from_vs <- nr_get_occ(vs)
      if (!is.null(occ_from_vs)) {
        occ_pts <- occ_from_vs
      }
    }

    if (is.null(suitable_env)) {
      suit_from_vs <- nr_get_suitable_all(vs)
      if (!is.null(suit_from_vs)) {
        suitable_env <- suit_from_vs
      }
    }
  }


  ## 1. Palette / color handling

  palettes <- list(
    default = list(
      bg           = "#F0F0F0FF",
      suitable_env = "#FED789FF",
      occ_fill     = "#B4BF3AFF",
      occ_stroke   = "#B4BF3AFF",
      dist         = "YlOrRd"
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


  ## 2. Occurrence points validation

  if (!is.null(occ_pts) && !all(c("x", "y") %in% names(occ_pts))) {
    stop("`occ_pts` must have columns named 'x' and 'y'.")
  }


  ## 3. Resolve show.in.plot modes (abbreviations + defaults)

  normalize_tokens <- function(x) {
    if (is.null(x)) return(NULL)
    if (!is.character(x)) {
      stop("`show.in.plot` must be NULL or a character vector.")
    }

    x <- tolower(x)

    map <- c(
      suit        = "suitability",
      suitability = "suitability",
      dist        = "distance",
      distance    = "distance",
      occ         = "occurrences",
      occs        = "occurrences",
      occurrences = "occurrences",
      bg          = "background",
      background  = "background",
      none        = "none"
    )

    unknown <- setdiff(x, names(map))

    if (length(unknown)) {
      stop(
        "Unknown entries in `show.in.plot`: ",
        paste(unique(unknown), collapse = ", "),
        ". Allowed values (abbrev or full) include: ",
        "'suit', 'suitability', 'dist', 'distance', ",
        "'occ', 'occurrences', 'bg', 'background', 'none'."
      )
    }
    unname(map[x])
  }

  tokens <- normalize_tokens(show.in.plot)

  # choose default if none specified
  if (is.null(tokens)) {
    has_suit_obj  <- !is.null(suitable_env)
    has_occ_obj   <- !is.null(occ_pts)

    if (has_suit_obj) {
      if (has_occ_obj) {
        tokens <- c("suitability", "occurrences")
      } else {
        tokens <- "suitability"
      }
    } else if (has_occ_obj) {
      tokens <- "occurrences"
    } else {
      tokens <- "none"
    }

    message(
      "Plot mode not specified; using default show.in.plot = c('",
      paste(tokens, collapse = "', '"), "').\n",
      "To change what is displayed, set `show.in.plot` explicitly, e.g. ",
      "show.in.plot = c('suit', 'dist', 'occ')."
    )
  } else {
    message(
      "Plot mode resolved as show.in.plot = c('",
      paste(tokens, collapse = "', '"), "')."
    )
  }

  # interpret tokens
  show_suit <- "suitability" %in% tokens
  show_dist <- "distance"    %in% tokens
  show_occ  <- "occurrences" %in% tokens

  # 'none' just means "no surfaces"; background always drawn
  surfaces_requested <- show_suit || show_dist


  ## 4. Standardize suitable_env into a data.frame (df prioritized)

  suitable_df  <- NULL
  has_suitable <- FALSE
  has_distance <- FALSE

  if (surfaces_requested) {

    if (is.null(suitable_env)) {
      stop(
        "Suitability/distance surfaces were requested via `show.in.plot`,\n",
        "but `suitable_env` could not be resolved. Either:\n",
        "  * pass suitable_env (e.g. output from get_suitable_env), or\n",
        "  * ensure `vs` contains a suitability object accessible via nr_get().\n",
        "If you only want background + occurrences, use show.in.plot = c('occ') or 'none'."
      )
    }

    # 4A. If the user passed a data.frame/matrix, use it directly (most efficient)
    if (is.data.frame(suitable_env) || is.matrix(suitable_env)) {

      suitable_df <- as.data.frame(suitable_env)

    } else {

      # 4B. Try to get a df via nr_get() first
      df_try <- suppressWarnings(nr_get_suitable_df(suitable_env))

      if (!is.null(df_try)) {
        suitable_df <- as.data.frame(df_try)
      } else {

        # 4C. Fall back to rasters: suitable + dist_sq
        r_suit <- NULL
        r_dist <- NULL

        if (show_suit) {
          r_suit <- suppressWarnings(nr_get_suitable(suitable_env))
        }
        if (show_dist) {
          r_dist <- suppressWarnings(nr_get_dist_sq(suitable_env))
        }

        if (inherits(r_suit, "Raster")) r_suit <- terra::rast(r_suit)
        if (inherits(r_dist, "Raster")) r_dist <- terra::rast(r_dist)

        if (show_suit && !is.null(r_suit) && !show_dist) {

          names(r_suit) <- "suitable"
          message(
            "Converting suitability raster to a data.frame via as.data.frame.nicheR().\n",
            "For large rasters this may take some time. For better performance, ",
            "consider passing a precomputed data.frame of suitable cells."
          )
          suitable_df <- as.data.frame.nicheR(r_suit)

        } else if (show_dist && !is.null(r_dist) && !show_suit) {

          names(r_dist) <- "dist_sq"
          message(
            "Converting distance raster to a data.frame via as.data.frame.nicheR().\n",
            "For large rasters this may take some time. For better performance, ",
            "consider passing a data.frame with x, y, and dist_sq."
          )
          suitable_df <- as.data.frame.nicheR(r_dist)

        } else if (show_suit && show_dist) {

          if (is.null(r_suit) || is.null(r_dist)) {
            stop(
              "show.in.plot includes both 'suitability' and 'distance', but both rasters\n",
              "could not be found. Ensure your object contains layers for 'suitable' and 'dist_sq',\n",
              "or provide a data.frame with columns x, y, and dist_sq."
            )
          }

          names(r_suit) <- "suitable"
          names(r_dist) <- "dist_sq"
          ras_stack <- c(r_suit, r_dist)
          message(
            "Converting suitability + distance rasters to a data.frame via as.data.frame.nicheR().\n",
            "For large rasters this may take some time. For better performance, ",
            "consider passing a data.frame with x, y, suitable, and dist_sq."
          )
          suitable_df <- as.data.frame.nicheR(ras_stack)
        }
      }
    }

    if (!is.null(suitable_df)) {
      if (!all(c("x", "y") %in% names(suitable_df))) {
        stop("`suitable_env` (after coercion) must contain 'x' and 'y' columns.")
      }
      has_suitable <- show_suit && nrow(suitable_df) > 0
      has_distance <- show_dist && ("dist_sq" %in% names(suitable_df))
      if (show_dist && !has_distance) {
        stop(
          "show.in.plot includes 'distance', but no 'dist_sq' column was found in suitable_env.\n",
          "Run get_suitable_env(..., distances = TRUE) or provide an object with 'dist_sq'."
        )
      }
    }
  }


  ## 5. Legend configuration (with distance note toggle)

  opts <- list(
    background_point = TRUE,
    suitable_point   = has_suitable,
    occurrence_point = !is.null(occ_pts) && show_occ,
    dist_note        = has_distance
  )

  legend_items <- data.frame(
    id    = c("background_point", "suitable_point", "occurrence_point", "dist_note"),
    type  = c("point",           "point",          "point",           "text"),
    label = c("Background environments",
              "Suitable environments",
              "Occurrences",
              "Distance surface: lighter = closer to centroid"),
    color = c(default_colors[["bg"]],
              default_colors[["suitable_env"]],
              default_colors[["occ_fill"]],
              NA_character_),
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


  ## 6. Compute map extent: env_bg -> suitable_env -> occ_pts

  x_lim <- y_lim <- NULL
  extent_source <- NULL

  if (!is.null(env_bg)) {
    env_bg_r <- NULL

    if (inherits(env_bg, "Raster")) {
      env_bg_r <- terra::rast(env_bg)
    } else if (inherits(env_bg, "SpatRaster")) {
      env_bg_r <- env_bg
    }

    if (!is.null(env_bg_r)) {
      ex <- terra::ext(env_bg_r)
      if (all(is.finite(as.vector(ex)))) {
        x_lim <- c(ex[1], ex[2])
        y_lim <- c(ex[3], ex[4])
        extent_source <- "env_bg"
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
          extent_source <- "env_bg_xy"
        }
      }
    }
  }

  if ((is.null(x_lim) || is.null(y_lim)) &&
      !is.null(suitable_df) && nrow(suitable_df) > 0) {

    coords <- suitable_df[, c("x", "y")]
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
      extent_source <- "suitable"
    }
  }

  if ((is.null(x_lim) || is.null(y_lim)) &&
      !is.null(occ_pts) && all(c("x", "y") %in% names(occ_pts))) {

    coords <- occ_pts[, c("x", "y")]
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
      extent_source <- "occ"
    }
  }

  if (!is.null(extent_source)) {
    msg_src <- switch(
      extent_source,
      "env_bg"    = "Extent derived from env_bg raster.",
      "env_bg_xy" = "Extent derived from env_bg x/y coordinates.",
      "suitable"  = "Extent derived from suitable_env coordinates.",
      "occ"       = "Extent derived from occurrence coordinates."
    )
    message(msg_src)
  } else {
    message("Could not infer extent from env_bg, suitable_env, or occ_pts; using global extent.")
  }


  ## 7. Basemap

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


  ## 8. Occurrence points (sf), with global toggle

  occ_pts_sp <- NULL
  if (!is.null(occ_pts)) {
    occ_pts_sp <- sf::st_as_sf(
      occ_pts[, c("x", "y")],
      coords = c("x", "y"),
      crs    = 4326
    )
  }

  add_occ <- function(p) {
    if (show_occ && !is.null(occ_pts_sp)) {
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


  ## 9. Build one or two map panels

  map_list <- list()

  if (has_suitable && show_suit) {
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

  if (has_distance && show_dist) {
    p_dist <- base_map +
      ggplot2::geom_tile(
        data = suitable_df,
        ggplot2::aes(x = x, y = y, fill = dist_sq)
      ) +
      ggplot2::scale_fill_distiller(
        name    = "Distance to centroid",
        palette = default_colors[["dist"]],
        guide   = "none"
        ) +
      ggplot2::ggtitle("Distance surface")

    p_dist <- add_occ(p_dist)

    map_list[["dist_sq"]] <- p_dist
  }

  # If no surfaces requested or available -> basemap + (optional) occs
  if (!surfaces_requested || length(map_list) == 0L) {
    p_base <- add_occ(base_map)
    return(ggpubr::ggarrange(p_base, legend_plot, widths = c(0.7, 0.3)))
  }


  ## 10. Arrange panels + legend

  if (length(map_list) == 1L) {
    main_panel <- map_list[[1]]

    ggpubr::ggarrange(
      main_panel,
      legend_plot,
      ncol    = 1,
      heights = c(0.7, 0.3)
    )

  } else {

    main_panel <- ggpubr::ggarrange(
      plotlist = map_list,
      nrow     = length(map_list),
      labels   = NULL
    )

    ggpubr::ggarrange(
      main_panel,
      legend_plot,
      ncol    = 1,
      heights = c(0.8, 0.2)
    )
  }
}
