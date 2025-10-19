#' @family plotting functions
#' @seealso [build_ellps()], [get_suitable_env()]
#' @import RColorBrewer ggplot2 dplyr
#' @export
plot_g_space <- function(env_bg,
                         n_bg = 10000,
                         niche = NULL,
                         show.suitable = FALSE,
                         show.distance = FALSE,
                         occ_pts = NULL,
                         show.occ.density = FALSE,
                         colors = NULL,
                         palette = "default"){


  # Color / palette handling
  default_colors <- list(
    bg       = "#FED789FF",
    suitable_env    = "#B4BF3AFF",
    occ_fill  = "black",
    occ_stroke= "black",
    dist = "YlOrRd"  # brewer palette name
  )

  if (!is.null(colors)) {
    # Merge user-provided named colors into defaults
    # valid names: background, suitable_fill, occurrence_fill, occurrence_stroke, distance_palette
    default_colors[names(colors)] <- colors
  }

  # Input validation & auto-fixes
  if (missing(env_bg) || is.null(env_bg)) {
    stop("`env_bg` must be provided (Raster*/SpatRaster) to set extent / background.")
  }

  # If both TRUE, turn distance off
  if (isTRUE(show.suitable) && isTRUE(show.distance)) {
    show.distance <- FALSE
    message("Both `show.suitable` and `show.distance` were TRUE. Setting `show.distance = FALSE` to avoid conflicts.")
  }

  # If either suitable or distance is requested, need a niche object
  if ((isTRUE(show.suitable) || isTRUE(show.distance)) && is.null(niche)) {
    stop("`niche` must not be NULL when `show.suitable` or `show.distance` is TRUE.")
  }

  # Occurrence points validation
  if (!is.null(occ_pts)) {
    if (!all(c("x", "y") %in% names(occ_pts))) {
      stop("`occ_pts` must have columns named 'x' and 'y'.")
    }
  }

  # Check options for  legend
  ## Toggle these on/off
  opts <- list(
    background_point = TRUE,
    suitable_point   = show.suitable,
    occurrence_point = !is.null(occ_pts)
  )

  ## 1) Legend spec (base R data.frame)
  legend_items <- data.frame(
    id       = c("background_point","suitable_point","occurrence_point"),
    type     = c("point","point","point"),
    label    = c("Background environments","Suitable evironments","Occurrence"),
    color    = c(default_colors[["bg"]], default_colors[["suitable_env"]], default_colors[["occ_fill"]]),
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
    x_text  <- 0.00
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
      coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 2.5), clip = "off") +
      theme_void() +
      theme(legend.position = "none")

    legend_plot <- legend_base +
      # points
      geom_point(
        data = filter(legend_items, type == "point"),
        aes(x = x_point, y = y, colour = color),
        size = 2
      ) +
      # labels
      geom_text(
        data = legend_items,
        aes(x = 0.01, y = y, label = label),
        hjust = 0
      ) +
      scale_colour_identity() +
      scale_linetype_identity()

    legend_plot
  }

  world <- map_data("world")

  return_plot <- ggplot( ) +
    geom_polygon(data = world,
                 aes(x = long, y = lat, group = group),
                 fill = default_colors[["bg"]]) +
    xlab("Longitude") + ylab("Latitude") +
    theme_bw()



  if(!is.null(occ_pts)){
    # Transform point into spatial points using sf
    occ_pts_sp <- sf::st_as_sf(occ_pts[,c("x", "y")],
                               coords = c("x", "y"))


  }


  if(isTRUE(show.suitable) || isTRUE(show.distance)){
    suitable_g_space <- get_suitable_env(niche = niche,
                                         env_bg = env_bg,
                                         out = "data.frame",
                                         distances = TRUE)

  }


  if(isTRUE(show.distance)){

    return_plot <- return_plot +
      geom_tile(data = suitable_g_space,
                aes(x = x, y = y, fill = dist_sq)) +
      scale_fill_distiller(name = "Distance to centroid",
                           palette = default_colors[["dist"]]) +
      theme(
        legend.position = "bottom"
      )
  }

  if(isTRUE(show.suitable)){

    return_plot <- return_plot +
      geom_tile(data = suitable_g_space,
                aes(x = x, y = y), fill = default_colors[["suitable_env"]])
  }

  if(!is.null(occ_pts)){

    return_plot <- return_plot +
      geom_sf(data = occ_pts_sp, aes(geometry = geometry),
              color = default_colors[["occ_fill"]],
              fill = default_colors[["occ_fill"]], pch = 21, size = 0.75)

  }


  ggpubr::ggarrange(return_plot, legend_plot,
                    widths = c(0.7, 0.3))


}


