# Title: Generate tab plot logic

# Description: E-space, G-space, and combined plots for the generate tab.
# Mirrors predict_tab_plot.R, with generated occurrences drawn over the
# surface they were sampled from.

# Date last updated: 08/06/2026


# Functions ---------------------------------------------------------------

generate_draw_espace_panel <- function(v1, v2, s, set = NULL, title = NULL){

  lims <- compute_lims(v1, v2, s)
  bg <- session_data$bg_df
  ttl <- if(!is.null(title)) title else paste(v1, "vs.", v2)

  if(!is.null(bg)){
    plot(bg[[v1]], bg[[v2]],
         col = s$bg_col,
         pch = s$pch_val,
         cex = s$cex_val,
         xlim = lims$xlim,
         ylim = lims$ylim,
         asp = if(!is.null(s$asp_espace) && s$asp_espace == "fixed") lims$asp else NA,
         xlab = v1,
         ylab = v2,
         main = ttl)
  } else {
    plot(NA, NA,
         xlim = lims$xlim,
         ylim = lims$ylim,
         asp = if(!is.null(s$asp_espace) && s$asp_espace == "fixed") lims$asp else NA,
         xlab = v1,
         ylab = v2,
         main = ttl)
  }

  # Occurrences drawn first so they sit under the background and ellipsoid
  draw_sets <- if(!is.null(set)) set else s$occ_set

  if(!is.null(draw_sets) && length(draw_sets) > 0){

    cols <- generate_set_colors(draw_sets, s)

    for(nm in draw_sets){

      occ <- s$occ_espace[[nm]]
      if(is.null(occ)) next
      if(!all(c(v1, v2) %in% names(occ))) next

      points(occ[[v1]], occ[[v2]],
             pch = s$occ_pch,
             col = cols[[nm]],
             cex = s$occ_cex)
    }
  }

  if(s$show_ell && s$has_ell){
    idx <- match(c(v1, v2), s$ell$var_names)
    add_ellipsoid(s$ell, dim = idx,
                  col_ell = s$ell_col,
                  lwd = s$ell_lwd,
                  lty = s$ell_lty)
  }

  if(s$show_centroid && !is.null(s$ell)){
    idx <- match(c(v1, v2), s$ell$var_names)
    points(s$ell$centroid[idx[1]], s$ell$centroid[idx[2]],
           pch = s$centroid_pch, col = s$centroid_col, cex = s$centroid_cex)
  }

}

generate_draw_espace_pairs <- function(vars, s){

  pairs <- t(combn(seq_along(vars), 2))
  n_pairs <- nrow(pairs)
  n_cols <- ceiling(sqrt(n_pairs))
  n_rows <- ceiling(n_pairs / n_cols)

  draw_sets <- s$occ_set
  show_key <- length(draw_sets) > 1

  if(show_key){
    cells <- c(seq_len(n_pairs), rep(0L, n_rows * n_cols - n_pairs))
    m <- matrix(cells, nrow = n_rows, ncol = n_cols, byrow = TRUE)
    m <- rbind(m, rep(n_pairs + 1L, n_cols))
    layout(m, heights = c(rep(1, n_rows), 0.28))
  } else {
    par(mfrow = c(n_rows, n_cols))
  }

  par(mar = c(4, 4, 2, 1))

  for(i in seq_len(n_pairs)){
    generate_draw_espace_panel(vars[pairs[i, 1]], vars[pairs[i, 2]], s)
  }

  if(show_key){

    cols <- generate_set_colors(draw_sets, s)
    labs <- vapply(draw_sets, generate_set_label, character(1), s = s)

    par(mar = c(0, 1, 0, 1))
    plot.new()
    legend("center",
           legend = labs,
           col = unname(cols),
           pch = s$occ_pch,
           pt.cex = 1.2,
           cex = 0.75,
           horiz = length(draw_sets) <= 3,
           bty = "n")
  }

  layout(1)
}

# One panel per occurrence set, same variable pair, so several sampling
# runs can be compared side by side
generate_draw_espace_sets <- function(v1, v2, s, sets){

  n <- length(sets)
  n_cols <- if(n <= 1) 1L else 2L
  n_rows <- ceiling(n / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 2, 1))

  for(nm in sets){
    n_pts <- if(!is.null(s$occ_espace[[nm]])) nrow(s$occ_espace[[nm]]) else 0L
    generate_draw_espace_panel(v1, v2, s, set = nm,
                               title = paste0(generate_set_label(nm, s), " (n = ", n_pts, ")"))
  }

  if(n > 1 && n %% 2 != 0) plot.new()
}

generate_draw_gspace_panel <- function(rast, s, title = NULL, set = NULL,
                                       col = NULL, legend = TRUE){

  ttl <- if(!is.null(title)) title else names(rast)[1]
  cols <- if(!is.null(col)) col else s$palette

  terra::plot(rast,
              col = cols,
              legend = legend,
              main = ttl,
              colNA = s$map_bg_col,
              axes = TRUE,
              xlab = "Longitude",
              ylab = "Latitude")

  nm <- if(!is.null(set)) set[1] else s$occ_set[1]

  if(!is.null(nm)){
    occ <- session_data$ellipsoid_occurrence_list[[s$ell$ell_id]][[nm]]
    if(!is.null(occ) && nrow(occ) > 0){
      points(occ$x, occ$y,
             pch = s$occ_pch,
             col = s$occ_col,
             cex = s$occ_cex)
    }
  }
}

# Distinct colours when several sets share a panel. Falls back to the
# single configured occurrence colour when only one is drawn.
generate_set_colors <- function(sets, s){

  if(length(sets) <= 1) return(setNames(s$occ_col, sets))

  setNames(pred_palette("Dark 3", n = length(sets)), sets)
}

# Reactives ---------------------------------------------------------------

generate_plot_vars <- reactive({

  ell <- session_data$current_ellipsoid

  if(!is.null(ell) && !is.null(ell$var_names)) return(ell$var_names)

  session_data$vars
})

# Occurrence set names for the current ellipsoid, one per prediction layer
# that was sampled from
generate_occ_sets <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(character(0))

  occ <- session_data$ellipsoid_occurrence_list[[ell$ell_id]]
  if(is.null(occ)) return(character(0))

  names(occ)
})

# Set key to display label. occ_set_labels() returns the inverse mapping,
# which suits a selectInput but not a lookup by key.
generate_occ_label_map <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(character(0))

  occ <- session_data$ellipsoid_occurrence_list[[ell$ell_id]]
  if(is.null(occ) || length(occ) == 0) return(character(0))

  labs <- occ_set_labels(occ)
  setNames(names(labs), unname(labs))
})

generate_set_label <- function(nm, s){
  lab <- s$occ_labels[nm]
  if(is.na(lab)) nm else unname(lab)
}

# Environmental values at the occurrence coordinates. generate_occ_for_ell()
# returns x and y only, so the values have to be extracted before the points
# can be drawn in environmental space. Cached per ellipsoid.
generate_occ_espace <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(NULL)

  occ <- session_data$ellipsoid_occurrence_list[[ell$ell_id]]
  if(is.null(occ) || length(occ) == 0) return(NULL)

  rast <- session_data$bg_raster
  vars <- ell$var_names

  lapply(occ, function(df){

    if(is.null(df) || nrow(df) == 0) return(NULL)

    # Sets with no coordinates are environmental values already. Detected by
    # structure rather than by the mode attribute, so this also covers the
    # non-spatial CSV case and any set written before mode was recorded.
    if(!all(c("x", "y") %in% names(df))) return(df)

    if(is.null(rast)) return(NULL)

    vals <- tryCatch(
      terra::extract(terra::subset(rast, vars), df[, c("x", "y")], ID = FALSE),
      error = function(e) NULL
    )

    if(is.null(vals)) return(NULL)

    cbind(df[, c("x", "y")], vals)
  })
})
# The raster an occurrence set was sampled from. Checks the biased list
# first, since a biased layer name never appears in the unbiased stack.
generate_source_raster <- function(ell_id, set){

  occ <- session_data$ellipsoid_occurrence_list[[ell_id]][[set]]
  layer <- occ_meta(occ, "layer", set)

  biased <- session_data$ellipsoid_prediction_list_biased[[ell_id]]
  if(inherits(biased, "SpatRaster") && layer %in% names(biased)){
    return(biased[[layer]])
  }

  pred <- session_data$ellipsoid_prediction_list[[ell_id]]
  if(inherits(pred, "SpatRaster") && layer %in% names(pred)){
    return(pred[[layer]])
  }

  NULL
}

generate_collect_plot_settings <- function(){

  has_ell <- !is.null(session_data$current_ellipsoid)

  list(
    has_ell = has_ell,
    ell = if(has_ell) session_data$current_ellipsoid else NULL,

    # No range lines on this tab
    show_lines = list(active = FALSE, ranges = NULL),

    show_ell = has_ell && get_input("generate_show_ellipsoid", TRUE),
    show_centroid = has_ell && get_input("generate_show_centroid", TRUE),

    occ_espace = generate_occ_espace(),
    occ_labels = generate_occ_label_map(),
    occ_set = get_input("generate_occ_set", NULL),

    occ_pch = as.numeric(get_input("generate_occ_pch", "16")),
    occ_cex = get_input("generate_occ_cex", 0.7),
    occ_col = get_input("generate_occ_col", "#c0392b"),

    palette = pred_palette(get_input("generate_palette", "viridis"),
                           get_input("generate_palette_rev", FALSE)),

    pch_val = {
      v <- get_input("generate_plot_pch", ".")
      if(v == ".") "." else as.numeric(v)
    },
    cex_val = get_input("generate_plot_cex", 0.3),
    bg_col = get_input("generate_plot_bg_col", "#B3B3B3"),

    map_bg_col = get_input("generate_plot_map_bg_col", "#F0F0F0"),

    ell_col = get_input("generate_plot_ell_col", "#000000"),
    ell_lwd = get_input("generate_plot_ell_lwd", 2),
    ell_lty = as.numeric(get_input("generate_plot_ell_lty", "1")),

    centroid_pch = as.numeric(get_input("generate_plot_centroid_pch", "8")),
    centroid_col = get_input("generate_plot_centroid_col", "#000000"),
    centroid_cex = get_input("generate_plot_centroid_cex", 1.5),

    zoom_mode_espace = get_input("generate_plot_zoom_mode_espace", "auto"),
    zoom_mode_combined = get_input("generate_plot_zoom_mode_combined", "auto"),
    asp_espace = get_input("generate_plot_asp_espace", "auto"),
    asp_combined = get_input("generate_plot_asp_combined", "auto"),

    export_cex = get_input("generate_export_cex", 1)

  )
}


# E-SPACE -----------------------------------------------------------------

output$generate_espace_plot_top_options_ui <- renderUI({

  vars <- generate_plot_vars()
  req(vars)

  # Which sets appear is chosen with the eye icons in the occurrence
  # summary, so nothing set-related belongs here
  fluidRow(
    column(width = 1, tags$span("Layout:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("generate_plot_espace_state",
                        label = NULL,
                        choiceNames = list(
                          tags$span("Pairs", class = "text-widget-inner"),
                          tags$span("2D", class = "text-widget-inner")
                        ),
                        choiceValues = c("generate_plot_pairs", "generate_plot_2d"),
                        selected = "generate_plot_pairs",
                        inline = TRUE)),

    conditionalPanel(
      "input.generate_plot_espace_state == 'generate_plot_2d'",
      column(width = 2, tags$span("Variable:", class = "text-widget-title")),
      column(width = 3,
             selectInput("generate_plot_2d_x", label = NULL, choices = vars,
                         selected = vars[1])),
      column(width = 3,
             selectInput("generate_plot_2d_y", label = NULL, choices = vars,
                         selected = vars[min(2, length(vars))]))
    )
  )
})

output$generate_espace_plot <- renderPlot({

  vars <- generate_plot_vars()
  req(vars)

  if(length(vars) < 2){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "E-space")
    text(0.5, 0.5, "Select at least two variables.", cex = 1, col = "grey50")
    return(invisible(NULL))
  }

  s <- generate_collect_plot_settings()
  req(s)

  show_sets <- generate_visible_sets()

  # Pairs view overlays every visible set on the same panels, colour-coded.
  # 2D view gives each its own panel.
  s$occ_set <- show_sets

  state <- if(!is.null(input$generate_plot_espace_state)){
    input$generate_plot_espace_state
  } else {
    "generate_plot_pairs"
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if(identical(state, "generate_plot_pairs")){

    generate_draw_espace_pairs(vars, s)

  } else {

    req(input$generate_plot_2d_x, input$generate_plot_2d_y)

    if(length(show_sets) > 0){

      generate_draw_espace_sets(input$generate_plot_2d_x,
                                input$generate_plot_2d_y,
                                s, show_sets)
    } else {

      par(mar = c(4, 4, 2, 1))
      generate_draw_espace_panel(input$generate_plot_2d_x,
                                 input$generate_plot_2d_y, s)
    }
  }
})

output$generate_espace_plot_bottom_options_ui <- renderUI({

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(
    column(width = 1,
           tags$span("Zoom:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("generate_plot_zoom_mode_espace", label = NULL,
                        choiceNames = list(
                          tags$span("Auto", class = "text-widget-inner"),
                          tags$span("Zoom in", class = "text-widget-inner")
                        ),
                        choiceValues = c("auto", "ellipsoid"),
                        selected = "auto",
                        inline = TRUE)),
    column(width = 4,
           tags$span(ell$ell_name, class = "text-center",
                     style = "font-size: 12px; color: #888; font-weight: 400;")),
    column(width = 1,
           tags$span("Aspect:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("generate_plot_asp_espace", label = NULL,
                        choiceNames = list(
                          tags$span("Auto", class = "text-widget-inner"),
                          tags$span("Fixed", class = "text-widget-inner")
                        ),
                        choiceValues = c("auto", "fixed"),
                        selected = "auto",
                        inline = TRUE))
  )
})


# G-SPACE -----------------------------------------------------------------

# Panel choice lives in the occurrence summary, so this row is empty. Kept
# as an output so the uiOutput in ui.R still resolves.
output$generate_gspace_plot_top_options_ui <- renderUI({
  NULL
})

output$generate_gspace_plot <- renderPlot({

  is_virtual <- identical(session_data$input_mode, "virtual")

  if(is_virtual || is.null(session_data$bg_raster)){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "G-space")
    text(0.5, 0.5,
         "Virtual mode on or no raster provided.\nG-space unavailable.",
         cex = 1, col = "grey50")
    return(invisible(NULL))
  }

  ell <- session_data$current_ellipsoid
  req(ell)

  sets <- generate_occ_sets()

  s <- generate_collect_plot_settings()

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if(length(sets) == 0){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "G-space")
    text(0.5, 0.5, "No occurrences generated for this ellipsoid.",
         cex = 1, col = "grey50")
    return(invisible(NULL))
  }

  show_sets <- generate_visible_sets()

  if(length(show_sets) == 0) show_sets <- sets[1]

  n <- length(show_sets)
  n_cols <- if(n <= 1) 1L else 2L
  n_rows <- ceiling(n / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2.5, 5))

  for(nm in show_sets){

    src <- generate_source_raster(ell$ell_id, nm)

    if(is.null(src)){
      plot.new()
      next
    }

    n_pts <- nrow(session_data$ellipsoid_occurrence_list[[ell$ell_id]][[nm]])

    generate_draw_gspace_panel(src, s,
                               title = paste0(generate_set_label(nm, s), " (n = ", n_pts, ")"),
                               set = nm)
  }

  if(n > 1 && n %% 2 != 0) plot.new()
})

output$generate_gspace_plot_bottom_options_ui <- renderUI({
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(class = "ell-row",
           column(width = 12,
                  tags$span(ell$ell_name, class = "text-center",
                            style = "font-size: 12px; color: #888; font-weight: 400;"))
  )

})


# COMBINED ----------------------------------------------------------------

output$generate_combined_plot_top_options_ui <- renderUI({

  ell_slot()

  req(session_data$bg_raster)

  vars <- generate_plot_vars()
  req(vars)

  fluidRow(
    column(width = 1, tags$span("Layout:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("generate_plot_combined_layout",
                        label = NULL,
                        choiceNames = list(
                          tags$span("Side by side", class = "text-widget-inner"),
                          tags$span("Stacked", class = "text-widget-inner")
                        ),
                        choiceValues = c("row", "col"),
                        selected = "row",
                        inline = TRUE)),
    column(width = 2, tags$span("Variables:", class = "text-widget-title")),
    column(width = 3,
           selectInput("generate_plot_combined_x", label = NULL,
                       choices = vars, selected = vars[1])),
    column(width = 3,
           selectInput("generate_plot_combined_y", label = NULL,
                       choices = vars, selected = vars[min(2, length(vars))]))
  )
})

output$generate_combined_plot <- renderPlot({

  is_virtual <- identical(session_data$input_mode, "virtual")

  if(is_virtual || is.null(session_data$bg_raster)){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "Combined")
    text(0.5, 0.5,
         "Virtual mode on or no raster provided.\nCombined view unavailable.",
         cex = 1, col = "grey50")
    return(invisible(NULL))
  }

  ell <- session_data$current_ellipsoid
  req(ell)

  vars <- generate_plot_vars()
  req(vars, length(vars) >= 2)
  req(input$generate_plot_combined_x, input$generate_plot_combined_y)

  sets <- generate_occ_sets()
  req(length(sets) > 0)

  s <- generate_collect_plot_settings()
  s$asp_espace <- s$asp_combined
  s$zoom_mode_espace <- s$zoom_mode_combined

  # Both halves show the same sets, chosen with the eye icons
  e_sets <- generate_visible_sets()
  g_sets <- e_sets

  if(length(e_sets) == 0) e_sets <- sets[1]
  if(length(g_sets) == 0) g_sets <- sets[1]

  n_e <- length(e_sets)
  n_g <- length(g_sets)

  e_block <- plot_panel_block(seq_len(n_e))
  g_block <- plot_panel_block(n_e + seq_len(n_g))

  layout_mode <- if(!is.null(input$generate_plot_combined_layout)){
    input$generate_plot_combined_layout
  } else {
    "row"
  }

  m <- if(identical(layout_mode, "row")){
    cbind(e_block, g_block)
  } else {
    rbind(e_block, g_block)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit({ layout(1); par(old_par) })

  layout(m)

  occ_list <- session_data$ellipsoid_occurrence_list[[ell$ell_id]]

  # E-space, one panel per selected set
  par(mar = c(4, 4, 2, 1))

  for(nm in e_sets){
    n_pts <- if(!is.null(occ_list[[nm]])) nrow(occ_list[[nm]]) else 0L
    generate_draw_espace_panel(input$generate_plot_combined_x,
                               input$generate_plot_combined_y,
                               s,
                               set = nm,
                               title = paste0(generate_set_label(nm, s), " (n = ", n_pts, ")"))
  }

  # G-space, one panel per selected set, each over the surface it came from
  par(mar = c(3, 3, 2.5, 5))

  for(nm in g_sets){

    src <- generate_source_raster(ell$ell_id, nm)

    if(is.null(src)){
      plot.new()
      next
    }

    n_pts <- if(!is.null(occ_list[[nm]])) nrow(occ_list[[nm]]) else 0L

    generate_draw_gspace_panel(src, s,
                               title = paste0(generate_set_label(nm, s), " (n = ", n_pts, ")"),
                               set = nm)
  }
})

output$generate_combined_plot_bottom_options_ui <- renderUI({

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)
  req(session_data$bg_raster)

  fluidRow(
    column(width = 1,
           tags$span("Zoom:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("generate_plot_zoom_mode_combined", label = NULL,
                        choiceNames = list(
                          tags$span("Auto", class = "text-widget-inner"),
                          tags$span("Zoom in", class = "text-widget-inner")
                        ),
                        choiceValues = c("auto", "ellipsoid"),
                        selected = "auto",
                        inline = TRUE)),
    column(width = 4,
           tags$span(ell$ell_name, class = "text-center",
                     style = "font-size: 12px; color: #888; font-weight: 400;")),
    column(width = 1,
           tags$span("Aspect:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("generate_plot_asp_combined", label = NULL,
                        choiceNames = list(
                          tags$span("Auto", class = "text-widget-inner"),
                          tags$span("Fixed", class = "text-widget-inner")
                        ),
                        choiceValues = c("auto", "fixed"),
                        selected = "auto",
                        inline = TRUE))
  )
})


# ADVANCED SETTINGS -------------------------------------------------------

output$generate_plot_settings_ui <- renderUI({

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  color_input <- function(id, label, value){
    column(width = 3,
           tags$span(label, class = "text-widget-title"),
           tags$div(
             style = "display: flex; align-items: center; gap: 8px;",
             tags$input(type = "color",
                        value = value,
                        oninput = sprintf("Shiny.setInputValue('%s', this.value)", id),
                        style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px; cursor: pointer;")
           ))
  }

  box(title = tagList(
    tags$span("Advanced Plot Settings", class = "text-section-header"),
    tags$span(icon("circle-info"),
              title = instructions$plot_settings,
              class = "tooltip-icon")),
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE,

    # Occurrence points
    fluidRow(
      column(width = 3,
             tags$span("Occurrence shape", class = "text-widget-title"),
             selectInput("generate_occ_pch", label = NULL,
                         choices = c("Filled circle" = "16",
                                     "Open circle" = "1",
                                     "Cross" = "4",
                                     "Filled triangle" = "17"),
                         selected = "16")),
      column(width = 3,
             tags$span("Occurrence size", class = "text-widget-title"),
             numericInput("generate_occ_cex", label = NULL, value = 0.7,
                          min = 0.1, max = 5, step = 0.1)),
      color_input("generate_occ_col", "Occurrence color", "#c0392b"),
      color_input("generate_plot_bg_col", "Background point color", "#B3B3B3")
    ),

    # Background points and palette
    fluidRow(
      column(width = 3,
             tags$span("Point shape (pch)", class = "text-widget-title"),
             selectInput("generate_plot_pch", label = NULL,
                         choices = c("Dot (.)" = ".",
                                     "Open circle" = "1",
                                     "Filled circle" = "16",
                                     "Square" = "15"),
                         selected = ".")),
      column(width = 3,
             tags$span("Point size (cex)", class = "text-widget-title"),
             numericInput("generate_plot_cex", label = NULL, value = 0.3,
                          min = 0.1, max = 5, step = 0.1)),
      column(width = 3,
             tags$span("Continuous palette", class = "text-widget-title"),
             selectInput("generate_palette", label = NULL,
                         choices = c("Viridis" = "viridis",
                                     "Plasma" = "Plasma",
                                     "Inferno" = "Inferno",
                                     "Rocket" = "Rocket",
                                     "Mako" = "Mako",
                                     "YlGnBu" = "YlGnBu"),
                         selected = "viridis")),
      column(width = 3,
             tags$span("Direction", class = "text-widget-title"),
             checkboxInput("generate_palette_rev", "Reverse palette",
                           value = FALSE))
    ),

    # Ellipsoid and centroid
    fluidRow(
      column(width = 3,
             checkboxInput("generate_show_ellipsoid", "Show ellipsoid",
                           value = TRUE)),
      color_input("generate_plot_ell_col", "Ellipsoid color", "#000000"),
      column(width = 3,
             tags$span("Ellipsoid line width", class = "text-widget-title"),
             numericInput("generate_plot_ell_lwd", label = NULL, value = 2,
                          min = 0.5, max = 6, step = 0.5)),
      column(width = 3,
             tags$span("Ellipsoid line type", class = "text-widget-title"),
             selectInput("generate_plot_ell_lty", label = NULL,
                         choices = c("Solid" = "1",
                                     "Dashed" = "2",
                                     "Dotted" = "3"),
                         selected = "1"))
    ),

    fluidRow(
      column(width = 3,
             checkboxInput("generate_show_centroid", "Show centroid",
                           value = TRUE)),
      column(width = 3,
             tags$span("Centroid shape (pch)", class = "text-widget-title"),
             selectInput("generate_plot_centroid_pch", label = NULL,
                         choices = c("Star" = "8",
                                     "Cross (X)" = "4",
                                     "Filled diamond" = "18"),
                         selected = "8")),
      color_input("generate_plot_centroid_col", "Centroid color", "#000000"),
      color_input("generate_plot_map_bg_col", "Map background (NA)", "#F0F0F0")
    ),

    br(),

    fluidRow(
      column(width = 12,
             class = "btn-spaced",
             actionButton("generate_open_export_modal",
                          tagList(icon("download"), "Export figure"),
                          class = "btn-default"))
    )
  )
})


# EXPORT ------------------------------------------------------------------

observeEvent(input$generate_open_export_modal, {

  showModal(modalDialog(
    title = "Export Figure",
    size = "m",

    fluidRow(
      column(width = 6,
             tags$span("File type", class = "text-widget-title"),
             radioButtons("generate_export_filetype", label = NULL,
                          choiceNames = list(
                            tags$span("PNG", class = "text-widget-inner"),
                            tags$span("PDF", class = "text-widget-inner"),
                            tags$span("SVG", class = "text-widget-inner")
                          ),
                          choiceValues = c("png", "pdf", "svg"),                          selected = if(!is.null(input$generate_export_filetype))
                            input$generate_export_filetype else "png",
                          inline = TRUE)),
      column(width = 6,
             tags$span("Unit", class = "text-widget-title"),
             radioButtons("generate_export_unit", label = NULL,
                          choiceNames = list(
                            tags$span("mm", class = "text-widget-inner"),
                            tags$span("inches", class = "text-widget-inner"),
                            tags$span("px", class = "text-widget-inner")
                          ),
                          choiceValues = c("mm", "in", "px"),
                          selected = if(!is.null(input$generate_export_unit))
                            input$generate_export_unit else "mm",
                          inline = TRUE))
    ),

    fluidRow(
      column(width = 4,
             tags$span("Text size (cex)", class = "text-widget-title"),
             numericInput("generate_export_cex", label = NULL,
                          value = if(!is.null(input$generate_export_cex))
                            input$generate_export_cex else 1,
                          min = 0.5, max = 3, step = 0.1))
    ),

    uiOutput("generate_export_settings_ui"),

    footer = tagList(
      modalButton("Cancel"),
      downloadButton("generate_confirm_export", "Export", class = "btn-save",
                     onclick = "Shiny.setInputValue('generate_export_done',
                     Math.random(), {priority: 'event'});")
    ),
    easyClose = TRUE
  ))
})

output$generate_export_settings_ui <- renderUI({

  ext <- if(!is.null(input$generate_export_filetype)){
    input$generate_export_filetype
  } else {
    "png"
  }

  unit <- if(!is.null(input$generate_export_unit)){
    input$generate_export_unit
  } else {
    "mm"
  }

  defaults <- switch(unit,
                     "mm" = list(val = 166, min = 50, max = 500, step = 1),
                     "in" = list(val = 6.54, min = 1, max = 20, step = 0.1),
                     "px" = list(val = 1961, min = 400, max = 6000, step = 100))

  tagList(
    p(tagList(icon("circle-info"),
              " Default is a standard full-page publication figure (166 x 166 mm)."),
      style = "font-size: 12px; color: #666; margin-bottom: 8px;"),

    fluidRow(
      column(width = if(ext == "png") 4 else 6,
             tags$span(paste0("Width (", unit, ")"), class = "text-widget-title"),
             numericInput("generate_export_width_val", label = NULL,
                          value = defaults$val, min = defaults$min,
                          max = defaults$max, step = defaults$step)),
      column(width = if(ext == "png") 4 else 6,
             tags$span(paste0("Height (", unit, ")"), class = "text-widget-title"),
             numericInput("generate_export_height_val", label = NULL,
                          value = defaults$val, min = defaults$min,
                          max = defaults$max, step = defaults$step)),
      if(ext == "png"){
        column(width = 4,
               tags$span("Resolution (dpi)", class = "text-widget-title"),
               numericInput("generate_export_res", label = NULL,
                            value = 300, min = 72, max = 600, step = 50))
      }
    ),

    if(ext != "png") p("PDF and SVG do not require a resolution setting.",
                       class = "text-instruction")
  )
})

output$generate_confirm_export <- downloadHandler(

  filename = function(){
    ext <- if(!is.null(input$generate_export_filetype)){
      input$generate_export_filetype
    } else {
      "png"
    }
    tab <- if(!is.null(input$generate_plot_tabs)){
      input$generate_plot_tabs
    } else {
      "generate_espace_plot_tab"
    }
    prefix <- switch(tab,
                     "generate_espace_plot_tab" = "generate_espace_plot",
                     "generate_gspace_plot_tab" = "generate_gspace_plot",
                     "generate_combined_plot_tab" = "generate_combined_plot",
                     "generate_plot")
    paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
  },

  content = function(file){

    ext <- if(!is.null(input$generate_export_filetype)){
      input$generate_export_filetype
    } else {
      "png"
    }
    tab <- if(!is.null(input$generate_plot_tabs)){
      input$generate_plot_tabs
    } else {
      "generate_espace_plot_tab"
    }

    plot_open_device(file, ext,
                     w_val = input$generate_export_width_val,
                     h_val = input$generate_export_height_val,
                     unit = input$generate_export_unit,
                     res = input$generate_export_res,
                     cex_val = input$generate_export_cex)
    on.exit(dev.off())

    s <- generate_collect_plot_settings()
    req(s)

    s$occ_set <- generate_visible_sets()

    vars <- generate_plot_vars()
    req(vars)

    ell <- session_data$current_ellipsoid
    req(ell)

    sets <- generate_occ_sets()

    switch(tab,

           "generate_espace_plot_tab" = {

             state <- if(!is.null(input$generate_plot_espace_state)){
               input$generate_plot_espace_state
             } else {
               "generate_plot_pairs"
             }

             if(identical(state, "generate_plot_pairs")){
               generate_draw_espace_pairs(vars, s)
             } else if(length(sets) > 0){
               show_sets <- generate_visible_sets()
               generate_draw_espace_sets(input$generate_plot_2d_x,
                                         input$generate_plot_2d_y, s, show_sets)
             } else {
               par(mar = c(4, 4, 2, 1))
               generate_draw_espace_panel(input$generate_plot_2d_x,
                                          input$generate_plot_2d_y, s)
             }
           },

           "generate_gspace_plot_tab" = {

             req(length(sets) > 0)

             show_sets <- generate_visible_sets()
             if(length(show_sets) == 0) show_sets <- unname(sets)[1]

             n <- length(show_sets)
             n_cols <- if(n <= 1) 1L else 2L
             par(mfrow = c(ceiling(n / n_cols), n_cols), mar = c(3, 3, 2.5, 5))

             for(nm in show_sets){
               src <- generate_source_raster(ell$ell_id, nm)
               if(is.null(src)){ plot.new(); next }
               n_pts <- nrow(session_data$ellipsoid_occurrence_list[[ell$ell_id]][[nm]])
               generate_draw_gspace_panel(src, s,
                                          title = paste0(generate_set_label(nm, s), " (n = ", n_pts, ")"),
                                          set = nm)
             }

             if(n > 1 && n %% 2 != 0) plot.new()
           },

           "generate_combined_plot_tab" = {

             req(length(sets) > 0)
             req(input$generate_plot_combined_x, input$generate_plot_combined_y)

             s$asp_espace <- s$asp_combined
             s$zoom_mode_espace <- s$zoom_mode_combined

             e_sets <- generate_visible_sets()
             g_sets <- e_sets
             if(length(e_sets) == 0) e_sets <- unname(sets)[1]
             if(length(g_sets) == 0) g_sets <- unname(sets)[1]

             e_block <- plot_panel_block(seq_len(length(e_sets)))
             g_block <- plot_panel_block(length(e_sets) + seq_len(length(g_sets)))

             lay <- if(!is.null(input$generate_plot_combined_layout)){
               input$generate_plot_combined_layout
             } else {
               "row"
             }

             layout(if(identical(lay, "row")) cbind(e_block, g_block)
                    else rbind(e_block, g_block))

             occ_list <- session_data$ellipsoid_occurrence_list[[ell$ell_id]]

             par(mar = c(4, 4, 2, 1))
             for(nm in e_sets){
               n_pts <- if(!is.null(occ_list[[nm]])) nrow(occ_list[[nm]]) else 0L
               generate_draw_espace_panel(input$generate_plot_combined_x,
                                          input$generate_plot_combined_y, s,
                                          set = nm,
                                          title = paste0(generate_set_label(nm, s), " (n = ", n_pts, ")"))
             }

             par(mar = c(3, 3, 2.5, 5))
             for(nm in g_sets){
               src <- generate_source_raster(ell$ell_id, nm)
               if(is.null(src)){ plot.new(); next }
               n_pts <- if(!is.null(occ_list[[nm]])) nrow(occ_list[[nm]]) else 0L
               generate_draw_gspace_panel(src, s,
                                          title = paste0(generate_set_label(nm, s), " (n = ", n_pts, ")"),
                                          set = nm)
             }

             layout(1)
           }
    )
  }
)

observeEvent(input$generate_export_done, {

  removeModal()

  showNotification("Figure exported successfully.",
                   type = "message", duration = 4)
})
