# Title: Predict tab plot logic
# Description: E-space, G-space, and combined plots for the predict tab.
#              Mirrors build_tab_plot.R, with no range lines and with
#              overlays driven by the stored predictions.
# Date last updated: 08/05/2026

# Functions -----------------------------------------------------------------

# Continuous palette for prediction layers. hcl.colors is base R since 3.6,
# so no extra dependency is needed.
pred_palette <- function(name = "viridis", reverse = FALSE, n = 100){
  if(!name %in% hcl.pals()) name <- "viridis"
  cols <- grDevices::hcl.colors(n, palette = name)
  if(isTRUE(reverse)) rev(cols) else cols
}

# Maps a numeric vector onto palette colours. Returns NA for NA values so
# points() skips them.
pred_colors <- function(vals, pal, rng = NULL){
  if(is.null(rng)) rng <- range(vals, na.rm = TRUE)
  if(!all(is.finite(rng)) || diff(rng) == 0) return(rep(pal[1], length(vals)))
  idx <- cut(vals, breaks = seq(rng[1], rng[2], length.out = length(pal) + 1),
             labels = FALSE, include.lowest = TRUE)
  pal[idx]
}

# Standalone colour bar drawn in its own panel. Called after the plot grid
# so every panel shares one legend.
pred_legend_panel <- function(rng, pal, label = ""){

  if(is.null(rng) || !all(is.finite(rng))) return(invisible(NULL))

  old_mar <- par("mar")
  on.exit(par(mar = old_mar))
  par(mar = c(2.5, 4, 0.5, 4))

  plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "")

  xs <- seq(0.15, 0.85, length.out = length(pal) + 1)
  rect(xs[-length(xs)], 0.45, xs[-1], 0.8, col = pal, border = NA)
  rect(0.15, 0.45, 0.85, 0.8, border = "#888", lwd = 0.5)

  ticks <- seq(0.15, 0.85, length.out = 5)
  vals <- seq(rng[1], rng[2], length.out = 5)
  segments(ticks, 0.45, ticks, 0.38, col = "#888", lwd = 0.5)
  text(ticks, 0.32, format(round(vals, 2)), cex = 0.65, adj = c(0.5, 1))

  if(nzchar(label)) text(0.5, 0.9, label, cex = 0.75, adj = c(0.5, 0))
}

# Value range for the overlay, shared by every panel in a grid
pred_layer_range <- function(vals_df, layer){
  if(is.null(vals_df) || is.null(layer)) return(NULL)
  if(!layer %in% names(vals_df)) return(NULL)
  rng <- range(vals_df[[layer]], na.rm = TRUE)
  if(!all(is.finite(rng))) return(NULL)
  rng
}

predict_draw_espace_panel <- function(v1, v2, s, layer = NULL, title = NULL){

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

  # Overlay. Values come from the stored prediction, never a fresh predict()
  vals_df <- s$espace_vals
  lyr <- if(!is.null(layer)) layer else s$espace_layer

  if(!is.null(bg) && !is.null(lyr)){

    if(identical(lyr, "binary")){

      # Binary suitability is the same regardless of which layer is shown
      base_lyr <- if("suitability_trunc" %in% names(vals_df)){
        "suitability_trunc"
      # } else if("suitability" %in% names(vals_df)){
      #   "suitability"
      } else {
        NULL
      }

      if(!is.null(base_lyr)){
        v <- vals_df[[base_lyr]]
        keep <- !is.na(v) & v > 0
        if(any(keep)){
          points(bg[[v1]][keep], bg[[v2]][keep],
                 col = s$suitable_col,
                 pch = s$suitable_pch,
                 cex = s$suitable_cex)
        }
      }

    } else if(lyr %in% names(vals_df)){

      v <- vals_df[[lyr]]
      keep <- !is.na(v)

      if(any(keep)){
        cols <- pred_colors(v[keep], s$palette, rng = s$layer_rng)
        points(bg[[v1]][keep], bg[[v2]][keep],
               col = cols,
               pch = s$suitable_pch,
               cex = s$suitable_cex)
      }
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

predict_draw_espace_pairs <- function(vars, s){

  pairs <- t(combn(seq_along(vars), 2))
  n_pairs <- nrow(pairs)
  n_cols <- ceiling(sqrt(n_pairs))
  n_rows <- ceiling(n_pairs / n_cols)

  show_bar <- !is.null(s$layer_rng)

  if(show_bar){
    # Panels numbered row-wise to match mfrow ordering. A 0 leaves a cell
    # empty and creates no region, so no padding calls are needed. The last
    # row is one short full-width region for the shared legend.
    cells <- c(seq_len(n_pairs), rep(0L, n_rows * n_cols - n_pairs))
    m <- matrix(cells, nrow = n_rows, ncol = n_cols, byrow = TRUE)
    m <- rbind(m, rep(n_pairs + 1L, n_cols))
    layout(m, heights = c(rep(1, n_rows), 0.3))
  } else {
    par(mfrow = c(n_rows, n_cols))
  }

  par(mar = c(4, 4, 2, 1))

  for(i in seq_len(n_pairs)){
    predict_draw_espace_panel(vars[pairs[i, 1]], vars[pairs[i, 2]], s)
  }

  if(show_bar) pred_legend_panel(s$layer_rng, s$palette, s$espace_layer)

  layout(1)
}

# One panel per predicted layer for the same variable pair
predict_draw_espace_layers <- function(v1, v2, s, layers){

  n <- length(layers)
  n_cols <- if(n <= 1) 1L else 2L
  n_rows <- ceiling(n / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 2, 1))

  for(lyr in layers){
    predict_draw_espace_panel(v1, v2, s, layer = lyr, title = lyr)
  }

  if(n %% 2 != 0 && n > 1) plot.new()
}


predict_draw_gspace_panel <- function(rast, s, title = NULL,
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
}

# Binary suitable area from a prediction stack. Uses the truncated layer
# when present, since that is what defines the ellipsoid boundary.
predict_binary_layer <- function(pred){

  if(is.null(pred) || !inherits(pred, "SpatRaster")) return(NULL)

  base_lyr <- if("suitability_trunc" %in% names(pred)){
    "suitability_trunc"
  } else if("suitability" %in% names(pred)){
    "suitability"
  } else {
    return(NULL)
  }

  binary <- pred[[base_lyr]] > 0
  names(binary) <- "suitable"
  binary
}


# Reactives ---------------------------------------------------------------

# On-the-fly suitability for G-space, used only until a real prediction
# exists for this ellipsoid. Returns NULL once one does, so the stored
# prediction is used instead and this never runs twice for the same work.
predict_raster_vis <- reactive({

  ell <- session_data$current_ellipsoid
  req(ell)
  req(session_data$bg_raster)

  if(!is.null(session_data$ellipsoid_prediction_list[[ell$ell_id]])) return(NULL)

  res <- tryCatch(
    predict(ell,
            newdata = terra::subset(session_data$bg_raster, ell$var_names),
            include_suitability = TRUE,
            include_mahalanobis = FALSE,
            suitability_truncated = TRUE,
            verbose = FALSE),
    error = function(e) NULL
  )

  if(is.null(res)) return(NULL)

  Reduce(c, res)
})

# Variables to plot. The current ellipsoid's own set, so a version built
# from different variables still plots correctly.
predict_plot_vars <- reactive({

  ell <- session_data$current_ellipsoid

  if(!is.null(ell) && !is.null(ell$var_names)) return(ell$var_names)

  session_data$vars
})

# Layer names available for the current ellipsoid's stored prediction
predict_pred_layers <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(character(0))

  pred <- session_data$ellipsoid_prediction_list[[ell$ell_id]]
  if(is.null(pred)) return(character(0))

  # keep_data = TRUE means predict() stacks the predictor layers in front of
  # the outputs, so they are stripped here exactly as in the data frame branch
  if(inherits(pred, "SpatRaster")){
    return(setdiff(names(pred), ell$var_names))
  }

  if(is.data.frame(pred)){
    return(setdiff(names(pred), c("x", "y", ell$var_names)))
  }

  character(0)
})

# Prediction values aligned to bg_df rows, for the E-space overlay.
# Extracted once per ellipsoid and cached, so changing a plot setting does
# not re-extract and predict() is never re-run here.
predict_espace_values <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(NULL)

  pred <- session_data$ellipsoid_prediction_list[[ell$ell_id]]
  if(is.null(pred)) return(NULL)

  bg <- session_data$bg_df
  if(is.null(bg)) return(NULL)

  # CSV-only mode: predict() returned a data frame already aligned to bg_df
  if(is.data.frame(pred)) return(pred)

  if(!inherits(pred, "SpatRaster")) return(NULL)
  if(!all(c("x", "y") %in% colnames(bg))) return(NULL)

  tryCatch(
    terra::extract(pred, bg[, c("x", "y")], ID = FALSE),
    error = function(e) NULL
  )
})

# Called at the top of every draw function. No range lines on this tab:
# show_lines is always inactive so compute_lims from build_tab_plot.R can
# be reused unchanged.
predict_collect_plot_settings <- function(){

  has_ell <- !is.null(session_data$current_ellipsoid)

  list(
    has_ell = has_ell,
    ell = if(has_ell) session_data$current_ellipsoid else NULL,

    show_lines = list(active = FALSE, ranges = NULL),

    show_ell = has_ell && get_input("predict_show_ellipsoid", TRUE),
    show_centroid = has_ell && get_input("predict_show_centroid", TRUE),

    espace_vals = predict_espace_values(),
    espace_layer = get_input("predict_espace_layer", NULL),
    layer_rng = NULL,
    palette = pred_palette(get_input("predict_palette", "viridis"),
                           get_input("predict_palette_rev", FALSE)),

    pch_val = {
      v <- get_input("predict_plot_pch", ".")
      if(v == ".") "." else as.numeric(v)
    },
    cex_val = get_input("predict_plot_cex", 0.3),
    bg_col = get_input("predict_plot_bg_col", "#B3B3B3"),

    suitable_pch = as.numeric(get_input("predict_plot_suitable_pch", "16")),
    suitable_cex = get_input("predict_plot_suitable_cex", 0.3),
    suitable_col = get_input("predict_plot_suitable_col", "#097a21"),
    unsuitable_col = get_input("predict_plot_unsuitable_col", "#D3D3D3"),
    map_bg_col = get_input("predict_plot_map_bg_col", "#F0F0F0"),

    ell_col = get_input("predict_plot_ell_col", "#000000"),
    ell_lwd = get_input("predict_plot_ell_lwd", 2),
    ell_lty = as.numeric(get_input("predict_plot_ell_lty", "1")),

    centroid_pch = as.numeric(get_input("predict_plot_centroid_pch", "8")),
    centroid_col = get_input("predict_plot_centroid_col", "#000000"),
    centroid_cex = get_input("predict_plot_centroid_cex", 1.5),

    zoom_mode_espace = get_input("predict_plot_zoom_mode_espace", "auto"),
    zoom_mode_combined = get_input("predict_plot_zoom_mode_combined", "auto"),

    asp_espace = get_input("predict_plot_asp_espace", "auto"),
    asp_combined = get_input("predict_plot_asp_combined", "auto")
  )
}


# SELECTOR AND WORKING SLOT -----------------------------------------------

# Picking a specific ellipsoid in the selector loads it into the working
# slot so the plots follow. "All versions" leaves the slot alone, so the
# library keeps control in that case.
observeEvent(input$predict_ellipsoid_selected, {

  sel <- input$predict_ellipsoid_selected
  req(sel)

  if(identical(sel, "all")) return()

  ell <- session_data$ellipsoid_list[[sel]]
  req(ell)

  if(identical(session_data$current_ellipsoid$ell_id, sel)) return()

  set_working_ellipsoid(ell, mode = "view")
})

# The reverse direction: clicking view in the library moves the selector,
# and the layer checkboxes follow whatever that ellipsoid already predicted.
observeEvent(ell_slot(), {

  ell <- session_data$current_ellipsoid
  req(ell)

  if(!ell$ell_id %in% names(session_data$ellipsoid_list)) return()

  if(!identical(input$predict_ellipsoid_selected, ell$ell_id)){
    updateSelectInput(session, "predict_ellipsoid_selected",
                      selected = ell$ell_id)
  }

  lyrs <- predict_pred_layers()

  updateCheckboxInput(session, "predict_suitability",
                      value = "suitability" %in% lyrs)
  updateCheckboxInput(session, "predict_suitability_trunc",
                      value = "suitability_trunc" %in% lyrs)
  updateCheckboxInput(session, "predict_mahalanobis",
                      value = "mahalanobis" %in% lyrs)
  updateCheckboxInput(session, "predict_mahalanobis_trunc",
                      value = "mahalanobis_trunc" %in% lyrs)

}, ignoreInit = TRUE)


# E-SPACE and G-SPACE -----------------------------------------------------------------

output$predict_espace_plot_top_options_ui <- renderUI({

  ell_slot()

  vars <- predict_plot_vars()
  req(vars)

  lyrs <- predict_pred_layers()
  has_pred <- length(lyrs) > 0

  # Overlay choices only exist once something has been predicted
  layer_choices <- if(has_pred){
    c(setNames("binary", "Binary suitability"),
      setNames(lyrs, lyrs),
      setNames("none", "No overlay"))
  } else {
    NULL
  }

  fluidRow(

    column(width = 1, tags$span("Layout:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("predict_plot_espace_state",
                        label = NULL,
                        choiceNames = list(
                          tags$span("Pairs", class = "text-widget-inner"),
                          tags$span("2D", class = "text-widget-inner")
                        ),
                        choiceValues = c("predict_plot_pairs", "predict_plot_2d"),
                        selected = "predict_plot_pairs",
                        inline = TRUE)),

    conditionalPanel(
      "input.predict_plot_espace_state == 'predict_plot_2d'",
      column(width = 2, tags$span("Variables:", class = "text-widget-title")),
      column(width = 3,
             selectInput("predict_plot_2d_x", label = NULL, choices = vars,
                         selected = vars[1])),
      column(width = 3,
             selectInput("predict_plot_2d_y", label = NULL, choices = vars,
                         selected = vars[min(2, length(vars))]))
    ),

    if(has_pred){
      # In 2D every predicted layer gets its own panel, so the selector
      # only applies to the pairs view
      conditionalPanel(
        "input.predict_plot_espace_state == 'predict_plot_pairs'",
        column(width = 2, tags$span("Prediction:", class = "text-widget-title")),
        column(width = 4,
               selectInput("predict_espace_layer",
                           label = NULL,
                           choices = layer_choices,
                           selected = "binary"))
      )
    }
  )
})

output$predict_espace_plot <- renderPlot({

  vars <- predict_plot_vars()
  req(vars)

  if(length(vars) < 2){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "E-space")
    text(0.5, 0.5, "Select at least two variables.", cex = 1, col = "grey50")
    return(invisible(NULL))
  }

  s <- predict_collect_plot_settings()
  req(s)

  if(identical(s$espace_layer, "none")) s$espace_layer <- NULL

  if(!is.null(s$espace_layer) && !identical(s$espace_layer, "binary")){
    s$layer_rng <- pred_layer_range(s$espace_vals, s$espace_layer)
  }

  lyrs <- predict_pred_layers()

  state <- if(!is.null(input$predict_plot_espace_state)){
    input$predict_plot_espace_state
  } else {
    "predict_plot_pairs"
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if(identical(state, "predict_plot_pairs")){

    predict_draw_espace_pairs(vars, s)

  } else {

    req(input$predict_plot_2d_x, input$predict_plot_2d_y)

    if(length(lyrs) > 0){
      # One panel per predicted layer for this variable pair
      predict_draw_espace_layers(input$predict_plot_2d_x,
                                 input$predict_plot_2d_y,
                                 s, lyrs)
    } else {
      par(mar = c(4, 4, 2, 1))
      predict_draw_espace_panel(input$predict_plot_2d_x,
                                input$predict_plot_2d_y, s)
    }
  }
})

output$predict_espace_plot_bottom_options_ui <- renderUI({
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(
    column(width = 1,
           tags$span("Zoom:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("predict_plot_zoom_mode_espace", label = NULL,
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
           radioButtons("predict_plot_asp_espace", label = NULL,
                        choiceNames = list(
                          tags$span("Auto", class = "text-widget-inner"),
                          tags$span("Fixed", class = "text-widget-inner")
                        ),
                        choiceValues = c("auto", "fixed"),
                        selected = "auto",
                        inline = TRUE))
  )
})

output$predict_gspace_plot_top_options_ui <- renderUI({

  ell_slot()

  req(session_data$bg_raster)
  req(session_data$current_ellipsoid)

  lyrs <- predict_pred_layers()

  # Nothing to switch between until a prediction exists: the map is the
  # on-the-fly suitable area either way
  req(length(lyrs) > 0)

  has_binary <- any(c("suitability_trunc", "Mahalanobis_trunc") %in% lyrs)

  choice_names <- list(tags$span("Prediction layers", class = "text-widget-inner"))
  choice_values <- "layers"

  if(has_binary){
    choice_names <- c(choice_names,
                      list(tags$span("Suitable area", class = "text-widget-inner")))
    choice_values <- c(choice_values, "binary")
  }

  fluidRow(
    column(width = 1, tags$span("Show:", class = "text-widget-title")),
    column(width = 6,
           radioButtons("predict_plot_gspace_show",
                        label = NULL,
                        choiceNames = choice_names,
                        choiceValues = choice_values,
                        selected = "layers",
                        inline = TRUE))
  )

})

output$predict_gspace_plot <- renderPlot({

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

  if(is.null(ell)){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "G-space")
    text(0.5, 0.5, "No ellipsoid selected.", cex = 1, col = "grey50")
    return(invisible(NULL))
  }

  s <- predict_collect_plot_settings()

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  pred <- session_data$ellipsoid_prediction_list[[ell$ell_id]]

  # Nothing predicted yet: preview the suitable area from the ellipsoid
  if(is.null(pred) || !inherits(pred, "SpatRaster")){

    vis <- predict_raster_vis()

    if(is.null(vis)){
      plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
           xlab = "", ylab = "", main = "G-space")
      text(0.5, 0.5, "Preview unavailable for this ellipsoid.",
           cex = 1, col = "grey50")
      return(invisible(NULL))
    }

    binary <- predict_binary_layer(vis)
    req(binary)

    par(mar = c(4, 4, 2, 4))
    predict_draw_gspace_panel(binary, s,
                              title = "Suitable area (preview)",
                              col = c(s$unsuitable_col, s$suitable_col),
                              legend = FALSE)
    return(invisible(NULL))
  }

  show_mode <- if(!is.null(input$predict_plot_gspace_show)){
    input$predict_plot_gspace_show
  } else {
    "layers"
  }

  if(identical(show_mode, "binary")){

    binary <- predict_binary_layer(pred)

    if(is.null(binary)){
      plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
           xlab = "", ylab = "", main = "G-space")
      text(0.5, 0.5, "No suitability layer in this prediction.",
           cex = 1, col = "grey50")
      return(invisible(NULL))
    }

    par(mar = c(4, 4, 2, 4))
    predict_draw_gspace_panel(binary, s,
                              title = "Suitable area",
                              col = c(s$unsuitable_col, s$suitable_col),
                              legend = FALSE)
    return(invisible(NULL))
  }

  # One panel per predicted layer, each with its own legend since the
  # scales differ between suitability and Mahalanobis
  lyrs <- predict_pred_layers()
  n <- length(lyrs)

  n_cols <- if(n <= 1) 1L else 2L
  n_rows <- ceiling(n / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 4))

  for(lyr in lyrs){
    predict_draw_gspace_panel(pred[[lyr]], s, title = lyr)
  }

  if(n > 1 && n %% 2 != 0) plot.new()
})

output$predict_gspace_plot_bottom_options_ui <- renderUI({
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(class = "ell-row",
           column(width = 12,
           tags$span(ell$ell_name, class = "text-center",
                     style = "font-size: 12px; color: #888; font-weight: 400;"))
  )

})

output$predict_combined_plot_top_options_ui <- renderUI({

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  vars <- predict_plot_vars()
  req(vars)

  lyrs <- predict_pred_layers()
  has_binary <- any(c("suitability", "suitability_trunc") %in% lyrs)

  # Each side picks its own layers, so E-space can show the binary while
  # the map shows a continuous one
  choices <- character(0)
  if(has_binary || length(lyrs) == 0) choices <- c("Binary suitability" = "binary")
  if(length(lyrs) > 0) choices <- c(choices, setNames(lyrs, lyrs))

  tagList(
    fluidRow(
      column(width = 1, tags$span("Layout:", class = "text-widget-title")),
      column(width = 3,
             radioButtons("predict_plot_combined_layout",
                          label = NULL,
                          choiceNames = list(
                            tags$span("Side by side", class = "text-widget-inner"),
                            tags$span("Stacked", class = "text-widget-inner")
                          ),
                          choiceValues = c("row", "col"),
                          selected = "row",
                          inline = FALSE)),

      column(width = 2, tags$span("Variables:", class = "text-widget-title")),
      column(width = 3,
             selectInput("predict_plot_combined_x", label = NULL,
                         choices = vars, selected = vars[1])),
      column(width = 3,
             selectInput("predict_plot_combined_y", label = NULL,
                         choices = vars, selected = vars[min(2, length(vars))]))
    ),

    fluidRow(
      column(width = 6,
             checkboxGroupInput("predict_combined_espace_layers",
                                label = tags$span("E-space shows:",
                                                  class = "text-widget-title"),
                                choices = choices,
                                selected = choices[1],
                                inline = TRUE)),
      column(width = 6,
             checkboxGroupInput("predict_combined_gspace_layers",
                                label = tags$span("G-space shows:",
                                                  class = "text-widget-title"),
                                choices = choices,
                                selected = choices[1],
                                inline = TRUE))
    ),
    tags$small("Up to four layers per side.", class = "text-muted-small")
  )
})

output$predict_combined_plot <- renderPlot({

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

  vars <- predict_plot_vars()
  req(vars, length(vars) >= 2)
  req(input$predict_plot_combined_x, input$predict_plot_combined_y)

  s <- predict_collect_plot_settings()
  s$asp_espace <- s$asp_combined
  s$zoom_mode_espace <- s$zoom_mode_combined

  pred <- session_data$ellipsoid_prediction_list[[ell$ell_id]]

  # Before a real prediction the only thing either side can show is the
  # on-the-fly suitable area
  if(is.null(pred) || !inherits(pred, "SpatRaster")){
    vis <- predict_raster_vis()
    e_layers <- "binary"
    g_layers <- "binary"
    g_source <- vis
  } else {
    e_layers <- head(input$predict_combined_espace_layers, 4)
    g_layers <- head(input$predict_combined_gspace_layers, 4)
    g_source <- pred
  }

  if(length(e_layers) == 0) e_layers <- "binary"
  if(length(g_layers) == 0) g_layers <- "binary"

  n_e <- length(e_layers)
  n_g <- length(g_layers)

  e_block <- plot_panel_block(seq_len(n_e))
  g_block <- plot_panel_block(n_e + seq_len(n_g))

  layout_mode <- if(!is.null(input$predict_plot_combined_layout)){
    input$predict_plot_combined_layout
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
  par(mar = c(4, 4, 2, 1))

  # E-space, one panel per selected layer
  for(lyr in e_layers){
    predict_draw_espace_panel(input$predict_plot_combined_x,
                              input$predict_plot_combined_y,
                              s,
                              layer = lyr,
                              title = if(identical(lyr, "binary")) "Suitable" else lyr)
  }

  # G-space, one panel per selected layer
  par(mar = c(3, 3, 2, 4))

  for(lyr in g_layers){

    if(identical(lyr, "binary")){
      binary <- predict_binary_layer(g_source)
      if(is.null(binary)){
        plot.new()
        next
      }
      predict_draw_gspace_panel(binary, s,
                                title = "Suitable area",
                                col = c(s$unsuitable_col, s$suitable_col),
                                legend = FALSE)
    } else if(!is.null(g_source) && lyr %in% names(g_source)){
      predict_draw_gspace_panel(g_source[[lyr]], s, title = lyr)
    } else {
      plot.new()
    }
  }
})

output$predict_combined_plot_bottom_options_ui <- renderUI({
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(
    column(width = 1,
           tags$span("Zoom:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("predict_plot_zoom_mode_combined", label = NULL,
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
           radioButtons("predict_plot_asp_combined", label = NULL,
                        choiceNames = list(
                          tags$span("Auto", class = "text-widget-inner"),
                          tags$span("Fixed", class = "text-widget-inner")
                        ),
                        choiceValues = c("auto", "fixed"),
                        selected = "auto",
                        inline = TRUE))
  )
})

output$predict_plot_settings_ui <- renderUI({

  ell_slot()

  has_ell <- !is.null(isolate(session_data$current_ellipsoid))
  has_raster <- !is.null(session_data$bg_raster)
  has_pred <- length(predict_pred_layers()) > 0

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

    # Background points
    fluidRow(
      column(width = 3,
             tags$span("Point shape (pch)", class = "text-widget-title"),
             selectInput("predict_plot_pch", label = NULL,
                         choices = c("Dot (.)" = ".",
                                     "Open circle" = "1",
                                     "Filled circle" = "16",
                                     "Square" = "15",
                                     "Triangle" = "17",
                                     "Cross" = "3"),
                         selected = ".")),
      column(width = 3,
             tags$span("Point size (cex)", class = "text-widget-title"),
             numericInput("predict_plot_cex", label = NULL, value = 0.3,
                          min = 0.1, max = 5, step = 0.1)),
      color_input("predict_plot_bg_col", "Background point color", "#B3B3B3"),
      column(width = 3,
             tags$span("Overlay point size", class = "text-widget-title"),
             numericInput("predict_plot_suitable_cex", label = NULL, value = 0.3,
                          min = 0.1, max = 5, step = 0.1))
    ),

    # Prediction palette
    if(has_pred) fluidRow(
      column(width = 3,
             tags$span("Continuous palette", class = "text-widget-title"),
             selectInput("predict_palette", label = NULL,
                         choices = c("Viridis" = "viridis",
                                     "Plasma" = "Plasma",
                                     "Inferno" = "Inferno",
                                     "Rocket" = "Rocket",
                                     "Mako" = "Mako",
                                     "YlGnBu" = "YlGnBu",
                                     "Heat" = "Heat"),
                         selected = "viridis")),
      column(width = 3,
             tags$span("Direction", class = "text-widget-title"),
             checkboxInput("predict_palette_rev", "Reverse palette",
                           value = FALSE)),
      column(width = 3,
             tags$span("Overlay shape (pch)", class = "text-widget-title"),
             selectInput("predict_plot_suitable_pch", label = NULL,
                         choices = c("Filled circle" = "16",
                                     "Square" = "15",
                                     "Dot (.)" = "."),
                         selected = "16")),
      color_input("predict_plot_suitable_col", "Binary suitable color", "#097a21")
    ),

    # Ellipsoid and centroid
    if(has_ell) tagList(
      fluidRow(
        column(width = 3,
               checkboxInput("predict_show_ellipsoid", "Show ellipsoid",
                             value = TRUE)),
        color_input("predict_plot_ell_col", "Ellipsoid color", "#000000"),
        column(width = 3,
               tags$span("Ellipsoid line width", class = "text-widget-title"),
               numericInput("predict_plot_ell_lwd", label = NULL, value = 2,
                            min = 0.5, max = 6, step = 0.5)),
        column(width = 3,
               tags$span("Ellipsoid line type", class = "text-widget-title"),
               selectInput("predict_plot_ell_lty", label = NULL,
                           choices = c("Solid" = "1",
                                       "Dashed" = "2",
                                       "Dotted" = "3"),
                           selected = "1"))
      ),

      fluidRow(
        column(width = 3,
               checkboxInput("predict_show_centroid", "Show centroid",
                             value = TRUE)),
        column(width = 3,
               tags$span("Centroid shape (pch)", class = "text-widget-title"),
               selectInput("predict_plot_centroid_pch", label = NULL,
                           choices = c("Cross (X)" = "4",
                                       "Star" = "8",
                                       "Filled diamond" = "18",
                                       "Filled circle" = "16"),
                           selected = "8")),
        color_input("predict_plot_centroid_col", "Centroid color", "#000000"),
        column(width = 3,
               tags$span("Centroid size (cex)", class = "text-widget-title"),
               numericInput("predict_plot_centroid_cex", label = NULL, value = 1.5,
                            min = 0.5, max = 5, step = 0.5))
      )
    ),

    # Map
    if(has_raster) fluidRow(
      color_input("predict_plot_unsuitable_col", "Unsuitable color", "#D3D3D3"),
      color_input("predict_plot_map_bg_col", "Map background (NA)", "#F0F0F0")
    ),

    br(),

    fluidRow(
      column(width = 12,
             class = "btn-spaced",
             actionButton("predict_open_export_modal",
                          tagList(icon("download"), "Export figure"),
                          class = "btn-default"))
    )
  )
})

# EXPORT ------------------------------------------------------------------

observeEvent(input$predict_open_export_modal, {
  showModal(modalDialog(
    title = "Export Figure",
    size  = "m",
    fluidRow(
      column(width = 6,
             tags$span("File type", class = "text-widget-title"),
             radioButtons("predict_export_filetype", label = NULL,
                          choiceNames = list(
                            tags$span("PNG", class = "text-widget-inner"),
                            tags$span("PDF", class = "text-widget-inner"),
                            tags$span("SVG", class = "text-widget-inner")
                          ),
                          choiceValues = c("png", "pdf", "svg"),                          selected = if(!is.null(input$predict_export_filetype))
                            input$predict_export_filetype else "png",
                          inline = TRUE)),
      column(width = 6,
             tags$span("Unit", class = "text-widget-title"),
             radioButtons("predict_export_unit", label = NULL,
                          choiceNames = list(
                            tags$span("mm", class = "text-widget-inner"),
                            tags$span("inches", class = "text-widget-inner"),
                            tags$span("px", class = "text-widget-inner")
                          ),
                          choiceValues = c("mm", "in", "px"),
                          selected = if(!is.null(input$predict_export_unit))
                            input$predict_export_unit else "mm",
                          inline = TRUE))
    ),

    # cex sits outside the renderUI so it never resets when type or unit change
    fluidRow(
      column(width = 4,
             tags$span("Text size (cex)", class = "text-widget-title"),
             numericInput("predict_export_cex", label = NULL,
                          value = if(!is.null(input$predict_export_cex))
                            input$predict_export_cex else 1,
                          min = 0.5, max = 3, step = 0.1))
    ),

    uiOutput("predict_export_settings_ui"),

    br(),

    footer = tagList(
      modalButton("Cancel"),
      downloadButton("predict_confirm_export", "Export", class = "btn-save",
                     onclick = "Shiny.setInputValue('predict_export_done',
                     Math.random(), {priority: 'event'});")
    ),
    easyClose = TRUE
  ))
})

output$predict_export_settings_ui <- renderUI({

  ext <- if(!is.null(input$predict_export_filetype)) input$predict_export_filetype else "png"
  unit <- if(!is.null(input$predict_export_unit)) input$predict_export_unit else "mm"

  defaults <- switch(unit,
                     "mm" = list(val = 166, min = 50, max = 500, step = 1),
                     "in" = list(val = 6.54, min = 1, max = 20, step = 0.1),
                     "px" = list(val = 1961, min = 400, max = 6000, step = 100))

  size_row <- fluidRow(
    column(width = if(ext == "png") 4 else 6,
           tags$span(paste0("Width (", unit, ")"), class = "text-widget-title"),
           numericInput("predict_export_width_val", label = NULL,
                        value = defaults$val, min = defaults$min,
                        max = defaults$max, step = defaults$step)),
    column(width = if(ext == "png") 4 else 6,
           tags$span(paste0("Height (", unit, ")"), class = "text-widget-title"),
           numericInput("predict_export_height_val", label = NULL,
                        value = defaults$val, min = defaults$min,
                        max = defaults$max, step = defaults$step)),
    if(ext == "png"){
      column(width = 4,
             tags$span("Resolution (dpi)", class = "text-widget-title"),
             numericInput("predict_export_res", label = NULL,
                          value = 300, min = 72, max = 600, step = 50))
    }
  )

  tagList(
    p(tagList(icon("circle-info"),
              " Default is a standard full-page publication figure (166 x 166 mm)."),
      style = "font-size: 12px; color: #666; margin-bottom: 8px;"),
    size_row,
    if(ext != "png") p("PDF and SVG do not require a resolution setting.",
                       class = "text-instruction")
  )
})

output$predict_confirm_export <- downloadHandler(
  filename = function(){
    ext <- if(!is.null(input$predict_export_filetype)) input$predict_export_filetype else "png"
    tab <- if(!is.null(input$predict_plot_tabs)) input$predict_plot_tabs else "predict_espace_plot_tab"
    prefix <- switch(tab,
                     "predict_espace_plot_tab" = "predict_espace_plot",
                     "predict_gspace_plot_tab" = "predict_gspace_plot",
                     "predict_combined_plot_tab" = "predict_combined_plot",
                     "predict_plot")
    paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
  },
  content = function(file){

    ext <- if(!is.null(input$predict_export_filetype)) input$predict_export_filetype else "png"
    tab <- if(!is.null(input$predict_plot_tabs)) input$predict_plot_tabs else "predict_espace_plot_tab"

    plot_open_device(file, ext,
                     w_val = input$predict_export_width_val,
                     h_val = input$predict_export_height_val,
                     unit = input$predict_export_unit,
                     res = input$predict_export_res,
                     cex_val = input$predict_export_cex)
    on.exit(dev.off())

    s <- predict_collect_plot_settings()
    req(s)

    if(identical(s$espace_layer, "none")) s$espace_layer <- NULL
    if(!is.null(s$espace_layer) && !identical(s$espace_layer, "binary")){
      s$layer_rng <- pred_layer_range(s$espace_vals, s$espace_layer)
    }

    vars <- predict_plot_vars()
    req(vars)

    ell <- session_data$current_ellipsoid
    pred <- if(!is.null(ell)) session_data$ellipsoid_prediction_list[[ell$ell_id]] else NULL

    switch(tab,

           "predict_espace_plot_tab" = {
             state <- if(!is.null(input$predict_plot_espace_state)){
               input$predict_plot_espace_state
             } else {
               "predict_plot_pairs"
             }
             lyrs <- predict_pred_layers()

             if(identical(state, "predict_plot_pairs")){
               predict_draw_espace_pairs(vars, s)
             } else if(length(lyrs) > 0){
               predict_draw_espace_layers(input$predict_plot_2d_x,
                                          input$predict_plot_2d_y, s, lyrs)
             } else {
               par(mar = c(4, 4, 2, 1))
               predict_draw_espace_panel(input$predict_plot_2d_x,
                                         input$predict_plot_2d_y, s)
             }
           },

           "predict_gspace_plot_tab" = {
             req(inherits(pred, "SpatRaster"))
             mode <- if(!is.null(input$predict_plot_gspace_show)){
               input$predict_plot_gspace_show
             } else {
               "layers"
             }

             if(identical(mode, "binary")){
               binary <- predict_binary_layer(pred)
               req(binary)
               par(mar = c(4, 4, 2, 4))
               predict_draw_gspace_panel(binary, s,
                                         title = "Suitable area",
                                         col = c(s$unsuitable_col, s$suitable_col),
                                         legend = FALSE)
             } else {
               lyrs <- predict_pred_layers()
               n_cols <- if(length(lyrs) <= 1) 1L else 2L
               par(mfrow = c(ceiling(length(lyrs) / n_cols), n_cols),
                   mar = c(3, 3, 2, 4))
               for(lyr in lyrs) predict_draw_gspace_panel(pred[[lyr]], s, title = lyr)
             }
           },

           "predict_combined_plot_tab" = {
             req(inherits(pred, "SpatRaster"))
             req(input$predict_plot_combined_x, input$predict_plot_combined_y)

             s$asp_espace <- s$asp_combined
             s$zoom_mode_espace <- s$zoom_mode_combined

             e_layers <- head(input$predict_combined_espace_layers, 4)
             g_layers <- head(input$predict_combined_gspace_layers, 4)
             if(length(e_layers) == 0) e_layers <- "binary"
             if(length(g_layers) == 0) g_layers <- "binary"

             e_block <- plot_panel_block(seq_len(length(e_layers)))
             g_block <- plot_panel_block(length(e_layers) + seq_len(length(g_layers)))

             lay <- if(!is.null(input$predict_plot_combined_layout)){
               input$predict_plot_combined_layout
             } else {
               "row"
             }

             layout(if(identical(lay, "row")) cbind(e_block, g_block)
                    else rbind(e_block, g_block))

             par(mar = c(4, 4, 2, 1))
             for(lyr in e_layers){
               predict_draw_espace_panel(input$predict_plot_combined_x,
                                         input$predict_plot_combined_y, s,
                                         layer = lyr,
                                         title = if(identical(lyr, "binary")) "Suitable" else lyr)
             }

             par(mar = c(3, 3, 2, 4))
             for(lyr in g_layers){
               if(identical(lyr, "binary")){
                 binary <- predict_binary_layer(pred)
                 if(is.null(binary)){ plot.new(); next }
                 predict_draw_gspace_panel(binary, s, title = "Suitable area",
                                           col = c(s$unsuitable_col, s$suitable_col),
                                           legend = FALSE)
               } else if(lyr %in% names(pred)){
                 predict_draw_gspace_panel(pred[[lyr]], s, title = lyr)
               } else {
                 plot.new()
               }
             }

             layout(1)
           }
    )
  }
)

observeEvent(input$predict_export_done, {

  removeModal()

  showNotification("Figure exported successfully.",
                   type = "message", duration = 4)
})
