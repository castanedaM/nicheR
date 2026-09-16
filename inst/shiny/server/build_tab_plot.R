# Title: Plot logic
# Description: Handle e-space, g-space, and combined plots
# Date last updated: 08/04/2026

# Functions -----------------------------------------------------------------

compute_lims <- function(v1, v2, s){

  ell_valid <- s$has_ell &&
    !is.null(s$ell$cov_matrix) &&
    all(is.finite(s$ell$cov_matrix))

  if(s$zoom_mode_espace == "ellipsoid" && ell_valid){
    idx <- match(c(v1, v2), s$ell$var_names)
    if(any(is.na(idx))) return(list(xlim = c(0, 1), ylim = c(0, 1), asp = NA))
    ell_pts <- ellipsoid_boundary_2d(s$ell, n_segments = 100, dim = idx)
    pad_x <- diff(range(ell_pts[, 1])) * 0.1
    pad_y <- diff(range(ell_pts[, 2])) * 0.1
    xlim <- range(ell_pts[, 1]) + c(-pad_x, pad_x)
    ylim <- range(ell_pts[, 2]) + c(-pad_y, pad_y)
    return(list(xlim = xlim, ylim = ylim, asp = diff(ylim) / diff(xlim)))
  }

  bg <- session_data$bg_df
  ranges <- s$show_lines$ranges

  # All data frames keyed by v1/v2 so rbind works correctly
  pts_xy <- if(!is.null(bg)){
    bg[, c(v1, v2)]
  } else if(!is.null(ranges)){
    df <- data.frame(x = c(ranges$mins[[v1]], ranges$maxs[[v1]]),
                     y = c(ranges$mins[[v2]], ranges$maxs[[v2]]))
    names(df) <- c(v1, v2)
    df
  } else {
    df <- data.frame(x = c(0, 1), y = c(0, 1))
    names(df) <- c(v1, v2)
    df
  }

  # Range line positions keyed by v1/v2 to match pts_xy
  range_pts <- if(!is.null(ranges) &&
                  !is.null(ranges$mins[[v1]]) &&
                  !is.null(ranges$maxs[[v1]])){
    df <- data.frame(x = c(ranges$mins[[v1]], ranges$maxs[[v1]]),
                     y = c(ranges$mins[[v2]], ranges$maxs[[v2]]))
    names(df) <- c(v1, v2)
    df
  } else {
    NULL
  }

  if(ell_valid){
    idx <- match(c(v1, v2), s$ell$var_names)
    if(!any(is.na(idx))){
      ell_pts <- ellipsoid_boundary_2d(s$ell, n_segments = 100, dim = idx)
      all_pts <- if(!is.null(range_pts)) rbind(pts_xy, range_pts) else pts_xy
      lims <- safe_lims(all_pts, ell_pts)
      return(c(lims, list(asp = NA)))
    }
  }

  # No ellipsoid: expand bg limits to include range lines
  all_x <- c(pts_xy[, 1], if(!is.null(range_pts)) range_pts[[v1]])
  all_y <- c(pts_xy[, 2], if(!is.null(range_pts)) range_pts[[v2]])

  list(xlim = range(all_x, na.rm = TRUE),
       ylim = range(all_y, na.rm = TRUE),
       asp = NA)
}

build_draw_espace_panel <- function(v1, v2, s){

  lims <- compute_lims(v1, v2, s)
  bg <- session_data$bg_df

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
         main = paste(v1, "vs.", v2))
  } else {
    plot(NA, NA,
         xlim = lims$xlim,
         ylim = lims$ylim,
         asp = if(!is.null(s$asp_espace) && s$asp_espace == "fixed") lims$asp else NA,
         xlab = v1,
         ylab = v2,
         main = paste(v1, "vs.", v2))
  }

  # Suitable area overlay using predict() on bg_df
  if(s$show_suitable_espace && !is.null(bg) && s$has_ell){
    pred_df <- tryCatch(
      predict(s$ell,
              newdata = bg,
              include_suitability = TRUE,
              include_mahalanobis = FALSE,
              suitability_truncated = TRUE,
              verbose = FALSE),
      error = function(e) NULL
    )

    if(!is.null(pred_df)){
      suitable <- pred_df[!is.na(pred_df$suitability_trunc) &
                            pred_df$suitability_trunc > 0, ]
      if(nrow(suitable) > 0){
        points(suitable[[v1]], suitable[[v2]],
               col = s$suitable_col,
               pch = s$suitable_pch,
               cex = s$suitable_cex)
      }
    }
  }

  if(isTRUE(s$show_lines$active)){
    ranges <- s$show_lines$ranges
    abline(v = ranges$mins[[v1]], col = s$xline_col, lwd = s$line_lwd, lty = 2)
    abline(v = ranges$maxs[[v1]], col = s$xline_col, lwd = s$line_lwd, lty = 2)
    abline(h = ranges$mins[[v2]], col = s$yline_col, lwd = s$line_lwd, lty = 2)
    abline(h = ranges$maxs[[v2]], col = s$yline_col, lwd = s$line_lwd, lty = 2)
  }

  if(s$show_ell){
    idx <- match(c(v1, v2), s$ell$var_names)
    add_ellipsoid(s$ell, dim = idx,
                  col_ell = s$ell_col,
                  lwd = s$ell_lwd,
                  lty = s$ell_lty)
  }

  if(s$show_centroid && !is.null(s$ell)){
    idx <- match(c(v1, v2), s$ell$var_names)
    c_pos <- if(!is.null(s$centroid_preview_val)){
      s$centroid_preview_val
    } else {
      s$ell$centroid
    }
    points(c_pos[idx[1]], c_pos[idx[2]],
           pch = s$centroid_pch, col = s$centroid_col, cex = s$centroid_cex)
  }
}

build_draw_espace_pairs <- function(vars, s){
  pairs <- t(combn(seq_along(vars), 2))
  n_pairs <- nrow(pairs)
  n_cols <- ceiling(sqrt(n_pairs))
  n_rows <- ceiling(n_pairs / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 2, 1))
  for(i in seq_len(n_pairs)) build_draw_espace_panel(vars[pairs[i, 1]], vars[pairs[i, 2]], s)
}

build_draw_gspace_panel <- function(rast, s, title = NULL, col = NULL){
  map_bg_col <- s$map_bg_col
  ttl  <- if(!is.null(title)) title else names(rast)[1]

  if(!is.null(col)){
    terra::plot(rast,
                col = col,
                legend = FALSE,
                main = ttl,
                colNA = map_bg_col,
                axes = TRUE,
                xlab = "Longitude",
                ylab = "Latitude")
  } else {
    terra::plot(rast,
                main = ttl,
                colNA = map_bg_col,
                axes = TRUE,
                xlab = "Longitude",
                ylab = "Latitude")
  }
}

build_draw_gspace_all <- function(vars, s){
  n_cols <- 2
  n_rows <- ceiling(length(vars) / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 2, 1))
  for(i in seq_len(length(vars)))
    build_draw_gspace_panel(session_data$bg_raster[[vars[i]]], s)
}

# Called at the top of every draw function.
# Returns a plain list so drawing functions are pure and testable.
collect_plot_settings <- function(){

  # Range lines. Follow the live inputs when there is no ellipsoid yet, or
  # when the user has edited the ranges of an existing one, so edits preview
  # before Rebuild commits them. Otherwise show the ranges the current
  # ellipsoid was built from, which do not move with the centroid.
  cur_ell <- session_data$current_ellipsoid

  live_ranges <- tryCatch(
    withCallingHandlers(
      {
        rp <- range_preview()
        if(is.null(rp)) NULL else list(mins = rp$mins, maxs = rp$maxs)
      },
      shiny.silent.error = function(e) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )

  ranges <- if(is.null(cur_ell)){
    live_ranges
  } else if(!is.null(live_ranges)){
    live_ranges
  } else if(!is.null(cur_ell$range_inputs)){
    list(mins = cur_ell$range_inputs$min,
         maxs = cur_ell$range_inputs$max)
  } else {
    NULL
  }

  has_ell <- !is.null(session_data$current_ellipsoid)

  list(
    has_ell = has_ell,
    ell = if(has_ell) session_data$current_ellipsoid else NULL,

    show_ell = has_ell && get_input("build_show_ellipsoid",TRUE),
    show_centroid = has_ell && get_input("build_show_centroid", TRUE),
    show_suitable_espace = has_ell && get_input("build_show_suitable_espace", TRUE),
    show_suitable_gspace = has_ell && get_input("build_show_suitable_gspace", TRUE),

    pch_val = {
      v <- get_input("build_plot_pch", ".")
      if(v == ".") "." else as.numeric(v)
    },
    cex_val = get_input("build_plot_cex", 0.3),
    bg_col = get_input("build_plot_bg_col", "#B3B3B3"),

    suitable_pch = as.numeric(get_input("build_plot_suitable_pch", "16")),
    suitable_cex = get_input("build_plot_suitable_cex", 0.3),
    suitable_col = get_input("build_plot_suitable_col", "#097a21"),
    unsuitable_col = get_input("build_plot_unsuitable_col", "#D3D3D3"),
    map_bg_col = get_input("build_plot_map_bg_col", "#F0F0F0"),

    show_lines = {
      show <- get_input("build_show_range_lines", TRUE)
      list(active = show && !is.null(ranges), ranges = ranges)
    },

    xline_col = get_input("build_plot_xline_col","#E10000"),
    yline_col = get_input("build_plot_yline_col","#0004D5"),
    line_lwd = get_input("build_plot_line_lwd", 2),

    ell_col = get_input("build_plot_ell_col", "#000000"),
    ell_lwd = get_input("build_plot_ell_lwd", 2),
    ell_lty = as.numeric(get_input("build_plot_ell_lty", "1")),

    centroid_pch = as.numeric(get_input("build_plot_centroid_pch", "8")),
    centroid_col = get_input("build_plot_centroid_col", "#000000"),
    centroid_cex = get_input("build_plot_centroid_cex", 1.5),

    zoom_mode_espace = get_input("build_plot_zoom_mode_espace","auto"),
    zoom_mode_combined = get_input("build_plot_zoom_mode_combined","auto"),

    asp_espace = get_input("build_plot_asp_espace", "auto"),
    asp_combined = get_input("build_plot_asp_combined", "auto"),
    export_cex = get_input("build_export_cex", 1),

    centroid_preview_val = tryCatch(
      withCallingHandlers(
        centroid_preview(),
        shiny.silent.error = function(e) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )
  )
}


# Reactives ---------------------------------------------------------------

# Drives the conditionalPanel that shows the range-line toggle
output$build_ellipsoid_exists <- reactive({
  !is.null(session_data$current_ellipsoid)
})
outputOptions(output, "build_ellipsoid_exists", suspendWhenHidden = FALSE)


# Prediction raster for G-space and Combined tabs only.
# E-space uses predict() on bg_df directly inside build_draw_espace_panel().
pred_raster_vis <- reactive({

  req(session_data$current_ellipsoid)
  req(session_data$bg_raster)

  ell <- session_data$current_ellipsoid
  vars <- ell$var_names

  tryCatch(
    predict(ell,
            newdata = session_data$bg_raster[[vars]],
            include_suitability = TRUE,
            include_mahalanobis = FALSE,
            suitability_truncated = TRUE,
            verbose = FALSE),
    error = function(e) NULL
  )

})


# Variables to plot, reactive to live selections before confirm.
# Before vars are confirmed: first 6 non-spatial columns from bg data.
# After confirm: session_data$vars.
plot_vars <- reactive({

  if(!is.null(session_data$vars)){
    return(session_data$vars)
  }

  if(identical(session_data$input_mode, "virtual")) return(NULL)

  all_vars <- available_vars()
  if(is.null(all_vars)) return(NULL)

  n_slots <- min(length(all_vars), MAX_DIMS)

  live_active <- vapply(seq_len(n_slots), function(i){
    val <- input[[paste0("build_var_active_", i)]]
    if(is.null(val)) TRUE else isTRUE(val)
  }, logical(1))

  live_select <- vapply(seq_len(n_slots), function(i){
    val <- input[[paste0("build_var_select_", i)]]
    if(is.null(val)) all_vars[i] else val
  }, character(1))

  selected <- unique(live_select[live_active])

  if(length(selected) == 0) return(head(all_vars, MAX_DIMS))

  selected
})


# Selector observers ------------------------------------------------------

# Selects the x and y based on the available variables, prevent form selection a
# 1:1
observeEvent({
  input$build_plot_espace_state
  input$build_plot_2d_x
  input$build_plot_2d_y
  session_data$vars
}, {
  vars <- plot_vars()
  req(vars)
  req(input$build_plot_espace_state == "build_plot_2d")
  update_axis_selectors("build_plot_2d_x", "build_plot_2d_y", vars)
}, ignoreInit = FALSE)

observeEvent({
  input$build_plot_combined_x
  input$build_plot_combined_y
  session_data$vars
}, {
  vars <- plot_vars()
  req(vars)
  update_axis_selectors("build_plot_combined_x", "build_plot_combined_y", vars)
}, ignoreInit = FALSE)



# Outputs -----------------------------------------------------------------

# E-Space
output$build_espace_plot_top_options_ui <- renderUI({

  vars <- session_data$vars
  req(vars)

  fluidRow(
    column(width = 12,
           column(width = 1, tags$span("Layout:", class = "text-widget-title")),
           column(width = 3,
                  radioButtons("build_plot_espace_state",
                               label = NULL,
                               choiceNames = list(
                                 tags$span("All pairs", class = "text-widget-inner"),
                                 tags$span("2D", class = "text-widget-inner")),
                               choiceValues = c("build_plot_pairs", "build_plot_2d"),
                               selected = "build_plot_pairs",
                               inline = TRUE)),

           conditionalPanel(
             "input.build_plot_espace_state == 'build_plot_2d'",
             column(width = 2, tags$span("Variables:", class = "text-widget-title")),
             column(width = 3,
                    selectInput("build_plot_2d_x",
                                label = NULL,
                                choices = character(0))),
             column(width = 3,
                    selectInput("build_plot_2d_y",
                                label = NULL,
                                choices = character(0)))
           )
    )
  )

})

output$build_espace_plot_bottom_options_ui <- renderUI({
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(
    column(width = 1,
           tags$span("Zoom:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("build_plot_zoom_mode_espace", label = NULL,
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
           radioButtons("build_plot_asp_espace", label = NULL,
                        choiceNames = list(
                          tags$span("Auto", class = "text-widget-inner"),
                          tags$span("Fixed", class = "text-widget-inner")
                        ),
                        choiceValues = c("auto", "fixed"),
                        selected = "auto",
                        inline = TRUE))
  )
})

output$build_espace_plot <- renderPlot({

  vars <- plot_vars()
  req(vars)

  s <- collect_plot_settings()

  req(s)

  state <- if(!is.null(input$build_plot_espace_state)) input$build_plot_espace_state else "build_plot_pairs"

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  switch(state,
         "build_plot_pairs" = build_draw_espace_pairs(vars, s),
         "build_plot_2d"= {
           req(input$build_plot_2d_x, input$build_plot_2d_y)
           par(mar = c(4, 4, 2, 1))
           build_draw_espace_panel(input$build_plot_2d_x, input$build_plot_2d_y, s)
         }
  )

})

# G-space
output$build_gspace_plot_top_options_ui <- renderUI({
  req(session_data$bg_raster)

  ell_slot()

  vars <- plot_vars()
  req(vars)

  has_ell <- !is.null(isolate(session_data$current_ellipsoid))

  show_choices_names <- if(has_ell){
    list(tags$span("Within ranges", class = "text-widget-inner"),
         tags$span("Suitable area", class = "text-widget-inner"))
  } else {
    list(tags$span("Within ranges", class = "text-widget-inner"))
  }
  show_choices_values <- if(has_ell) c("range", "suitable") else "range"


  fluidRow(
    column(width = 1, tags$span("Show:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("build_plot_gspace_show",
                        label = NULL,
                        choiceNames = show_choices_names,
                        choiceValues = show_choices_values,
                        selected = if(has_ell) "suitable" else "range",
                        inline = FALSE)),

    # Layer controls only apply to the range view, since the suitable area
    # is a single map for all variables
    conditionalPanel(
      "input.build_plot_gspace_show == 'range'",
      column(width = 1, tags$span("Layer:", class = "text-widget-title")),
      column(width = 1,
             radioButtons("build_plot_gspace_state",
                          label = NULL,
                          choiceNames = list(
                            tags$span("All", class = "text-widget-inner"),
                            tags$span("One", class = "text-widget-inner")),
                          choiceValues = c("build_plot_all", "build_plot_one"),
                          selected = "build_plot_all",
                          inline = FALSE))
    ),

    conditionalPanel(
      "input.build_plot_gspace_show == 'range' && input.build_plot_gspace_state == 'build_plot_one'",
      column(width = 2, tags$span("Variable:", class = "text-widget-title")),
      column(width = 4,
             selectInput("build_plot_gspace_lyr",
                         label = NULL,
                         choices = vars))
    )
  )
})

output$build_gspace_plot <- renderPlot({

  vars <- plot_vars()
  req(vars)

  is_virtual <- identical(session_data$input_mode, "virtual_mode")
  if(is_virtual || is.null(session_data$bg_raster)){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "G-space")
    text(0.5, 0.5, "Virtual mode on or no raster provided.\nG-space unavailable.",
         cex = 1, col = "grey50")
    return(invisible(NULL))
  }

  s <- collect_plot_settings()
  rast <- session_data$bg_raster
  lyr <- if(!is.null(input$build_plot_gspace_lyr)) input$build_plot_gspace_lyr else vars[1]
  state <- if(!is.null(input$build_plot_gspace_state)) input$build_plot_gspace_state else "build_plot_all"

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  has_ell <- s$has_ell

  # Before an ellipsoid exists, show which cells fall inside the current
  # range inputs. This updates live as the user edits the ranges.
  binarize_by_range <- function(v){

    rng <- s$show_lines$ranges

    if(is.null(rng) || is.null(rng$mins[[v]]) || is.null(rng$maxs[[v]])){
      return(rast[[v]])
    }

    binary <- terra::classify(rast[[v]],
                              rcl = matrix(c(-Inf, rng$mins[[v]], NA,
                                             rng$mins[[v]], rng$maxs[[v]], 1,
                                             rng$maxs[[v]], Inf, NA),
                                           ncol = 3, byrow = TRUE),
                              include.lowest = TRUE, right = FALSE)
    names(binary) <- v
    binary
  }


  show_mode <- if(!is.null(input$build_plot_gspace_show)){
    input$build_plot_gspace_show
  } else if(s$has_ell){
    "suitable"
  } else {
    "range"
  }

  if(identical(show_mode, "suitable") && s$has_ell){

    pred <- tryCatch(pred_raster_vis(), error = function(e) NULL)

    if(is.null(pred) || !("suitability_trunc" %in% names(pred))){
      plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
           xlab = "", ylab = "", main = "G-space")
      text(0.5, 0.5, "Prediction unavailable for this ellipsoid.",
           cex = 1, col = "grey50")
      return(invisible(NULL))
    }

    binary <- pred[["suitability_trunc"]] > 0
    names(binary) <- "suitable"

    par(mar = c(4, 4, 2, 4))
    build_draw_gspace_panel(binary, s,
                            title = "Suitable area",
                            col = c(s$unsuitable_col, s$suitable_col))

  } else {

    within_col <- c("#D3D3D3", "#E07B39")
    outside_col <- s$map_bg_col
    has_range <- !is.null(s$show_lines$ranges)

    if(state == "build_plot_all"){
      n_cols <- 2L
      n_rows <- ceiling(length(vars) / n_cols)
      par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 3))
      for(v in vars){
        if(has_range){
          build_draw_gspace_panel(binarize_by_range(v), s,
                                  title = paste0(v, " (within ranges)"),
                                  col = c(outside_col, within_col))
        } else {
          build_draw_gspace_panel(rast[[v]], s, title = v)
        }
      }
      if(length(vars) %% 2 != 0) plot.new()
    } else {
      par(mar = c(4, 4, 2, 4))
      if(has_range){
        build_draw_gspace_panel(binarize_by_range(lyr), s,
                                title = paste0(lyr, " (within ranges)"),
                                col = c(outside_col, within_col))
      } else {
        build_draw_gspace_panel(rast[[lyr]], s, title = lyr)
      }
    }

  }
})

output$build_gspace_plot_bottom_options_ui <- renderUI({
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(class = "ell-row",
    column(width = 12,
           tags$span(ell$ell_name, class = "text-center",
                     style = "font-size: 12px; color: #888; font-weight: 400;"))
  )

})


# Combined
output$build_combined_plot <- renderPlot({

  is_virtual <- identical(session_data$input_mode, "virtual_mode")
  if(is_virtual || is.null(session_data$bg_raster)){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "Combined")
    text(0.5, 0.5, "Virtual mode on or raster not provided.\nG-space unavailable.",
         cex = 1, col = "grey50")
    return(invisible(NULL))
  }

  req(input$build_plot_combined_x, input$build_plot_combined_y)

  vars <- plot_vars()
  req(vars)

  s <- collect_plot_settings()
  s$asp_espace <- s$asp_combined
  s$zoom_mode_espace <- s$zoom_mode_combined

  layout <- if(!is.null(input$build_plot_combined_layout)) input$build_plot_combined_layout else "col"
  mfrow  <- if(layout == "col") c(2, 1) else c(1, 2)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = mfrow, mar = c(4, 4, 2, 1))

  build_draw_espace_panel(input$build_plot_combined_x, input$build_plot_combined_y, s)

  # G-space: binary suitability if ellipsoid exists, else first raster layer
  gspace_rast <- if(s$show_suitable_gspace && s$has_ell){
    pred <- tryCatch(pred_raster_vis(), error = function(e) NULL)
    if(!is.null(pred)) pred[["suitability_trunc"]] > 0 else session_data$bg_raster[[1]]
  } else {
    session_data$bg_raster[[1]]
  }

  build_draw_gspace_panel(gspace_rast, s,
                          title = "G-space",
                          col = if(s$show_suitable_gspace && s$has_ell)
                            c(s$unsuitable_col, s$suitable_col) else NULL)
})

output$build_combined_plot_top_options_ui <- renderUI({

  vars <- plot_vars()
  req(vars)

  fluidRow(
    column(width = 1, tags$span("Layout:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("build_plot_combined_layout",
                        label = NULL,
                        choiceNames = list(
                          tags$span("Stacked", class = "text-widget-inner"),
                          tags$span("Side by side", class = "text-widget-inner")),
                        choiceValues = c("col", "row"),
                        selected = "col",
                        inline = FALSE)),
    column(width = 2, tags$span("Variables:", class = "text-widget-title")),
    column(width = 3,
           selectInput("build_plot_combined_x", label = NULL, choices = character(0))),
    column(width = 3,
           selectInput("build_plot_combined_y", label = NULL, choices = character(0)))
  )
})

output$build_combined_plot_bottom_options_ui <- renderUI({
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(
    column(width = 1,
           tags$span("Zoom:", class = "text-widget-title")),
    column(width = 3,
           radioButtons("build_plot_zoom_mode_combined", label = NULL,
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
           radioButtons("build_plot_asp_combined", label = NULL,
                        choiceNames = list(
                          tags$span("Auto", class = "text-widget-inner"),
                          tags$span("Fixed", class = "text-widget-inner")
                        ),
                        choiceValues = c("auto", "fixed"),
                        selected = "auto",
                        inline = TRUE))
  )
})

# Ellipsoid library, this shows all version of the base elliposid created
output$build_ellipsoid_info_ui <- renderUI({

  # Live on purpose: this is a running summary and should track every
  # slider move, unlike the option panels which must not be torn down.
  ell <- session_data$current_ellipsoid
  req(ell)

  ref_ell <- ell_reset_target(ell)
  vars <- ell$var_names
  n_vars <- length(vars)

  # Volume change relative to base
  vol_current <- ell$volume
  vol_base <- if(!is.null(ref_ell) && !is.null(ref_ell$volume)){
    ref_ell$volume
  } else {
    vol_current
  }

  vol_pct <- if(!is.null(ref_ell) && vol_base > 0){
    round((vol_current - vol_base) / vol_base * 100, 1)
  } else {
    0
  }

  vol_icon<- if(vol_pct > 0) icon("arrow-trend-up") else if(vol_pct < 0) icon("arrow-trend-down") else icon("minus")
  vol_color <- if(vol_pct > 0) "#097a21" else if(vol_pct < 0) "#e74c3c" else "#888"

  # Covariance pairs that differ from zero
  pairs <- t(combn(vars, 2))
  pair_names <- apply(pairs, 1, function(p) paste(p, collapse = " / "))
  cov_vals <- apply(pairs, 1, function(p){
    round(ell$cov_matrix[p[1], p[2]], 4)
  })

  nonzero <- cov_vals != 0

  # Covariance table rows
  cov_rows <- if(any(nonzero)){
    lapply(which(nonzero), function(i){
      val <- cov_vals[i]
      color <- if(val > 0) "#097a21" else "#e74c3c"
      icn <- if(val > 0) icon("arrow-up") else icon("arrow-down")
      tags$tr(
        tags$td(pair_names[i],
                style = "font-size: 12px; color: #666; padding: 3px 6px;"),
        tags$td(style = paste0("color:", color, "; padding: 3px 6px; font-size: 12px;"),
                icn, " ", format(val, nsmall = 4))
      )
    })
  } else {
    list(tags$tr(tags$td(
      colspan = "2",
      style = "font-size: 12px; color: #aaa; padding: 3px 6px;",
      "All covariances at zero"
    )))
  }

  # Centroid rows
  centroid_rows <- lapply(vars, function(v){

    cur_val<- round(ell$centroid[v], 3)
    ref_val <- if(!is.null(ref_ell)) round(ref_ell$centroid[v], 3) else cur_val
    delta <- round(cur_val - ref_val, 3)
    color <- if(delta > 0) "#097a21" else if(delta < 0) "#e74c3c" else "#aaa"
    delta_str <- if(delta > 0) paste0("+", delta) else if(delta < 0) as.character(delta) else "no change"

    tags$tr(
      tags$td(v,
              style = "font-size: 12px; color: #666; padding: 3px 6px;"),
      tags$td(cur_val,
              style = "font-size: 12px; color: #555; padding: 3px 6px;"),
      tags$td(delta_str,
              style = paste0("color:", color, "; padding: 3px 6px; font-size: 12px;"))
    )
  })

  box(title = tagList(
    tags$span("Ellipsoid summary", class = "text-section-header"),
    tags$span(paste0(" — ", ell$ell_name),
              style = "font-size: 12px; color: #888; font-weight: 400; margin-left: 4px;")),
    width = 12,
    collapsible = TRUE,
    collapsed = FALSE,

    fluidRow(
      column(width = 4,
             tags$div(
               tags$span("Dimensions", class = "text-widget-title"),
               tags$p(paste0(n_vars, "D (", paste(vars, collapse = ", "), ")"),
                      style = "font-size: 12px; color: #555; margin: 2px 0 12px;"),

               tags$span("Confidence level", class = "text-widget-title"),
               tags$p(paste0(round(ell$cl * 100, 1), "%"),
                      style = "font-size: 12px; color: #555; margin: 2px 0 12px;"),

               tags$span("Volume", class = "text-widget-title"),
               tags$p(
                 tags$span(format(round(vol_current, 2), big.mark = ","),
                           style = "font-size: 12px; color: #555;"),
                 # tags$span(
                 #   style = paste0("font-size: 11px; color:", vol_color,
                 #                  "; margin-left: 6px;"),
                 #   vol_icon, " ", abs(vol_pct), "% vs original"
                 # ),
                 style = "margin: 2px 0 12px;"
               )
             )
      ),

      # Covariance summary
      column(width = 4,
             tags$span("Covariance adjustments", class = "text-widget-title"),
             tags$table(
               style = "width: 100%; margin-top: 4px;",
               tags$tbody(cov_rows)
             )
      ),

      # Centroid
      column(width = 4,
             tags$span("Centroid", class = "text-widget-title"),
             tags$table(
               style = "width: 100%; margin-top: 4px;",
               tags$tbody(centroid_rows)
             )
      )
    )
  )
})

# Advanced plot settings UI
output$build_plot_settings_ui <- renderUI({

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)

  has_ell <- !is.null(ell)
  has_raster <- !is.null(session_data$bg_raster)

  box(title = tagList(
    tags$span("Advanced Plot Settings", class = "text-section-header"),
    tags$span(icon("circle-info"),
              title = instructions$plot_settings,
              class = "tooltip-icon")),
    width= 12,
    collapsible = TRUE,
    collapsed = TRUE,

    # Background points
    fluidRow(
      column(width = 4,
             tags$span("Point shape (pch)", class = "text-widget-title"),
             selectInput("build_plot_pch", label = NULL,
                         choices= c("Dot (.)"= ".",
                                    "Open circle" = "1",
                                    "Filled circle" = "16",
                                    "Square"= "15",
                                    "Triangle"= "17",
                                    "Cross" = "3"),
                         selected = ".")),
      column(width = 4,
             tags$span("Point size (cex)", class = "text-widget-title"),
             numericInput("build_plot_cex", label = NULL, value = 0.3,
                          min = 0.1, max = 5, step = 0.1)),
      column(width = 4,
             tags$span("Background point color", class = "text-widget-title"),
             tags$div(
               style = "display: flex; align-items: center; gap: 8px;",
               tags$input(type = "color",
                          value = "#B3B3B3",
                          oninput = "Shiny.setInputValue('build_plot_bg_col', this.value)",
                          style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px;
 cursor: pointer;")
             )
      )
    ),

    # Range lines
    conditionalPanel(
      condition = "input.build_range_method_choice != null && input.build_range_method_choice != '' || output.build_ellipsoid_exists",
      fluidRow(
        column(width = 3,
               checkboxInput("build_show_range_lines", "Show range lines", value = TRUE)),
        column(width = 3,
               tags$span("X-line color", class = "text-widget-title"),
               tags$div(style = "display: flex; align-items: center; gap: 8px;",
                        tags$input(type = "color",
                                   value = "#E10000",
                                   oninput = "Shiny.setInputValue('build_plot_xline_col', this.value)",
                                   style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px;
 cursor: pointer;"))
        ),
        column(width = 3,
               tags$span("Y-line color", class = "text-widget-title"),
               tags$div(
                 style = "display: flex; align-items: center; gap: 8px;",
                 tags$input(type = "color",
                            value = "#0004D5",
                            oninput = "Shiny.setInputValue('build_plot_yline_col', this.value)",
                            style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px;
 cursor: pointer;")
               )
        ),
        column(width = 3,
               tags$span("Line width", class = "text-widget-title"),
               numericInput("build_plot_line_lwd", label = NULL, value = 2,
                            min = 0.5, max = 6, step = 0.5))
      )
    ),

    # Ellipsoid, centroid, suitable area (only when ellipsoid built)
    if(has_ell) tagList(
      fluidRow(
        column(width = 3,
               checkboxInput("build_show_ellipsoid", "Show ellipsoid", value = TRUE)),
        column(width = 3,
               tags$span("Ellipsoid color", class = "text-widget-title"),
               tags$div(
                 style = "display: flex; align-items: center; gap: 8px;",
                 tags$input(type = "color",
                            value = "#000000",
                            oninput = "Shiny.setInputValue('build_plot_ell_col', this.value)",
                            style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px;
 cursor: pointer;")
               )
        ),
        column(width = 3,
               tags$span("Ellipsoid line width", class = "text-widget-title"),
               numericInput("build_plot_ell_lwd", label = NULL, value = 2,
                            min = 0.5, max = 6, step = 0.5)),
        column(width = 3,
               tags$span("Ellipsoid line type", class = "text-widget-title"),
               selectInput("build_plot_ell_lty", label = NULL,
                           choices= c("Solid"= "1",
                                      "Dashed" = "2",
                                      "Dotted" = "3"),
                           selected = "1"))
      ),

      fluidRow(
        column(width = 3,
               checkboxInput("build_show_centroid", "Show centroid", value = TRUE)),
        column(width = 3,
               tags$span("Centroid shape (pch)", class = "text-widget-title"),
               selectInput("build_plot_centroid_pch", label = NULL,
                           choices = c("Cross (X)"= "4",
                                       "Star" = "8",
                                       "Filled diamond" = "18",
                                       "Filled circle"= "16"),
                           selected = "8")),
        column(width = 3,
               tags$span("Centroid color", class = "text-widget-title"),
               tags$div(
                 style = "display: flex; align-items: center; gap: 8px;",
                 tags$input(type = "color",
                            value = "#000000",
                            oninput = "Shiny.setInputValue('build_plot_centroid_col', this.value)",
                            style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px;
 cursor: pointer;")
               )
        ),
        column(width = 3,
               tags$span("Centroid size (cex)", class = "text-widget-title"),
               numericInput("build_plot_centroid_cex", label = NULL, value = 1.5,
                            min = 0.5, max = 5, step = 0.5))
      ),

      fluidRow(
        column(width = 3,
               checkboxInput("build_show_suitable_espace",
                             "Show suitable area (E-space)", value = TRUE)),
        column(width = 3,
               checkboxInput("build_show_suitable_gspace",
                             "Show suitable area (G-space)", value = TRUE)),
        column(width = 3,
               tags$span("Suitable area color", class = "text-widget-title"),
               tags$div(
                 style = "display: flex; align-items: center; gap: 8px;",
                 tags$input(type = "color",
                            value = "#097a21",
                            oninput = "Shiny.setInputValue('build_plot_suitable_col', this.value)",
                            style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px;
 cursor: pointer;")
               )

        ),
        column(width = 3,
               tags$span("Unsuitable area color", class = "text-widget-title"),
               tags$div(
                 style = "display: flex; align-items: center; gap: 8px;",
                 tags$input(type = "color",
                            value = "#D3D3D3",
                            oninput = "Shiny.setInputValue('build_plot_unsuitable_col', this.value)",
                            style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px;
 cursor: pointer;")
               )
        )
      ),

      fluidRow(
        column(width = 3,
               tags$span("Suitable point shape (pch)", class = "text-widget-title"),
               selectInput("build_plot_suitable_pch", label = NULL,
                           choices= c("Dot (.)"= ".",
                                      "Open circle" = "1",
                                      "Filled circle" = "16",
                                      "Square"= "15",
                                      "Triangle"= "17",
                                      "Cross" = "3"),
                           selected = "16")),
        column(width = 3,
               tags$span("Suitable point size (cex)", class = "text-widget-title"),
               numericInput("build_plot_suitable_cex", label = NULL, value = 0.3,
                            min = 0.1, max = 5, step = 0.1))
      )
    ),

    # Map background color
    if(has_raster) fluidRow(
      column(width = 4,
             tags$span("Map background color", class = "text-widget-title"),
             tags$div(
               style = "display: flex; align-items: center; gap: 8px;",
               tags$input(type = "color",
                          value = "#F0F0F0",
                          oninput = "Shiny.setInputValue('build_plot_map_bg_col', this.value)",
                          style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px;
 cursor: pointer;")
             )
      )
    ),

    br(),

    # Export button and settings sit below the tabBox, outside all panels
    fluidRow(
      column(width = 4,
             actionButton("build_open_export_modal",
                          tagList(icon("download"), "Export Figure"),
                          class = "btn-default"))
    )
  )
})

# Export Logic ------------------------------------------------------------

plot_open_device <- function(file, ext, w_val, h_val, unit, res, cex_val){

  if(is.null(unit)) unit <- "mm"
  if(is.null(w_val)) w_val <- 166
  if(is.null(h_val)) h_val <- 166
  if(is.null(res)) res <- 300
  if(is.null(cex_val)) cex_val <- 1

  to_inches <- function(val){
    switch(unit, "mm" = val / 25.4, "in" = val, "px" = val / res)
  }
  to_px <- function(val){
    switch(unit,
           "mm" = round(val / 25.4 * res),
           "in" = round(val * res),
           "px" = round(val))
  }

  if(ext == "png"){
    png(file, width = to_px(w_val), height = to_px(h_val), res = res)
  } else if(ext == "pdf"){
    pdf(file, width = to_inches(w_val), height = to_inches(h_val))
  } else if(ext == "svg"){
    svg(file, width = to_inches(w_val), height = to_inches(h_val))
  }

  # Set every cex parameter so individual plot calls do not override them
  par(cex = cex_val,
      cex.axis = cex_val,
      cex.lab = cex_val,
      cex.main = cex_val * 1.1,
      cex.sub = cex_val * 0.9)
}

output$build_export_settings_ui <- renderUI({

  ext  <- if(!is.null(input$build_export_filetype)) input$build_export_filetype else "png"
  unit <- if(!is.null(input$build_export_unit))     input$build_export_unit     else "mm"

  defaults <- switch(unit,
                     "mm" = list(val = 166,  min = 50,  max = 500,  step = 1),
                     "in" = list(val = 6.54, min = 1,   max = 20,   step = 0.1),
                     "px" = list(val = 1961, min = 400, max = 6000, step = 100)
  )

  if(ext == "png"){
    tagList(
      p(tagList(icon("circle-info"),
                " Default is a standard full-page publication figure (166 x 166 mm at 300 dpi)."),
        style = "font-size: 12px; color: #666; margin-bottom: 8px;"),
      fluidRow(
        column(width = 4,
               tags$span(paste0("Width (", unit, ")"), class = "text-widget-title"),
               numericInput("build_export_width_val", label = NULL,
                            value = defaults$val,
                            min   = defaults$min,
                            max   = defaults$max,
                            step  = defaults$step)),
        column(width = 4,
               tags$span(paste0("Height (", unit, ")"), class = "text-widget-title"),
               numericInput("build_export_height_val", label = NULL,
                            value = defaults$val,
                            min   = defaults$min,
                            max   = defaults$max,
                            step  = defaults$step)),
        column(width = 4,
               tags$span("Resolution (dpi)", class = "text-widget-title"),
               numericInput("build_export_res", label = NULL,
                            value = 300, min = 72, max = 600, step = 50))
      )
    )
  } else {
    tagList(
      p(tagList(icon("circle-info"),
                " Default is a standard full-page publication figure (166 x 166 mm)."),
        style = "font-size: 12px; color: #666; margin-bottom: 8px;"),
      fluidRow(
        column(width = 6,
               tags$span("Width (mm)", class = "text-widget-title"),
               numericInput("build_export_width_val", label = NULL,
                            value = 166, min = 50, max = 500, step = 1)),
        column(width = 6,
               tags$span("Height (mm)", class = "text-widget-title"),
               numericInput("build_export_height_val", label = NULL,
                            value = 166, min = 50, max = 500, step = 1))
      ),
      p("PDF and SVG do not require a resolution setting.",
        class = "text-instruction")
    )
  }
})

observeEvent(input$build_open_export_modal, {
  showModal(modalDialog(
    title = "Export Figure",
    size  = "m",
    fluidRow(
      column(width = 6,
             tags$span("File type", class = "text-widget-title"),
             radioButtons("build_export_filetype", label = NULL,
                          choiceNames = list(
                            tags$span("PNG", class = "text-widget-inner"),
                            tags$span("PDF", class = "text-widget-inner"),
                            tags$span("SVG", class = "text-widget-inner")
                          ),
                          choiceValues = c("png", "pdf", "svg"),
                          selected = if(!is.null(input$build_export_filetype))
                            input$build_export_filetype else "png",
                          inline   = TRUE)),
      column(width = 6,
             tags$span("Unit", class = "text-widget-title"),
             radioButtons("build_export_unit", label = NULL,
                          choiceNames  = list(
                            tags$span("mm", class = "text-widget-inner"),
                            tags$span("inches", class = "text-widget-inner"),
                            tags$span("px", class = "text-widget-inner")
                          ),
                          choiceValues = c("mm", "in", "px"),
                          selected = if(!is.null(input$build_export_unit))
                            input$build_export_unit else "mm",
                          inline = TRUE))
    ),

    # cex outside renderUI so it never resets when file type or unit changes
    fluidRow(
      column(width = 4,
             tags$span("Text size (cex)", class = "text-widget-title"),
             numericInput("build_export_cex", label = NULL,
                          value = if(!is.null(input$build_export_cex)) input$build_export_cex else 1,
                          min   = 0.5, max = 3, step = 0.1)),
      column(width = 8,
             br(),
             uiOutput("build_export_cex_msg_ui"))
    ),

    uiOutput("build_export_settings_ui"),

    footer = tagList(
      modalButton("Cancel"),
      downloadButton("build_confirm_export", "Export", class = "btn-save",
                     onclick = "Shiny.setInputValue('build_export_done',
                     Math.random(), {priority: 'event'});")
    ),
    easyClose = TRUE
  ))
})

output$build_confirm_export <- downloadHandler(
  filename = function(){
    ext <- if(!is.null(input$build_export_filetype)) input$build_export_filetype else "png"
    tab <- if(!is.null(input$build_plot_tabs)) input$build_plot_tabs else "build_espace_plot_tab"
    prefix <- switch(tab,
                     "build_espace_plot_tab" = "build_espace_plot",
                     "build_gspace_plot_tab" = "build_gspace_plot",
                     "build_combined_plot_tab" = "build_combined_plot",
                     "build_espace_plot")
    paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)

  },
  content = function(file){
    ext <- if(!is.null(input$build_export_filetype)) input$build_export_filetype else "png"
    tab <- if(!is.null(input$build_plot_tabs)) input$build_plot_tabs else "build_espace_plot_tab"
    s <- collect_plot_settings()
    req(s)

    plot_open_device(file, ext,
                     w_val = input$build_export_width_val,
                     h_val = input$build_export_height_val,
                     unit = input$build_export_unit,
                     res = input$build_export_res,
                     cex_val = input$build_export_cex)
    on.exit(dev.off())

    switch(tab,
           "build_espace_plot_tab" = {
             vars <- plot_vars()
             req(vars)
             state <- if(!is.null(input$build_plot_espace_state)) input$build_plot_espace_state else "build_plot_pairs"
             switch(state,
                    "build_plot_pairs" = build_draw_espace_pairs(vars, s),
                    "build_plot_2d" = {
                      req(input$build_plot_2d_x, input$build_plot_2d_y)
                      par(mar = c(4, 4, 2, 1))
                      build_draw_espace_panel(input$build_plot_2d_x, input$build_plot_2d_y, s)
                    }
             )
           },
           "build_gspace_plot_tab" = {
             req(session_data$bg_raster)
             vars <- plot_vars()
             req(vars)

             lyr <- if(!is.null(input$build_plot_gspace_lyr)){
               input$build_plot_gspace_lyr
             } else {
               vars[1]
             }

             par(mar = c(4, 4, 2, 4))
             build_draw_gspace_panel(session_data$bg_raster[[lyr]], s, title = lyr)
           },
           "build_combined_plot_tab" = {
             req(input$build_plot_combined_x, input$build_plot_combined_y)
             layout <- if(!is.null(input$build_plot_combined_layout)) input$build_plot_combined_layout else "col"
             mfrow <- if(layout == "col") c(2, 1) else c(1, 2)
             par(mfrow = mfrow, mar = c(4, 4, 2, 1))
             build_draw_espace_panel(input$build_plot_combined_x, input$build_plot_combined_y, s)
             gspace_rast <- if(s$show_suitable_gspace && s$has_ell){
               pred <- tryCatch(pred_raster_vis(), error = function(e) NULL)
               if(!is.null(pred)) pred[["suitability_trunc"]] > 0
               else session_data$bg_raster[[1]]
             } else {
               session_data$bg_raster[[1]]
             }

             build_draw_gspace_panel(gspace_rast, s,
                                     title = "G-space",
                                     col = if(s$show_suitable_gspace && s$has_ell)
                                       c(s$unsuitable_col, s$suitable_col) else NULL)
           }
    )
  }
)


observeEvent(input$build_export_done, {

  removeModal()

  showNotification("Figure exported successfully.",
                   type = "message", duration = 4)
})



