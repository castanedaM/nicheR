# Title: Bias Tab Plot server

# Description: Plots for the bias tab. Raw and processed bias layers, the
# composite surface, and a side-by-side comparison of an unbiased
# prediction against every biased version of that layer.

# Date Last Updated: 08/05/2026


# Settings ----------------------------------------------------------------

# Collected once per draw. Defaults hold until bias_plot_settings_ui exists.
bias_collect_plot_settings <- function(){
  list(
    palette = pred_palette(get_input("bias_palette", "viridis"),
                           get_input("bias_palette_rev", FALSE)),
    map_bg_col = get_input("bias_plot_map_bg_col", "#F0F0F0")
  )
}


# BIAS INPUT LAYERS -------------------------------------------------------

# Raw layers on the left, processed on the right once bias is prepared
output$bias_layers_plot <- renderPlot({

  req(session_data$bias_raster)

  has_prepared <- !is.null(session_data$prepared_bias) &&
    !is.null(session_data$prepared_bias$processed_layers)

  s <- bias_collect_plot_settings()

  rast <- session_data$bias_raster
  n_raw <- terra::nlyr(rast)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if(has_prepared){

    processed <- session_data$prepared_bias$processed_layers
    n_proc <- terra::nlyr(processed)
    n_rows <- max(n_raw, n_proc)

    par(mfrow = c(n_rows, 2), mar = c(3, 3, 2, 4))

    for(i in seq_len(n_rows)){

      if(i <= n_raw){
        terra::plot(rast[[i]],
                    main = paste0("Raw: ", names(rast[[i]])),
                    axes = TRUE,
                    col = s$palette,
                    colNA = s$map_bg_col)
      } else {
        plot.new()
      }

      if(i <= n_proc){
        terra::plot(processed[[i]],
                    main = names(processed[[i]]),
                    axes = TRUE,
                    col = s$palette,
                    colNA = s$map_bg_col)
      } else {
        plot.new()
      }
    }

  } else {

    par(mfrow = c(n_raw, 1), mar = c(3, 3, 2, 4))

    for(i in seq_len(n_raw)){
      terra::plot(rast[[i]],
                  main = paste0("Raw: ", names(rast[[i]])),
                  axes = TRUE,
                  col = s$palette,
                  colNA = s$map_bg_col)
    }
  }
})


# COMPOSITE SURFACE -------------------------------------------------------

output$bias_composite_plot <- renderPlot({

  req(session_data$prepared_bias)
  req(session_data$prepared_bias$composite_surface)

  s <- bias_collect_plot_settings()

  composite <- session_data$prepared_bias$composite_surface

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(3, 3, 2, 4))

  terra::plot(composite,
              main = paste0("Composite: ",
                            session_data$prepared_bias$combination_formula),
              axes = TRUE,
              col = s$palette,
              colNA = s$map_bg_col)
})


# PREDICTION AND BIASED G-SPACE -------------------------------------------

output$bias_gspace_plot_layer_select_ui <- renderUI({

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  pred_result <- session_data$ellipsoid_prediction_list[[ell$ell_id]]
  req(inherits(pred_result, "SpatRaster"))

  layer_choices <- names(pred_result)

  # Keep the current layer across re-renders when the new ellipsoid also
  # has it, so switching versions does not reset the comparison
  keep <- if(!is.null(input$bias_gspace_layer) &&
             input$bias_gspace_layer %in% layer_choices){
    input$bias_gspace_layer
  } else if("suitability_trunc" %in% layer_choices){
    "suitability_trunc"
  } else {
    layer_choices[1]
  }

  fluidRow(
    column(width = 3, tags$span("Prediction layer:",
                                class = "text-widget-title")),
    column(width = 6,
           selectInput("bias_gspace_layer",
                       label = NULL,
                       choices = layer_choices,
                       selected = keep))
  )
})

output$bias_gspace_plot_bottom_options_ui <- renderUI({
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  fluidRow(class = "ell-row",
           column(width = 12,
                  tags$span(ell$ell_name, class = "text-center",
                            style = "font-size: 12px; color: #888; font-weight: 400;"))
  )

})

# The apply-bias selector drives the working slot, the same way the Predict
# tab works.
observeEvent(input$bias_ellipsoid_selected, {

  sel <- input$bias_ellipsoid_selected
  req(sel)

  if(identical(sel, "all")) return()

  ell <- session_data$ellipsoid_list[[sel]]
  req(ell)

  if(identical(session_data$current_ellipsoid$ell_id, sel)) return()

  set_working_ellipsoid(ell, mode = "view")
})

output$bias_gspace_plot <- renderPlot({

  req(session_data$bg_raster)

  ell <- session_data$current_ellipsoid
  req(ell)

  layer <- input$bias_gspace_layer
  req(layer)

  s <- bias_collect_plot_settings()

  id <- ell$ell_id
  pred_result <- session_data$ellipsoid_prediction_list[[id]]
  bias_result <- session_data$ellipsoid_prediction_list_biased[[id]]

  has_pred <- inherits(pred_result, "SpatRaster")
  has_bias <- inherits(bias_result, "SpatRaster")

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if(!has_pred || !layer %in% names(pred_result)){
    plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.5, "No prediction for this ellipsoid.",
         cex = 1.1, col = "grey50")
    return(invisible(NULL))
  }

  # No trailing underscore, so this matches whether or not apply_bias()
  # appends the direction to the layer name
  matched_bias <- if(has_bias){
    names(bias_result)[startsWith(names(bias_result), paste0(layer, "_biased"))]
  } else {
    character(0)
  }

  if(length(matched_bias) == 0){
    par(mar = c(4, 3, 2.5, 5))
    terra::plot(pred_result[[layer]],
                main = layer,
                axes = TRUE,
                col = s$palette,
                colNA = s$map_bg_col)
    mtext(paste0("No biased layer yet for '", layer, "'"),
          side = 1, line = 2.5, cex = 0.8, col = "grey50")
    return(invisible(NULL))
  }

  # Unbiased first, then every biased version of that layer. Two columns so
  # each map stays large. Legends are per panel on purpose: multiplying by
  # bias compresses the range, so a shared scale would hide the effect.
  panels <- c(layer, matched_bias)
  n <- length(panels)

  n_cols <- if(n == 1) 1L else 2L
  n_rows <- ceiling(n / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2.5, 5), oma = c(0, 0, 2, 0))

  for(i in seq_along(panels)){

    nm <- panels[i]
    r <- if(i == 1) pred_result[[nm]] else bias_result[[nm]]
    ttl <- if(i == 1) paste0(nm, " (unbiased)") else nm

    terra::plot(r,
                main = ttl,
                axes = TRUE,
                col = s$palette,
                colNA = s$map_bg_col)
  }

  if(n > 1 && n %% 2 != 0) plot.new()

})

output$bias_plot_settings_ui <- renderUI({

  req(session_data$bias_raster)

  box(title = tagList(
    tags$span("Advanced Plot Settings", class = "text-section-header"),
    tags$span(icon("circle-info"),
              title = instructions$bias_plot_settings,
              class = "tooltip-icon")),
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE,

    fluidRow(
      column(width = 4,
             tags$span("Continuous palette", class = "text-widget-title"),
             selectInput("bias_palette", label = NULL,
                         choices = c("Viridis" = "viridis",
                                     "Plasma" = "Plasma",
                                     "Inferno" = "Inferno",
                                     "Rocket" = "Rocket",
                                     "Mako" = "Mako",
                                     "YlGnBu" = "YlGnBu",
                                     "Heat" = "Heat"),
                         selected = "viridis")),
      column(width = 4,
             tags$span("Direction", class = "text-widget-title"),
             checkboxInput("bias_palette_rev", "Reverse palette",
                           value = FALSE)),
      column(width = 4,
             tags$span("Map background (NA)", class = "text-widget-title"),
             tags$div(
               style = "display: flex; align-items: center; gap: 8px;",
               tags$input(type = "color",
                          value = "#F0F0F0",
                          oninput = "Shiny.setInputValue('bias_plot_map_bg_col', this.value)",
                          style = "width: 40px; height: 32px; padding: 2px;
 border: 1px solid #ccc; border-radius: 4px; cursor: pointer;")
             ))
    ),

    fluidRow(
      column(width = 12,
             class = "btn-spaced",
             actionButton("bias_open_export_modal",
                          tagList(icon("download"), "Export figure"),
                          class = "btn-default"))
    )
  )
})

observeEvent(input$bias_open_export_modal, {

  showModal(modalDialog(
    title = "Export Figure",
    size = "m",

    fluidRow(
      column(width = 6,
             tags$span("File type", class = "text-widget-title"),
             radioButtons("bias_export_filetype", label = NULL,
                          choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                          selected = if(!is.null(input$bias_export_filetype))
                            input$bias_export_filetype else "png",
                          inline = TRUE)),
      column(width = 6,
             tags$span("Unit", class = "text-widget-title"),
             radioButtons("bias_export_unit", label = NULL,
                          choiceNames = list(
                            tags$span("mm", class = "text-widget-inner"),
                            tags$span("inches", class = "text-widget-inner"),
                            tags$span("px", class = "text-widget-inner")
                          ),
                          choiceValues = c("mm", "in", "px"),
                          selected = if(!is.null(input$bias_export_unit))
                            input$bias_export_unit else "mm",
                          inline = TRUE))
    ),

    # cex sits outside the renderUI so it does not reset when type or unit change
    fluidRow(
      column(width = 4,
             tags$span("Text size (cex)", class = "text-widget-title"),
             numericInput("bias_export_cex", label = NULL,
                          value = if(!is.null(input$bias_export_cex))
                            input$bias_export_cex else 1,
                          min = 0.5, max = 3, step = 0.1))
    ),

    uiOutput("bias_export_settings_ui"),

    footer = tagList(
      modalButton("Cancel"),
      downloadButton("bias_confirm_export", "Export", class = "btn-save",
                     onclick = "Shiny.setInputValue('bias_export_done',
                     Math.random(), {priority: 'event'});")
    ),
    easyClose = TRUE
  ))
})

output$bias_export_settings_ui <- renderUI({

  ext <- if(!is.null(input$bias_export_filetype)) input$bias_export_filetype else "png"
  unit <- if(!is.null(input$bias_export_unit)) input$bias_export_unit else "mm"

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
             numericInput("bias_export_width_val", label = NULL,
                          value = defaults$val, min = defaults$min,
                          max = defaults$max, step = defaults$step)),
      column(width = if(ext == "png") 4 else 6,
             tags$span(paste0("Height (", unit, ")"), class = "text-widget-title"),
             numericInput("bias_export_height_val", label = NULL,
                          value = defaults$val, min = defaults$min,
                          max = defaults$max, step = defaults$step)),
      if(ext == "png"){
        column(width = 4,
               tags$span("Resolution (dpi)", class = "text-widget-title"),
               numericInput("bias_export_res", label = NULL,
                            value = 300, min = 72, max = 600, step = 50))
      }
    ),

    if(ext != "png") p("PDF and SVG do not require a resolution setting.",
                       class = "text-instruction")
  )
})

output$bias_confirm_export <- downloadHandler(

  filename = function(){
    ext <- if(!is.null(input$bias_export_filetype)) input$bias_export_filetype else "png"
    tab <- if(!is.null(input$bias_plot_tabs)) input$bias_plot_tabs else "bias_input_layers_plot_tab"
    prefix <- switch(tab,
                     "bias_input_layers_plot_tab" = "bias_layers_plot",
                     "bias_composite_plot_tab" = "bias_composite_plot",
                     "bias_gspace_plot_tab" = "bias_gspace_plot",
                     "bias_plot")
    paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
  },

  content = function(file){

    ext <- if(!is.null(input$bias_export_filetype)) input$bias_export_filetype else "png"
    tab <- if(!is.null(input$bias_plot_tabs)) input$bias_plot_tabs else "bias_input_layers_plot_tab"

    s <- bias_collect_plot_settings()

    plot_open_device(file, ext,
                     w_val = input$bias_export_width_val,
                     h_val = input$bias_export_height_val,
                     unit = input$bias_export_unit,
                     res = input$bias_export_res,
                     cex_val = input$bias_export_cex)
    on.exit(dev.off())

    switch(tab,

           "bias_input_layers_plot_tab" = {

             req(session_data$bias_raster)

             rast <- session_data$bias_raster
             n_raw <- terra::nlyr(rast)

             has_prepared <- !is.null(session_data$prepared_bias) &&
               !is.null(session_data$prepared_bias$processed_layers)

             if(has_prepared){
               processed <- session_data$prepared_bias$processed_layers
               n_proc <- terra::nlyr(processed)
               n_rows <- max(n_raw, n_proc)
               par(mfrow = c(n_rows, 2), mar = c(3, 3, 2, 4))

               for(i in seq_len(n_rows)){
                 if(i <= n_raw){
                   terra::plot(rast[[i]], main = paste0("Raw: ", names(rast[[i]])),
                               axes = TRUE, col = s$palette, colNA = s$map_bg_col)
                 } else {
                   plot.new()
                 }
                 if(i <= n_proc){
                   terra::plot(processed[[i]], main = names(processed[[i]]),
                               axes = TRUE, col = s$palette, colNA = s$map_bg_col)
                 } else {
                   plot.new()
                 }
               }
             } else {
               par(mfrow = c(n_raw, 1), mar = c(3, 3, 2, 4))
               for(i in seq_len(n_raw)){
                 terra::plot(rast[[i]], main = paste0("Raw: ", names(rast[[i]])),
                             axes = TRUE, col = s$palette, colNA = s$map_bg_col)
               }
             }
           },

           "bias_composite_plot_tab" = {

             req(session_data$prepared_bias)
             req(session_data$prepared_bias$composite_surface)

             par(mar = c(3, 3, 2, 4))
             terra::plot(session_data$prepared_bias$composite_surface,
                         main = paste0("Composite: ",
                                       session_data$prepared_bias$combination_formula),
                         axes = TRUE, col = s$palette, colNA = s$map_bg_col)
           },

           "bias_gspace_plot_tab" = {

             ell <- session_data$current_ellipsoid
             req(ell)

             layer <- input$bias_gspace_layer
             req(layer)

             pred_result <- session_data$ellipsoid_prediction_list[[ell$ell_id]]
             bias_result <- session_data$ellipsoid_prediction_list_biased[[ell$ell_id]]

             req(inherits(pred_result, "SpatRaster"))
             req(layer %in% names(pred_result))

             matched_bias <- if(inherits(bias_result, "SpatRaster")){
               names(bias_result)[startsWith(names(bias_result),
                                             paste0(layer, "_biased"))]
             } else {
               character(0)
             }

             panels <- c(layer, matched_bias)
             n <- length(panels)

             n_cols <- if(n == 1) 1L else 2L
             n_rows <- ceiling(n / n_cols)

             par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2.5, 5),
                 oma = c(0, 0, 2, 0))

             for(i in seq_along(panels)){
               nm <- panels[i]
               r <- if(i == 1) pred_result[[nm]] else bias_result[[nm]]
               ttl <- if(i == 1) paste0(nm, " (unbiased)") else nm
               terra::plot(r, main = ttl, axes = TRUE,
                           col = s$palette, colNA = s$map_bg_col)
             }

             if(n > 1 && n %% 2 != 0) plot.new()

             mtext(ell$ell_name, outer = TRUE, cex = 1, font = 2)
           }
    )
  }
)


observeEvent(input$bias_export_done, {

  removeModal()

  showNotification("Figure exported successfully.",
                   type = "message", duration = 4)
})
