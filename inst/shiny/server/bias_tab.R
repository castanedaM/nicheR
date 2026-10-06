# Title: Bias Tab server

# Description: Server for the bias tab. Uploads one or more sampling bias
# rasters, prepares them into a composite surface, and applies that surface
# to saved predictions. The whole tab is optional and can be skipped.

# Date Last Updated: 10/06/2026


# BIASED RECORDS ----------------------------------------------------------

# Removes every record set that was drawn from a biased layer. Called
# wherever the biased predictions are cleared, since those sets would
# otherwise point at surfaces that no longer exist. Returns how many sets
# were removed.
drop_biased_records <- function(){

  occ <- session_data$ellipsoid_records_list
  if(length(occ) == 0) return(0L)

  n_removed <- 0L

  for(id in names(occ)){

    is_biased <- vapply(occ[[id]], function(df){
      isTRUE(attr(df, "biased", exact = TRUE))
    }, logical(1))

    n_removed <- n_removed + sum(is_biased)
    occ[[id]] <- occ[[id]][!is_biased]
  }

  # An ellipsoid with no sets left goes back to "not generated"
  occ <- occ[vapply(occ, length, integer(1)) > 0]

  session_data$ellipsoid_records_list <- occ

  n_removed
}


# SKIP AND CONTINUE -------------------------------------------------------

# Sits above the inputs so a user redirected here from Predict can move on
# without providing a bias layer
output$bias_skip_ui <- renderUI({

  div(class = "action-btn-row",
      actionButton(inputId = "bias_skip_btn",
                   label = tagList("Skip bias", icon("forward")),
                   class = "btn-back")
  )
})

observeEvent(input$bias_skip_btn, {

  showModal(modalDialog(
    title = "Skip sampling bias?",
    p(instructions$bias_skip, class = "text-instruction"),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("bias_confirm_skip_btn",
                   "Yes, skip",
                   class = "btn-continue")
    ),
    easyClose = TRUE
  ))
})

observeEvent(input$bias_confirm_skip_btn, {
  removeModal()
  updateTabItems(session, "sidebar_menu", selected = "generate_tab")
})

output$bias_next_step_ui <- renderUI({

  req(length(session_data$ellipsoid_prediction_list_biased) > 0)

  div(class = "action-btn-row",
      actionButton(inputId = "bias_next_step_btn",
                   label = tagList("Continue",
                                   icon("arrow-right")),
                   class = "btn-save")
  )
})

observeEvent(input$bias_next_step_btn, {
  updateTabItems(session, "sidebar_menu", selected = "generate_tab")
})


# INPUT BIAS LAYERS -------------------------------------------------------

# Loads the selected bias file(s) once, so the preview and the Upload
# button do not parse the same files twice
bias_raster_upload <- reactive({

  files <- input$bias_raster_file
  if(is.null(files)) return(NULL)

  if(nrow(files) > 10){
    showNotification("Maximum 10 raster files allowed.",
                     type = "warning", duration = 4)
    return(NULL)
  }

  rasters <- lapply(seq_len(nrow(files)), function(i){
    ext <- tolower(tools::file_ext(files$name[i]))
    tryCatch(
      load_raster_file(files$datapath[i], ext),
      error = function(e){
        showNotification(paste0("Could not load ", files$name[i], ": ", e$message),
                         type = "error", duration = 4)
        NULL
      }
    )
  })

  rasters <- Filter(Negate(is.null), rasters)

  if(length(rasters) == 0){
    showNotification("No raster files could be loaded.",
                     type = "error", duration = 4)
    return(NULL)
  }

  if(length(rasters) == 1) return(rasters[[1]])

  tryCatch(
    do.call(c, rasters),
    error = function(e){
      showNotification(paste0("Could not stack bias rasters: ", e$message,
                              " Check that all files have matching resolution, ",
                              "extent, and CRS."),
                       type = "error", duration = 6)
      NULL
    }
  )
})

output$bias_raster_print <- renderPrint({
  rast <- bias_raster_upload()
  req(rast)
  print(rast)
})

observeEvent(input$bias_upload_btn, {

  rast <- bias_raster_upload()

  if(is.null(rast)){
    showNotification("No bias file selected.", type = "warning", duration = 4)
    return()
  }

  session_data$bias_raster <- rast
  session_data$bias_source <- input$bias_raster_file$name

  session_data$prepared_bias <- NULL
  session_data$bias_settings <- NULL
  session_data$ellipsoid_prediction_list_biased <- list()
  drop_biased_records()

  showNotification(paste0(terra::nlyr(rast), " bias layer(s) loaded successfully."),
                   type = "message", duration = 4)
})

observeEvent(input$bias_example_btn, {

  rast <- tryCatch(
    terra::rast(system.file("extdata", "ma_biases.tif", package = "nicheR")),
    error = function(e){
      showNotification(paste("Could not load example bias:", e$message),
                       type = "error", duration = 4)
      NULL
    }
  )

  req(rast)

  session_data$bias_raster <- rast
  session_data$bias_source <- "example"
  session_data$prepared_bias <- NULL
  session_data$bias_settings <- NULL
  session_data$ellipsoid_prediction_list_biased <- list()
  drop_biased_records()

  showNotification("Example bias raster loaded.", type = "message", duration = 4)
})

observeEvent(input$bias_edit_upload_link, {

  showModal(modalDialog(
    title = "Edit bias raster?",
    p(instructions$bias_edit_upload, class = "text-instruction"),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("bias_confirm_edit_upload_btn",
                   "Yes, edit",
                   class = "btn-cancel")
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$bias_confirm_edit_upload_btn, {

  removeModal()

  session_data$bias_raster <- NULL
  session_data$bias_source <- NULL
  session_data$prepared_bias <- NULL
  session_data$bias_settings <- NULL
  session_data$ellipsoid_prediction_list_biased <- list()
  n_removed <- drop_biased_records()
  shinyjs::reset("bias_raster_file")

  showNotification(paste0("Bias raster cleared. Upload a new file.",
                          if(n_removed > 0){
                            paste0(" ", n_removed, " record set(s) from ",
                                   "biased layers were removed.")
                          }),
                   type = "message", duration = 4)
})

output$bias_upload_ui <- renderUI({

  if(identical(session_data$input_mode, "virtual")){
    return(column(12,
                  p(instructions$bias_virtual_unavailable, class = "text-instruction")))
  }

  # Bias is inherently geographic, so it needs a raster study area
  if(is.null(session_data$bg_raster)){
    return(
      box(title = tags$span("Bias inputs", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          p(instructions$bias_no_gspace, class = "text-instruction"))
    )
  }

  has_pred <- length(session_data$ellipsoid_prediction_list) > 0
  has_bias <- !is.null(session_data$bias_raster)

  if(!has_pred){
    return(
      box(title = tags$span("Bias inputs", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          p(instructions$bias_needs_prediction, class = "text-instruction"))
    )
  }

  if(has_bias){
    return(
      box(title = tags$span("Bias inputs", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          p(paste0(terra::nlyr(session_data$bias_raster),
                   " bias layer(s) loaded."),
            class = "text-instruction"),
          fluidRow(
            column(width = 12, class = "btn-spaced",
                   actionLink("bias_edit_upload_link",
                              label = tagList(icon("pen"), "Edit bias raster")))
          )
      )
    )
  }

  box(title = tags$span("Bias input layers", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$bias_input, class = "text-instruction"),

      fluidRow(
        column(width = 12,
               fileInput(inputId = "bias_raster_file",
                         label = tags$span("Sampling bias layer/s (raster)",
                                           class = "text-widget-title"),
                         multiple = TRUE,
                         accept = c(".tif", ".tiff", ".rds")))
      ),

      verbatimTextOutput("bias_raster_print"),

      fluidRow(
        column(width = 12,
               div(class = "action-btn-row",
                   actionButton("bias_upload_btn",
                                "Upload bias",
                                class = "btn-continue"),
                   if(identical(session_data$input_mode, "example")){
                     actionButton("bias_example_btn",
                                  "Example bias",
                                  class = "btn-back")
                   }
               ))
      )
  )
})


# PREPARE BIAS ------------------------------------------------------------

output$bias_effect_directions_ui <- renderUI({

  req(session_data$bias_raster)

  rast <- session_data$bias_raster
  lyrs <- names(rast)

  stats <- lapply(lyrs, function(nm){
    vals <- terra::values(rast[[nm]], na.rm = TRUE)
    list(mean = round(mean(vals, na.rm = TRUE), 3),
         sd = round(sd(vals, na.rm = TRUE), 3))
  })

  header <- fluidRow(
    class = "ell-row",
    column(width = 4,
           tags$span("Layer", class = "text-widget-title")),
    column(width = 3,
           tags$div(class = "tooltip-label-row",
                    tags$span("Mean (SD)", class = "text-widget-title"),
                    tags$span(icon("circle-info"),
                              title = instructions$bias_layer_stats_tooltip,
                              class = "tooltip-icon"))),
    column(width = 5,
           tags$div(class = "tooltip-label-row",
                    tags$span("Effect direction", class = "text-widget-title"),
                    tags$span(icon("circle-info"),
                              title = instructions$bias_direction_tooltip,
                              class = "tooltip-icon")))
  )

  rows <- lapply(seq_along(lyrs), function(i){
    nm <- lyrs[i]
    fluidRow(
      class = "ell-row",
      column(width = 4,
             tags$span(nm, class = "text-widget-inner")),
      column(width = 3,
             tags$span(paste0(stats[[i]]$mean, " (", stats[[i]]$sd, ")"),
                       class = "text-widget-inner")),
      column(width = 5,
             radioButtons(inputId = paste0("bias_dir_", i),
                          label = NULL,
                          choiceNames = list(
                            tags$span("Direct", class = "text-widget-inner"),
                            tags$span("Inverse", class = "text-widget-inner")),
                          choiceValues = c("direct", "inverse"),
                          selected = "direct",
                          inline = TRUE))
    )
  })

  tagList(header, tagList(rows))
})

observeEvent(input$bias_prepare_btn, {

  req(session_data$bias_raster)

  lyrs <- names(session_data$bias_raster)

  # Indexed rather than keyed by layer name, so a name with a space or a
  # dash cannot produce an unreachable input id
  effect_direction <- vapply(seq_along(lyrs), function(i){
    val <- input[[paste0("bias_dir_", i)]]
    if(is.null(val)) "direct" else val
  }, character(1))

  names(effect_direction) <- lyrs

  mask_na <- isTRUE(input$bias_mask_na == "TRUE")

  result <- tryCatch(
    prepare_bias(bias_surface = session_data$bias_raster,
                 effect_direction = effect_direction,
                 include_composite = TRUE,
                 include_processed_layers = TRUE,
                 mask_na = mask_na,
                 verbose = FALSE),
    error = function(e){
      showNotification(paste("Bias preparation failed:", e$message),
                       type = "error", duration = 5)
      NULL
    }
  )

  req(result)

  session_data$prepared_bias <- result

  # mask_na leaves no trace in the returned object, and the per layer
  # directions survive only as a string that has to be parsed back apart, so
  # both are recorded here for the session report
  session_data$bias_settings <- list(effect_direction = effect_direction,
                                     mask_na = mask_na)
  session_data$ellipsoid_prediction_list_biased <- list()
  drop_biased_records()

  showNotification("Bias prepared successfully.", type = "message", duration = 4)
})

observeEvent(input$bias_edit_prepare_link, {

  showModal(modalDialog(
    title = "Edit bias preparation?",
    p(instructions$bias_edit_prepare, class = "text-instruction"),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("bias_confirm_edit_prepare_btn",
                   "Yes, edit",
                   class = "btn-cancel")
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$bias_confirm_edit_prepare_btn, {

  removeModal()

  session_data$prepared_bias <- NULL
  session_data$bias_settings <- NULL
  session_data$ellipsoid_prediction_list_biased <- list()
  n_removed <- drop_biased_records()

  showNotification(paste0("Bias preparation cleared. Adjust settings and re-prepare.",
                          if(n_removed > 0){
                            paste0(" ", n_removed, " record set(s) from ",
                                   "biased layers were removed.")
                          }),
                   type = "message", duration = 4)
})

output$bias_prepare_ui <- renderUI({

  req(session_data$bias_raster)

  if(!is.null(session_data$prepared_bias)){

    result <- session_data$prepared_bias

    return(
      box(title = tags$span("Prepare bias raster", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          p(paste0("Formula: ", result$combination_formula),
            class = "text-instruction"),
          fluidRow(
            column(width = 12, class = "btn-spaced",
                   actionLink("bias_edit_prepare_link",
                              label = tagList(icon("pen"), "Edit bias preparation")))
          )
      )
    )
  }

  box(title = tags$span("Prepare bias raster", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$bias_prepare, class = "text-instruction"),

      uiOutput("bias_effect_directions_ui"),

      fluidRow(
        column(width = 12,
               tags$div(class = "tooltip-label-row",
                        tags$span("NA overlap handling", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$bias_mask_na_tooltip,
                                  class = "tooltip-icon")),
               radioButtons(inputId = "bias_mask_na",
                            label = NULL,
                            choiceNames = list(
                              tags$span("Union", class = "text-widget-inner"),
                              tags$span("Intersection", class = "text-widget-inner")),
                            choiceValues = c("FALSE", "TRUE"),
                            selected = "FALSE",
                            inline = TRUE))
      ),

      fluidRow(
        column(width = 12, class = "action-btn-row",
               actionButton("bias_prepare_btn",
                            "Prepare bias",
                            class = "btn-continue"))
      )
  )
})


# APPLY BIAS --------------------------------------------------------------

# TRUE while the apply form is open on top of an existing set of results
bias_show_apply_form <- reactiveVal(FALSE)

output$bias_ellipsoid_selector_ui <- renderUI({

  req(length(session_data$ellipsoid_prediction_list) > 0)

  pred_ids <- names(session_data$ellipsoid_prediction_list)
  versions <- session_data$ellipsoid_list

  ell_choices <- c(
    "All versions" = "all",
    setNames(pred_ids,
             vapply(pred_ids, function(id){
               ell <- versions[[id]]
               if(!is.null(ell) && !is.null(ell$ell_name)) ell$ell_name else id
             }, character(1)))
  )

  sel <- isolate(input$bias_ellipsoid_selected)

  keep <- if(!is.null(sel) && sel %in% ell_choices) sel else "all"

  selectInput(inputId = "bias_ellipsoid_selected",
              label = tagList(
                tags$span("Ellipsoid version", class = "text-widget-title"),
                tags$span(icon("circle-info"),
                          title = instructions$bias_ellipsoid_select_tooltip,
                          class = "tooltip-icon")
              ),
              choices = ell_choices,
              selected = keep)
})

output$bias_apply_ui <- renderUI({

  req(session_data$prepared_bias)
  req(length(session_data$ellipsoid_prediction_list) > 0)

  has_applied <- length(session_data$ellipsoid_prediction_list_biased) > 0

  if(has_applied && !isTRUE(bias_show_apply_form())){

    n_layers <- sum(vapply(session_data$ellipsoid_prediction_list_biased,
                           function(rast) terra::nlyr(rast),
                           numeric(1)))

    return(
      box(title = tags$span("Apply bias", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          # The way back into the form is Add new bias layer, at the bottom of
          # the ellipsoid library
          p(paste0(n_layers, " biased layer(s) across ",
                   length(session_data$ellipsoid_prediction_list_biased),
                   " ellipsoid(s). To add more, use Add new bias layer in ",
                   "the ellipsoid library."),
            class = "text-instruction")
      )
    )
  }

  box(title = tags$span("Apply bias", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$bias_apply, class = "text-instruction"),

      uiOutput("bias_ellipsoid_selector_ui"),

      uiOutput("bias_layer_selector_ui"),

      fluidRow(
        column(width = 12,
               tags$div(class = "tooltip-label-row",
                        tags$span("Prediction layer effect direction", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$bias_apply_direction_tooltip,
                                  class = "tooltip-icon")),
               radioButtons("bias_effect_direction",
                            label = NULL,
                            choiceNames = list(
                              tags$span("Direct", class = "text-widget-inner"),
                              tags$span("Inverse", class = "text-widget-inner")),
                            choiceValues = c("direct", "inverse"),
                            selected = "direct",
                            inline = TRUE))
      ),

      fluidRow(
        column(width = 12, class = "action-btn-row",
               if(has_applied){
                 actionButton("bias_cancel_apply_btn",
                              "Done",
                              class = "btn-back")
               },
               actionButton("bias_apply_btn",
                            label = "Apply bias",
                            class = "btn-continue"))
      )
  )
})
# Layers offered depend on which version is selected, so this is a separate
# output. Reading the version input inside bias_apply_ui would rebuild the
# whole box, and the direction radio with it, on every version change.
output$bias_layer_selector_ui <- renderUI({

  pred_list <- session_data$ellipsoid_prediction_list
  req(length(pred_list) > 0)

  sel_ell <- input$bias_ellipsoid_selected
  if(is.null(sel_ell)) sel_ell <- "all"

  ids <- if(identical(sel_ell, "all")) names(pred_list) else sel_ell

  # The union across the selected versions rather than the first one's layers,
  # since versions can be predicted with different layer sets. A layer a given
  # version does not carry is skipped for that version at apply time.
  layers <- unique(unlist(lapply(ids, function(id){
    p <- pred_list[[id]]
    ell <- session_data$ellipsoid_list[[id]]
    if(!inherits(p, "SpatRaster") || is.null(ell)) return(character(0))
    report_pred_layer_names(p, ell)
  })))

  if(length(layers) == 0) return(NULL)

  choices <- c("All prediction layers" = "all_pred", setNames(layers, layers))

  # Isolated, since this is read only to keep a still-valid choice across a
  # re-render and not to make the output depend on the input it creates
  prev <- isolate(input$bias_prediction_layer)

  keep <- if(!is.null(prev) && prev %in% choices){
    prev
  } else if("suitability_trunc" %in% layers){
    "suitability_trunc"
  } else if("suitability" %in% layers){
    "suitability"
  } else {
    layers[1]
  }

  fluidRow(
    column(width = 12,
           tags$div(class = "tooltip-label-row",
                    tags$span("Prediction layer", class = "text-widget-title"),
                    tags$span(icon("circle-info"),
                              title = instructions$bias_layer_select_tooltip,
                              class = "tooltip-icon")),
           selectInput("bias_prediction_layer",
                       label = NULL,
                       choices = choices,
                       selected = keep))
  )
})

# Add new bias layer, in the ellipsoid library. Reopens the Apply bias box as
# the form and scrolls up to it, since the box sits above the library. The
# short wait lets the box redraw before the page moves.
observeEvent(input$bias_add_layer_btn, {

  bias_show_apply_form(TRUE)

  shinyjs::runjs("
    setTimeout(function(){
      var box = document.getElementById('bias_apply_ui');
      if(box) box.scrollIntoView({behavior: 'smooth', block: 'start'});
    }, 150);
  ")
})

observeEvent(input$bias_cancel_apply_btn, {
  bias_show_apply_form(FALSE)
})

observeEvent(input$bias_apply_btn, {

  req(session_data$prepared_bias)
  req(length(session_data$ellipsoid_prediction_list) > 0)
  req(input$bias_prediction_layer)
  req(input$bias_ellipsoid_selected)

  pred_list <- session_data$ellipsoid_prediction_list

  selected_ids <- if(identical(input$bias_ellipsoid_selected, "all")){
    names(pred_list)
  } else {
    input$bias_ellipsoid_selected
  }

  direction <- if(!is.null(input$bias_effect_direction)){
    input$bias_effect_direction
  } else {
    "direct"
  }

  is_all_pred <- identical(input$bias_prediction_layer, "all_pred")

  current_biased <- session_data$ellipsoid_prediction_list_biased
  if(is.null(current_biased)) current_biased <- list()

  n_success <- 0L
  n_skipped <- 0L

  for(id in selected_ids){

    pred <- pred_list[[id]]

    if(is.null(pred) || !inherits(pred, "SpatRaster")){
      showNotification(paste0("No raster prediction found for ", id, ". Skipping."),
                       type = "warning", duration = 4)
      next
    }

    # Layers read per ellipsoid rather than from the first one, since
    # versions can be predicted with different layer sets. Only what predict()
    # added counts: the stored prediction is made with keep_data = TRUE, so
    # names(pred) also holds the environmental variables, and those are not
    # something to bias.
    target_layers <- if(is_all_pred){
      report_pred_layer_names(pred, session_data$ellipsoid_list[[id]])
    } else {
      input$bias_prediction_layer
    }

    for(layer in target_layers){

      if(!layer %in% names(pred)){
        showNotification(paste0("Layer '", layer, "' not found for ", id, ". Skipping."),
                         type = "warning", duration = 4)
        next
      }

      prev_layers <- if(id %in% names(current_biased)){
        terra::nlyr(current_biased[[id]])
      } else {
        0L
      }

      current_biased <- apply_bias_to_list(
        biased_list = current_biased,
        ell_id = id,
        pred_rast = pred,
        layer = layer,
        prepared_bias = session_data$prepared_bias,
        direction = direction
      )

      new_layers <- if(id %in% names(current_biased)){
        terra::nlyr(current_biased[[id]])
      } else {
        0L
      }

      if(new_layers > prev_layers){
        n_success <- n_success + 1L
      } else {
        n_skipped <- n_skipped + 1L
      }
    }
  }

  if(n_success == 0L && n_skipped > 0L){
    showNotification("All selected combinations already exist. No new layers added.",
                     type = "warning", duration = 4)
    return()
  }

  if(n_success == 0L){
    showNotification("Bias application failed for all selected combinations.",
                     type = "error", duration = 4)
    return()
  }

  session_data$ellipsoid_prediction_list_biased <- current_biased
  bias_show_apply_form(FALSE)

  msg <- paste0(n_success, " new biased layer(s) added.")
  if(n_skipped > 0L){
    msg <- paste0(msg, " ", n_skipped, " already existed and were skipped.")
  }

  showNotification(msg, type = "message", duration = 4)

  # Jump to the comparison panel so the result is visible immediately
  updateTabsetPanel(session, "bias_plot_tabs",
                    selected = "bias_gspace_plot_tab")
})


# ELLIPSOID LIBRARY -------------------------------------------------------

output$bias_ellipsoid_library_ui <- renderUI({

  cur_ell <- session_data$current_ellipsoid
  versions <- session_data$ellipsoid_list
  ids <- names(versions)

  req(!is.null(cur_ell) || length(ids) > 0)

  biased <- names(session_data$ellipsoid_prediction_list_biased)

  working_row <- if(!is.null(cur_ell)){
    fluidRow(
      class = "ell-row",
      style = "background: #f0f7f0; border-radius: 4px; margin-bottom: 6px; padding: 4px 0;",
      column(width = 5,
             tags$span(icon("eye"),
                       tags$span(paste0(" ", cur_ell$ell_name),
                                 class = "text-widget-inner",
                                 style = "color: #097a21; font-weight: 500;"))),
      column(width = 4,
             tags$span(ell_lineage_label(cur_ell),
                       style = "font-size: 11px; color: #aaa;")),
      column(width = 3,
             tags$span("View-only", style = "font-size: 11px; color: #aaa;"))
    )
  }

  # The composite and the prepared layers exist once bias is prepared, with
  # or without anything applied, so that is what turns the download on
  has_prepared <- !is.null(session_data$prepared_bias)

  rows <- lapply(ids, function(id){

    ell <- versions[[id]]

    n_biased <- if(id %in% biased){
      terra::nlyr(session_data$ellipsoid_prediction_list_biased[[id]])
    } else {
      0L
    }

    fluidRow(
      class = "ell-row",
      style = "padding: 2px 0;",
      column(width = 5,
             tags$span(ell$ell_name, class = "text-widget-inner"),
             tags$br(),
             tags$span(id, style = "font-size: 10px; color: #bbb;")),
      column(width = 4,
             tags$span(ell_lineage_label(ell),
                       style = "font-size: 11px; color: #aaa;"),
             tags$br(),
             tags$span(if(n_biased > 0){
               paste0(n_biased, " biased layer(s)")
             } else {
               "no bias applied"
             },
             style = paste0("font-size: 10px; color: ",
                            if(n_biased > 0) "#097a21;" else "#bbb;"))),
      column(width = 3,
             class = "ell-actions",
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('bias_ell_view', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("View ", ell$ell_name, " (read-only)"),
                    icon("eye")),
             # Nothing to download until bias has been prepared
             if(has_prepared){
               tags$a(href = "#",
                      onclick = sprintf("Shiny.setInputValue('bias_ell_download', '%s', {priority: 'event'}); return false;", id),
                      title = paste0("Download bias layers for ", ell$ell_name),
                      icon("download"))
             } else {
               tags$span(icon("download"),
                         title = "Prepare bias first",
                         style = "color: #ddd; cursor: not-allowed;")
             },
             tags$a(href = "#",
                    class = "ell-action-danger",
                    onclick = sprintf("Shiny.setInputValue('bias_ell_delete', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("Delete ", ell$ell_name),
                    icon("trash-can")))
    )
  })

  # Sits under the last ellipsoid, like Create new ellipsoid on the Build
  # tab. Shown only once bias has been applied and the Apply bias box has
  # collapsed, since until then that box is already open as the form.
  has_applied <- length(session_data$ellipsoid_prediction_list_biased) > 0

  add_layer_row <- if(has_applied && !is.null(session_data$prepared_bias) &&
                      !isTRUE(bias_show_apply_form())){
    tagList(
      tags$hr(style = "margin: 8px 0;"),
      fluidRow(
        column(width = 12, class = "btn-spaced",
               actionLink("bias_add_layer_btn",
                          label = tagList(icon("plus"),
                                          "Add new bias layer")))
      )
    )
  }

  box(title = tags$span("Ellipsoid library", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$bias_library, class = "text-instruction"),

      if(!is.null(cur_ell)){
        tagList(working_row, tags$hr(style = "margin: 8px 0;"))
      },

      if(length(ids) > 0){
        tagList(
          fluidRow(
            class = "ell-row",
            style = "padding: 2px 0;",
            column(width = 5, tags$span("Name", class = "text-widget-title")),
            column(width = 4, tags$span("Details", class = "text-widget-title")),
            column(width = 3, tags$span("Actions", class = "text-widget-title"))
          ),
          tagList(rows)
        )
      } else {
        p(instructions$bias_library_empty, class = "text-muted-small")
      },

      add_layer_row
  )
})

# View, read-only
observeEvent(input$bias_ell_view, {

  ell <- session_data$ellipsoid_list[[input$bias_ell_view]]
  req(ell)

  set_working_ellipsoid(ell, mode = "view")

  if(!identical(input$bias_ellipsoid_selected, ell$ell_id) &&
     ell$ell_id %in% names(session_data$ellipsoid_prediction_list)){
    updateSelectInput(session, "bias_ellipsoid_selected",
                      selected = ell$ell_id)
  }

  showNotification(paste0("Viewing ", ell$ell_name, "."),
                   type = "message", duration = 3)
})

# Delete, asks first
observeEvent(input$bias_ell_delete, {

  ell <- session_data$ellipsoid_list[[input$bias_ell_delete]]
  req(ell)

  session_data$pending_ell_delete <- input$bias_ell_delete

  n_children <- sum(vapply(session_data$ellipsoid_list, function(e){
    identical(e$parent_id, ell$ell_id)
  }, logical(1)))

  showModal(modalDialog(
    title = paste0("Delete ", ell$ell_name, "?"),
    p(instructions$bias_delete_ell, class = "text-instruction"),
    if(n_children > 0){
      p(paste0(n_children, " ellipsoid(s) were copied from this one. ",
               "They will be kept, but will no longer have a parent."),
        class = "text-muted-small")
    },
    footer = tagList(
      modalButton("Cancel"),
      actionButton("bias_confirm_ell_delete_btn",
                   "Yes, delete",
                   class = "btn-cancel")
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$bias_confirm_ell_delete_btn, {

  id <- session_data$pending_ell_delete
  req(id)

  nm <- session_data$ellipsoid_list[[id]]$ell_name

  removeModal()

  session_data$ellipsoid_list[[id]] <- NULL
  session_data$ellipsoid_prediction_list[[id]] <- NULL
  session_data$prediction_settings[[id]] <- NULL
  session_data$ellipsoid_prediction_list_biased[[id]] <- NULL
  session_data$ellipsoid_records_list[[id]] <- NULL
  session_data$pending_ell_delete <- NULL

  # Copies of the deleted ellipsoid are kept and become roots
  session_data$ellipsoid_list <- lapply(session_data$ellipsoid_list, function(e){
    if(identical(e$parent_id, id)) e$parent_id <- NULL
    e
  })

  cur <- session_data$current_ellipsoid

  if(identical(cur$ell_id, id)){
    clear_working_ellipsoid()
    showNotification(paste0(nm, " deleted. Go back to Build to create a new ellipsoid."),
                     type = "message", duration = 4)
    return()
  }

  if(identical(cur$parent_id, id)){
    cur$parent_id <- NULL
    session_data$current_ellipsoid <- cur
  }

  showNotification(paste0(nm, " deleted."),
                   type = "message", duration = 3)
})


# BIAS DOWNLOAD -----------------------------------------------------------

# Which ellipsoid the download modal on screen belongs to. Kept out of
# session_data so it is not written into a saved session.
bias_dl_target <- reactiveVal(NULL)

# The composite first, then each prepared layer, as one stack. These are the
# same for every ellipsoid.
bias_dl_surfaces <- function(){

  pb <- session_data$prepared_bias
  if(is.null(pb)) return(NULL)

  comp <- pb$composite_surface
  proc <- pb$processed_layers

  if(!inherits(comp, "SpatRaster")) return(NULL)
  if(!inherits(proc, "SpatRaster")) return(comp)

  c(comp, proc)
}

# The results of Apply bias for one ellipsoid. Layers made from an
# environmental variable are left out, a session saved before the Apply step
# was fixed can still carry them.
bias_dl_applied <- function(id){

  r <- session_data$ellipsoid_prediction_list_biased[[id]]
  ell <- session_data$ellipsoid_list[[id]]

  if(!inherits(r, "SpatRaster") || is.null(ell)) return(NULL)

  src <- sub("_biased_(direct|inverse)$", "", names(r))
  keep <- names(r)[!src %in% c("x", "y", ell$var_names)]

  if(length(keep) == 0) return(NULL)

  r[[keep]]
}

# Whether the bias surfaces and the biased predictions sit on one grid.
# Bias is prepared on its own grid and resampled to the prediction when it
# is applied, so the two can differ, and then they cannot share a file.
bias_dl_same_grid <- function(surfaces, applied){
  if(is.null(surfaces) || is.null(applied)) return(TRUE)
  isTRUE(terra::compareGeom(surfaces, applied, stopOnError = FALSE))
}

# What is ticked in the modal, checked against what exists now. The applied
# checkboxes are not drawn for an ellipsoid with no biased layers, and Shiny
# keeps the last value of an input that is no longer drawn.
bias_dl_selection <- function(id){

  surfaces <- bias_dl_surfaces()
  applied <- bias_dl_applied(id)

  surf_sel <- intersect(input$bias_dl_surfaces, names(surfaces))
  app_sel <- intersect(input$bias_dl_applied, names(applied))

  mixed <- length(surf_sel) > 0 && length(app_sel) > 0

  list(surfaces = surfaces,
       applied = applied,
       surf_sel = surf_sel,
       app_sel = app_sel,
       ok = (length(surf_sel) + length(app_sel)) > 0 &&
         (!mixed || bias_dl_same_grid(surfaces, applied)))
}

bias_export_name <- function(ell){
  paste0(ell$ell_id, "_bias_", format(Sys.Date(), "%Y%m%d"))
}

observeEvent(input$bias_ell_download, {

  id <- input$bias_ell_download
  ell <- session_data$ellipsoid_list[[id]]
  surfaces <- bias_dl_surfaces()
  req(ell, surfaces)

  bias_dl_target(id)

  applied <- bias_dl_applied(id)
  same_grid <- bias_dl_same_grid(surfaces, applied)

  surf_names <- names(surfaces)

  # The composite is always the first layer of the stack
  surf_labels <- c(paste0("Composite (", surf_names[1], ")"), surf_names[-1])

  # Everything is ticked when it can all go in one file, as in the prediction
  # download. On different grids only the applied results start ticked.
  surf_selected <- if(same_grid) surf_names else character(0)

  showModal(modalDialog(
    title = paste0("Download bias for ", ell$ell_name),
    p(instructions$bias_download, class = "text-instruction"),

    checkboxGroupInput("bias_dl_surfaces",
                       label = tags$span("Composite and prepared layers",
                                         class = "text-widget-title"),
                       choiceNames = as.list(surf_labels),
                       choiceValues = surf_names,
                       selected = surf_selected),

    if(!is.null(applied)){
      checkboxGroupInput("bias_dl_applied",
                         label = tags$span(paste0("Applied to ", ell$ell_name),
                                           class = "text-widget-title"),
                         choices = names(applied),
                         selected = names(applied))
    } else {
      p("No bias has been applied to this ellipsoid yet.",
        class = "text-muted-small")
    },

    if(!same_grid){
      p(instructions$bias_download_grids, class = "text-muted-small")
    },

    radioButtons("bias_dl_format",
                 label = tags$span("Format", class = "text-widget-title"),
                 choices = c("SpatRaster (.tif)" = "tif",
                             "Data frame (.csv)" = "csv"),
                 selected = "tif",
                 inline = TRUE),

    tags$small("Data frames always include the x and y coordinates.",
               class = "text-muted-small"),
    br(), br(),

    footer = tagList(
      modalButton("Close"),
      downloadButton("bias_dl_btn", "Download", class = "btn-continue")
    ),
    easyClose = FALSE
  ))
})

# The button is off while nothing is ticked, and while the ticks mix layers
# from two different grids
observe({
  id <- bias_dl_target()
  req(id)
  shinyjs::toggleState("bias_dl_btn", condition = bias_dl_selection(id)$ok)
})

output$bias_dl_btn <- downloadHandler(

  filename = function(){
    ell <- session_data$ellipsoid_list[[bias_dl_target()]]
    nm <- gsub("[^A-Za-z0-9._-]", "", bias_export_name(ell))
    ext <- if(identical(input$bias_dl_format, "tif")) ".tif" else ".csv"
    paste0(substr(nm, 1, 60), ext)
  },

  content = function(file){

    sel <- bias_dl_selection(bias_dl_target())

    if(!isTRUE(sel$ok)){
      stop("No layers selected, or the selected layers are on different grids.")
    }

    # Subsetting only, nothing is prepared or applied again
    parts <- list()

    if(length(sel$surf_sel) > 0){
      parts <- c(parts, list(sel$surfaces[[sel$surf_sel]]))
    }

    if(length(sel$app_sel) > 0){
      parts <- c(parts, list(sel$applied[[sel$app_sel]]))
    }

    out <- if(length(parts) == 1) parts[[1]] else c(parts[[1]], parts[[2]])

    if(identical(input$bias_dl_format, "tif")){
      terra::writeRaster(out, file, filetype = "GTiff", overwrite = TRUE)
    } else {
      # na.rm = NA drops only cells that are NA in every layer
      df <- terra::as.data.frame(out, xy = TRUE, na.rm = NA)
      utils::write.csv(df, file, row.names = FALSE)
    }
  }
)
