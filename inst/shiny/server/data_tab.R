# Title: Data tab server logic
# Description: Handles file upload, validation, and variable selection
# Date last updated: 09/29/2026

# Start session button in about
observeEvent(input$about_start_session_btn, {
  updateTabItems(session, "sidebar_menu",
                 selected = "build_tab")
})

# Reactives ---------------------------------------------------------------

# Variable names available for selection, with coordinate columns removed.
# Used by the selector UI, the mutual-exclusion observer, and Confirm.
available_vars <- reactive({
  nms <- if(!is.null(session_data$bg_raster)){
    names(session_data$bg_raster)
  } else if(!is.null(session_data$bg_df)){
    colnames(session_data$bg_df)
  } else {
    return(NULL)
  }

  is_coord <- grepl(X_COL_PATTERN, nms, ignore.case = TRUE) |
    grepl(Y_COL_PATTERN, nms, ignore.case = TRUE)

  nms[!is_coord]
})


# Loads the selected raster file(s) once so the preview and the Upload
# button do not parse the same files twice.
build_raster_upload <- reactive({
  files <- input$build_raster_file
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
      showNotification(paste0("Could not stack rasters: ", e$message,
                              " Check that all files have matching resolution, ",
                              "extent, and CRS."),
                       type = "error", duration = 6)
      NULL
    }
  )
})


# Observer Events ---------------------------------------------------------

observeEvent(input$build_clear_files, {
  session_data$bg_raster <- NULL
  session_data$bg_df <- NULL
  session_data$file_type <- NULL
  session_data$vars <- NULL

  shinyjs::reset("build_raster_file")
  shinyjs::reset("build_df_file")

  shinyjs::show("build_raster_file_box")
  shinyjs::show("build_df_file_box")

  removeNotification("raster_preview_msg")
  removeNotification("df_preview_msg")
})

observeEvent(input$build_raster_file, {
  req(input$build_raster_file)

  session_data$file_type <- "raster"
  shinyjs::hide("build_df_file_box")

  showNotification("Raster file(s) selected, preview loading...",
                   id = "raster_preview_msg",
                   type = "message",
                   duration = NULL)
})

observeEvent(input$build_df_file, {
  req(input$build_df_file)

  session_data$file_type <- "df"
  shinyjs::hide("build_raster_file_box")

  showNotification("CSV file selected, preview loading...",
                   id = "df_preview_msg",
                   type = "message",
                   duration = NULL)
})

observeEvent(input$build_data_upload_btn, {

  session_data$input_mode <- "bg_layers"
  session_data$bg_raster <- NULL

  if(is.null(session_data$file_type)){
    showNotification("No files selected. Please choose a file before uploading.",
                     type = "warning", duration = 4)
    return()
  }

  # Raster: g-space always available
  if(identical(session_data$file_type, "raster")){
    rast <- build_raster_upload()
    if(is.null(rast)) return()

    session_data$bg_raster <- rast
    session_data$bg_df <- terra::as.data.frame(rast, xy = TRUE, na.rm = TRUE)

    showNotification("Raster loaded successfully.", type = "message", duration = 4)
    updateTabsetPanel(session, "build_tabs", selected = "build_range_tab")
    return()
  }


  # CSV: g-space only if coordinates are present and form a regular grid
  if(identical(session_data$file_type, "df")){
    ext <- tolower(tools::file_ext(input$build_df_file$name))

    df <- tryCatch(
      load_df_file(input$build_df_file$datapath, ext),
      error = function(e){
        showNotification(paste("Could not load file:", e$message),
                         type = "error", duration = 4)
        NULL
      }
    )
    if(is.null(df)) return()

    session_data$bg_df <- df

    x_col <- names(df)[grepl(X_COL_PATTERN, names(df), ignore.case = TRUE)][1]
    y_col <- names(df)[grepl(Y_COL_PATTERN, names(df), ignore.case = TRUE)][1]

    if(is.na(x_col) || is.na(y_col)){
      showNotification(instructions$build_no_spatial_cols, type = "warning", duration = 8)
    } else {
      rast <- tryCatch(
        terra::rast(df[, c(x_col, y_col, setdiff(names(df), c(x_col, y_col)))],
                    type = "xyz"),
        error = function(e) NULL
      )

      if(is.null(rast)){
        showNotification(instructions$build_irregular_grid, type = "warning", duration = 8)
      } else {
        session_data$bg_raster <- rast
      }
    }

    showNotification("File loaded successfully.", type = "message", duration = 4)
    updateTabsetPanel(session, "build_tabs", selected = "build_range_tab")
    return()
  }

})

observeEvent(input$build_continue_virtual_btn, {
  session_data$input_mode <- "virtual"

  updateTabsetPanel(session, "build_tabs", selected = "build_range_tab")
})

observeEvent(input$build_continue_example_btn, {
  session_data$input_mode <- "example"

  session_data$bg_raster <- terra::rast(system.file("extdata", "ma_bios.tif",
                                                    package = "nicheR"))
  session_data$bg_df <- terra::as.data.frame(session_data$bg_raster,
                                             xy = TRUE, na.rm = TRUE)

  updateTabsetPanel(session, "build_tabs", selected = "build_range_tab")
})

observeEvent({
  lapply(seq_len(MAX_DIMS), function(i) input[[paste0("build_var_select_", i)]])
  lapply(seq_len(MAX_DIMS), function(i) input[[paste0("build_var_active_", i)]])
}, {

  vars <- available_vars()
  if(is.null(vars)) return()

  req(vars)

  n_slots <- min(length(vars), MAX_DIMS)

  active <- vapply(seq_len(n_slots), function(i)
    isTRUE(input[[paste0("build_var_active_", i)]]), logical(1))

  current <- vapply(seq_len(n_slots), function(i){
    val <- input[[paste0("build_var_select_", i)]]
    if(is.null(val)) vars[i] else val
  }, character(1))

  # only active slots block choices from others
  active_selections <- current[active]

  lapply(seq_len(n_slots), function(i){
    others <- if(active[i]) setdiff(active_selections, current[i]) else active_selections
    available <- c(current[i], setdiff(vars, others))
    updateSelectInput(session,
                      inputId = paste0("build_var_select_", i),
                      choices = available,
                      selected = current[i])
  })
}, ignoreNULL = TRUE, ignoreInit = TRUE)

observeEvent(input$build_confirm_variables_btn, {

  if(identical(session_data$input_mode, "virtual")){

    req(input$build_virtual_n_dims)
    n_dims <- input$build_virtual_n_dims

    vars <- vapply(seq_len(n_dims), function(i){
      val <- input[[paste0("build_virtual_var_name_", i)]]
      if(is.null(val) || !nzchar(val)) paste0("var", i) else val
    }, character(1))

    if(length(unique(vars)) != length(vars)){
      showNotification("Variable names must be unique.", type = "error", duration = 4)
      return()
    }

    session_data$vars <- vars

    showNotification(paste("Defined variables:", paste(vars, collapse = ", ")),
                     type = "message", duration = 4)

  } else {

    all_vars <- available_vars()
    if(is.null(all_vars)) return()

    n_slots  <- min(length(all_vars), MAX_DIMS)

    active <- vapply(seq_len(n_slots), function(i)
      isTRUE(input[[paste0("build_var_active_", i)]]), logical(1))

    vars <- vapply(seq_len(n_slots),
                   function(i){
                     val <- input[[paste0("build_var_select_", i)]]
                     if(is.null(val)) all_vars[i] else val
                   }, character(1))[active]

    session_data$vars <- vars

    showNotification(paste("Selected variables:", paste(vars, collapse = ", ")),
                     type = "message", duration = 4)
  }

  # Clear everything downstream that depended on the old variable set
  session_data$ellipsoid_list <- list()
  session_data$current_ellipsoid <- NULL

})

observeEvent(input$build_edit_variables_link, {

  # Nothing built yet, so there is nothing to lose. An unsaved working
  # ellipsoid counts as something built, it is not in the library yet.
  if(length(session_data$ellipsoid_list) == 0 &&
     is.null(session_data$current_ellipsoid)){
    session_data$vars <- NULL

    # Ranges
    session_data$session_range <- NULL
    session_data$df_range <- NULL

    # Ellipsoids
    session_data$current_ellipsoid <- NULL

    return()
  }

  showModal(modalDialog(
    title = "Edit variables?",
    p(instructions$build_edit_variables_tooltip, class = "text-instruction"),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("build_confirm_edit_variables_btn",
                   "Yes, edit variables",
                   class = "btn-warning")
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$build_confirm_edit_variables_btn, {
  removeModal()

  session_data$vars <- NULL

  session_data$session_range <- NULL
  session_data$df_range <- NULL

  session_data$ellipsoid_list <- list()
  session_data$current_ellipsoid <- NULL

  session_data$ellipsoid_prediction_list <-list()

  session_data$bias_raster <- NULL
  session_data$prepared_bias <- NULL
  session_data$ellipsoid_prediction_list_biased <- list()


  session_data$sampling_mask <- NULL
  session_data$ellipsoid_occurrence_list <- list()

})

# Reset logic if user changes from one input to the other of the input data
observeEvent(input$build_data_input_type_choice, {

  session_data$input_mode <- NULL
  session_data$file_type <- NULL

  session_data$bg_raster <- NULL
  session_data$bg_df <- NULL
  session_data$vars <- NULL

  session_data$session_range <- NULL
  session_data$df_range <- NULL

  session_data$ellipsoid_list <-list()
  session_data$current_ellipsoid <- NULL

  session_data$ellipsoid_prediction_list <-list()

  session_data$bias_raster <- NULL
  session_data$prepared_bias <- NULL
  session_data$ellipsoid_prediction_list_biased <- list()


  session_data$sampling_mask <- NULL
  session_data$ellipsoid_occurrence_list <- list()

  updateRadioButtons(session, "build_range_method_choice", selected = character(0))

}, ignoreInit = TRUE)

# Render Outputs ----------------------------------------------------------

output$build_data_input_type_ui <- renderUI({
  req(input$build_data_input_type_choice)

  switch(input$build_data_input_type_choice,
         "bg_layers" = {

           tagList(box(title = tags$span("Load Background Layers", class = "text-section-header"),
                       width = 12,
                       p(instructions$build_data_upload, class = "text-instruction"),
                       fluidRow(
                         column(
                           width = 12,
                           div(id = "build_raster_file_box",
                               fileInput(inputId = "build_raster_file",
                                         label = tagList(
                                           tags$span("Background layers (raster)", class = "text-widget-title"),
                                           tags$span(icon("circle-info"),
                                                     title = instructions$build_raster_file_tooltip,
                                                     class = "tooltip-icon")),
                                         multiple = TRUE,
                                         accept = c(".tif", ".tiff", ".rds"))),

                           div(id = "build_df_file_box",
                               fileInput(inputId = "build_df_file",
                                         label = tagList(
                                           tags$span("Background layers (CSV)", class = "text-widget-title"),
                                           tags$span(icon("circle-info"),
                                                     title = instructions$build_df_file_tooltip,
                                                     class = "tooltip-icon")),
                                         multiple = FALSE,
                                         accept = c("text/csv",
                                                    "text/comma-separated-values",
                                                    "text/plain",
                                                    ".csv", ".rds"))),

                           actionLink("build_clear_files",
                                      label = tagList(icon("rotate-left"), "Clear files")),

                           br(), br(),

                           verbatimTextOutput("build_raster_print"),
                           tableOutput("build_df_header")
                         )
                       ),

                       fluidRow(
                         column(
                           width = 12,
                           div(class = "action-btn-row",
                               actionButton(inputId = "build_data_upload_btn",
                                            label = "Upload",
                                            class = "btn-continue"))
                         )
                       )
           ))
         },

         "prev_session" = {

           tagList(box(title = tags$span("Load Previous Session", class = "text-section-header"),
                       width = 12,
                       p(instructions$build_prev_session, class = "text-instruction"),
                       fileInput(inputId = "build_session_file",
                                 label = tags$span("Session file (.rds)", class = "text-widget-title"),
                                 multiple = FALSE,
                                 accept = c(".rds")),
                       fluidRow(
                         column(
                           width = 12,
                           div(class = "action-btn-row",
                               actionButton(inputId = "build_load_session_btn",
                                            label = "Load Session",
                                            class = "btn-continue"))
                         )
                       )
           ))

         },
         "virtual_mode" = {

           tagList(box(title = tags$span("Virtual Mode", class = "text-section-header"),
                       width = 12,
                       p(instructions$build_virtual_mode, class = "text-instruction"),
                       fluidRow(
                         column(
                           width = 12,
                           div(class = "action-btn-row",
                               actionButton(inputId = "build_continue_virtual_btn",
                                            label = "Continue",
                                            class = "btn-continue"))
                         )
                       )
           ))

         },

         "example_data" = {

           tagList(
             box(title = tags$span("Example Data", class = "text-section-header"),
                 width = 12,
                 p(instructions$build_example_data, class = "text-instruction"),
                 fluidRow(
                   column(
                     width = 12,
                     div(class = "action-btn-row",
                         actionButton(inputId = "build_continue_example_btn",
                                      label = "Continue",
                                      class = "btn-continue"))
                   )
                 )
             )
           )

         }

  )
})

output$build_raster_print <- renderPrint({
  req(identical(session_data$file_type, "raster"))

  rast <- build_raster_upload()
  removeNotification("raster_preview_msg")
  req(rast)
  print(rast)
})

output$build_df_header <- renderTable({
  req(identical(session_data$file_type, "df"))
  req(input$build_df_file)

  ext <- tolower(tools::file_ext(input$build_df_file$name))
  result <- tryCatch(
    load_df_file(input$build_df_file$datapath, ext),
    error = function(e) stop(safeError(e))
  )

  removeNotification("df_preview_msg")
  head(result)
})

output$build_variable_selector_ui <- renderUI({

  # If variables are already confirmed, show a collapsed summary box
  if(!is.null(session_data$vars)){
    return(
      box(title = tags$span("Variables", class = "text-section-header"),
          id = "build_vars_box",
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          p(paste("Selected variables:", paste(session_data$vars, collapse = ", ")),
            class = "text-instruction"),
          # Always shown. Showing it only when the library has entries would
          # make this box redraw on every ellipsoid change, and a redraw
          # collapses it again right after Create new ellipsoid opens it.
          p("New ellipsoids use these variables. Changing them clears the library.",
            class = "text-muted-small"),
          fluidRow(
            column(12,
                   actionLink("build_edit_variables_link",
                              label = tagList(icon("pen"), "Edit variables")))
          )
      )
    )
  }

  # Virtual mode: ask for number of dimensions, then name each one
  if(identical(session_data$input_mode, "virtual")){

    n_dims <- if(!is.null(input$build_virtual_n_dims)) input$build_virtual_n_dims else 2

    name_rows <- lapply(seq_len(n_dims), function(i){
      fluidRow(
        column(width = 4, class = "var-label",
               tags$span(paste("Variable", i), class = "text-widget-inner")),
        column(width = 8,
               textInput(inputId = paste0("build_virtual_var_name_", i),
                         label = NULL,
                         value = paste0("var", i)))
      )
    })

    return(
      box(title = tags$span("Define Variables", class = "text-section-header"),
          id = "build_vars_box",
          width = 12,
          p(instructions$build_virtual_variables, class = "text-instruction"),
          fluidRow(
            column(width = 4, tags$span("Number of dimensions",
                                        class = "text-widget-title")),
            column(width = 4,
                   numericInput(inputId = "build_virtual_n_dims",
                                label = NULL,
                                value = n_dims,
                                min = 2, max = MAX_DIMS, step = 1))
          ),
          name_rows,
          fluidRow(
            column(12,
                   div(class = "action-btn-row",
                       actionButton("build_confirm_variables_btn",
                                    "Confirm",
                                    class = "btn-continue")))
          )
      )
    )
  }

  # Default: bg_layers / example_data / prev_session
  vars <- available_vars()
  if(is.null(vars)) return()

  n_slots <- min(length(vars), MAX_DIMS)
  all_choices <- vars

  var_slots <- lapply(seq_len(n_slots), function(i){
    fluidRow(
      column(width = 2, class = "checkbox-align",
             checkboxInput(paste0("build_var_active_", i),
                           label = NULL,
                           value = TRUE)),
      column(width = 10,
             conditionalPanel(paste0("input.build_var_active_", i, " == true"),
                              selectInput(paste0("build_var_select_", i),
                                          paste("Variable", i),
                                          all_choices,
                                          selected = vars[i])),
             conditionalPanel(paste0("input.build_var_active_", i, " == false"),
                              tags$div(class = "selector-disabled",
                                       selectInput(paste0("build_var_select_", i, "_ghost"),
                                                   paste("Variable", i),
                                                   all_choices,
                                                   selected = vars[i])))
      )
    )
  })

  box(title = tags$span("Select Variables", class = "text-section-header"),
      id = "build_vars_box",
      width = 12,
      p(instructions$build_variable_settings, class = "text-instruction"),
      var_slots,
      fluidRow(
        column(12,
               div(class = "action-btn-row",
                   actionButton("build_confirm_variables_btn",
                                "Confirm",
                                class = "btn-continue")))
      )
  )
})
