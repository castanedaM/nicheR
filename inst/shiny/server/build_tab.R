# Title: Build Tab Server Script

# Description: This script includes the code for the range, confidence level,
# covariance, and centroid mover logic, plus the ellipsoid library and the
# export of single ellipsoids. Functions from nicheR: build_ellipsoid(),
# ellipsoid_calculator(), update_ellipsoid_covariance(),
# update_ellipsoid_centroid(), ranges_from_data(), ranges_from_stats(),
# save_nicheR()

# Author: Mariana Castaneda-Guzman

# Date last updated: 09/15/2026


# DEBUG -------------------------------------------------------------------

# Which ellipsoid the covariance sliders currently on screen belong to.
# Set when they render, checked before any slider value is applied.
cov_owner <- reactiveVal(NULL)

# Debug message for covariance

# Which ellipsoid the centroid sliders currently on screen belong to
centroid_owner <- reactiveVal(NULL)

# Which ellipsoid the confidence level box currently on screen belongs to
cl_owner <- reactiveVal(NULL)

# Pairwise covariances as a flat named vector, in the same order as the
# sliders.
cov_upper <- function(ell){
  if(is.null(ell)) return(NULL)
  pairs <- t(combn(ell$var_names, 2))
  setNames(sapply(seq_len(nrow(pairs)),
                  function(i) ell$cov_matrix[pairs[i, 1], pairs[i, 2]]),
           apply(pairs, 1, paste, collapse = "-"))
}

# REACTIVE VALUES ---------------------------------------------------------

covariance_set <- reactiveVal(FALSE)
centroid_set <- reactiveVal(FALSE)
ell_mode <- reactiveVal("edit")

# Bumped whenever a different ellipsoid enters the working slot. The panels
# below depend on this instead of current_ellipsoid, so editing a slider
# does not tear down the panel holding that slider.
ell_slot <- reactiveVal(0L)

# Puts an ellipsoid into the working slot. Lineage lives on the ellipsoid
# itself as parent_id, so nothing else needs to be tracked here.
# Bumping ell_slot is what makes the isolated panels rebuild.
set_working_ellipsoid <- function(ell, mode = "edit"){

  session_data$current_ellipsoid <- ell
  ell_mode(mode)

  ell_slot(ell_slot() + 1L)

  covariance_set(FALSE)
  centroid_set(FALSE)
  range_dirty(FALSE)
}

# Empties the working slot. Separate from the above so every caller goes
# through one of the two and ell_slot can never be missed.
clear_working_ellipsoid <- function(){

  session_data$current_ellipsoid <- NULL
  ell_mode("edit")

  ell_slot(ell_slot() + 1L)

  covariance_set(FALSE)
  centroid_set(FALSE)
  range_dirty(FALSE)
}

carry_ell_meta <- function(new_ell, old_ell){
  new_ell$ell_id <- old_ell$ell_id
  new_ell$ell_name <- old_ell$ell_name
  new_ell$parent_id <- old_ell$parent_id
  new_ell$range_method <- old_ell$range_method
  new_ell$range_inputs <- old_ell$range_inputs
  new_ell
}

# TRUE once the user edits a range input after an ellipsoid is in the slot.
# Range lines follow the live inputs while this is TRUE, so edits preview
# without being committed until Rebuild.
range_dirty <- reactiveVal(FALSE)

# RANGES ------------------------------------------------------------------

output$build_range_method_choice_ui <- renderUI({

  req(session_data$vars)

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  is_view <- identical(ell_mode(), "view")

  # View mode, ranges and cl are read-only
  if(is_view && !is.null(ell)){

    range_rows <- lapply(ell$var_names, function(v){
      fluidRow(
        column(width = 4, tags$span(v, class = "text-widget-inner")),
        column(width = 4, tags$span(round(ell$ranges["min", v], 2),
                                    class = "text-widget-inner text-center")),
        column(width = 4, tags$span(round(ell$ranges["max", v], 2),
                                    class = "text-widget-inner text-center"))
      )
    })

    return(
      box(title = tags$span("Ranges", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          p(instructions$build_view_only, class = "text-instruction"),
          fluidRow(
            column(width = 4, tags$span("Variable", class = "text-widget-title text-center")),
            column(width = 4, tags$span("Min", class = "text-widget-title text-center")),
            column(width = 4, tags$span("Max", class = "text-widget-title text-center"))
          ),
          range_rows,
          tags$hr(style = "margin: 8px 0;"),
          fluidRow(
            column(width = 5, tags$span("Ellipsoid confidence level",
                                        class = "text-widget-title")),
            column(width = 4, tags$span(ell$cl,
                                        class = "text-widget-inner text-center"))
          )
      )
    )
  }

  # cl belongs to the ellipsoid, not to a range method, so it lives here
  # rather than inside build_range_method_ui. That panel only exists once a
  # method is picked, and the radio comes back unselected every time a new
  # ellipsoid enters the slot.
  cl_owner(if(is.null(ell)) NULL else ell$ell_id)

  cl_row <- fluidRow(
    column(width = 5,
           tags$div(class = "tooltip-label-row",
                    tags$span("Ellipsoid confidence level",
                              class = "text-widget-title"),
                    tags$span(icon("circle-info"),
                              title = instructions$build_cl_ell_tooltip,
                              class = "tooltip-icon"))),
    column(width = 4,
           numericInput(inputId = "build_cl_ell",
                        label = NULL,
                        value = if(is.null(ell)) 0.95 else ell$cl,
                        min = 0.5, max = 0.999,
                        step = 0.01))
  )

  box(title = tags$span("Ranges", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = !is.null(ell),
      p(instructions$build_range_choice, class = "text-instruction"),
      tags$hr(style = "margin: 8px 0;"),
      cl_row,
      br(),
      radioButtons("build_range_method_choice",
                   label = tagList(tags$span("Select how to define the ranges for your ellipsoid ",
                                             class = "text-widget-title"),
                                   tags$span(icon("circle-info"),
                                             title = paste0("Manual: ",
                                                            instructions$build_range_manual_tooltip, "\n",
                                                            "Stats: ",
                                                            instructions$build_range_stats_tooltip, "\n",
                                                            "Data: ",
                                                            instructions$build_range_data_tooltip),
                                             class = "tooltip-icon")),
                   choiceNames = list(tags$span("Manual", class = "text-widget-inner"),
                                      tags$span("From Data", class = "text-widget-inner"),
                                      tags$span("From Stats", class = "text-widget-inner")),
                   choiceValues = c("man", "df", "stats"),
                   selected = character(0)),
      uiOutput("build_range_method_ui")
  )
})

output$build_range_method_ui <- renderUI({

  req(input$build_range_method_choice, session_data$vars)

  vars <- session_data$vars

  # No background data in virtual mode, defaults fall back to a standard scale
  has_bg_df <- !is.null(session_data$bg_df)

  # Isolated on purpose. Reading current_ellipsoid reactively rebuilds this
  # whole panel on every covariance and centroid slider move, which resets
  # the range inputs to their defaults mid-edit. The parent already tears
  # this panel down whenever the slot changes, so the label stays correct.
  build_label <- if(is.null(isolate(session_data$current_ellipsoid))){
    "Initialize Ellipsoid"
  } else {
    "Rebuild Ellipsoid"
  }

  build_btn <- fluidRow(
    column(width = 12,
           div(class = "action-btn-row",
               actionButton("build_init_ell_btn",
                            build_label,
                            class = "btn-continue"))
    )
  )

  reset_btn <- fluidRow(
    column(width = 12,
           actionLink("build_reset_ranges_link",
                      label = tagList(icon("rotate-left"), "Reset to defaults")))
  )

  dft_values <- build_range_defaults(vars, has_bg_df)

  switch(input$build_range_method_choice,

         "man" = {
           header <- fluidRow(
             column(width = 12,
                    p(instructions$build_range_manual,
                      class = "text-instruction")),
             column(width = 4,
                    tags$span("Variable", class = "text-widget-title text-center")),
             column(width = 4,
                    tags$span("Min", class = "text-widget-title text-center")),
             column(width = 4,
                    tags$span("Max", class = "text-widget-title text-center"))
           )

           var_rows <- lapply(vars, function(v){
             fluidRow(
               column(width = 4,
                      class = "var-label",
                      tags$span(v, class = "text-widget-inner")),
               column(width = 4,
                      numericInput(inputId = paste0("build_min_", v),
                                   label = NULL,
                                   value = dft_values$min[[v]],
                                   step = 0.5)),
               column(width = 4,
                      numericInput(inputId = paste0("build_max_", v),
                                   label = NULL,
                                   value = dft_values$max[[v]],
                                   step = 0.5))
             )
           })

           column(width = 12, header, var_rows, reset_btn, br(), build_btn)
         },

         "df" = {

           upload_row <- fluidRow(
             column(width = 12,
                    p(instructions$build_range_data, class = "text-instruction")),
             column(width = 12,
                    fileInput(inputId = "build_df_range_file",
                              label = tags$span("Choose CSV file with range data",
                                                class = "text-widget-title"),
                              multiple = FALSE,
                              accept = c("text/csv",
                                         "text/comma-separated-values",
                                         "text/plain",
                                         ".csv", ".rds")))
           )

           # The file is parsed by its own observer, this only reads the result
           obs_ranges <- NULL
           shared_vars <- character(0)

           if(!is.null(session_data$df_range)){
             shared_vars <- intersect(vars, colnames(session_data$df_range))
             obs_ranges <- lapply(shared_vars, function(v){
               list(min = round(min(session_data$df_range[, v], na.rm = TRUE), 2),
                    max = round(max(session_data$df_range[, v], na.rm = TRUE), 2))
             })
             names(obs_ranges) <- shared_vars
           }

           header1 <- fluidRow(
             column(width = 4,
                    tags$span("Variable", class = "text-widget-title text-center")),
             column(width = 4,
                    tags$span("Observed min", class = "text-widget-title text-center")),
             column(width = 4,
                    tags$span("Observed max", class = "text-widget-title text-center"))
           )

           header2 <- fluidRow(
             column(width = 4,
                    tags$span("Variable", class = "text-widget-title text-center")),
             column(width = 4,
                    tags$div(class = "tooltip-label-row",
                             tags$span("Expand min (%)",
                                       class = "text-widget-title text-center"),
                             tags$span(icon("circle-info"),
                                       title = instructions$build_expand_range_tooltip,
                                       class = "tooltip-icon"))),
             column(width = 4,
                    tags$div(class = "tooltip-label-row",
                             tags$span("Expand max (%)",
                                       class = "text-widget-title text-center"),
                             tags$span(icon("circle-info"),
                                       title = instructions$build_expand_range_tooltip,
                                       class = "tooltip-icon")))
           )

           rows1 <- lapply(shared_vars, function(v){
             fluidRow(
               column(width = 4, class = "var-label",
                      tags$span(v, class = "text-widget-inner")),
               column(width = 4,
                      tags$span(format(obs_ranges[[v]]$min, nsmall = 2),
                                class = "text-widget-inner text-center")),
               column(width = 4,
                      tags$span(format(obs_ranges[[v]]$max, nsmall = 2),
                                class = "text-widget-inner text-center")))
           })

           rows2 <- lapply(shared_vars, function(v){
             fluidRow(
               column(width = 4, class = "var-label",
                      tags$span(v, class = "text-widget-inner")),
               column(width = 4,
                      numericInput(inputId = paste0("build_expand_min_df_", v),
                                   label = NULL,
                                   value = dft_values$expand_min[[v]],
                                   step = 5)),
               column(width = 4,
                      numericInput(inputId = paste0("build_expand_max_df_", v),
                                   label = NULL,
                                   value = dft_values$expand_max[[v]],
                                   step = 5))
             )
           })

           range_rows <- if(length(shared_vars) > 0){
             tagList(header1, rows1, br(), header2, rows2)
           } else {
             NULL
           }

           column(width = 12, upload_row, range_rows, reset_btn, br(), build_btn)
         },

         "stats" = {

           # Stats needs its own cl. This one is the normal interval that
           # becomes the ranges, it is not the ellipsoid cutoff, and editing
           # it is a range edit like any other.
           cl_row <- fluidRow(
             column(width = 5,
                    tags$div(class = "tooltip-label-row",
                             tags$span("Range confidence level for range",
                                       class = "text-widget-title text-center"),
                             tags$span(icon("circle-info"),
                                       title = instructions$build_cl_interval_tooltip,
                                       class = "tooltip-icon"))),
             column(width = 4,
                    numericInput(inputId = "build_range_cl_stats",
                                 label = NULL,
                                 value = 0.95,
                                 min = 0, max = 1,
                                 step = 0.05))
           )

           header1 <- fluidRow(
             column(width = 12,
                    p(instructions$build_range_stats, class = "text-instruction")),
             column(width = 4,
                    tags$span("Variable", class = "text-widget-title text-center")),
             column(width = 4,
                    tags$span("Mean", class = "text-widget-title text-center")),
             column(width = 4,
                    tags$span("SD", class = "text-widget-title text-center"))
           )

           header2 <- fluidRow(
             column(width = 4,
                    tags$span("Variable", class = "text-widget-title text-center")),
             column(width = 4,
                    tags$div(class = "tooltip-label-row",
                             tags$span("Expand min (%)",
                                       class = "text-widget-title text-center"),
                             tags$span(icon("circle-info"),
                                       title = instructions$build_expand_range_tooltip,
                                       class = "tooltip-icon"))),
             column(width = 4,
                    tags$div(class = "tooltip-label-row",
                             tags$span("Expand max (%)",
                                       class = "text-widget-title text-center"),
                             tags$span(icon("circle-info"),
                                       title = instructions$build_expand_range_tooltip,
                                       class = "tooltip-icon")))
           )

           var_rows1 <- lapply(vars, function(v){
             fluidRow(
               column(width = 4, class = "var-label",
                      tags$span(v, class = "text-widget-inner")),
               column(width = 4,
                      numericInput(inputId = paste0("build_mean_", v),
                                   label = NULL,
                                   value = dft_values$mean[[v]],
                                   step = 0.5)),
               column(width = 4,
                      numericInput(inputId = paste0("build_sd_", v),
                                   label = NULL,
                                   value = dft_values$sd[[v]],
                                   step = 0.5))
             )
           })

           var_rows2 <- lapply(vars, function(v){
             fluidRow(
               column(width = 4, class = "var-label",
                      tags$span(v, class = "text-widget-inner")),
               column(width = 4,
                      numericInput(inputId = paste0("build_expand_min_stats_", v),
                                   label = NULL,
                                   value = dft_values$expand_min[[v]],
                                   step = 5)),
               column(width = 4,
                      numericInput(inputId = paste0("build_expand_max_stats_", v),
                                   label = NULL,
                                   value = dft_values$expand_max[[v]],
                                   step = 5))
             )
           })

           column(width = 12, header1, var_rows1, header2, var_rows2,
                  cl_row, reset_btn, br(), build_btn)
         }
  )
})

# Parses the uploaded range CSV once, outside any render function
observeEvent(input$build_df_range_file, {

  req(session_data$vars)

  ext <- tolower(tools::file_ext(input$build_df_range_file$name))

  df_range <- tryCatch(
    load_df_file(input$build_df_range_file$datapath, ext),
    error = function(e) NULL
  )

  if(is.null(df_range)){
    session_data$df_range <- NULL
    showNotification(instructions$build_range_file_unreadable,
                     type = "error", duration = 6)
    return()
  }

  shared_vars <- intersect(session_data$vars, colnames(df_range))

  if(length(shared_vars) != length(session_data$vars)){
    session_data$df_range <- NULL
    showNotification(instructions$build_range_file_mismatch,
                     type = "warning", duration = 8)
    return()
  }

  session_data$df_range <- df_range
})

# Collects the current range inputs and turns them into min/max bounds.
# Returns both the computed bounds and the raw inputs, so the raw inputs
# can be stored on the ellipsoid and restored when it is loaded again.
range_preview <- reactive({

  if(is.null(input$build_range_method_choice) || is.null(session_data$vars)){
    return(NULL)
  }

  vars <- session_data$vars
  zeros <- setNames(as.list(rep(0, length(vars))), vars)

  switch(input$build_range_method_choice,

         "man" = {
           mins <- setNames(sapply(vars, function(v) input[[paste0("build_min_", v)]]), vars)
           maxs <- setNames(sapply(vars, function(v) input[[paste0("build_max_", v)]]), vars)

           # Inputs may not have rendered yet, bail out instead of comparing
           # incompatible types
           if(!is.numeric(mins) || !is.numeric(maxs)) return(NULL)
           if(length(mins) != length(vars) || length(maxs) != length(vars)) return(NULL)
           if(any(is.na(mins)) || any(is.na(maxs))) return(NULL)
           if(any(maxs <= mins)) return(NULL)

           list(mins = mins,
                maxs = maxs,
                inputs = list(method = "man",
                              min = as.list(mins),
                              max = as.list(maxs),
                              mean = zeros,
                              sd = zeros,
                              expand_min = zeros,
                              expand_max = zeros))
         },

         "df" = {

           df_range <- session_data$df_range
           if(is.null(df_range)) return(NULL)

           expand_min <- setNames(sapply(vars, function(v){
             val <- input[[paste0("build_expand_min_df_", v)]]
             if(is.null(val) || is.na(val)) 0 else val
           }), vars)

           expand_max <- setNames(sapply(vars, function(v){
             val <- input[[paste0("build_expand_max_df_", v)]]
             if(is.null(val) || is.na(val)) 0 else val
           }), vars)

           ranges_df <- ranges_from_data(data = df_range[, vars, drop = FALSE],
                                         expand_min = as.list(expand_min),
                                         expand_max = as.list(expand_max))

           list(mins = ranges_df["min", ],
                maxs = ranges_df["max", ],
                inputs = list(method = "df",
                              min = as.list(ranges_df["min", ]),
                              max = as.list(ranges_df["max", ]),
                              mean = zeros,
                              sd = zeros,
                              expand_min = as.list(expand_min),
                              expand_max = as.list(expand_max)))
         },

         "stats" = {

           means <- setNames(sapply(vars, function(v) input[[paste0("build_mean_", v)]]), vars)
           sds <- setNames(sapply(vars, function(v) input[[paste0("build_sd_", v)]]), vars)

           if(!is.numeric(means) || !is.numeric(sds)) return(NULL)
           if(any(is.na(means)) || any(is.na(sds))) return(NULL)
           if(any(sds <= 0)) return(NULL)

           expand_min <- setNames(sapply(vars, function(v){
             val <- input[[paste0("build_expand_min_stats_", v)]]
             if(is.null(val) || is.na(val)) 0 else val
           }), vars)

           expand_max <- setNames(sapply(vars, function(v){
             val <- input[[paste0("build_expand_max_stats_", v)]]
             if(is.null(val) || is.na(val)) 0 else val
           }), vars)

           cl <- if(!is.null(input$build_range_cl_stats)) input$build_range_cl_stats else 0.95

           range_stats <- ranges_from_stats(mean = means, sd = sds, cl = cl,
                                            expand_min = as.list(expand_min),
                                            expand_max = as.list(expand_max))

           list(mins = range_stats["min", ],
                maxs = range_stats["max", ],
                inputs = list(method = "stats",
                              min = as.list(range_stats["min", ]),
                              max = as.list(range_stats["max", ]),
                              mean = as.list(means),
                              sd = as.list(sds),
                              expand_min = as.list(expand_min),
                              expand_max = as.list(expand_max),
                              cl = cl))
         }
  )
})

# Marks the range inputs as edited so the plot previews them. Ignores the
# initial render, which fires once as the widgets are created.
observeEvent({
  req(session_data$vars)
  vars <- session_data$vars
  list(
    input$build_range_method_choice,
    input$build_range_cl_stats,
    lapply(vars, function(v) input[[paste0("build_min_", v)]]),
    lapply(vars, function(v) input[[paste0("build_max_", v)]]),
    lapply(vars, function(v) input[[paste0("build_mean_", v)]]),
    lapply(vars, function(v) input[[paste0("build_sd_", v)]]),
    lapply(vars, function(v) input[[paste0("build_expand_min_df_", v)]]),
    lapply(vars, function(v) input[[paste0("build_expand_max_df_", v)]]),
    lapply(vars, function(v) input[[paste0("build_expand_min_stats_", v)]]),
    lapply(vars, function(v) input[[paste0("build_expand_max_stats_", v)]])
  )
}, {
  if(!is.null(session_data$current_ellipsoid)) range_dirty(TRUE)
}, ignoreInit = TRUE)

observeEvent(input$build_reset_ranges_link, {

  req(input$build_range_method_choice, session_data$vars)

  vars <- session_data$vars

  dft_values <- build_range_defaults(vars, !is.null(session_data$bg_df))

  switch(input$build_range_method_choice,

         "man" = {
           lapply(vars, function(v){
             updateNumericInput(session, paste0("build_min_", v),
                                value = dft_values$min[[v]])
             updateNumericInput(session, paste0("build_max_", v),
                                value = dft_values$max[[v]])
           })
         },

         "df" = {
           lapply(vars, function(v){
             updateNumericInput(session, paste0("build_expand_min_df_", v),
                                value = dft_values$expand_min[[v]])
             updateNumericInput(session, paste0("build_expand_max_df_", v),
                                value = dft_values$expand_max[[v]])
           })
         },

         "stats" = {
           lapply(vars, function(v){
             updateNumericInput(session, paste0("build_mean_", v),
                                value = dft_values$mean[[v]])
             updateNumericInput(session, paste0("build_sd_", v),
                                value = dft_values$sd[[v]])
             updateNumericInput(session, paste0("build_expand_min_stats_", v),
                                value = dft_values$expand_min[[v]])
             updateNumericInput(session, paste0("build_expand_max_stats_", v),
                                value = dft_values$expand_max[[v]])
           })
         }
  )

  # Only the interval cl is reset here. The ellipsoid cl is not a range
  # input, so Reset to defaults leaves it alone.
  updateNumericInput(session, "build_range_cl_stats", value = 0.95)
})


# CONFIDENCE LEVEL --------------------------------------------------------

# cl only rescales the boundary. Sigma and the centroid go straight back to
# ellipsoid_calculator(), which is why this is applied live and never marks
# the ranges dirty.
observeEvent(input$build_cl_ell, {

  ell <- session_data$current_ellipsoid
  req(ell)
  req(!identical(ell_mode(), "view"))

  # The box on screen still belongs to a previous ellipsoid
  if(!identical(cl_owner(), ell$ell_id)){
    return()
  }

  cl <- input$build_cl_ell

  if(is.null(cl) || is.na(cl) || cl <= 0 || cl >= 1){
    return()
  }

  if(abs(cl - ell$cl) < 1e-8){
    return()
  }

  updated_ell <- tryCatch(
    ellipsoid_calculator(cov_matrix = ell$cov_matrix,
                         centroid = ell$centroid,
                         cl = cl,
                         verbose = FALSE),
    error = function(e){
      showNotification(paste("Confidence level update failed:", e$message),
                       type = "error", duration = 5)
      NULL
    }
  )

  req(updated_ell)

  # ellipsoid_calculator() does not carry ranges, build_ellipsoid() assigns
  # them afterwards. sigma is (max - min) / 6 and does not depend on cl, so
  # the ranges this ellipsoid was built from are still the right ones.
  updated_ell$ranges <- ell$ranges

  session_data$current_ellipsoid <- carry_ell_meta(updated_ell, ell)

}, ignoreNULL = TRUE, ignoreInit = TRUE)


# ELLIPSOID BUILD ---------------------------------------------------------

observeEvent(input$build_init_ell_btn, {

  req(session_data$vars)

  ranges <- range_preview()

  if(is.null(ranges)){
    showNotification(instructions$build_range_invalid,
                     type = "error", duration = 5)
    return()
  }

  vars <- session_data$vars
  range_df <- as.data.frame(rbind(unlist(ranges$mins), unlist(ranges$maxs)),
                            row.names = c("min", "max"))
  colnames(range_df) <- vars

  # The cl box is always on screen, so whatever it shows is what gets built.
  # A rebuild keeps the cl the user last set instead of dropping back to a
  # default, because the box is redrawn from the working ellipsoid.
  cl <- input$build_cl_ell

  if(is.null(cl) || is.na(cl) || cl <= 0 || cl >= 1){
    cl <- 0.95
  }

  new_ell <- tryCatch(
    build_ellipsoid(range = range_df, cl = cl, verbose = FALSE),
    error = function(e){
      showNotification(paste("Error building ellipsoid:", e$message),
                       type = "error", duration = 6)
      NULL
    }
  )

  req(new_ell)

  # Keep the identity if we are rebuilding a version that is already saved,
  # so the library shows Update rather than Save. Otherwise tag a new one.
  cur <- session_data$current_ellipsoid
  is_saved <- !is.null(cur) && cur$ell_id %in% names(session_data$ellipsoid_list)

  if(is_saved){
    new_ell$ell_id <- cur$ell_id
    new_ell$ell_name <- cur$ell_name
  } else {
    new_ell <- tag_ellipsoid(new_ell)
  }

  # Store the range inputs so this ellipsoid can be reloaded and edited later
  new_ell$range_method <- ranges$inputs$method
  new_ell$range_inputs <- ranges$inputs

  session_data$session_range <- ranges$inputs

  # New ranges mean new geometry, so any inherited lineage no longer
  # describes this ellipsoid. It becomes its own root.
  new_ell$parent_id <- NULL

  set_working_ellipsoid(new_ell, mode = "edit")

  showNotification("Ellipsoid built successfully.",
                   type = "message", duration = 4)
})


# ELLIPSOID LIBRARY -------------------------------------------------------

output$build_ellipsoid_library_ui <- renderUI({

  cur_ell <- session_data$current_ellipsoid
  versions <- session_data$ellipsoid_list
  ids <- names(versions)

  req(!is.null(cur_ell) || length(ids) > 0)

  is_view <- identical(ell_mode(), "view")
  is_saved <- !is.null(cur_ell) && !is.null(cur_ell$ell_id) &&
    cur_ell$ell_id %in% ids

  # Working slot, the ellipsoid every plot and every downstream tab uses
  working_row <- if(!is.null(cur_ell)){
    fluidRow(
      class = "ell-row",
      style = "background: #f0f7f0; border-radius: 4px; margin-bottom: 6px; padding: 4px 0;",
      column(width = 4,
             tags$span(icon(if(is_view) "eye" else "pen"),
                       tags$span(paste0(" ", cur_ell$ell_name),
                                 class = "text-widget-inner",
                                 style = "color: #097a21; font-weight: 500;"),
                       tags$span(if(is_view) "(viewing)" else "(editing)",
                                 style = "font-size: 11px; color: #aaa;"))),
      column(width = 4,
             tags$span(ell_lineage_label(cur_ell),
                       style = "font-size: 11px; color: #aaa;")),
      column(width = 4,
             if(is_view){
               actionButton("build_exit_view_btn",
                            "Edit",
                            class = "btn-continue")
             } else if(is_saved){
               actionButton("build_update_ell_btn",
                            "Update",
                            class = "btn-continue")
             } else {
               actionButton("build_save_ell_btn",
                            "Save",
                            class = "btn-continue")
             })
    )
  }

  rows <- lapply(ids, function(id){

    ell <- versions[[id]]

    fluidRow(
      class = "ell-row",
      style = "padding: 2px 0;",
      column(width = 4,
             tags$span(ell$ell_name, class = "text-widget-inner"),
             tags$br(),
             tags$span(id, style = "font-size: 10px; color: #bbb;")),
      column(width = 4,
             tags$span(ell_lineage_label(ell),
                       style = "font-size: 11px; color: #aaa;")),
      column(width = 4,
             class = "ell-actions",
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('build_ell_view', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("View ", ell$ell_name, " (read-only)"),
                    icon("eye")),
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('build_ell_edit', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("Edit ", ell$ell_name),
                    icon("pen-to-square")),
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('build_ell_add', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("New copy from ", ell$ell_name),
                    icon("plus")),
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('build_ell_export', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("Export ", ell$ell_name),
                    icon("download")),
             tags$a(href = "#",
                    class = "ell-action-danger",
                    onclick = sprintf("Shiny.setInputValue('build_ell_delete', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("Delete ", ell$ell_name),
                    icon("trash-can")))
    )
  })

  box(title = tags$span("Ellipsoid library", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$build_library, class = "text-instruction"),

      if(!is.null(cur_ell)){
        tagList(working_row, tags$hr(style = "margin: 8px 0;"))
      },

      if(length(ids) > 0){
        tagList(
          fluidRow(
            class = "ell-row",
            style = "padding: 2px 0;",
            column(width = 4, tags$span("Name", class = "text-widget-title text-center")),
            column(width = 4, tags$span("Details", class = "text-widget-title text-center")),
            column(width = 4, tags$span("Actions", class = "text-widget-title text-center"))
          ),
          tagList(rows)
        )
      } else {
        p(instructions$build_library_empty, class = "text-muted-small")
      }
  )
})

# View, read-only
observeEvent(input$build_ell_view, {

  ell <- session_data$ellipsoid_list[[input$build_ell_view]]
  req(ell)

  set_working_ellipsoid(ell, mode = "view")

  showNotification(paste0("Viewing ", ell$ell_name, ". Editing is locked."),
                   type = "message", duration = 3)
})

observeEvent(input$build_exit_view_btn, {

  req(session_data$current_ellipsoid)

  ell_mode("edit")
  ell_slot(ell_slot() + 1L)

  showNotification("Editing unlocked.", type = "message", duration = 3)
})

# Edit, loads a saved version back into the working slot
observeEvent(input$build_ell_edit, {

  ell <- session_data$ellipsoid_list[[input$build_ell_edit]]
  req(ell)

  # Range panel is restored before the slot changes, so the re-render
  # triggered by ell_slot already sees the right values
  if(!is.null(ell$range_inputs)){
    session_data$session_range <- ell$range_inputs
  }

  set_working_ellipsoid(ell, mode = "edit")

  showNotification(paste0(ell$ell_name, " loaded for editing."),
                   type = "message", duration = 3)
})

# Add, creates a new unsaved copy from any existing version
observeEvent(input$build_ell_add, {

  ell <- session_data$ellipsoid_list[[input$build_ell_add]]
  req(ell)

  new_ell <- tag_ellipsoid(ell)
  new_ell$parent_id <- ell$ell_id

  if(!is.null(ell$range_inputs)){
    session_data$session_range <- ell$range_inputs
  }

  set_working_ellipsoid(new_ell, mode = "edit")

  showNotification(paste0("New copy created from ", ell$ell_name, "."),
                   type = "message", duration = 3)
})

# Delete, asks first
observeEvent(input$build_ell_delete, {

  ell <- session_data$ellipsoid_list[[input$build_ell_delete]]
  req(ell)

  session_data$pending_ell_delete <- input$build_ell_delete

  n_children <- sum(vapply(session_data$ellipsoid_list, function(e){
    identical(e$parent_id, ell$ell_id)
  }, logical(1)))

  showModal(modalDialog(
    title = paste0("Delete ", ell$ell_name, "?"),
    p(instructions$build_delete_ell, class = "text-instruction"),
    if(n_children > 0){
      p(paste0(n_children, " ellipsoid(s) were copied from this one. ",
               "They will be kept, but will no longer have a parent."),
        class = "text-muted-small")
    },
    footer = tagList(
      modalButton("Cancel"),
      actionButton("build_confirm_ell_delete_btn",
                   "Yes, delete",
                   class = "btn-cancel")
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$build_confirm_ell_delete_btn, {

  id <- session_data$pending_ell_delete
  req(id)

  nm <- session_data$ellipsoid_list[[id]]$ell_name

  removeModal()

  session_data$ellipsoid_list[[id]] <- NULL
  session_data$ellipsoid_prediction_list[[id]] <- NULL
  session_data$ellipsoid_prediction_list_biased[[id]] <- NULL
  session_data$pending_ell_delete <- NULL

  # Copies of the deleted ellipsoid, captured before reparenting so the
  # message reports only what this delete changed
  orphaned <- names(session_data$ellipsoid_list)[
    vapply(session_data$ellipsoid_list,
           function(e) identical(e$parent_id, id), logical(1))]

  # orphans to say root
  session_data$ellipsoid_list <- lapply(session_data$ellipsoid_list, function(e){
    if(identical(e$parent_id, id)) e$parent_id <- NULL
    e
  })

  cur <- session_data$current_ellipsoid


  if(identical(cur$ell_id, id)){
    clear_working_ellipsoid()
    showNotification(paste0(nm, " deleted. Set ranges to build a new ellipsoid."),
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

# Save, names a new version
observeEvent(input$build_save_ell_btn, {

  req(session_data$current_ellipsoid)

  default_name <- session_data$current_ellipsoid$ell_name

  showModal(modalDialog(
    title = "Save ellipsoid version",
    p(instructions$build_save_ell_modal, class = "text-instruction"),
    textInput("build_ell_save_name",
              label = NULL,
              value = default_name,
              placeholder = default_name),
    tags$small("Max 30 characters.", class = "text-muted-small"),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("build_confirm_save_ell_btn",
                   "Save",
                   class = "btn-save")
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$build_confirm_save_ell_btn, {

  req(session_data$current_ellipsoid)

  raw_name <- if(!is.null(input$build_ell_save_name) &&
                 nzchar(trimws(input$build_ell_save_name))){
    input$build_ell_save_name
  } else {
    session_data$current_ellipsoid$ell_name
  }

  clean_name <- gsub("\\s+", "_", trimws(raw_name))
  clean_name <- substr(clean_name, 1, 30)

  ell_to_save <- session_data$current_ellipsoid
  ell_to_save$ell_name <- clean_name

  session_data$ellipsoid_list[[ell_to_save$ell_id]] <- ell_to_save
  session_data$current_ellipsoid <- ell_to_save

  removeModal()

  showModal(modalDialog(
    title = paste0(clean_name, " saved."),
    footer = tagList(
      div(class = "action-btn-row",
          actionButton("build_next_done_btn",
                       label = "Done",
                       class = "btn-save"))
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$build_next_done_btn, {
  removeModal()
})

# Update, overwrites a version already in the library
observeEvent(input$build_update_ell_btn, {

  ell <- session_data$current_ellipsoid
  req(ell)

  session_data$ellipsoid_list[[ell$ell_id]] <- ell

  showNotification(paste0(ell$ell_name, " updated."),
                   type = "message", duration = 4)
})

# Button and UI to move to the next step
output$build_next_step_ui <- renderUI({

  req(length(session_data$ellipsoid_list) > 0)

  div(class = "action-btn-row",
      actionButton(inputId = "build_next_step_btn",
                   label = tagList("Continue",
                                   icon("arrow-right")),
                   class = "btn-save")
  )
})

observeEvent(input$build_next_step_btn, {
  target <- if(identical(session_data$input_mode, "virtual")){
    "generate_tab"
  } else {
    "predict_tab"
  }

  updateTabItems(session, "sidebar_menu", selected = target)
})


# ELLIPSOID EXPORT --------------------------------------------------------

# Default file name is the internal id plus today, so an exported file can
# always be traced back to the library entry it came from. The user can
# overwrite it in the modal.
ell_export_name <- function(ell){
  paste0(ell$ell_id, "_", format(Sys.Date(), "%Y%m%d"))
}

# The working ellipsoid may not be in the library yet, so fall back to it
ell_export_target <- function(id){
  ell <- session_data$ellipsoid_list[[id]]
  cur <- session_data$current_ellipsoid
  if(is.null(ell) && !is.null(cur) && identical(cur$ell_id, id)){
    ell <- cur
  }
  ell
}

observeEvent(input$build_ell_export, {

  ell <- ell_export_target(input$build_ell_export)
  req(ell)

  session_data$pending_ell_export <- input$build_ell_export

  showModal(modalDialog(
    title = paste0("Export ", ell$ell_name),
    p(instructions$build_export_ell, class = "text-instruction"),
    textInput("build_ell_export_name",
              label = NULL,
              value = ell_export_name(ell),
              placeholder = ell_export_name(ell)),
    tags$small("Saved as .rds. Letters, numbers, dashes and underscores only.",
               class = "text-muted-small"),
    footer = tagList(
      modalButton("Close"),
      downloadButton("build_ell_export_dl", "Download", class = "btn-save")
    ),
    easyClose = FALSE
  ))
})

output$build_ell_export_dl <- downloadHandler(

  filename = function(){
    ell <- ell_export_target(session_data$pending_ell_export)
    nm <- input$build_ell_export_name
    if(is.null(nm) || !nzchar(trimws(nm))){
      nm <- ell_export_name(ell)
    }
    nm <- gsub("\\s+", "_", trimws(nm))
    nm <- gsub("[^A-Za-z0-9._-]", "", nm)
    nm <- sub("\\.rds$", "", nm, ignore.case = TRUE)
    paste0(substr(nm, 1, 50), ".rds")
  },

  content = function(file){

    ell <- ell_export_target(session_data$pending_ell_export)

    if(is.null(ell)){
      stop("No ellipsoid available to export.")
    }

    # save_nicheR() appends .rds when the path does not end in it, and the
    # temp path Shiny provides has no extension. Writing straight to file
    # would leave the download empty, so write beside it and rename.
    tmp <- paste0(file, ".rds")
    save_nicheR(object = ell, file = tmp, overwrite = TRUE)
    file.rename(tmp, file)
  }
)


# COVARIANCE --------------------------------------------------------------

output$build_covariance_ui <- renderUI({

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  # View mode, show the matrix without controls
  if(identical(ell_mode(), "view")){

    pairs <- t(combn(ell$var_names, 2))
    pair_names <- apply(pairs, 1, function(p) paste(p, collapse = "-"))

    cov_rows <- lapply(seq_along(pair_names), function(i){
      fluidRow(
        column(width = 6,
               tags$span(pair_names[i], class = "text-widget-inner")),
        column(width = 6,
               tags$span(round(ell$cov_matrix[pairs[i, 1], pairs[i, 2]], 4),
                         class = "text-widget-inner text-center"))
      )
    })

    return(
      box(title = tags$span("Covariances", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          p(instructions$build_view_only, class = "text-instruction"),
          cov_rows
      )
    )
  }

  if(isTRUE(covariance_set())){
    return(
      box(title = tags$span("Covariances", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          p(instructions$build_cov_set, class = "text-instruction"),
          fluidRow(
            column(width = 12,
                   actionLink("build_edit_cov_link",
                              label = tagList(icon("pen"), "Edit covariances")))
          )
      )
    )
  }

  box(title = tags$span("Covariances", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$build_covariance, class = "text-instruction"),
      uiOutput("build_covariance_sliders_ui"),
      br(),
      fluidRow(
        column(width = 12,
               div(class = "action-btn-row",
                   actionButton("build_set_cov_btn",
                                "Set Covariances",
                                class = "btn-continue")))
      )
  )
})

output$build_covariance_sliders_ui <- renderUI({

  # Depend on the slot so the sliders redraw when the working ellipsoid
  # changes. Without this the panel renders once and never again.
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  cov_owner(ell$ell_id)

  vars <- ell$var_names
  pairs <- t(combn(vars, 2))
  pair_names <- apply(pairs, 1, function(p) paste(p, collapse = "-"))

  lims <- ell$cov_limits
  rownames(lims) <- pair_names

  sliders <- lapply(seq_along(pair_names), function(i){
    pn <- pair_names[i]
    min_val <- lims[pn, "min"]
    max_val <- lims[pn, "max"]
    step <- round((max_val - min_val) / 100, 4)
    cur_cov <- ell$cov_matrix[pairs[i, 1], pairs[i, 2]]

    fluidRow(
      column(width = 10,
             sliderInput(inputId = paste0("build_cov_", i),
                         label = pn,
                         min = round(min_val, 2),
                         max = round(max_val, 2),
                         value = round(cur_cov, 4),
                         step = step)),
      column(width = 2,
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('build_cov_reset_pair', %d, {priority: 'event'}); return false;", i),
                    tags$span(icon("rotate-left"),
                              title = instructions$build_covariance_reset_tooltip,
                              class = "tooltip-icon")))
    )
  })

  reset_all_btn <- fluidRow(
    column(width = 12,
           tags$a(href = "#",
                  onclick = "Shiny.setInputValue('build_cov_reset_all', Math.random(), {priority: 'event'}); return false;",
                  tagList(icon("rotate-left"), "Reset all to zero")))
  )

  tagList(sliders, reset_all_btn)
})

# Applies slider changes to the working ellipsoid
observeEvent({
  req(session_data$vars)
  n_pairs <- ncol(combn(length(session_data$vars), 2))
  lapply(seq_len(n_pairs), function(i) input[[paste0("build_cov_", i)]])
}, {

  ell <- session_data$current_ellipsoid
  req(ell)
  req(!identical(session_data$ell_mode, "view"))

  # The sliders on screen still belong to a previous ellipsoid. Their
  # values describe that one, not this one, so applying them would
  # overwrite the copy we just made.
  if(!identical(cov_owner(), ell$ell_id)){
    return()
  }

  pairs <- t(combn(ell$var_names, 2))
  pair_names <- apply(pairs, 1, function(p) paste(p, collapse = "-"))

  cov_list <- lapply(seq_along(pair_names),
                     function(i) input[[paste0("build_cov_", i)]])

  if(any(vapply(cov_list, is.null, logical(1)))){
    return()
  }

  # Sliders for a newly loaded ellipsoid may not have rendered yet
  if(any(vapply(cov_list, is.null, logical(1)))) return()

  cov_vals <- setNames(unlist(cov_list), pair_names)

  current_cov <- setNames(sapply(seq_len(nrow(pairs)), function(i){
    ell$cov_matrix[pairs[i, 1], pairs[i, 2]]
  }), pair_names)


  if(all(abs(cov_vals - current_cov) < 1e-8)){
    return()
  }

  if(all(abs(cov_vals - current_cov) < 1e-8)) return()

  updated_ell <- tryCatch(
    update_ellipsoid_covariance(object = ell,
                                covariance = cov_vals,
                                verbose = FALSE),
    error = function(e){
      showNotification(paste("Covariance update failed:", e$message),
                       type = "error", duration = 5)
      NULL
    }
  )

  req(updated_ell)

  session_data$current_ellipsoid <- carry_ell_meta(updated_ell, ell)

  # Rotating one pair narrows what the others can take
  remaining <- updated_ell$cov_limits_remaining

  if(!is.null(remaining)){
    if(is.null(rownames(remaining))) rownames(remaining) <- pair_names

    lapply(seq_along(pair_names), function(i){
      pn <- pair_names[i]
      if(!pn %in% rownames(remaining)) return(NULL)

      cur_val <- input[[paste0("build_cov_", i)]]
      if(is.null(cur_val)) return(NULL)

      new_min <- remaining[pn, "min"]
      new_max <- remaining[pn, "max"]
      clamped <- max(new_min, min(new_max, cur_val))

      updateSliderInput(session,
                        inputId = paste0("build_cov_", i),
                        min = round(new_min, 2),
                        max = round(new_max, 2),
                        value = round(clamped, 2))
    })
  }
}, ignoreNULL = TRUE, ignoreInit = TRUE)

# Resets one covariance pair to zero, index arrives as the input value
observeEvent(input$build_cov_reset_pair, {

  ell <- session_data$current_ellipsoid
  req(ell)

  i <- as.integer(input$build_cov_reset_pair)
  pair_names <- apply(t(combn(ell$var_names, 2)), 1,
                      function(p) paste(p, collapse = "-"))
  req(i >= 1, i <= length(pair_names))

  lims <- ell$cov_limits
  rownames(lims) <- pair_names

  # Restore the full limits as well as the value, since earlier rotations
  # narrowed this slider and the panel is isolated
  updateSliderInput(session,
                    inputId = paste0("build_cov_", i),
                    min = round(lims[pair_names[i], "min"], 2),
                    max = round(lims[pair_names[i], "max"], 2),
                    value = 0)
})

# Resets every covariance pair to zero
observeEvent(input$build_cov_reset_all, {

  ell <- session_data$current_ellipsoid
  req(ell)

  pair_names <- apply(t(combn(ell$var_names, 2)), 1,
                      function(p) paste(p, collapse = "-"))

  lims <- ell$cov_limits
  rownames(lims) <- pair_names

  lapply(seq_along(pair_names), function(i){
    updateSliderInput(session,
                      inputId = paste0("build_cov_", i),
                      min = round(lims[pair_names[i], "min"], 2),
                      max = round(lims[pair_names[i], "max"], 2),
                      value = 0)
  })
})

observeEvent(input$build_set_cov_btn, {
  covariance_set(TRUE)
})

observeEvent(input$build_edit_cov_link, {
  covariance_set(FALSE)
})

# CENTROID ----------------------------------------------------------------

output$build_centroid_mover_ui <- renderUI({

  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  # View mode, show the centroid without controls
  if(identical(ell_mode(), "view")){

    centroid_rows <- lapply(ell$var_names, function(v){
      fluidRow(
        column(width = 6, tags$span(v, class = "text-widget-inner")),
        column(width = 6, tags$span(round(ell$centroid[v], 3),
                                    class = "text-widget-inner text-center"))
      )
    })

    return(
      box(title = tags$span("Centroid Mover", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          p(instructions$build_view_only, class = "text-instruction"),
          centroid_rows
      )
    )
  }

  if(isTRUE(centroid_set())){
    return(
      box(title = tags$span("Centroid Mover", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          p(instructions$build_centroid_set, class = "text-instruction"),
          fluidRow(
            column(width = 12,
                   actionLink("build_edit_centroid_link",
                              label = tagList(icon("pen"), "Edit centroid")))
          )
      )
    )
  }

  box(title = tags$span("Centroid Mover", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = !isTRUE(covariance_set()),
      p(instructions$build_centroid_mover, class = "text-instruction"),
      uiOutput("build_centroid_sliders_ui"),
      fluidRow(
        column(width = 12,
               div(class = "action-btn-row",
                   actionButton("build_set_centroid_btn",
                                "Set Centroid",
                                class = "btn-continue"))
        )
      )
  )
})

output$build_centroid_sliders_ui <- renderUI({

  # Depend on the slot so the sliders redraw when the working ellipsoid
  # changes. Without this the panel renders once and never again.
  ell_slot()

  ell <- isolate(session_data$current_ellipsoid)
  req(ell)

  centroid_owner(ell$ell_id)

  vars <- ell$var_names
  centroid <- ell$centroid

  # The reset target has to sit inside the slider range or it gets clamped
  target <- isolate(ell_reset_target(ell))
  target_centroid <- if(!is.null(target)) target$centroid else centroid

  sliders <- lapply(seq_along(vars), function(j){

    v <- vars[j]

    # Background data gives the widest sensible range. In virtual mode there
    # is none, so fall back to the ellipsoid's own spread along this variable.
    if(!is.null(session_data$bg_df)){
      spread <- 3 * sd(session_data$bg_df[, v], na.rm = TRUE)
      min_val <- round(min(session_data$bg_df[, v], na.rm = TRUE) - spread, 2)
      max_val <- round(max(session_data$bg_df[, v], na.rm = TRUE) + spread, 2)
    } else {
      spread <- 3 * sqrt(ell$cov_matrix[v, v])
      min_val <- round(centroid[v] - spread, 2)
      max_val <- round(centroid[v] + spread, 2)
    }

    min_val <- min(min_val, round(centroid[v], 2), round(target_centroid[v], 2))
    max_val <- max(max_val, round(centroid[v], 2), round(target_centroid[v], 2))

    step <- (max_val - min_val) / 100
    step <- signif(step, 2)

    fluidRow(
      column(width = 12,
             sliderInput(inputId = paste0("build_centroid_", j),
                         label = v,
                         min = min_val,
                         max = max_val,
                         value = round(centroid[v], 2),
                         step = step))
    )
  })

  reset_all_btn <- fluidRow(
    column(width = 12,
           tags$a(href = "#",
                  onclick = "Shiny.setInputValue('build_centroid_reset_all', Math.random(), {priority: 'event'}); return false;",
                  tagList(icon("rotate-left"), ell_reset_label(ell))))
  )

  tagList(sliders, br(), reset_all_btn, br())
})

# Applies slider changes to the working ellipsoid
observeEvent({
  req(session_data$vars)
  lapply(seq_along(session_data$vars),
         function(j) input[[paste0("build_centroid_", j)]])
}, {

  ell <- session_data$current_ellipsoid
  req(ell)
  req(!identical(ell_mode(), "view"))

  # The sliders on screen still belong to a previous ellipsoid, so their
  # values describe that one and must not be applied to this one
  if(!identical(centroid_owner(), ell$ell_id)){
    return()
  }

  vars <- ell$var_names

  centroid_list <- lapply(seq_along(vars),
                          function(j) input[[paste0("build_centroid_", j)]])

  if(any(vapply(centroid_list, is.null, logical(1)))){
    return()
  }

  centroid_vals <- setNames(unlist(centroid_list), vars)
  current_centroid <- setNames(sapply(vars, function(v) ell$centroid[v]), vars)

  # The sliders are drawn with the centroid rounded to two decimals, so on
  # first render they report values that differ from the stored centroid by up
  # to half a step. Comparing against the rounded value makes the untouched
  # case exact instead of relying on a tolerance.
  changed <- centroid_vals != round(current_centroid, 2)

  if(!any(changed)){
    return()
  }

  # Only the variables the user moved are taken from the sliders. Taking all of
  # them would round the untouched ones to two decimals on every move.
  new_centroid <- current_centroid
  new_centroid[changed] <- centroid_vals[changed]

  updated_ell <- tryCatch(
    update_ellipsoid_centroid(ell,
                              new_centroid = new_centroid,
                              verbose = FALSE),
    error = function(e){
      showNotification(paste("Centroid update failed:", e$message),
                       type = "error", duration = 5)
      NULL
    }
  )

  req(updated_ell)

  # range_inputs is deliberately left alone. It records the ranges this
  # ellipsoid was built from, which is what a reset with no parent
  # rebuilds against.
  session_data$current_ellipsoid <- carry_ell_meta(updated_ell, ell)

}, ignoreNULL = TRUE, ignoreInit = TRUE)

observeEvent(input$build_set_centroid_btn, {
  centroid_set(TRUE)
})

observeEvent(input$build_edit_centroid_link, {
  centroid_set(FALSE)
}, ignoreInit = TRUE)

# Resets the centroid to the parent ellipsoid, or to what this ellipsoid's
# own ranges produce when it has no parent
observeEvent(input$build_centroid_reset_all, {

  ell <- session_data$current_ellipsoid
  req(ell)

  target <- ell_reset_target(ell)

  if(is.null(target)){
    showNotification("No reset target available for this ellipsoid.",
                     type = "warning", duration = 4)
    return()
  }

  lapply(seq_along(ell$var_names), function(j){
    v <- ell$var_names[j]
    updateSliderInput(session,
                      inputId = paste0("build_centroid_", j),
                      value = round(target$centroid[[v]], 2))
  })
})
