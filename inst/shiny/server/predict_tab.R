# Title: Predict Tab server

# Description: Server for the predict tab. Projects one or more saved
# ellipsoids onto the background data, and lists the saved ellipsoids in a
# view-only library.

# Date Last Updated: 10/06/2026


# PREDICTION LAYERS -------------------------------------------------------

# Default layer checkboxes, the same as the old static UI
PREDICT_LAYER_DEFAULTS <- c(predict_suitability = TRUE,
                            predict_suitability_trunc = FALSE,
                            predict_mahalanobis = TRUE,
                            predict_mahalanobis_trunc = FALSE)

# Values the checkboxes are drawn with. View sets this, and the stamp makes
# the panel redraw even when the values match the last ones set, since the
# user may have changed the boxes by hand in between.
predict_layer_state <- reactiveVal(list(values = PREDICT_LAYER_DEFAULTS,
                                        stamp = 0L))

output$predict_layers_ui <- renderUI({

  vals <- predict_layer_state()$values

  layer_box <- function(id, label){
    checkboxInput(id,
                  label = tags$span(label, class = "text-widget-inner"),
                  value = isTRUE(vals[[id]]))
  }

  tagList(
    fluidRow(
      column(width = 8,
             tags$span("Prediction layers to include",
                       class = "text-widget-title"))
    ),
    fluidRow(
      column(width = 6, layer_box("predict_suitability", "Suitability")),
      column(width = 6, layer_box("predict_suitability_trunc", "Suitability (truncated)"))
    ),
    fluidRow(
      column(width = 6, layer_box("predict_mahalanobis", "Mahalanobis")),
      column(width = 6, layer_box("predict_mahalanobis_trunc", "Mahalanobis (truncated)"))
    )
  )
})


# ADVANCED SETTINGS -------------------------------------------------------

# Truncation level for the ellipsoid picked in the selector, NULL for All
# versions. If its current prediction used an adjustment, that level is shown,
# so the box always matches the stored layers. Otherwise it is the
# ellipsoid's own cl. With several ellipsoids there is no single default.
predict_selected_cl <- reactive({
  sel <- input$predict_ellipsoid_selected
  if(is.null(sel) || identical(sel, "all")) return(NULL)
  ell <- session_data$ellipsoid_list[[sel]]
  if(is.null(ell)) return(NULL)
  used <- session_data$prediction_settings[[sel]]$adjust_truncation_level
  if(!is.null(used)) used else ell$cl
})

# Drawn once. Only the value follows the selection, through the observer
# below, so the box does not snap shut every time the selection changes.
output$predict_advanced_settings_ui <- renderUI({

  req(!identical(session_data$input_mode, "virtual"))

  cl <- isolate(predict_selected_cl())

  box(title = tags$span("Advanced prediction settings",
                        class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = TRUE,
      p(instructions$predict_trunc_default, class = "text-muted-small"),
      fluidRow(
        column(width = 8,
               tagList(tags$span("Truncation level adjustment",
                                 class = "text-widget-title"),
                       tags$span(icon("circle-info"),
                                 title = instructions$predict_adjust_trunc_tooltip,
                                 class = "tooltip-icon"))),
        column(width = 4,
               numericInput(inputId = "predict_adjust_trunc",
                            label = NULL,
                            value = if(is.null(cl)) NA else cl,
                            min = 0.0001,
                            max = 0.99999,
                            step = 0.05))
      )
  )
})

# Picking another ellipsoid, or changing its cl on Build, resets the box to
# that ellipsoid's own cl. All versions empties it.
observeEvent(predict_selected_cl(), {

  cl <- predict_selected_cl()

  if(is.null(cl)){
    shinyjs::runjs("$('#predict_adjust_trunc').val('').trigger('change');")
  } else {
    updateNumericInput(session, "predict_adjust_trunc", value = cl)
  }
}, ignoreNULL = FALSE, ignoreInit = TRUE)


output$predict_next_step_ui <- renderUI({

  req(length(session_data$ellipsoid_prediction_list) > 0)

  div(class = "action-btn-row",
      actionButton(inputId = "predict_next_step_btn",
                   label = tagList("Continue",
                                   icon("arrow-right")),
                   class = "btn-save")
  )

})

observeEvent(input$predict_next_step_btn, {
  updateTabItems(session, "sidebar_menu", selected = "bias_tab")
})

output$predict_ellipsoid_selector_ui <- renderUI({

  req(length(session_data$ellipsoid_list) > 0)

  if(identical(session_data$input_mode, "virtual")){
    return(p(instructions$predict_virtual_unavailable, class = "text-instruction"))
  }

  versions <- session_data$ellipsoid_list

  ell_choices <- c("All versions" = "all",
                   setNames(names(versions),
                            vapply(versions,
                                   function(ell) ell$ell_name,
                                   character(1))))

  # Keep the current choice across re-renders, otherwise saving an
  # ellipsoid on the Build tab resets this back to All versions
  keep <- if(!is.null(input$predict_ellipsoid_selected) &&
             input$predict_ellipsoid_selected %in% ell_choices){
    input$predict_ellipsoid_selected
  } else {
    "all"
  }

  selectInput(inputId = "predict_ellipsoid_selected",
              label = tagList(
                tags$span("Ellipsoid Version", class = "text-widget-title"),
                tags$span(icon("circle-info"),
                          title = instructions$predict_ellipsoid_select_tooltip,
                          class = "tooltip-icon")
              ),
              choices  = ell_choices,
              selected = keep)
})

observeEvent(input$predict_run_btn, {

  req(length(session_data$ellipsoid_list) > 0)
  req(input$predict_ellipsoid_selected)

  versions <- session_data$ellipsoid_list

  has_raster <- !is.null(session_data$bg_raster)

  if(!has_raster && is.null(session_data$bg_df)){
    showNotification("No background data available for prediction.",
                     type = "error", duration = 4)
    return()
  }

  layers <- c(isTRUE(input$predict_suitability),
              isTRUE(input$predict_suitability_trunc),
              isTRUE(input$predict_mahalanobis),
              isTRUE(input$predict_mahalanobis_trunc))

  if(!any(layers)){
    showNotification(instructions$predict_no_layers,
                     type = "warning", duration = 5)
    return()
  }

  trunc_val <- input$predict_adjust_trunc

  # Guarded for length because an empty numeric input gives numeric(0), which
  # makes the && chain return NA and if() throw
  has_trunc <- length(trunc_val) == 1 && is.finite(trunc_val)

  # Blank, or equal to the ellipsoid's own cl, means no adjustment. Only a
  # value the user actually changed reaches predict(), so a prediction is
  # never truncated at a level nobody chose.
  trunc_for <- function(ell){
    if(!has_trunc) return(NULL)
    if(!is.null(ell$cl) && abs(trunc_val - ell$cl) < 1e-8) return(NULL)
    trunc_val
  }

  predict_one <- function(ell){

    newdata <- if(has_raster){
      terra::subset(session_data$bg_raster, ell$var_names)
    } else {
      session_data$bg_df[, ell$var_names, drop = FALSE]
    }

    tryCatch(
      predict(ell,
              newdata = newdata,
              adjust_truncation_level = trunc_for(ell),
              include_suitability = isTRUE(input$predict_suitability),
              suitability_truncated = isTRUE(input$predict_suitability_trunc),
              include_mahalanobis = isTRUE(input$predict_mahalanobis),
              mahalanobis_truncated = isTRUE(input$predict_mahalanobis_trunc),
              keep_data = TRUE,
              verbose = FALSE),
      error = function(e){
        showNotification(paste0(ell$ell_name, " prediction failed: ", e$message),
                         type = "error", duration = 5)
        NULL
      }
    )
  }

  ids <- if(input$predict_ellipsoid_selected == "all"){
    names(versions)
  } else {
    input$predict_ellipsoid_selected
  }

  n_ok <- 0L

  # Assigned per id rather than replacing the whole list, so one failed
  # ellipsoid does not discard predictions that already succeeded
  for(id in ids){

    ell <- versions[[id]]
    if(is.null(ell)) next

    pred <- predict_one(ell)
    if(is.null(pred)) next

    # Always a single stacked SpatRaster, which is the shape every
    # downstream tab checks for
    session_data$ellipsoid_prediction_list[[id]] <- pred

    # predict() arguments that cannot be recovered from the returned layers.
    # Two predictions at different truncation levels are identical in name and
    # shape, so the report has no way to tell them apart without this. The
    # layer requests are stored for the same reason: predict() also returns
    # the untruncated layers, and the truncated Mahalanobis comes back under
    # the plain Mahalanobis name, so the layers do not show what was asked.
    session_data$prediction_settings[[id]] <- list(
      adjust_truncation_level = trunc_for(ell),
      keep_data = TRUE,
      layers_requested = c(
        predict_suitability = isTRUE(input$predict_suitability),
        predict_suitability_trunc = isTRUE(input$predict_suitability_trunc),
        predict_mahalanobis = isTRUE(input$predict_mahalanobis),
        predict_mahalanobis_trunc = isTRUE(input$predict_mahalanobis_trunc)))

    # A new prediction invalidates anything derived from the old one
    session_data$ellipsoid_prediction_list_biased[[id]] <- NULL
    session_data$ellipsoid_records_list[[id]] <- NULL

    n_ok <- n_ok + 1L
  }

  if(n_ok == 0L){
    showNotification("No predictions were completed.",
                     type = "error", duration = 4)
    return()
  }

  msg <- if(length(ids) > 1){
    paste0("Prediction completed for ", n_ok, " of ", length(ids), " ellipsoids.")
  } else {
    paste0(versions[[ids]]$ell_name, ": prediction completed.")
  }

  showNotification(msg, type = "message", duration = 4)
})


# ELLIPSOID LIBRARY -------------------------------------------------------

output$predict_ellipsoid_library_ui <- renderUI({

  cur_ell <- session_data$current_ellipsoid
  versions <- session_data$ellipsoid_list
  ids <- names(versions)

  req(!is.null(cur_ell) || length(ids) > 0)

  predicted <- names(session_data$ellipsoid_prediction_list)

  # Working slot, the ellipsoid the plots on this tab use. Read-only here,
  # editing happens on the Build tab.
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
             tags$span("Settings loaded", style = "font-size: 11px; color: #aaa;"))
    )
  }

  rows <- lapply(ids, function(id){

    ell <- versions[[id]]
    has_pred <- id %in% predicted

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
             tags$span(if(has_pred) "predicted" else "not predicted",
                       style = paste0("font-size: 10px; color: ",
                                      if(has_pred) "#097a21;" else "#bbb;"))),
      column(width = 3,
             class = "ell-actions",
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('predict_ell_view', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("View ", ell$ell_name, " (read-only)"),
                    icon("eye")),
             # Only a predicted ellipsoid has anything to download
             if(has_pred){
               tags$a(href = "#",
                      onclick = sprintf("Shiny.setInputValue('predict_ell_download', '%s', {priority: 'event'}); return false;", id),
                      title = paste0("Download prediction for ", ell$ell_name),
                      icon("download"))
             } else {
               tags$span(icon("download"),
                         title = "Predict this ellipsoid first",
                         style = "color: #ddd; cursor: not-allowed;")
             },
             tags$a(href = "#",
                    class = "ell-action-danger",
                    onclick = sprintf("Shiny.setInputValue('predict_ell_delete', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("Delete ", ell$ell_name),
                    icon("trash-can")))
    )
  })

  box(title = tags$span("Ellipsoid library", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$predict_library, class = "text-instruction"),

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
        p(instructions$predict_library_empty, class = "text-muted-small")
      }
  )
})

# Guess at which layer checkboxes a stored prediction came from, used only for
# predictions made before layers_requested was stored. It cannot tell a
# truncated Mahalanobis apart, since it comes back as plain Mahalanobis.
pred_layer_flags <- function(pred){
  nms <- tolower(names(pred))
  c(predict_suitability = "suitability" %in% nms,
    predict_suitability_trunc = any(grepl("^suitability_?trunc", nms)),
    predict_mahalanobis = "mahalanobis" %in% nms,
    predict_mahalanobis_trunc = any(grepl("^mahalanobis_?trunc", nms)))
}

# View, read-only. Loading the ellipsoid bumps ell_slot, and the slot
# observer in predict_tab_plot.R moves the selector and sets the layer
# checkboxes from what was requested. The truncation box follows through
# predict_selected_cl(). Syncing here as well made two observers fight.
observeEvent(input$predict_ell_view, {

  id <- input$predict_ell_view
  ell <- session_data$ellipsoid_list[[id]]
  req(ell)

  set_working_ellipsoid(ell, mode = "view")

  pred <- session_data$ellipsoid_prediction_list[[id]]

  showNotification(paste0("Viewing ", ell$ell_name,
                          if(is.null(pred)) ", not predicted yet." else "."),
                   type = "message", duration = 3)
})

# Delete, asks first
observeEvent(input$predict_ell_delete, {

  ell <- session_data$ellipsoid_list[[input$predict_ell_delete]]
  req(ell)

  session_data$pending_ell_delete <- input$predict_ell_delete

  n_children <- sum(vapply(session_data$ellipsoid_list, function(e){
    identical(e$parent_id, ell$ell_id)
  }, logical(1)))

  showModal(modalDialog(
    title = paste0("Delete ", ell$ell_name, "?"),
    p(instructions$predict_delete_ell, class = "text-instruction"),
    if(n_children > 0){
      p(paste0(n_children, " ellipsoid(s) were copied from this one. ",
               "They will be kept, but will no longer have a parent."),
        class = "text-muted-small")
    },
    footer = tagList(
      modalButton("Cancel"),
      actionButton("predict_confirm_ell_delete_btn",
                   "Yes, delete",
                   class = "btn-cancel")
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$predict_confirm_ell_delete_btn, {

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


# PREDICTION DOWNLOAD -----------------------------------------------------

# Coordinate columns, matched with the same patterns the data tab uses
pred_xy_cols <- function(nms){
  nms[grepl(X_COL_PATTERN, nms, ignore.case = TRUE) |
        grepl(Y_COL_PATTERN, nms, ignore.case = TRUE)]
}

# Prediction layers, meaning everything predict() added. The environmental
# variables and coordinates that keep_data carries along are left out.
pred_layer_names <- function(pred, ell){
  nms <- names(pred)
  setdiff(nms, c(ell$var_names, pred_xy_cols(nms)))
}

# The prediction as a data frame with x and y always first. A raster
# prediction gets its cell coordinates. A data frame prediction came from
# bg_df without coordinates, so they are taken back from bg_df by row.
pred_as_df <- function(pred){

  if(inherits(pred, "SpatRaster")){
    # na.rm = NA drops only cells that are NA in every layer, so the
    # truncated layers keep their NAs inside the background
    return(terra::as.data.frame(pred, xy = TRUE, na.rm = NA))
  }

  df <- as.data.frame(pred)
  if(length(pred_xy_cols(names(df))) >= 2) return(df)

  bg <- session_data$bg_df
  x_col <- names(bg)[grepl(X_COL_PATTERN, names(bg), ignore.case = TRUE)][1]
  y_col <- names(bg)[grepl(Y_COL_PATTERN, names(bg), ignore.case = TRUE)][1]

  if(!is.na(x_col) && !is.na(y_col) && nrow(bg) == nrow(df)){
    df <- cbind(bg[, c(x_col, y_col), drop = FALSE], df)
  }

  df
}

# Columns or layers to write, in the prediction's own order
pred_keep_names <- function(pred, ell, layers, include_env){
  want <- c(if(isTRUE(include_env)) ell$var_names, layers)
  names(pred)[names(pred) %in% want]
}

pred_export_name <- function(ell){
  paste0(ell$ell_id, "_prediction_", format(Sys.Date(), "%Y%m%d"))
}

observeEvent(input$predict_ell_download, {

  id <- input$predict_ell_download
  ell <- session_data$ellipsoid_list[[id]]
  pred <- session_data$ellipsoid_prediction_list[[id]]
  req(ell, pred)

  session_data$pending_pred_download <- id

  layers <- pred_layer_names(pred, ell)
  is_rast <- inherits(pred, "SpatRaster")
  has_env <- any(ell$var_names %in% names(pred))

  # Only a raster prediction can be written back as a raster
  format_choices <- if(is_rast){
    c("SpatRaster (.tif)" = "tif", "Data frame (.csv)" = "csv")
  } else {
    c("Data frame (.csv)" = "csv")
  }

  xy_note <- if(!is_rast && length(pred_xy_cols(names(pred_as_df(pred)))) < 2){
    p("The background data has no coordinate columns, so the file will not include x and y.",
      class = "text-muted-small")
  }

  showModal(modalDialog(
    title = paste0("Download prediction for ", ell$ell_name),
    p(instructions$predict_download, class = "text-instruction"),

    checkboxGroupInput("predict_dl_layers",
                       label = tags$span("Prediction layers",
                                         class = "text-widget-title"),
                       choices = layers,
                       selected = layers),

    if(has_env){
      checkboxInput("predict_dl_env",
                    label = tags$span(paste0("Include environmental variables (",
                                             paste(ell$var_names, collapse = ", "), ")"),
                                      class = "text-widget-inner"),
                    value = FALSE)
    },

    radioButtons("predict_dl_format",
                 label = tags$span("Format", class = "text-widget-title"),
                 choices = format_choices,
                 selected = format_choices[[1]],
                 inline = TRUE),

    tags$small("Data frames always include the x and y coordinates.",
               class = "text-muted-small"),
    xy_note,
    br(), br(),

    footer = tagList(
      modalButton("Close"),
      downloadButton("predict_dl_btn", "Download", class = "btn-continue")
    ),
    easyClose = FALSE
  ))
})

# Nothing selected means nothing to write, so the button is off until
# something is
observe({
  req(session_data$pending_pred_download)
  shinyjs::toggleState("predict_dl_btn",
                       condition = length(input$predict_dl_layers) > 0 ||
                         isTRUE(input$predict_dl_env))
})

output$predict_dl_btn <- downloadHandler(

  filename = function(){
    ell <- session_data$ellipsoid_list[[session_data$pending_pred_download]]
    nm <- pred_export_name(ell)
    if(is.null(nm) || !nzchar(trimws(nm))) nm <- pred_export_name(ell)
    nm <- gsub("\\s+", "_", trimws(nm))
    nm <- gsub("[^A-Za-z0-9._-]", "", nm)
    nm <- sub("\\.(tif|tiff|csv)$", "", nm, ignore.case = TRUE)
    ext <- if(identical(input$predict_dl_format, "tif")) ".tif" else ".csv"
    paste0(substr(nm, 1, 60), ext)
  },

  content = function(file){

    id <- session_data$pending_pred_download
    ell <- session_data$ellipsoid_list[[id]]
    pred <- session_data$ellipsoid_prediction_list[[id]]

    if(is.null(ell) || is.null(pred)){
      stop("No prediction available to download.")
    }

    keep <- pred_keep_names(pred, ell, input$predict_dl_layers,
                            input$predict_dl_env)

    if(length(keep) == 0){
      stop("No layers selected.")
    }

    # Subsetting only, no new prediction
    if(identical(input$predict_dl_format, "tif") && inherits(pred, "SpatRaster")){
      terra::writeRaster(pred[[keep]], file, filetype = "GTiff", overwrite = TRUE)
    } else {
      df <- pred_as_df(pred)
      xy <- pred_xy_cols(names(df))
      utils::write.csv(df[, c(xy, keep), drop = FALSE], file, row.names = FALSE)
    }
  }
)
