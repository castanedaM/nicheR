# Title: Generate Tab server

# Description: Server for the generate tab. Samples virtual data points
# from a prediction surface, biased or unbiased, for one or more
# saved ellipsoids.

# Date Last Updated: 10/06/2026


# POINT SET VISIBILITY ----------------------------------------------------

# Point sets currently drawn in the plots. Keys are prefixed with the
# ellipsoid id, since two ellipsoids can produce sets with identical
# parameter signatures.
occ_visible <- reactiveVal(character(0))

# Builds a visibility key. Kept as a function so the separator lives in
# one place.
occ_vis_key <- function(ell_id, set) paste0(ell_id, "::", set)

# Panels stay readable up to a 2x2 grid
OCC_MAX_VISIBLE <- 4L

# Sets the user has toggled visible for the current ellipsoid, in the order
# they were toggled. Falls back to the first set so the plots are never
# blank when data points exist.
generate_visible_sets <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(character(0))

  occ <- session_data$ellipsoid_point_sets_list[[ell$ell_id]]
  if(is.null(occ) || length(occ) == 0) return(character(0))

  prefix <- paste0(ell$ell_id, "::")

  keys <- grep(prefix, occ_visible(), value = TRUE, fixed = TRUE)
  sets <- sub(prefix, "", keys, fixed = TRUE)

  sets <- sets[sets %in% names(occ)]

  if(length(sets) == 0) return(names(occ)[1])

  sets
})

observeEvent(input$generate_toggle_set, {

  idx <- generate_occ_index()
  req(idx)

  i <- as.integer(input$generate_toggle_set)
  req(i >= 1, i <= nrow(idx))

  r <- idx[i, ]
  key <- occ_vis_key(r$ell_id, r$set)

  cur <- occ_visible()

  if(key %in% cur){
    occ_visible(setdiff(cur, key))
    return()
  }

  # Count only the sets visible for this ellipsoid, so other ellipsoids'
  # selections do not eat the budget
  same_ell <- grep(paste0("^", r$ell_id, "::"), cur, value = TRUE)

  if(length(same_ell) >= OCC_MAX_VISIBLE){
    showNotification(instructions$generate_max_visible,
                     type = "warning", duration = 4)
    return()
  }

  occ_visible(c(cur, key))
})


# CONTROLS ----------------------------------------------------------------

# TRUE while the form is open on top of an existing result, so the user can
# regenerate without losing what is already there until they confirm.
generate_show_form <- reactiveVal(FALSE)

# Collapsed state shown once sets exist, in either mode. Factored out so
# the raster and virtual branches cannot drift.
generate_summary_box <- function(){

  occ <- session_data$ellipsoid_point_sets_list

  n_sets <- sum(vapply(occ, length, integer(1)))
  n_pts <- sum(vapply(occ, function(ell_res){
    sum(vapply(ell_res, function(df){
      if(is.null(df)) 0L else nrow(df)
    }, integer(1)))
  }, integer(1)))

  box(title = tags$span("4. Generate data", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = TRUE,
      p(paste0(n_pts, " data point(s) across ", n_sets, " set(s) from ",
               length(occ), " ellipsoid(s)."),
        class = "text-instruction"),
      fluidRow(
        column(width = 12, class = "btn-spaced",
               actionLink("generate_edit_link",
                          label = tagList(icon("pen"), "Generate again")))
      )
  )
}

output$generate_controls_ui <- renderUI({

  is_virtual <- identical(session_data$input_mode, "virtual")

  if(is_virtual){
    return(generate_virtual_controls())
  }

  has_pred <- length(session_data$ellipsoid_prediction_list) > 0

  if(!has_pred){
    return(
      box(title = tags$span("4. Generate data", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          p(instructions$generate_needs_prediction, class = "text-instruction"))
    )
  }

  occ <- session_data$ellipsoid_point_sets_list
  has_occ <- length(occ) > 0

  if(has_occ && !isTRUE(generate_show_form())){
    return(generate_summary_box())
  }

  box(title = tags$span("4. Generate data", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$generate_intro, class = "text-instruction"),

      uiOutput("generate_ellipsoid_selector_ui"),

      fluidRow(
        column(width = 6,
               tags$span("Number of points", class = "text-widget-title"),
               numericInput("generate_n_occ",
                            label = NULL,
                            value = 100,
                            min = 1,
                            max = 1000000,
                            step = 10)),
        column(width = 6,
               tags$div(class = "tooltip-label-row",
                        tags$span("Random seed", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$generate_seed_tooltip,
                                  class = "tooltip-icon")),
               numericInput("generate_seed",
                            label = NULL,
                            value = 123,
                            min = 1,
                            step = 1))
      ),

      # Surfaces come first. The sampling strategy is drawn with them, in
      # generate_surface_ui, so it sits under the unbiased surfaces it
      # applies to.
      fluidRow(
        column(width = 12,
               tags$div(class = "tooltip-label-row",
                        tags$span("Surface", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$generate_surface_tooltip,
                                  class = "tooltip-icon")),
               uiOutput("generate_surface_ui"))
      ),

      # Shared by both kinds of surface: sample_data() and
      # sample_biased_data() both take strict. One setting per run.
      fluidRow(
        column(width = 12,
               tags$div(class = "tooltip-label-row",
                        tags$span("Strict filtering", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$generate_strict_tooltip,
                                  class = "tooltip-icon")),
               radioButtons("generate_strict",
                            label = NULL,
                            choiceNames = list(
                              tags$span("True", class = "text-widget-inner"),
                              tags$span("False", class = "text-widget-inner")),
                            choiceValues = c("TRUE", "FALSE"),
                            selected = "TRUE",
                            inline = TRUE))
      ),

      # The mask is a raster, so the box is left out of a table session
      if(!is.null(session_data$bg_raster)){
        box(title = tagList(
          tags$span("Advanced settings", class = "text-section-header"),
          tags$span(icon("circle-info"),
                    title = instructions$generate_advanced_tooltip,
                    class = "tooltip-icon")),
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,

          fluidRow(
            column(width = 6,
                   fileInput("generate_mask_file",
                             label = tags$span("Sampling mask (optional)",
                                               class = "text-widget-title"),
                             multiple = FALSE,
                             accept = c(".tif", ".tiff", ".rds")))
          ),
          p(instructions$generate_mask, class = "text-instruction")
        )
      },

      fluidRow(
        column(width = 12,
               div(class = "action-btn-row",
                   actionButton("generate_run_btn",
                                label = "Generate",
                                class = "btn-continue"))
        )
      )
  )
})

output$generate_ellipsoid_selector_ui <- renderUI({

  is_virtual <- identical(session_data$input_mode, "virtual")

  # Virtual mode samples straight from the ellipsoid, so every saved
  # version is available. Raster mode can only sample what was predicted.
  ids <- if(is_virtual){
    names(session_data$ellipsoid_list)
  } else {
    names(session_data$ellipsoid_prediction_list)
  }

  req(length(ids) > 0)

  versions <- session_data$ellipsoid_list

  ell_choices <- c(
    "All versions" = "all",
    setNames(ids, vapply(ids, function(id){
      ell <- versions[[id]]
      if(!is.null(ell) && !is.null(ell$ell_name)) ell$ell_name else id
    }, character(1)))
  )

  keep <- if(!is.null(input$generate_ellipsoid_selected) &&
             input$generate_ellipsoid_selected %in% ell_choices){
    input$generate_ellipsoid_selected
  } else {
    "all"
  }

  selectInput(inputId = "generate_ellipsoid_selected",
              label = tagList(
                tags$span("Ellipsoid version", class = "text-widget-title"),
                tags$span(icon("circle-info"),
                          title = instructions$generate_ellipsoid_select_tooltip,
                          class = "tooltip-icon")
              ),
              choices = ell_choices,
              selected = keep)
})

# Virtual mode samples from the ellipsoid directly, so there is no
# prediction, no layer, no mask, and no sampling strategy.
generate_virtual_controls <- function(){

  ell <- session_data$current_ellipsoid

  if(is.null(ell)){
    return(
      box(title = tags$span("Generate data", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          p(instructions$generate_needs_ellipsoid, class = "text-instruction"))
    )
  }

  occ <- session_data$ellipsoid_point_sets_list
  has_occ <- length(occ) > 0

  if(has_occ && !isTRUE(generate_show_form())){
    return(generate_summary_box())
  }

  box(title = tags$span("4. Generate data", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$generate_virtual_intro, class = "text-instruction"),

      uiOutput("generate_ellipsoid_selector_ui"),

      fluidRow(
        column(width = 6,
               tags$span("Number of points", class = "text-widget-title"),
               numericInput("generate_n_occ", label = NULL,
                            value = 100, min = 1, max = 1000000, step = 10)),
        column(width = 6,
               tags$div(class = "tooltip-label-row",
                        tags$span("Random seed", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$generate_seed_tooltip,
                                  class = "tooltip-icon")),
               numericInput("generate_seed", label = NULL,
                            value = 123, min = 1, step = 1))
      ),

      fluidRow(
        column(width = 12,
               tags$div(class = "tooltip-label-row",
                        tags$span("Truncate to the ellipsoid",
                                  class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$generate_truncate_tooltip,
                                  class = "tooltip-icon")),
               radioButtons("generate_truncate", label = NULL,
                            choiceNames = list(
                              tags$span("True", class = "text-widget-inner"),
                              tags$span("False", class = "text-widget-inner")),
                            choiceValues = c("TRUE", "FALSE"),
                            selected = "TRUE",
                            inline = TRUE))
      ),

      fluidRow(
        column(width = 12,
               tags$div(class = "tooltip-label-row",
                        tags$span("Point distribution", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$generate_virtual_sampling_tooltip,
                                  class = "tooltip-icon")),
               uiOutput("generate_virtual_sampling_ui"))
      ),

      fluidRow(
        column(width = 12,
               div(class = "action-btn-row",
                   actionButton("generate_run_btn",
                                label = "Generate",
                                class = "btn-continue")))
      )
  )
}

# Edge and uniform need truncation, so they are removed rather than
# offered and then rejected
output$generate_virtual_sampling_ui <- renderUI({
  truncate <- !identical(input$generate_truncate, "FALSE")

  sampling_labels <- if(truncate) c("Centroid", "Edge", "Uniform") else "Centroid"
  sampling_values <- if(truncate) c("centroid", "edge", "uniform") else "centroid"

  keep <- if(!is.null(input$generate_virtual_sampling) &&
             input$generate_virtual_sampling %in% sampling_values){
    input$generate_virtual_sampling
  } else {
    "centroid"
  }

  tagList(
    radioButtons("generate_virtual_sampling", label = NULL,
                 choiceNames = lapply(sampling_labels,
                                      function(x) tags$span(x, class = "text-widget-inner")),
                 choiceValues = sampling_values,
                 selected = keep, inline = TRUE),
    if(!truncate){
      p(instructions$generate_virtual_sampling_needs_truncate,
        class = "text-muted-small")
    }
  )
})

# Picking a specific version loads it into the working slot so the plots
# follow. "All versions" leaves the slot alone, so the library keeps control.
observeEvent(input$generate_ellipsoid_selected, {

  sel <- input$generate_ellipsoid_selected
  req(sel)

  if(identical(sel, "all")) return()

  ell <- session_data$ellipsoid_list[[sel]]
  req(ell)

  if(identical(session_data$current_ellipsoid$ell_id, sel)) return()

  set_working_ellipsoid(ell, mode = "view")
})

output$generate_surface_ui <- renderUI({

  req(length(session_data$ellipsoid_prediction_list) > 0)

  pred_list <- session_data$ellipsoid_prediction_list
  bias_list <- session_data$ellipsoid_prediction_list_biased

  sel <- input$generate_ellipsoid_selected

  ids <- if(!is.null(sel) && !identical(sel, "all")) sel else names(pred_list)

  versions <- session_data$ellipsoid_list

  # Only prediction layers are surfaces. predict() runs with keep_data = TRUE,
  # so the stored raster also carries the environmental variables, and those
  # are not something to sample from. pred_layer_names() is shared with the
  # Predict tab's download, so both tabs agree on what a layer is.
  # A table session stores each prediction as a data frame, and its columns
  # are filtered the same way as the layers of a raster.
  layer_names <- function(lst){
    unique(unlist(lapply(ids, function(id){
      r <- lst[[id]]
      if(inherits(r, "SpatRaster") || is.data.frame(r)){
        pred_layer_names(r, versions[[id]])
      } else {
        character(0)
      }
    })))
  }

  unbiased_lyrs <- layer_names(pred_list)

  # Biased layers made from an environmental variable are left out. "All
  # prediction layers" on the Bias tab used to bias the variables kept in the
  # prediction too, and a session saved before that was fixed still has them.
  bias_lyrs <- unique(unlist(lapply(ids, function(id){
    r <- bias_list[[id]]
    if(!inherits(r, "SpatRaster")) return(character(0))
    nms <- pred_layer_names(r, versions[[id]])
    src <- sub("_(centroid|edge|uniform)_biased$", "", nms)
    nms[!src %in% c("x", "y", versions[[id]]$var_names)]
  })))

  bias_lyrs <- setdiff(bias_lyrs, unbiased_lyrs)

  req(length(c(unbiased_lyrs, bias_lyrs)) > 0)

  has_bias <- length(bias_lyrs) > 0

  # Current ticks are read with isolate(), so ticking a box does not redraw
  # the whole block. It is redrawn only when the ellipsoid or the available
  # layers change, and the ticks that still apply are kept.
  cur_unbiased <- isolate(input$generate_surface)
  cur_biased <- isolate(input$generate_surface_biased)
  cur_sampling <- isolate(input$generate_sampling)

  keep_biased <- intersect(cur_biased, bias_lyrs)

  # Falls back to a default only when nothing is ticked on either side, so
  # a biased-only selection is not given an unbiased surface on redraw
  keep_unbiased <- if(any(cur_unbiased %in% unbiased_lyrs)){
    intersect(cur_unbiased, unbiased_lyrs)
  } else if(length(keep_biased) > 0){
    character(0)
  } else if("suitability_trunc" %in% unbiased_lyrs){
    "suitability_trunc"
  } else {
    unbiased_lyrs[1]
  }

  keep_sampling <- if(is.null(cur_sampling)) "centroid" else cur_sampling

  note_style <- "font-size: 10px; color: #aaa; margin: 4px 0 8px;"

  surface_boxes <- function(id, layers, selected, inline, col = ""){

    if(length(layers) == 0){
      return(tags$p("None available.", style = note_style))
    }

    checkboxGroupInput(id,
                       label = NULL,
                       choiceNames = lapply(layers, function(nm){
                         tags$span(nm, style = col, class = "text-widget-inner")
                       }),
                       choiceValues = layers,
                       selected = selected,
                       inline = inline)
  }

  # Applies to the unbiased surfaces only. Checkboxes, so several strategies
  # can be run with the same number of points and seed.
  sampling_block <- tagList(
    tags$div(class = "tooltip-label-row",
             tags$span("Sampling strategy", class = "text-widget-title"),
             tags$span(icon("circle-info"),
                       title = instructions$generate_sampling_tooltip,
                       class = "tooltip-icon")),
    checkboxGroupInput("generate_sampling",
                       label = NULL,
                       choiceNames = list(
                         tags$span("Centroid", class = "text-widget-inner"),
                         tags$span("Edge", class = "text-widget-inner"),
                         tags$span("Uniform", class = "text-widget-inner")
                       ),
                       choiceValues = c("centroid", "edge", "uniform"),
                       selected = keep_sampling,
                       inline = TRUE),
    uiOutput("generate_method_msg_ui")
  )

  # No biased layer: one column, as before
  if(!has_bias){
    return(tagList(
      surface_boxes("generate_surface", unbiased_lyrs, keep_unbiased, inline = TRUE),
      sampling_block
    ))
  }

  # Biased layers get their own column, split off by a vertical line. The
  # sampling strategy stays on the unbiased side, so it does not read as if
  # a biased surface were sampled toward the centroid or the edge.
  fluidRow(
    column(width = 6,
           tags$span("Unbiased", class = "text-widget-title"),
           surface_boxes("generate_surface", unbiased_lyrs, keep_unbiased,
                         inline = FALSE),
           sampling_block),
    column(width = 6,
           style = "border-left: 1px solid #ddd;",
           tags$span("Biased", class = "text-widget-title",
                     style = "color: #c47c16;"),
           surface_boxes("generate_surface_biased", bias_lyrs, keep_biased,
                         inline = FALSE, col = "color: #c47c16;"),
           uiOutput("generate_method_biased_msg_ui"))
  )
})

# Shown only while a biased surface is ticked, the same way the unbiased
# message comes and goes
output$generate_method_biased_msg_ui <- renderUI({

  req(length(input$generate_surface_biased) > 0)

  tags$p(icon("circle-info"), " ",
         paste0("Method detected: ", occ_sampling_biased, "."),
         style = "font-size: 10px; color: #aaa; margin: 4px 0 8px;")
})

# Mirrors how generate_occ_for_ell picks a method for the unbiased surfaces,
# so the user sees what will happen before pressing Generate
output$generate_method_msg_ui <- renderUI({

  req(input$generate_surface)
  req(length(input$generate_surface) > 0)

  methods <- unique(vapply(input$generate_surface, function(layer){
    if(grepl("mahalanobis", layer, ignore.case = TRUE)){
      "mahalanobis"
    } else {
      "suitability"
    }
  }, character(1)))

  fluidRow(
    column(width = 12,
           tags$p(icon("circle-info"), " ",
                  paste0("Method(s) detected: ", paste(methods, collapse = ", "), "."),
                  style = "font-size: 10px; color: #aaa; margin: 4px 0 8px;"))
  )
})


# GENERATE ----------------------------------------------------------------

observeEvent(input$generate_run_btn, {

  # VIRTUAL MODE ----------------------------------------------------------
  # Samples straight from the ellipsoid, so there is no prediction, no
  # layer, no mask, and no sampling strategy. Handled first because the
  # req() calls below reference inputs that do not exist here.
  if(identical(session_data$input_mode, "virtual")){

    req(input$generate_n_occ)
    req(input$generate_ellipsoid_selected)

    versions <- session_data$ellipsoid_list
    req(length(versions) > 0)

    selected_ids <- if(identical(input$generate_ellipsoid_selected, "all")){
      names(versions)
    } else {
      input$generate_ellipsoid_selected
    }

    n <- as.integer(input$generate_n_occ)
    truncate <- !identical(input$generate_truncate, "FALSE")
    sampling <- if(!is.null(input$generate_virtual_sampling)){
      input$generate_virtual_sampling
    } else {
      "centroid"
    }

    seed <- if(!is.null(input$generate_seed) && is.finite(input$generate_seed)){
      as.integer(input$generate_seed)
    } else {
      123L
    }

    # Same signature for every ellipsoid, which is correct: sets are keyed
    # within an ellipsoid, so two ellipsoids sharing parameters share a key.
    # Sampling is the value passed to virtual_data().
    attrs <- list(mode = "virtual", layer = "virtual", n_occ = n,
                  seed = seed,
                  sampling = sampling,
                  strict = NA, mask = "none",
                  truncate = truncate)

    key <- occ_set_signature(attrs)

    n_success <- 0L

    for(id in selected_ids){

      ell <- versions[[id]]
      if(is.null(ell)) next

      df <- generate_virtual_for_ell(ell, n, truncate, sampling, seed)
      if(is.null(df) || nrow(df) == 0) next

      for(f in occ_set_fields) attr(df, f) <- attrs[[f]]
      attr(df, "created") <- format(Sys.time(), "%Y-%m-%d %H:%M")

      if(is.null(session_data$ellipsoid_point_sets_list[[id]])){
        session_data$ellipsoid_point_sets_list[[id]] <- list()
      }

      session_data$ellipsoid_point_sets_list[[id]][[key]] <- df

      # New sets are shown by default, up to the panel limit
      vis <- occ_visible()
      vis_key <- occ_vis_key(id, key)
      same_ell <- grep(paste0("^", id, "::"), vis, value = TRUE)

      if(!vis_key %in% vis && length(same_ell) < OCC_MAX_VISIBLE){
        occ_visible(c(vis, vis_key))
      }

      n_success <- n_success + 1L
    }

    if(n_success == 0L){
      showNotification("Virtual generation failed for all selected ellipsoids.",
                       type = "error", duration = 5)
      return()
    }

    generate_show_form(FALSE)

    msg <- if(length(selected_ids) > 1){
      paste0("Virtual points generated for ", n_success, " of ",
             length(selected_ids), " ellipsoids.")
    } else {
      paste0(n, " virtual points generated.")
    }

    showNotification(msg, type = "message", duration = 4)
    return()
  }


  # RASTER MODE -----------------------------------------------------------

  req(input$generate_ellipsoid_selected)
  req(input$generate_n_occ)

  pred_list <- session_data$ellipsoid_prediction_list
  bias_list <- session_data$ellipsoid_prediction_list_biased

  selected_ids <- if(identical(input$generate_ellipsoid_selected, "all")){
    names(pred_list)
  } else {
    input$generate_ellipsoid_selected
  }

  # Biased layers that exist for the selected ellipsoids. The biased
  # checkboxes are only drawn when there is a biased layer, and Shiny keeps
  # the last value of an input that is no longer drawn, so the ticks are
  # checked against what is there now.
  bias_available <- unique(unlist(lapply(selected_ids, function(id){
    r <- bias_list[[id]]
    if(inherits(r, "SpatRaster")) names(r) else character(0)
  })))

  unbiased_layers <- input$generate_surface
  biased_layers <- intersect(input$generate_surface_biased, bias_available)

  if(length(unbiased_layers) == 0 && length(biased_layers) == 0){
    showNotification(instructions$generate_no_surface,
                     type = "warning", duration = 5)
    return()
  }

  # The strategy only matters when an unbiased surface is ticked
  sampling_set <- input$generate_sampling

  if(length(unbiased_layers) > 0 && length(sampling_set) == 0){
    showNotification("Select at least one sampling strategy for the unbiased surfaces.",
                     type = "warning", duration = 5)
    return()
  }

  # Default to TRUE rather than NULL if the radio has not reported yet
  strict <- !identical(input$generate_strict, "FALSE")

  n_occ <- as.integer(input$generate_n_occ)

  # One run per sampling strategy for the unbiased surfaces, plus one for
  # the biased surfaces, which have no strategy since the surface values are
  # the weights. Strict applies to all of them.
  runs <- list()

  if(length(unbiased_layers) > 0){
    for(samp in sampling_set){
      runs[[length(runs) + 1]] <- list(layers = unbiased_layers, sampling = samp)
    }
  }

  if(length(biased_layers) > 0){
    runs[[length(runs) + 1]] <- list(layers = biased_layers, sampling = NA)
  }

  seed <- if(!is.null(input$generate_seed) && is.finite(input$generate_seed)){
    as.integer(input$generate_seed)
  } else {
    123L
  }

  # A mask is a raster, so it only applies when the session has one. In a
  # table session the box is not drawn, but an input keeps its last value.
  use_mask <- !is.null(input$generate_mask_file) &&
    !is.null(session_data$bg_raster)

  sampling_mask <- if(use_mask){
    ext <- tolower(tools::file_ext(input$generate_mask_file$name))
    tryCatch(
      load_raster_file(input$generate_mask_file$datapath, ext),
      error = function(e){
        showNotification(paste("Could not load sampling mask:", e$message),
                         type = "error", duration = 5)
        NULL
      }
    )
  } else {
    NULL
  }

  session_data$sampling_mask <- sampling_mask

  n_success <- 0L
  n_attempted <- 0L

  mask_name <- if(use_mask){
    input$generate_mask_file$name
  } else {
    "none"
  }

  # Assigned per id rather than replacing the whole list, so one failed
  # ellipsoid does not discard data points that already succeeded
  for(id in selected_ids){

    for(run in runs){

      res <- generate_occ_for_ell(ell_id = id,
                                  pred_list = pred_list,
                                  biased_list = bias_list,
                                  layers = run$layers,
                                  n_occ = n_occ,
                                  sampling = run$sampling,
                                  strict = strict,
                                  sampling_mask = sampling_mask,
                                  seed = seed,
                                  var_names = session_data$ellipsoid_list[[id]]$var_names)

      n_attempted <- n_attempted + length(run$layers)

      for(layer in names(res)){

        df <- res[[layer]]
        if(is.null(df) || nrow(df) == 0) next

        # Set by generate_occ_for_ell, which knows which function drew it
        biased <- isTRUE(attr(df, "biased", exact = TRUE))

        # Metadata travels with the set so the summary can report how it was
        # made, and so the source raster can still be found from the layer.
        # Biased layers have no sampling strategy, since the surface values
        # are the weights, so they are labeled as bias weighted instead.
        # mode "raster" means drawn from a prediction surface, whether that
        # surface is a raster or, in a table session, a data frame.
        attrs <- list(layer = layer,
                      n_occ = n_occ,
                      seed = seed,
                      sampling = if(biased) occ_sampling_biased else run$sampling,
                      strict = strict,
                      mask = mask_name,
                      mode = "raster",
                      truncate = NA)

        for(f in occ_set_fields) attr(df, f) <- attrs[[f]]

        attr(df, "created") <- format(Sys.time(), "%Y-%m-%d %H:%M")

        key <- occ_set_signature(attrs)

        if(is.null(session_data$ellipsoid_point_sets_list[[id]])){
          session_data$ellipsoid_point_sets_list[[id]] <- list()
        }

        session_data$ellipsoid_point_sets_list[[id]][[key]] <- df

        # New sets are shown by default, up to the panel limit
        vis <- occ_visible()
        vis_key <- occ_vis_key(id, key)
        same_ell <- grep(paste0("^", id, "::"), vis, value = TRUE)

        if(!vis_key %in% vis && length(same_ell) < OCC_MAX_VISIBLE){
          occ_visible(c(vis, vis_key))
        }

        n_success <- n_success + 1L
      }
    }
  }

  if(n_success == 0L){
    showNotification("Data generation failed for all selected combinations.",
                     type = "error", duration = 5)
    return()
  }

  generate_show_form(FALSE)

  n_skipped <- n_attempted - n_success

  msg <- paste0(n_success, " point set(s) generated.")
  if(n_skipped > 0L) msg <- paste0(msg, " ", n_skipped, " skipped.")

  showNotification(msg, type = "message", duration = 4)
})

observeEvent(input$generate_edit_link, {
  generate_show_form(TRUE)
})

# POINT SETS --------------------------------------------------------------

# Flat index of every point set, so the summary can render rows and
# the download handlers can be created against stable numeric ids.
generate_occ_index <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(NULL)

  occ <- session_data$ellipsoid_point_sets_list[[ell$ell_id]]
  if(is.null(occ) || length(occ) == 0) return(NULL)

  lab_of <- generate_occ_label_map()

  rows <- lapply(names(occ), function(nm){

    df <- occ[[nm]]
    if(is.null(df)) return(NULL)

    data.frame(
      ell_id = ell$ell_id,
      ell_name = ell$ell_name,
      set = nm,
      label = lab_of[[nm]],
      layer = occ_meta(df, "layer", nm),
      n = nrow(df),
      seed = occ_meta(df, "seed"),
      n_occ = occ_meta(df, "n_occ"),
      sampling = occ_sampling(df),
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if(length(rows) == 0) return(NULL)

  do.call(rbind, rows)
})


# One flattening function, used by all three download scopes. The metadata
# is read from the set itself, so the columns cannot differ between scopes.
# idx_rows needs ell_id, ell_name, and set.
generate_occ_table <- function(idx_rows){

  occ <- session_data$ellipsoid_point_sets_list

  parts <- lapply(seq_len(nrow(idx_rows)), function(j){
    r <- idx_rows[j, ]
    df <- occ[[r$ell_id]][[r$set]]
    if(is.null(df) || nrow(df) == 0) return(NULL)

    out <- data.frame(ell_id = r$ell_id,
                      ell_name = r$ell_name,
                      layer = occ_meta(df, "layer", r$set),
                      n = occ_meta(df, "n_occ"),
                      seed = occ_meta(df, "seed"),
                      sampling = occ_sampling(df),
                      strict = occ_meta(df, "strict"),
                      truncate = occ_meta(df, "truncate"),
                      df,
                      stringsAsFactors = FALSE)

    # Each mode keeps the column that applies to it. Raster sets have strict,
    # and carry truncation in the layer name. Virtual sets have no layer name
    # to carry truncation, and no strict.
    if(identical(occ_meta(df, "mode"), "virtual")){
      out$strict <- NULL
    } else {
      out$truncate <- NULL
    }

    out
  })

  do.call(rbind, parts)
}

# Handlers are rebuilt whenever the index changes, so the numeric ids in the
# summary always point at the right set
observe({

  idx <- generate_occ_index()
  req(idx)

  lapply(seq_len(nrow(idx)), function(i){
    local({
      my_i <- i
      output[[paste0("generate_dl_set_", my_i)]] <- downloadHandler(
        filename = function(){
          r <- generate_occ_index()[my_i, ]
          nm <- gsub("[^A-Za-z0-9_-]", "_",
                     paste0(r$ell_name, "_", r$layer, "_n", r$n_occ, "_seed", r$seed))
          paste0(nm, ".csv")
        },
        content = function(file){
          r <- generate_occ_index()[my_i, , drop = FALSE]
          write.csv(generate_occ_table(r), file, row.names = FALSE)
        }
      )
    })
  })
})

output$generate_dl_ell <- downloadHandler(
  filename = function(){
    idx <- generate_occ_index()
    nm <- gsub("[^A-Za-z0-9_-]", "_", idx$ell_name[1])
    paste0(nm, "_data_points.csv")
  },
  content = function(file){
    idx <- generate_occ_index()
    req(idx)
    write.csv(generate_occ_table(idx), file, row.names = FALSE)
  }
)

output$generate_dl_all <- downloadHandler(

  filename = function(){
    paste0("nicheR_data_points_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file){

    occ <- session_data$ellipsoid_point_sets_list
    req(length(occ) > 0)

    versions <- session_data$ellipsoid_list

    # Every set of every ellipsoid, flattened by the same function the other
    # two downloads use
    idx <- do.call(rbind, lapply(names(occ), function(id){
      if(length(occ[[id]]) == 0) return(NULL)
      data.frame(ell_id = id,
                 ell_name = if(!is.null(versions[[id]])) versions[[id]]$ell_name else id,
                 set = names(occ[[id]]),
                 stringsAsFactors = FALSE)
    }))
    req(!is.null(idx))

    out <- generate_occ_table(idx)
    req(!is.null(out))

    write.csv(out, file, row.names = FALSE)
  }
)

observeEvent(input$generate_delete_set, {

  idx <- generate_occ_index()
  req(idx)

  i <- as.integer(input$generate_delete_set)
  req(i >= 1, i <= nrow(idx))

  r <- idx[i, ]

  session_data$ellipsoid_point_sets_list[[r$ell_id]][[r$set]] <- NULL

  # Drop the ellipsoid entry entirely once its last set is gone, so the
  # library status returns to "not generated"
  if(length(session_data$ellipsoid_point_sets_list[[r$ell_id]]) == 0){
    session_data$ellipsoid_point_sets_list[[r$ell_id]] <- NULL
  }

  occ_visible(setdiff(occ_visible(), occ_vis_key(r$ell_id, r$set)))

  showNotification(paste0("Set removed: ", r$layer, " (seed ", r$seed, ")."),
                   type = "message", duration = 3)
})

observeEvent(input$generate_new_set_link, {
  generate_show_form(TRUE)
})

# POINT SET SUMMARY -------------------------------------------------------

output$generate_point_sets_summary_ui <- renderUI({

  ell <- session_data$current_ellipsoid

  idx <- generate_occ_index()

  if(is.null(idx)){
    return(
      box(title = tags$span("Point sets", class = "text-section-header"),
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          p(instructions$generate_summary_empty, class = "text-muted-small"))
    )
  }

  biased <- session_data$ellipsoid_prediction_list_biased[[ell$ell_id]]

  is_biased <- function(layer){
    inherits(biased, "SpatRaster") && layer %in% names(biased)
  }

  rows <- lapply(seq_len(nrow(idx)), function(i){

    r <- idx[i, ]
    shown <- occ_vis_key(r$ell_id, r$set) %in% occ_visible()

    fluidRow(
      class = "ell-row",
      style = "padding: 2px 0;",
      column(width = 5,
             tags$span(r$label, class = "text-widget-inner",
                       style = if(is_biased(r$layer)) "color: #c47c16;" else "")),
      column(width = 4,
             tags$span(paste0(r$n, " points"), class = "text-widget-inner"),
             if(!is.na(r$n_occ) && r$n != r$n_occ){
               tagList(tags$br(),
                       tags$span(paste0("requested ", r$n_occ),
                                 style = "font-size: 10px; color: #bbb;"))
             }),
      column(width = 3,
             class = "ell-actions",
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('generate_toggle_set', %d, {priority: 'event'}); return false;", i),
                    title = if(shown) "Hide from plots" else "Show in plots",
                    style = if(shown) "color: #097a21;" else "color: #ccc;",
                    icon(if(shown) "eye" else "eye-slash")),
             downloadLink(paste0("generate_dl_set_", i),
                          label = icon("download")),
             tags$a(href = "#",
                    class = "ell-action-danger",
                    onclick = sprintf("Shiny.setInputValue('generate_delete_set', %d, {priority: 'event'}); return false;", i),
                    title = "Delete this set",
                    icon("trash-can")))
    )
  })

  box(title = tags$span(paste0("Point sets: ", idx$ell_name[1]),
                        class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$generate_summary, class = "text-instruction"),

      fluidRow(
        column(width = 4, class = "btn-spaced",
               actionLink("generate_new_set_link",
                          label = tagList(icon("circle-plus"), "New set"))),
        column(width = 4, class = "btn-spaced",
               downloadLink("generate_dl_ell",
                            label = tagList(icon("download"), "This ellipsoid"))),
        column(width = 4, class = "btn-spaced",
               downloadLink("generate_dl_all",
                            label = tagList(icon("download"), "All ellipsoids")))
      ),

      br(),

      fluidRow(
        class = "ell-row",
        style = "padding: 4px 0; border-top: 1px solid #eee;",
        column(width = 5, tags$span("Set", class = "text-widget-title")),
        column(width = 4, tags$span("Points", class = "text-widget-title")),
        column(width = 3, tags$span("Actions", class = "text-widget-title"))
      ),

      tagList(rows)
  )
})

# ELLIPSOID LIBRARY -------------------------------------------------------

output$generate_ellipsoid_library_ui <- renderUI({

  cur_ell <- session_data$current_ellipsoid
  versions <- session_data$ellipsoid_list
  ids <- names(versions)

  req(!is.null(cur_ell) || length(ids) > 0)

  predicted <- names(session_data$ellipsoid_prediction_list)
  biased <- names(session_data$ellipsoid_prediction_list_biased)
  generated <- names(session_data$ellipsoid_point_sets_list)

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
             tags$span("View-only", style = "font-size: 11px; color: #aaa;"))
    )
  }

  rows <- lapply(ids, function(id){

    ell <- versions[[id]]

    status <- if(id %in% generated){
      list(txt = "data points generated", col = "#097a21")
    } else if(id %in% biased){
      list(txt = "biased, not generated", col = "#aaa")
    } else if(id %in% predicted){
      list(txt = "predicted, not generated", col = "#aaa")
    } else {
      list(txt = "not predicted", col = "#bbb")
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
             tags$span(status$txt,
                       style = paste0("font-size: 10px; color: ", status$col, ";"))),
      column(width = 3,
             class = "ell-actions",
             tags$a(href = "#",
                    onclick = sprintf("Shiny.setInputValue('generate_ell_view', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("View ", ell$ell_name, " (read-only)"),
                    icon("eye")),
             tags$a(href = "#",
                    class = "ell-action-danger",
                    onclick = sprintf("Shiny.setInputValue('generate_ell_delete', '%s', {priority: 'event'}); return false;", id),
                    title = paste0("Delete ", ell$ell_name),
                    icon("trash-can")))
    )
  })

  box(title = tags$span("Ellipsoid library", class = "text-section-header"),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      p(instructions$generate_library, class = "text-instruction"),

      if(!is.null(cur_ell)){
        tagList(working_row, tags$hr(style = "margin: 8px 0;"))
      },

      br(),

      if(length(ids) > 0){
        tagList(
          fluidRow(
            class = "ell-row",
            style = "padding: 2px 0;",
            column(width = 5, tags$span("Name", class = "text-widget-title")),
            column(width = 4, tags$span("Built from", class = "text-widget-title")),
            column(width = 3, tags$span("Actions", class = "text-widget-title"))
          ),
          tagList(rows)
        )
      } else {
        p(instructions$generate_library_empty, class = "text-muted-small")
      }
  )
})

# View, read-only. Also points the ellipsoid selector at it, as the Predict
# tab does, when it has a prediction to sample from.
observeEvent(input$generate_ell_view, {

  id <- input$generate_ell_view
  ell <- session_data$ellipsoid_list[[id]]
  req(ell)

  set_working_ellipsoid(ell, mode = "view")

  selectable <- if(identical(session_data$input_mode, "virtual")){
    names(session_data$ellipsoid_list)
  } else {
    names(session_data$ellipsoid_prediction_list)
  }

  if(id %in% selectable){
    updateSelectInput(session, "generate_ellipsoid_selected", selected = id)
  }

  showNotification(paste0("Viewing ", ell$ell_name, "."),
                   type = "message", duration = 3)
})

# Delete, asks first
observeEvent(input$generate_ell_delete, {

  ell <- session_data$ellipsoid_list[[input$generate_ell_delete]]
  req(ell)

  session_data$pending_ell_delete <- input$generate_ell_delete

  n_children <- sum(vapply(session_data$ellipsoid_list, function(e){
    identical(e$parent_id, ell$ell_id)
  }, logical(1)))

  showModal(modalDialog(
    title = paste0("Delete ", ell$ell_name, "?"),
    p(instructions$generate_delete_ell, class = "text-instruction"),
    if(n_children > 0){
      p(paste0(n_children, " ellipsoid(s) were copied from this one. ",
               "They will be kept, but will no longer have a parent."),
        class = "text-muted-small")
    },
    footer = tagList(
      modalButton("Cancel"),
      actionButton("generate_confirm_ell_delete_btn",
                   "Yes, delete",
                   class = "btn-cancel")
    ),
    easyClose = FALSE
  ))
})

observeEvent(input$generate_confirm_ell_delete_btn, {

  id <- session_data$pending_ell_delete
  req(id)

  nm <- session_data$ellipsoid_list[[id]]$ell_name

  removeModal()

  session_data$ellipsoid_list[[id]] <- NULL
  session_data$ellipsoid_prediction_list[[id]] <- NULL
  session_data$prediction_settings[[id]] <- NULL
  session_data$ellipsoid_prediction_list_biased[[id]] <- NULL
  session_data$ellipsoid_point_sets_list[[id]] <- NULL
  session_data$pending_ell_delete <- NULL

  # Copies of the deleted ellipsoid are kept and become roots
  session_data$ellipsoid_list <- lapply(session_data$ellipsoid_list, function(e){
    if(identical(e$parent_id, id)) e$parent_id <- NULL
    e
  })

  cur <- session_data$current_ellipsoid

  occ_visible(grep(paste0("^", id, "::"), occ_visible(),
                   value = TRUE, invert = TRUE))

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
