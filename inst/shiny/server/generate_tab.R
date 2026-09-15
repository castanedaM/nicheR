# Title: Generate Tab server

# Description: Server for the generate tab. Samples virtual occurrence
# points from a prediction surface, biased or unbiased, for one or more
# saved ellipsoids.

# Date Last Updated: 08/06/2026


# OCCURRENCE SET VISIBILITY -----------------------------------------------

# Occurrence sets currently drawn in the plots. Keys are prefixed with the
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
# blank when occurrences exist.
generate_visible_sets <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(character(0))

  occ <- session_data$ellipsoid_occurrence_list[[ell$ell_id]]
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

  occ <- session_data$ellipsoid_occurrence_list

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
      p(paste0(n_pts, " occurrence(s) across ", n_sets, " set(s) from ",
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

  occ <- session_data$ellipsoid_occurrence_list
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

      fluidRow(
        column(width = 12,
               tags$div(class = "tooltip-label-row",
                        tags$span("Sampling strategy", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$generate_sampling_tooltip,
                                  class = "tooltip-icon")),
               radioButtons("generate_sampling",
                            label = NULL,
                            choiceNames = list(
                              tags$span("Centroid", class = "text-widget-inner"),
                              tags$span("Edge", class = "text-widget-inner"),
                              tags$span("Random", class = "text-widget-inner")
                            ),
                            choiceValues = c("centroid", "edge", "random"),
                            selected = "centroid",
                            inline = TRUE))
      ),

      fluidRow(
        column(width = 12,
               tags$div(class = "tooltip-label-row",
                        tags$span("Surface", class = "text-widget-title"),
                        tags$span(icon("circle-info"),
                                  title = instructions$generate_surface_tooltip,
                                  class = "tooltip-icon")),
               uiOutput("generate_surface_ui"))
      ),

      uiOutput("generate_method_msg_ui"),

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
      ),

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

  occ <- session_data$ellipsoid_occurrence_list
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
                                  title = instructions$generate_effect_tooltip,
                                  class = "tooltip-icon")),
               uiOutput("generate_effect_ui"))
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

# Inverse and uniform need truncation, so they are removed rather than
# offered and then rejected
output$generate_effect_ui <- renderUI({
  truncate <- !identical(input$generate_truncate, "FALSE")

  effect_labels <- if(truncate) c("Centroid", "Edge", "Random") else "Centroid"
  effect_values <- if(truncate) c("direct", "inverse", "uniform") else "direct"

  keep <- if(!is.null(input$generate_effect) &&
             input$generate_effect %in% effect_values){
    input$generate_effect
  } else {
    "direct"
  }

  tagList(
    radioButtons("generate_effect", label = NULL,
                 choiceNames = lapply(effect_labels,
                                      function(x) tags$span(x, class = "text-widget-inner")),
                 choiceValues = effect_values,
                 selected = keep, inline = TRUE),
    if(!truncate){
      p(instructions$generate_effect_needs_truncate, class = "text-muted-small")
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

  layer_names <- function(lst){
    unique(unlist(lapply(ids, function(id){
      r <- lst[[id]]
      if(inherits(r, "SpatRaster")) names(r) else character(0)
    })))
  }

  unbiased_lyrs <- layer_names(pred_list)
  bias_lyrs <- layer_names(bias_list)

  all_layers <- unique(c(unbiased_lyrs, bias_lyrs))
  req(length(all_layers) > 0)

  keep <- if(!is.null(input$generate_surface) &&
             any(input$generate_surface %in% all_layers)){
    intersect(input$generate_surface, all_layers)
  } else if("suitability_trunc" %in% unbiased_lyrs){
    "suitability_trunc"
  } else if(length(unbiased_lyrs) > 0){
    unbiased_lyrs[1]
  } else {
    all_layers[1]
  }

  checkboxGroupInput("generate_surface",
                     label = NULL,
                     choiceNames = lapply(all_layers, function(nm){
                       is_biased <- nm %in% bias_lyrs && !nm %in% unbiased_lyrs
                       tags$span(nm,
                                 style = if(is_biased) "color: #c47c16;" else "",
                                 class = "text-widget-inner")
                     }),
                     choiceValues = all_layers,
                     selected = keep,
                     inline = TRUE)
})

# Mirrors how generate_occ_for_ell picks a method, so the user sees what
# will happen before pressing Generate
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
    effect <- if(!is.null(input$generate_effect)) input$generate_effect else "direct"

    seed <- if(!is.null(input$generate_seed) && is.finite(input$generate_seed)){
      as.integer(input$generate_seed)
    } else {
      123L
    }

    # Same signature for every ellipsoid, which is correct: sets are keyed
    # within an ellipsoid, so two ellipsoids sharing parameters share a key
    attrs <- list(mode = "virtual", layer = "virtual", n_occ = n,
                  seed = seed, sampling = NA, strict = NA, mask = "none",
                  truncate = truncate, effect = effect)

    key <- occ_set_signature(attrs)

    n_success <- 0L

    for(id in selected_ids){

      ell <- versions[[id]]
      if(is.null(ell)) next

      df <- generate_virtual_for_ell(ell, n, truncate, effect, seed)
      if(is.null(df) || nrow(df) == 0) next

      for(f in occ_set_fields) attr(df, f) <- attrs[[f]]
      attr(df, "created") <- format(Sys.time(), "%Y-%m-%d %H:%M")

      if(is.null(session_data$ellipsoid_occurrence_list[[id]])){
        session_data$ellipsoid_occurrence_list[[id]] <- list()
      }

      session_data$ellipsoid_occurrence_list[[id]][[key]] <- df

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
  req(input$generate_sampling)

  if(is.null(input$generate_surface) || length(input$generate_surface) == 0){
    showNotification(instructions$generate_no_surface,
                     type = "warning", duration = 5)
    return()
  }

  pred_list <- session_data$ellipsoid_prediction_list
  bias_list <- session_data$ellipsoid_prediction_list_biased

  selected_ids <- if(identical(input$generate_ellipsoid_selected, "all")){
    names(pred_list)
  } else {
    input$generate_ellipsoid_selected
  }

  target_layers <- input$generate_surface
  n_occ <- as.integer(input$generate_n_occ)
  sampling <- input$generate_sampling

  # Default to TRUE rather than NULL if the radio has not reported yet
  strict <- !identical(input$generate_strict, "FALSE")

  seed <- if(!is.null(input$generate_seed) && is.finite(input$generate_seed)){
    as.integer(input$generate_seed)
  } else {
    123L
  }

  sampling_mask <- if(!is.null(input$generate_mask_file)){
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

  mask_name <- if(!is.null(input$generate_mask_file)){
    input$generate_mask_file$name
  } else {
    "none"
  }

  # Assigned per id rather than replacing the whole list, so one failed
  # ellipsoid does not discard occurrences that already succeeded
  for(id in selected_ids){

    res <- generate_occ_for_ell(ell_id = id,
                                pred_list = pred_list,
                                biased_list = bias_list,
                                layers = target_layers,
                                n_occ = n_occ,
                                sampling = sampling,
                                strict = strict,
                                sampling_mask = sampling_mask,
                                seed = seed)

    # which layer names for this ellipsoid came from the biased prediction
    biased_layers <- if(!is.null(bias_list[[id]])){
      names(bias_list[[id]])
    } else {
      character(0)
    }

    n_attempted <- n_attempted + length(target_layers)

    for(layer in names(res)){

      df <- res[[layer]]
      if(is.null(df) || nrow(df) == 0) next

      biased <- layer %in% biased_layers

      # Metadata travels with the set so the summary can report how it was
      # made, and so the source raster can still be found from the layer
      attrs <- list(layer = layer,
                    n_occ = n_occ,
                    seed = seed,
                    sampling = if(biased) NA else sampling,
                    strict = strict,
                    mask = mask_name,
                    mode = "raster",
                    truncate = NA,
                    effect = NA)

      for(f in occ_set_fields) attr(df, f) <- attrs[[f]]

      attr(df, "created") <- format(Sys.time(), "%Y-%m-%d %H:%M")

      key <- occ_set_signature(attrs)

      if(is.null(session_data$ellipsoid_occurrence_list[[id]])){
        session_data$ellipsoid_occurrence_list[[id]] <- list()
      }

      session_data$ellipsoid_occurrence_list[[id]][[key]] <- df
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

  if(n_success == 0L){
    showNotification("Data generation failed for all selected combinations.",
                     type = "error", duration = 5)
    return()
  }

  generate_show_form(FALSE)

  n_skipped <- n_attempted - n_success

  msg <- paste0(n_success, " occurrence set(s) generated.")
  if(n_skipped > 0L) msg <- paste0(msg, " ", n_skipped, " skipped.")

  showNotification(msg, type = "message", duration = 4)
})

observeEvent(input$generate_edit_link, {
  generate_show_form(TRUE)
})

# OCCURRENCE SETS ---------------------------------------------------------

# Flat index of every occurrence set, so the summary can render rows and
# the download handlers can be created against stable numeric ids.
generate_occ_index <- reactive({

  ell <- session_data$current_ellipsoid
  if(is.null(ell)) return(NULL)

  occ <- session_data$ellipsoid_occurrence_list[[ell$ell_id]]
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
      sampling = occ_meta(df, "sampling", ""),
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if(length(rows) == 0) return(NULL)

  do.call(rbind, rows)
})


# One flattening function, used by all three download scopes
generate_occ_table <- function(idx_rows){

  occ <- session_data$ellipsoid_occurrence_list

  parts <- lapply(seq_len(nrow(idx_rows)), function(j){
    r <- idx_rows[j, ]
    df <- occ[[r$ell_id]][[r$set]]
    if(is.null(df) || nrow(df) == 0) return(NULL)
    data.frame(ell_id = r$ell_id,
               ell_name = r$ell_name,
               layer = r$layer,
               seed = r$seed,
               sampling = r$sampling,
               df,
               stringsAsFactors = FALSE)
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
    paste0(nm, "_occurrences.csv")
  },
  content = function(file){
    idx <- generate_occ_index()
    req(idx)
    write.csv(generate_occ_table(idx), file, row.names = FALSE)
  }
)

output$generate_dl_all <- downloadHandler(

  filename = function(){
    paste0("nicheR_occurrences_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file){

    occ <- session_data$ellipsoid_occurrence_list
    req(length(occ) > 0)

    versions <- session_data$ellipsoid_list

    parts <- lapply(names(occ), function(id){
      lapply(names(occ[[id]]), function(nm){
        df <- occ[[id]][[nm]]
        if(is.null(df) || nrow(df) == 0) return(NULL)
        data.frame(ell_id = id,
                   ell_name = if(!is.null(versions[[id]])) versions[[id]]$ell_name else id,
                   layer = occ_meta(df, "layer", nm),
                   n_occ = occ_meta(df, "n_occ"),
                   seed = occ_meta(df, "seed"),
                   sampling = occ_meta(df, "sampling", ""),
                   df,
                   stringsAsFactors = FALSE)
      })
    })

    out <- do.call(rbind, unlist(parts, recursive = FALSE))
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

  session_data$ellipsoid_occurrence_list[[r$ell_id]][[r$set]] <- NULL

  # Drop the ellipsoid entry entirely once its last set is gone, so the
  # library status returns to "not generated"
  if(length(session_data$ellipsoid_occurrence_list[[r$ell_id]]) == 0){
    session_data$ellipsoid_occurrence_list[[r$ell_id]] <- NULL
  }

  occ_visible(setdiff(occ_visible(), occ_vis_key(r$ell_id, r$set)))

  showNotification(paste0("Set removed: ", r$layer, " (seed ", r$seed, ")."),
                   type = "message", duration = 3)
})

observeEvent(input$generate_new_set_link, {
  generate_show_form(TRUE)
})

# OCCURRENCE SUMMARY ------------------------------------------------------

output$generate_occurrence_summary_ui <- renderUI({

  ell <- session_data$current_ellipsoid

  idx <- generate_occ_index()

  if(is.null(idx)){
    return(
      box(title = tags$span("Occurrence sets", class = "text-section-header"),
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

  box(title = tags$span(paste0("Occurrence sets: ", idx$ell_name[1]),
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
  generated <- names(session_data$ellipsoid_occurrence_list)

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
      list(txt = "occurrences generated", col = "#097a21")
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

# View, read-only
observeEvent(input$generate_ell_view, {

  ell <- session_data$ellipsoid_list[[input$generate_ell_view]]
  req(ell)

  set_working_ellipsoid(ell, mode = "view")

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
  session_data$ellipsoid_prediction_list_biased[[id]] <- NULL
  session_data$ellipsoid_occurrence_list[[id]] <- NULL
  session_data$pending_ell_delete <- NULL

  # Copies of the deleted ellipsoid, captured before reparenting so the
  # message reports only what this delete changed
  orphaned <- names(session_data$ellipsoid_list)[
    vapply(session_data$ellipsoid_list,
           function(e) identical(e$parent_id, id), logical(1))]

  session_data$ellipsoid_list <- lapply(session_data$ellipsoid_list, function(e){
    if(identical(e$parent_id, id)) e$parent_id <- NULL
    e
  })

  dbg("DELETE ", id, "  reparented to root: ",
      if(length(orphaned) == 0) "none" else paste(orphaned, collapse = ", "))

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
