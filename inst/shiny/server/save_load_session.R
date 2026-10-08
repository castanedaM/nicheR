# Title: Save and Load Session Logic

# Description: Serialises the whole session to an .rds file and restores it.
# terra objects are C++ pointers and cannot be saved directly, so every
# raster is wrapped on the way out and unwrapped on the way in.

# Date last updated: 10/06/2026


# Bumped when the saved structure changes in a way older files cannot
# satisfy. Load warns rather than failing silently on a mismatch.
SESSION_VERSION <- 2L

# Every session_data field that holds rasters, and how deep they sit.
# "single" is one raster, "list" is a named list of rasters.
SESSION_RASTER_FIELDS <- list(
  bg_raster = "single",
  bias_raster = "single",
  sampling_mask = "single",
  ellipsoid_prediction_list = "list",
  ellipsoid_prediction_list_biased = "list"
)


# HELPERS -----------------------------------------------------------------

# Wraps or unwraps every raster in the session list. `fun` is terra::wrap
# or terra::unwrap. Failures leave the object untouched rather than
# dropping it, so a partial restore is still inspectable.
session_convert_rasters <- function(session_list, fun){

  for(nm in names(SESSION_RASTER_FIELDS)){

    obj <- session_list[[nm]]
    if(is.null(obj)) next

    if(identical(SESSION_RASTER_FIELDS[[nm]], "single")){
      session_list[[nm]] <- tryCatch(fun(obj), error = function(e) obj)
      next
    }

    if(length(obj) == 0) next

    session_list[[nm]] <- lapply(obj, function(r){
      tryCatch(fun(r), error = function(e) r)
    })
  }

  # prepared_bias holds two rasters inside a plain list
  pb <- session_list$prepared_bias

  if(!is.null(pb)){
    session_list$prepared_bias <- list(
      combination_formula = pb$combination_formula,
      composite_surface = tryCatch(fun(pb$composite_surface),
                                   error = function(e) pb$composite_surface),
      processed_layers = tryCatch(fun(pb$processed_layers),
                                  error = function(e) pb$processed_layers)
    )
  }

  session_list
}

# Fields that describe the current UI state rather than the user's work.
# Excluded from the save so a reloaded session starts clean.
SESSION_SKIP_FIELDS <- c("pending_ell_delete", "file_type")


# SAVE --------------------------------------------------------------------

output$save_session_btn <- downloadHandler(

  filename = function(){
    paste0("nicheR_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
  },

  content = function(file){

    session_list <- reactiveValuesToList(session_data)

    for(nm in SESSION_SKIP_FIELDS) session_list[[nm]] <- NULL

    # The ellipsoid id counter lives outside session_data. Without it, ids
    # restart at E1 after a reload and collide with restored ellipsoids.
    session_list$.ell_id_counter <- ell_id_counter()
    session_list$.session_version <- SESSION_VERSION
    session_list$.saved_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    session_list <- session_convert_rasters(session_list, terra::wrap)

    tryCatch(
      saveRDS(session_list, file),
      error = function(e){
        showNotification(paste("Failed to save session:", e$message),
                         type = "error", duration = 5)
      }
    )
  }
)


# LOAD --------------------------------------------------------------------

observeEvent(input$build_load_session_btn, {

  req(input$build_session_file)
  req(is.data.frame(input$build_session_file))
  req(file.exists(input$build_session_file$datapath))

  session_list <- tryCatch(
    readRDS(input$build_session_file$datapath),
    error = function(e){
      showNotification(paste("Could not read session file:", e$message),
                       type = "error", duration = 5)
      NULL
    }
  )

  req(session_list)

  if(!is.list(session_list)){
    showNotification(instructions$build_session_invalid,
                     type = "error", duration = 5)
    return()
  }

  saved_version <- session_list$.session_version

  if(is.null(saved_version) || saved_version < SESSION_VERSION){
    showNotification(instructions$build_session_old_version,
                     type = "warning", duration = 8)
  }

  session_list <- session_convert_rasters(session_list, terra::unwrap)

  # Restore the id counter before anything else, so any ellipsoid built
  # after loading gets a fresh id rather than reusing a restored one
  restored_counter <- session_list$.ell_id_counter

  if(!is.null(restored_counter) && is.finite(restored_counter)){
    ell_id_counter(as.integer(restored_counter))
  } else {
    # Older file, or a corrupted counter: derive from the highest id present
    ids <- names(session_list$ellipsoid_list)
    nums <- suppressWarnings(as.integer(sub("^E([0-9]+)_.*$", "\\1", ids)))
    nums <- nums[is.finite(nums)]
    ell_id_counter(if(length(nums) > 0) max(nums) else 0L)
  }

  # Private fields are not part of session_data
  session_list$.ell_id_counter <- NULL
  session_list$.session_version <- NULL
  session_list$.saved_at <- NULL

  # Sessions saved before this list was renamed carry it under an old name
  for(old in c("ellipsoid_occurrence_list", "ellipsoid_records_list")){
    if(!is.null(session_list[[old]])){
      session_list$ellipsoid_point_sets_list <- session_list[[old]]
      session_list[[old]] <- NULL
    }
  }

  # Read before anything is restored, to tell below whether vars changed
  vars_before <- session_data$vars

  for(nm in names(session_list)){
    session_data[[nm]] <- session_list[[nm]]
  }

  # State that lives outside session_data has to be reset by hand
  session_data$pending_ell_delete <- NULL
  session_data$file_type <- NULL

  ell_mode("edit")
  covariance_set(FALSE)
  centroid_set(FALSE)

  # Nothing in the restored list is the working slot, so pick one: the
  # first saved ellipsoid, or nothing if the session had none
  ids <- names(session_data$ellipsoid_list)

  if(length(ids) > 0){
    session_data$current_ellipsoid <- session_data$ellipsoid_list[[ids[1]]]
    session_data$vars <- session_data$current_ellipsoid$var_names
    session_data$session_range <- session_data$current_ellipsoid$range_inputs
  } else {
    session_data$current_ellipsoid <- NULL
  }

  # A change of variables empties the working slot through an observer in
  # build_tab.R. That observer runs after this handler, so it would undo the
  # slot filled above. The flag lets this one change through. It is only set
  # when vars really changed, since otherwise the observer does not run and
  # the flag would be left on for the next change.
  if(!identical(vars_before, session_data$vars)){
    vars_reset_skip(TRUE)
  }

  # Forces every isolated panel to rebuild against the restored state
  ell_slot(ell_slot() + 1L)

  showNotification(instructions$build_session_loaded,
                   type = "message", duration = 4)

  # Land on the furthest step the restored session actually reached
  if(length(session_data$ellipsoid_point_sets_list) > 0){
    updateTabItems(session, "sidebar_menu", selected = "generate_tab")
  } else if(length(session_data$ellipsoid_prediction_list_biased) > 0){
    updateTabItems(session, "sidebar_menu", selected = "bias_tab")
  } else if(length(session_data$ellipsoid_prediction_list) > 0){
    updateTabItems(session, "sidebar_menu", selected = "predict_tab")
  } else if(length(ids) > 0){
    updateTabItems(session, "sidebar_menu", selected = "build_tab")
    updateTabsetPanel(session, "build_tabs", selected = "build_range_tab")
  } else {
    updateTabItems(session, "sidebar_menu", selected = "build_tab")
  }

})
