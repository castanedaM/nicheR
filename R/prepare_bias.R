#' Prepare Sampling Bias Surfaces
#'
#' @export
prepare_bias <- function(bias_surface,
                         effect_direction = c("direct", "inverse"),
                         template_layer = NULL,
                         include_composite = TRUE,
                         include_processed_layers = FALSE,
                         mask_na = FALSE,
                         verbose = TRUE){

  gc()
  verbose_message(verbose, "Starting: prepare_bias()\n")

  # Basic Input checks --------------------------------------------------------

  if(missing(bias_surface) || is.null(bias_surface)){
    stop("'bias_surface' must be a SpatRaster or list of SpatRasters.")
  }

  if(!is.logical(mask_na) || length(mask_na) != 1L){
    stop("'mask_na' must be TRUE or FALSE.")
  }

  # Normalize bias_surface → list of single-layer rasters
  if(inherits(bias_surface, "SpatRaster")){
    verbose_message(verbose, "Step: splitting SpatRaster into layers...\n")
    bias_list <- terra::as.list(bias_surface)

    # If user did not provide template, use first layer
    if(is.null(template_layer)){
      template_layer <- bias_list[[1]]

      verbose_message(verbose,
                      "Step: bias_surface is a SpatRaster stack. Using first layer as template surface...\n")
    }

  }else if(is.list(bias_surface) &&
           all(vapply(bias_surface, inherits, logical(1), "SpatRaster"))){
    verbose_message(verbose, "Step: flattening list of SpatRasters...\n")
    bias_list <- unlist(lapply(bias_surface, terra::as.list), recursive = FALSE)
  }else{
    stop("'bias_surface' must be SpatRaster or list of SpatRasters.")
  }

  if(length(bias_list) == 0){
    stop("No usable bias layers provided.")
  }

  # 1. Determine template raster ------------------------------------------

  if(!is.null(template_layer) && inherits(template_layer, "SpatRaster") && !inherits(bias_surface, "SpatRaster")){

    verbose_message(verbose,
                    "Step: Using user-provided template_layer to align (crop/mask/resample) bias layers as needed...\n")
  }

  # choose template (if user didn't provide one)
  if(is.null(template_layer)){

    # pick finest resolution (smallest cell area)
    res_area <- vapply(bias_list, function(r){
      rr <- terra::res(r)
      rr[1] * rr[2]
    }, numeric(1))

    template_idx <- which.min(res_area)
    template_layer <- bias_list[[template_idx]]

    verbose_message(verbose,
                    paste0(
                      "Step: No template_layer provided. Using finest-resolution bias layer as template (index ",
                      template_idx, ").\n"
                    )
    )
  }

  # if union requested, expand template extent to union extent
  if(!isTRUE(mask_na)){

    union_ext <- Reduce(function(e, r) terra::union(e, terra::ext(r)),
                        x = bias_list,
                        init = terra::ext(bias_list[[1]]))

    template_layer <- terra::extend(template_layer, union_ext)

    verbose_message(verbose,
                    "Step: mask_na = FALSE. Expanding template to union extent (keeping finest resolution).\n"
    )
  }

  # Extent mismatch warning (do not stop)
  template_ext <- terra::ext(template_layer)

  extent_mismatch <- vapply(bias_list, function(r){
    !isTRUE(terra::ext(r) == template_ext)
  }, logical(1))

  if(any(extent_mismatch)){
    warning("Some bias layers have extents that do not match the template. They will be aligned during processing.", call. = FALSE)
  }

  # 2. Prepare effect direction ------------------------------------------

  effect_direction <- match.arg(effect_direction,
                                choices = c("direct", "inverse"),
                                several.ok = TRUE)

  if(length(effect_direction) == 1L){
    effect_direction <- rep(effect_direction, length(bias_list))
  }

  if(length(effect_direction) != length(bias_list)){
    stop("'effect_direction' must be length 1 or the same length as the number of bias layers.")
  }

  # Output logic checks ------------------------------------------------------

  if(!is.logical(include_composite) || length(include_composite) != 1L){
    stop("'include_composite' must be TRUE or FALSE.")
  }

  if(!is.logical(include_processed_layers) || length(include_processed_layers) != 1L){
    stop("'include_processed_layers' must be TRUE or FALSE.")
  }

  if(!include_composite && !include_processed_layers){
    include_composite <- TRUE
    verbose_message(verbose,
                    "Step: both 'include_composite' and 'include_processed_layers' were FALSE. Defaulting to include_composite = TRUE...\n"
    )
  }

  # 3. Process layers ------------------------------------------------------

  verbose_message(verbose,
                  "Step: standarizing (min/max) and applying direction of effect to ",
                  length(bias_list), " bias layer/s...\n"
  )

  directional_bias_list <- vector("list", length(bias_list))
  formula_entries <- character(length(bias_list))

  for(i in seq_along(bias_list)){

    raw <- bias_list[[i]]
    this_dir <- effect_direction[i]

    nm <- names(raw)
    nm <- if(is.null(nm) || nm == "" || nm == "lyr.1"){
      paste0("bias_", i)
    }else{
      nm[1]
    }

    same_grid <- terra::compareGeom(raw,
                                    template_layer,
                                    stopOnError = FALSE)

    if(!isTRUE(same_grid)){
      verbose_message(verbose, "Step: resampling bias layer ", i, " to match template surface...\n")
      aligned <- terra::resample(raw,
                                 template_layer,
                                 method = "near")
    }else{
      aligned <- raw
    }

    aligned <- terra::crop(aligned,
                           template_layer,
                           mask = TRUE)

    mm <- terra::minmax(aligned)
    min_val <- mm[1]
    max_val <- mm[2]

    if(!is.finite(min_val) || !is.finite(max_val) || (max_val - min_val) == 0){
      scaled <- aligned
      scaled[!is.na(scaled)] <- 1
    }else{
      scaled <- (aligned - min_val)/(max_val - min_val)
    }

    resampled_flag <- !isTRUE(same_grid)
    resample_tag <- if(resampled_flag) "_resampled" else ""

    if(this_dir == "inverse"){
      formula_entries[i] <- paste0("(1-", nm, ")")
      directional_bias_list[[i]] <- 1 - scaled
      names(directional_bias_list[[i]]) <- paste0("standarized_", nm, "_inverse", resample_tag)
    }else{
      formula_entries[i] <- nm
      directional_bias_list[[i]] <- scaled
      names(directional_bias_list[[i]]) <- paste0("standarized_", nm, "_direct", resample_tag)
    }
  }

  directional_bias_stack <- if(length(directional_bias_list) == 1L){
    directional_bias_list[[1]]
  }else{
    do.call(c, directional_bias_list)
  }

  # 4. Combine layers ----------------------------------------------------

  out_rast <- list()

  if(isTRUE(include_composite)){

    verbose_message(
      verbose,
      paste0(
        "Step: building standarized (min/max) directional composite bias surface (mask_na = ",
        mask_na, ")...\n"
      )
    )

    if(terra::nlyr(directional_bias_stack) > 1){

      if(isTRUE(mask_na)){
        # INTERSECTION: any NA in a pixel makes the composite NA
        composite_raster <- terra::app(directional_bias_stack, fun = function(x){
          if(any(is.na(x))) NA_real_ else prod(x)
        })
        out_rast$combination_formula <- paste(formula_entries, collapse = " * ")

      } else {
        # UNION: ignore NA and combine what exists
        composite_raster <- terra::app(directional_bias_stack, fun = function(x){
          if(all(is.na(x))){
            NA_real_
          } else if(sum(!is.na(x)) == 1L){
            x[!is.na(x)]
          } else {
            prod(x[!is.na(x)])
          }
        })
        out_rast$combination_formula <- paste(formula_entries, collapse = " * ")
      }

    } else {
      composite_raster <- directional_bias_stack
      out_rast$combination_formula <- formula_entries[1]
    }

    names(composite_raster) <- "standarized_composite_bias_surface"
    out_rast$composite_surface <- composite_raster

    if(isTRUE(include_processed_layers)){
      out_rast$processed_layers <- directional_bias_stack
    }

  } else {
    out_rast$processed_layers <- directional_bias_stack
    out_rast$combination_formula <- paste(formula_entries, collapse = " * ")
  }

  # 5. Build result --------------------------------------------------------

  verbose_message(verbose, "Done: prepare_bias()\n")
  gc()

  out_rast
}
