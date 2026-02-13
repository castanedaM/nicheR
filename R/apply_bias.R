#' Construct and Combine Sampling Bias Surfaces
#'
#' Builds standardized, optionally directional, and pooled sampling bias
#' surfaces from one or more user-provided rasters. Bias layers are aligned to
#' a common template, scaled to \eqn{[0, 1]}, optionally inverted, and combined
#' multiplicatively to produce a pooled bias surface.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Split \code{bias_surface} into single-layer rasters.
#'   \item Select a template raster: \code{suitable_env} (if provided) or the
#'   first bias layer.
#'   \item Optionally crop each bias layer to the template extent
#'   (\code{truncated = TRUE}).
#'   \item Resample each layer to match the template grid.
#'   \item Standardize each layer to \eqn{[0, 1]} using min and max values.
#'   \item Optionally invert layers when \code{bias_dir = -1}.
#'   \item Stack standardized directional layers.
#'   \item Optionally pool layers by multiplication to form a single bias
#'   surface.
#' }
#'
#' Cropping and masking are done using \code{terra::crop()} and
#' \code{terra::mask()}. Resampling uses nearest-neighbor interpolation
#' (\code{method = "near"}) to preserve sharp edges common in bias layers.
#'
#' @param bias_surface A \code{terra::SpatRaster} or a list of
#'   \code{terra::SpatRaster} objects. Multi-layer rasters are split internally
#'   into single-layer bias surfaces.
#' @param bias_dir Numeric vector of \code{1} (use layer as-is) or \code{-1}
#'   (invert via \eqn{1 - x}). A single value is recycled to match the number of
#'   bias layers.
#' @param suitable_env Optional. A \code{terra::SpatRaster} used as a
#'   mask/template raster for cropping, resampling, and masking.
#' @param out_bias Character. Controls returned outputs:
#'   \itemize{
#'     \item \code{"biased"} – pooled bias surface only
#'     \item \code{"standardized"} – directional standardized layers only
#'     \item \code{"both"} – return both
#'   }
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#' @param truncated Logical. If \code{TRUE} (default), each bias layer is cropped
#'   to the template extent before resampling. If \code{FALSE}, layers are only
#'   resampled and masked.
#'
#' @return
#' A list of class \code{"nicheR_bias_surface"} with:
#' \itemize{
#'   \item \code{pooled_bias_sp} – pooled bias raster (if requested)
#'   \item \code{directional_bias_stack} – stack of standardized and directional
#'   layers (if requested)
#'   \item \code{combination_formula} – character string showing how layers were
#'   combined (e.g., \code{"roads * (1-pop_density) * protected"})
#' }
#'
#' @seealso \code{\link{sample_data}}
#'
#' @export
apply_bias <- function(bias_surface,
                       bias_dir = 1,
                       suitable_env = NULL,
                       out_bias = c("biased", "standardized", "both"),
                       verbose = TRUE,
                       truncated = TRUE){
  # TO DO:
  # prepare(scaling, bias_direction) ->
  # apply(prepared_bias, suitable_layers, effect_direction) -> crop and mask for the layers
  # bias_sampling(n_occ, layer, seed = 1, sampling_mask = NULL)


# TO DO: Change bias
  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)

  gc()

  out_bias <- match.arg(out_bias)

  verbose_message("Starting: apply_bias()\n")

  # Basic Input checks --------------------------------------------------------

  if(missing(bias_surface) || is.null(bias_surface)){
    stop("'bias_surface' must be a SpatRaster or list of SpatRasters.")
  }

  if(!is.logical(truncated) || length(truncated) != 1L){
    stop("'truncated' must be TRUE or FALSE.")
  }

  # Normalize bias_surface → list of single-layer rasters
  if(inherits(bias_surface, "SpatRaster")){

    verbose_message("Step: splitting SpatRaster into layers...\n")
    bias_list <- terra::as.list(bias_surface)

  }else if(is.list(bias_surface) &&
           all(vapply(bias_surface, inherits, logical(1), "SpatRaster"))){

    verbose_message("Step: flattening list of SpatRasters...\n")
    bias_list <- unlist(lapply(bias_surface, terra::as.list), recursive = FALSE)

  }else{
    stop("'bias_surface' must be SpatRaster or list of SpatRasters.")
  }

  if(length(bias_list) == 0){
    stop("No usable bias layers provided.")
  }

  # 1. Determine template raster ------------------------------------------

  if(!is.null(suitable_env) && inherits(suitable_env, "SpatRaster") && isTRUE(truncated)){
    mask_ras <- suitable_env
    verbose_message("Step: Using suitable_env as mask/template...\n")
  }else{
    mask_ras <- bias_list[[1]]
    verbose_message("Step: Using first bias layer as mask/template...\n")
  }

  # 2. Prepare bias_dir ----------------------------------------------------

  if(length(bias_dir) == 1){
    bias_dir <- rep(bias_dir, length(bias_list))
  }

  if(!all(bias_dir %in% c(1, -1))){
    stop("bias_dir must contain only 1 or -1.")
  }

  directional_bias_list <- vector("list", length(bias_list))
  formula_entries <- character(length(bias_list))

  verbose_message("Step: processing ", length(bias_list), " bias layers...\n")

  # 3. Process layers ------------------------------------------------------

  for(i in seq_along(bias_list)){

    raw <- bias_list[[i]]
    this_dir <- bias_dir[i]

    nm <- names(raw)
    nm <- if(is.null(nm) || nm == "") paste0("bias_", i) else nm[1]

    # Crop before resampling
    if(isTRUE(truncated)){
      raw <- terra::crop(raw, mask_ras, mask = TRUE)
    }

    # Resample if needed
    needs_resample <- (!identical(terra::res(raw), terra::res(mask_ras)) ||
                         !identical(terra::ext(raw), terra::ext(mask_ras)))

    aligned <- if(needs_resample){
      terra::resample(raw, mask_ras, method = "near")
    }else{
      raw
    }

    # Standardize
    mm <- terra::minmax(aligned)
    min_val <- mm[1]
    max_val <- mm[2]

    if(!is.finite(min_val) || !is.finite(max_val)){
      scaled <- aligned
      scaled[!is.na(scaled)] <- 1
    }else if((max_val - min_val) == 0){
      scaled <- aligned
      scaled[!is.na(scaled)] <- 1
    }else{
      scaled <- (aligned - min_val)/(max_val - min_val)
    }

    # Apply direction
    directional <- if(this_dir == -1){
      formula_entries[i] <- paste0("(1-", nm, ")")
      1 - scaled
    }else{
      formula_entries[i] <- nm
      scaled
    }

    names(directional) <- nm
    directional_bias_list[[i]] <- directional
  }

  # Stack layers ----------------------------------------------------------

  directional_bias_stack <- if(length(directional_bias_list) == 1){
    directional_bias_list[[1]]
  }else{
    do.call(c, directional_bias_list)
  }

  # 4. Combine layers ------------------------------------------------------

  res <- list()

  if(out_bias %in% c("biased", "both")){

    verbose_message("Step: pooling directional layers...\n")

    pooled_bias_sp <- if(terra::nlyr(directional_bias_stack) > 1){
      terra::app(directional_bias_stack, fun = function(x){
        if(all(is.na(x))) NA_real_ else prod(x, na.rm = TRUE)
      })
    }else{
      directional_bias_stack
    }

    names(pooled_bias_sp) <- "pooled_bias"

    pooled_bias_sp <- terra::mask(pooled_bias_sp, mask_ras)

    res$pooled_bias_sp <- pooled_bias_sp
  }

  directional_bias_stack <- terra::mask(directional_bias_stack, mask_ras)

  if(out_bias %in% c("standardized", "both")){
    res$directional_bias_stack <- directional_bias_stack
  }

  # 5. Build result --------------------------------------------------------

  res$combination_formula <- paste(formula_entries, collapse = " * ")

  class(res) <- "nicheR_bias_surface"

  verbose_message("Done: apply_bias()\n")

  gc()

  res
}

