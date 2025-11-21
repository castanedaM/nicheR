#' Construct and Combine Sampling Bias Surfaces (Optimized, nr_get-integrated)
#'
#' This function ingests one or more user-provided bias layers and produces
#' standardized, optionally directional, and pooled sampling-bias surfaces.
#' It is designed for use within the NicheR virtual-species workflow and is
#' fully compatible with any NicheR object structure through the internal
#' accessor helper \code{nr_get()}.
#'
#' The function performs several optimized steps:
#' \enumerate{
#'   \item Crops each bias layer to the extent of the suitable environment
#'         (if provided) — a major performance improvement.
#'   \item Resamples each layer to match the resolution, extent, and grid of
#'         the template raster.
#'   \item Standardizes each raster to the \eqn{[0, 1]} range using
#'         \code{terra::minmax()}.
#'   \item Optionally inverts layers when \code{bias_dir = -1}.
#'   \item Combines all directional layers multiplicatively to form a pooled
#'         bias surface.
#'   \item Optionally masks outputs to the suitable environment extent.
#' }
#'
#' This procedure produces bias surfaces suitable for passed directly to
#' \code{\link{get_sample_occ}}, or for standalone inspection within the
#' NicheR workflow.
#'
#' @section Integration with NicheR:
#' Unlike older versions, this function does not try to decode nested object
#' structures manually. Instead, it uses \code{nr_get()} internally to extract
#' a template raster from any of the following:
#' \itemize{
#'   \item suitability objects returned by \code{\link{get_suitable_env}},
#'   \item full \code{NicheR_species} objects from
#'         \code{\link{create_virtual_species}},
#'   \item lists containing SpatRasters,
#'   \item direct SpatRaster input.
#' }
#' This ensures compatibility across all internal and user-facing NicheR
#' functions.
#'
#' @param bias_surface A \code{SpatRaster} or a list of \code{SpatRaster}
#'        objects. Multi-layer rasters are split internally into individual
#'        single-layer bias surfaces.
#'
#' @param bias_dir Numeric vector of \code{1} (use layer as-is) or \code{-1}
#'        (invert layer via \eqn{1 - x}). A single value is recycled to match
#'        the number of bias layers.
#'
#' @param suitable_env Optional. A suitability object, \code{NicheR_species}
#'        object, or any structure containing a suitable environment raster.
#'        When \code{out.bias} includes \code{"biased"}, a suitable raster is
#'        required and is extracted automatically using \code{nr_get()}.
#'
#' @param out.bias Character. Controls which bias outputs are returned:
#'        \itemize{
#'          \item \code{"standardized"} – only directional standardized layers,
#'          \item \code{"biased"} – only the pooled bias surface,
#'          \item \code{"both"} – return both types of output.
#'        }
#'
#' @param verbose Logical. If \code{TRUE}, print progress information.
#'
#' @param crop_to_suitable Logical. If \code{TRUE} (default), each bias layer
#'        is cropped to match the extent of the suitable environment before
#'        resampling. This provides substantial speed improvements.
#'
#' @return A list of class \code{nicheR_bias_surface} with elements:
#' \itemize{
#'   \item \code{pooled_bias_sp} – a single \code{SpatRaster} (if requested),
#'   \item \code{directional_bias_stack} – multilayer raster stack of all
#'          standardized + directional layers (if requested),
#'   \item \code{combination_formula} – a character string showing how layers
#'          were multiplied, e.g. \code{"roads * (1-pop_density) * protected"}.
#' }
#'
#' @details
#' The pooled bias raster is always masked to the suitable environment when
#' \code{suitable_env} is provided. Standardized layers are also masked unless
#' \code{crop_to_suitable = FALSE}.
#'
#' All raster alignment uses \code{terra::resample()} with a nearest-neighbor
#' method to preserve sharp edges typical of bias layers (e.g., distances,
#' road networks).
#'
#' @seealso
#' \code{\link{nr_get}} – universal accessor for NicheR structures
#' \code{\link{create_virtual_species}} – full virtual species pipeline
#' \code{\link{get_sample_occ}} – uses pooled bias surfaces for sampling
#'
#' @export
set_bias_surface <- function(bias_surface,
                                 bias_dir = 1,
                                 suitable_env = NULL,
                                 out.bias = c("biased", "standardized", "both"),
                                 verbose = TRUE,
                                 crop_to_suitable = TRUE) {

  gc()
  out.bias <- match.arg(out.bias)

  if (isTRUE(verbose)) message("Starting set_bias_surface()...")


  # 0. Input checks

  if (missing(bias_surface) || is.null(bias_surface))
    stop("'bias_surface' must be a SpatRaster or list of SpatRasters.")

  # Normalize bias_surface → list of single-layer rasters
  if (inherits(bias_surface, "SpatRaster")) {
    if (verbose) message("bias_surface is SpatRaster → splitting into layers...")
    bias_list <- terra::as.list(bias_surface)

  } else if (is.list(bias_surface) &&
             all(vapply(bias_surface, inherits, logical(1), "SpatRaster"))) {

    if (verbose) message("bias_surface is list → flattening into layers...")
    bias_list <- unlist(lapply(bias_surface, terra::as.list), recursive = FALSE)

  } else {
    stop("'bias_surface' must be SpatRaster or list of SpatRasters.")
  }

  if (length(bias_list) == 0)
    stop("No usable bias layers provided.")


  # 1. Determine mask / template raster using nr_get()

  mask_ras <- NULL

  if (out.bias %in% c("biased", "both")) {
    mask_ras <- nr_get(suitable_env, "suitable")
    if (is.null(mask_ras)) {
      stop("out.bias = '", out.bias, "' requires 'suitable_env' containing a suitable raster.")
    }
    if (verbose) message("Using suitable_env (via nr_get) as template/mask.")
  }

  template_raster <- if (!is.null(mask_ras)) mask_ras else bias_list[[1]]


  # 2. Prepare bias_dir

  if (length(bias_dir) == 1) bias_dir <- rep(bias_dir, length(bias_list))
  if (!all(bias_dir %in% c(1, -1)))
    stop("bias_dir must contain only 1 or -1.")

  directional_bias_list <- vector("list", length(bias_list))
  formula_entries       <- character(length(bias_list))

  if (verbose) message("Processing ", length(bias_list), " bias layers...")


  # 3. Process each layer: crop → resample → scale → invert

  for (i in seq_along(bias_list)) {

    raw <- bias_list[[i]]
    this_dir <- bias_dir[i]
    nm       <- names(raw)
    nm       <- if (is.null(nm) || nm == "") paste0("bias_", i) else nm[1]

    # Crop before resampling (huge speed-up)
    if (crop_to_suitable && !is.null(mask_ras)) {
      if (!identical(terra::ext(raw), terra::ext(template_raster))) {
        if (verbose) message("  • Cropping ", nm, " to suitable_env extent...")
        raw <- terra::crop(raw, template_raster)
      }
    }

    # Resample if needed
    needs_resample <- (!identical(terra::res(raw), terra::res(template_raster)) ||
                         !identical(terra::ext(raw), terra::ext(template_raster)))

    aligned <- if (needs_resample) {
      if (verbose) message("  • Resampling ", nm, "...")
      terra::resample(raw, template_raster, method = "near")
    } else raw

    # Standardize using terra::minmax()
    mm <- terra::minmax(aligned)
    min_val <- mm[1]; max_val <- mm[2]

    if (!is.finite(min_val) || !is.finite(max_val)) {
      scaled <- aligned; scaled[!is.na(scaled)] <- 1
    } else if ((max_val - min_val) == 0) {
      scaled <- aligned; scaled[!is.na(scaled)] <- 1
    } else {
      scaled <- (aligned - min_val) / (max_val - min_val)
    }

    # Apply direction (invert)
    directional <- if (this_dir == -1) {
      formula_entries[i] <- paste0("(1-", nm, ")")
      1 - scaled
    } else {
      formula_entries[i] <- nm
      scaled
    }

    names(directional) <- nm
    directional_bias_list[[i]] <- directional
  }

  # Stack layers
  directional_bias_stack <- if (length(directional_bias_list) == 1)
    directional_bias_list[[1]]
  else
    do.call(c, directional_bias_list)


  # 4. Combine layers (product)

  pooled_bias_sp <- NULL

  if (out.bias %in% c("biased", "both")) {
    if (verbose) message("Pooling directional layers...")

    pooled_bias_sp <- if (terra::nlyr(directional_bias_stack) > 1) {
      terra::app(directional_bias_stack, fun = function(x) {
        if (all(is.na(x))) NA_real_ else prod(x, na.rm = TRUE)
      })
    } else directional_bias_stack

    names(pooled_bias_sp) <- "pooled_bias"

    # Mask to suitable_env
    if (!is.null(mask_ras))
      pooled_bias_sp <- terra::mask(pooled_bias_sp, mask_ras)
  }

  # Mask the standardized stack as well
  if (!is.null(mask_ras))
    directional_bias_stack <- terra::mask(directional_bias_stack, mask_ras)


  # 5. Build result

  res <- list(
    pooled_bias_sp         = if (out.bias %in% c("biased", "both")) pooled_bias_sp else NULL,
    directional_bias_stack = if (out.bias %in% c("standardized", "both")) directional_bias_stack else NULL,
    combination_formula    = paste(formula_entries, collapse = " * ")
  )

  class(res) <- "nicheR_bias_surface"

  if (verbose) message("Completed set_bias_surface().")
  gc()
  return(res)
}
