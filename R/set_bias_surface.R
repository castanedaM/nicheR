#' Apply one or more bias surfaces to a suitability raster
#'
#' This function combines one or more bias layers (sampling probability
#' surfaces) and applies them to a suitability raster. Each bias layer is:
#' \itemize{
#'   \item aligned to the suitability raster,
#'   \item scaled to \eqn{[0, 1]},
#'   \item optionally reversed via \eqn{1 - x},
#'   \item multiplied with other bias layers to form a pooled bias surface.
#' }
#' The pooled bias is then multiplied with the ecological suitability raster to
#' obtain a biased suitability surface.
#'
#' All operations are performed with \pkg{terra} objects; no data.frame
#' conversion is used internally to keep memory usage manageable.
#'
#' @param suitable_env A required \code{terra::SpatRaster} representing the
#'   ecological suitability surface (values typically 0–1).
#' @param bias_surface One or more bias layers which may be:
#'   \itemize{
#'     \item a \code{terra::SpatRaster} (single layer or multi-layer),
#'     \item a \code{terra::SpatVector} (rasterized to match
#'           \code{suitable_env}),
#'     \item a character path to a raster or vector file,
#'     \item a \code{list} of any of the above.
#'   }
#' @param bias_dir Numeric values of \code{1} or \code{-1} controlling
#'   directionality. A length 1 value is recycled across all layers. A value
#'   of \code{-1} applies \code{1 - layer} after scaling to \eqn{[0, 1]}.
#' @param out Character, one of \code{"default"} or \code{"both"}, determining
#'   whether to return only the pooled bias and final suitability
#'   (\code{"default"}) or also the individual standardized, direction-adjusted
#'   bias layers (\code{"both"}).
#' @param verbose Logical. If \code{TRUE}, prints progress messages before and
#'   after major processing steps.
#'
#' @return A list of class \code{"biased_suitable_env"} containing:
#'   \item{pooled_bias_sp}{Pooled bias raster (scaled 0–1; \code{SpatRaster}).}
#'   \item{final_suitability_sp}{Suitability * pooled bias (\code{SpatRaster}).}
#'   \item{directional_bias_stack}{(If \code{out = "both"}) stack of all
#'         standardized, direction-corrected bias layers (\code{SpatRaster}).}
#'   \item{combination_formula}{String describing how layers were combined
#'         (e.g., \code{"Bias layers were combined: bias_1 * (1-bias_2)"}).}
#'   \item{is_niche_applied}{Logical, always \code{TRUE} for this function
#'         (reserved for potential future use).}
#'
#' @family niche functions
#' @importFrom terra rast vect resample values rasterize app nlyr as.list res ext xmin xmax ymin ymax
#' @export
set_bias_surface <- function(suitable_env,
                             bias_surface,
                             bias_dir = 1,
                             out = c("default", "both"),
                             verbose = TRUE) {

  gc()

  out <- match.arg(out)

  # Basic checks ---------------------------------------------------------------

  if (missing(suitable_env) || is.null(suitable_env)) {
    stop("'suitable_env' is required and must be a terra::SpatRaster.")
  }
  if (!inherits(suitable_env, "SpatRaster")) {
    stop("'suitable_env' must be a terra::SpatRaster.")
  }
  if (missing(bias_surface) || is.null(bias_surface)) {
    stop("'bias_surface' is required (SpatRaster, SpatVector, path, or list).")
  }

  # Normalize bias_surface into a list -----------------------------------------

  if (inherits(bias_surface, "list")) {
    bias_list <- bias_surface
  } else if (inherits(bias_surface, "SpatRaster")) {
    if (terra::nlyr(bias_surface) > 1) {
      bias_list <- terra::as.list(bias_surface)
    } else {
      bias_list <- list(bias_surface)
    }
  } else if (inherits(bias_surface, c("SpatVector", "character"))) {
    bias_list <- list(bias_surface)
  } else {
    stop("'bias_surface' must be a SpatRaster, SpatVector, path, or a list of these.")
  }

  if (length(bias_list) == 0) {
    stop("No bias layers provided.")
  }

  # Template is always the suitability raster ----------------------------------

  template_raster  <- suitable_env
  is_niche_applied <- TRUE

  # Validate bias_dir ----------------------------------------------------------

  if (length(bias_dir) == 1) {
    bias_dir <- rep(bias_dir, length(bias_list))
  } else if (length(bias_dir) != length(bias_list)) {
    stop("The length of 'bias_dir' (", length(bias_dir),
         ") must match the number of bias layers (", length(bias_list),
         ") or be 1.")
  }
  if (!all(bias_dir %in% c(1, -1))) {
    stop("'bias_dir' must only contain values of 1 or -1.")
  }

  # Process, align, and standardize each bias layer ----------------------------

  if (verbose) {
    message("Processing and standardizing ", length(bias_list),
            " bias layer(s) to match the suitability raster...")
  }

  directional_bias_list <- vector("list", length(bias_list))
  dir_message_parts     <- character(length(bias_list))

  for (i in seq_along(bias_list)) {

    bias_obj <- bias_list[[i]]
    this_dir <- bias_dir[i]

    # Try to derive a readable name
    layer_name <- names(bias_list[[i]])
    if (is.null(layer_name) || length(layer_name) == 0 || layer_name == "") {
      layer_name <- paste0("bias_", i)
    }

    # Load from path if needed
    if (inherits(bias_obj, "character")) {
      bias_obj <- tryCatch(
        terra::rast(bias_obj),
        error = function(e) {
          tryCatch(
            terra::vect(bias_obj),
            error = function(e2) {
              stop("Error loading bias layer ", i, ": ", e2$message)
            }
          )
        }
      )
    }

    # Convert SpatVector → SpatRaster using template
    if (inherits(bias_obj, "SpatVector")) {
      flds <- names(bias_obj)
      field_to_use <- if (length(flds) == 0) 1 else flds[1]
      if (length(flds) == 0) {
        warning("Bias layer ", i,
                " (SpatVector) has no attributes. Rasterizing as binary (1 = presence).",
                call. = FALSE)
      }
      bias_raster_raw <- terra::rasterize(
        bias_obj,
        template_raster,
        fun        = "max",
        background = 0,
        field      = field_to_use
      )
    } else if (inherits(bias_obj, "SpatRaster")) {
      bias_raster_raw <- bias_obj
    } else {
      stop("Bias layer ", i, " must be a SpatRaster, SpatVector, or a valid path.")
    }

    # Align resolution + extent to template
    if (!identical(terra::res(bias_raster_raw), terra::res(template_raster)) ||
        !identical(terra::ext(bias_raster_raw), terra::ext(template_raster))) {

      if (verbose) {
        message("  - Aligning bias layer ", i, " to suitability grid (resample).")
      }
      bias_raster_aligned <- terra::resample(bias_raster_raw, template_raster, method = "near")
    } else {
      bias_raster_aligned <- bias_raster_raw
    }

    # Normalize to [0, 1]
    vals      <- terra::values(bias_raster_aligned)
    min_val   <- min(vals, na.rm = TRUE)
    max_val   <- max(vals, na.rm = TRUE)
    range_val <- max_val - min_val

    if (range_val == 0) {
      warning("Bias layer ", i,
              " contains a single unique value. Setting all non-NA values to 1.",
              call. = FALSE)
      scaled_layer <- bias_raster_aligned
      scaled_layer[!is.na(scaled_layer)] <- 1
    } else {
      scaled_layer <- (bias_raster_aligned - min_val) / range_val
    }

    # Apply directionality
    if (this_dir == -1) {
      directional_layer    <- 1 - scaled_layer
      dir_message_parts[i] <- paste0("(1-", layer_name, ")")
    } else {
      directional_layer    <- scaled_layer
      dir_message_parts[i] <- layer_name
    }

    names(directional_layer)      <- layer_name
    directional_bias_list[[i]]    <- directional_layer
  }

  if (verbose) {
    message("Finished processing bias layers.")
  }

  # Stack directional layers ---------------------------------------------------

  directional_bias_stack <- if (length(directional_bias_list) == 1) {
    directional_bias_list[[1]]
  } else {
    do.call(c, directional_bias_list)
  }

  # Pool bias layers -----------------------------------------------------------

  if (verbose) {
    message("Pooling ", terra::nlyr(directional_bias_stack),
            " bias layer(s) into a single surface...")
  }

  if (terra::nlyr(directional_bias_stack) > 1) {
    pooled_bias_sp <- terra::app(
      directional_bias_stack,
      fun = function(x) {
        prod(x, na.rm = TRUE)
      }
    )
  } else {
    pooled_bias_sp <- directional_bias_stack
  }

  names(pooled_bias_sp) <- "pooled_bias"

  combination_formula <- paste0(
    "Bias layers were combined: ",
    paste(dir_message_parts, collapse = " * ")
  )

  if (verbose) {
    message("Finished pooling bias layers.")
  }

  # Final suitability ----------------------------------------------------------

  if (verbose) {
    message("Applying pooled bias to suitability raster...")
  }

  final_suitability_sp <- suitable_env * pooled_bias_sp
  names(final_suitability_sp) <- "final_suitability"

  if (verbose) {
    message("Finished computing biased suitability surface.")
  }

  # Build result object --------------------------------------------------------

  res <- list(
    pooled_bias_sp       = pooled_bias_sp,
    final_suitability_sp = final_suitability_sp,
    combination_formula  = combination_formula,
    is_niche_applied     = is_niche_applied
  )

  if (identical(tolower(out), "both")) {
    res$directional_bias_stack <- directional_bias_stack
  }

  class(res) <- "biased_suitable_env"

  gc()

  if (verbose) {
    message("set_bias_surface() completed successfully.")
  }

  return(res)
}

#' @export
print.biased_suitable_env <- function(x, ...) {
  if (!is.null(x$combination_formula)) {
    cat(x$combination_formula, "\n")
  }

  if (isTRUE(x$is_niche_applied)) {
    cat("Ecological suitability provided: final suitability = Suitability * Pooled Bias.\n")
  } else {
    cat("Ecological suitability is NULL (this should not occur for set_bias_surface()).\n")
  }

  invisible(x)
}
