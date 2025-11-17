#' Combine one or more bias surfaces into a pooled bias raster
#'
#' This helper function combines one or more sampling bias layers into a single
#' pooled bias surface. Each bias layer is:
#' \itemize{
#'   \item aligned to a common template grid,
#'   \item scaled to \eqn{[0, 1]},
#'   \item optionally reversed via \eqn{1 - x},
#'   \item multiplied with other bias layers to form a pooled bias surface.
#' }
#'
#' The resulting pooled bias raster can then be used elsewhere in NicheR to
#' modulate sampling or suitability, but this function itself only handles
#' combination of bias layers.
#'
#' @param bias_surface One or more bias layers which may be:
#'   \itemize{
#'     \item a \code{terra::SpatRaster} (single layer or multi-layer),
#'     \item a \code{terra::SpatVector} (rasterized to a template grid),
#'     \item a character path to a raster or vector file,
#'     \item a \code{list} of any of the above.
#'   }
#' @param bias_dir Numeric values of \code{1} or \code{-1} controlling
#'   directionality. A length 1 value is recycled across all layers. A value
#'   of \code{-1} applies \code{1 - layer} after scaling to \eqn{[0, 1]}.
#' @param template Optional \code{terra::SpatRaster} used as the template grid
#'   (resolution, extent, and CRS) to which all bias layers will be aligned.
#'   If \code{NULL}, the template is inferred from the first layer in
#'   \code{bias_surface} (see Details).
#' @param res Optional numeric resolution used to build a template raster when
#'   no \code{template} is provided and only vector inputs are available, or
#'   when \code{ext} is provided. Ignored if \code{template} is supplied.
#' @param ext Optional spatial extent defining the template grid when no
#'   \code{template} is provided. Can be:
#'   \itemize{
#'     \item a \code{terra::SpatRaster} or \code{terra::SpatVector}, in which
#'           case \code{terra::ext()} is used, or
#'     \item a numeric vector \code{c(xmin, xmax, ymin, ymax)}.
#'   }
#'   Must be used together with \code{res} if no \code{template} is given.
#' @param out Character, one of \code{"default"} or \code{"both"}, determining
#'   whether to return only the pooled bias surface (\code{"default"}) or also
#'   the individual standardized, direction-adjusted bias layers
#'   (\code{"both"}).
#' @param verbose Logical. If \code{TRUE}, prints progress messages before and
#'   after major processing steps.
#'
#' @return A list of class \code{"nicheR_bias_surface"} containing:
#'   \item{pooled_bias_sp}{Pooled bias raster (scaled 0–1; \code{SpatRaster}).}
#'   \item{directional_bias_stack}{(If \code{out = "both"}) stack of all
#'         standardized, direction-corrected bias layers (\code{SpatRaster}).}
#'   \item{combination_formula}{String describing how layers were combined
#'         (e.g., \code{"Bias layers were combined: bias_1 * (1-bias_2)"}).}
#'
#' @details
#' Template selection follows this order:
#' \enumerate{
#'   \item If \code{template} is provided, it is used directly.
#'   \item Otherwise, if \code{ext} and \code{res} are provided, a new template
#'         raster is created from these.
#'   \item Otherwise, the first element of \code{bias_surface} is used:
#'     \itemize{
#'       \item If it is a \code{SpatRaster}, that raster becomes the template.
#'       \item If it is a \code{SpatVector}, a template is created from its
#'             extent and \code{res}. In this case, \code{res} must not be
#'             \code{NULL}.
#'     }
#' }
#'
#' All subsequent bias layers are resampled or rasterized to match this
#' template grid.
#'
#' @importFrom terra rast vect resample values rasterize app nlyr as.list ext res xmin xmax ymin ymax
#' @export
set_bias_surface <- function(bias_surface,
                             bias_dir = 1,
                             template = NULL,
                             res = NULL,
                             ext = NULL,
                             out = c("default", "both"),
                             verbose = TRUE) {

  gc()

  out <- match.arg(out)

  # Basic checks ---------------------------------------------------------------

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

  # ---------------------------------------------------------------------------
  # Determine template raster
  # ---------------------------------------------------------------------------

  if (!is.null(template)) {
    if (!inherits(template, "SpatRaster")) {
      stop("'template' must be a terra::SpatRaster if provided.")
    }
    template_raster <- template

    if (verbose) {
      message("Using user-supplied template raster.")
    }

  } else if (!is.null(ext) && !is.null(res)) {

    # ext can be SpatRaster, SpatVector, or numeric vector
    if (inherits(ext, c("SpatRaster", "SpatVector"))) {
      ext_obj <- terra::ext(ext)
    } else if (is.numeric(ext) && length(ext) == 4) {
      ext_obj <- terra::ext(ext[1], ext[2], ext[3], ext[4])
    } else {
      stop("'ext' must be a SpatRaster, SpatVector, or numeric c(xmin, xmax, ymin, ymax).")
    }

    template_raster <- terra::rast(ext = ext_obj, resolution = res)

    if (verbose) {
      message("Created template raster from provided extent and resolution.")
    }

  } else {
    # Infer template from first bias layer
    first_bias_obj <- bias_list[[1]]

    if (inherits(first_bias_obj, "character")) {
      first_bias_obj <- tryCatch(
        terra::rast(first_bias_obj),
        error = function(e) {
          tryCatch(
            terra::vect(first_bias_obj),
            error = function(e2) {
              stop("Could not load first bias layer to infer template: ", e2$message)
            }
          )
        }
      )
    }

    if (inherits(first_bias_obj, "SpatRaster")) {
      template_raster <- first_bias_obj
      if (verbose) {
        message("Using first bias raster as template.")
      }
    } else if (inherits(first_bias_obj, "SpatVector")) {
      if (is.null(res)) {
        stop("First bias layer is a SpatVector and no 'res' was provided. ",
             "Please supply 'res' (resolution) or an explicit 'template'.")
      }
      ext_obj <- terra::ext(first_bias_obj)
      template_raster <- terra::rast(ext = ext_obj, resolution = res)
      if (verbose) {
        message("First bias layer is a vector; created template from its extent and 'res'.")
      }
    } else {
      stop("Could not determine template from first bias layer.")
    }
  }

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
            " bias layer(s) on the template grid...")
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
        background = NA,
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
        message("  - Aligning bias layer ", i, " to template grid (resample).")
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
        if (all(is.na(x))) {
          NA_real_           # keep fully-missing cells as NA
        } else {
          prod(x, na.rm = TRUE)
        }
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

  # Build result object --------------------------------------------------------

  res <- list(
    pooled_bias_sp      = pooled_bias_sp,
    combination_formula = combination_formula
  )

  if (identical(tolower(out), "both")) {
    res$directional_bias_stack <- directional_bias_stack
  }

  class(res) <- "nicheR_bias_surface"

  gc()

  if (verbose) {
    message("set_bias_surface() completed successfully.")
  }

  return(res)
}

#' @export
print.nicheR_bias_surface <- function(x, ...) {
  if (!is.null(x$combination_formula)) {
    cat(x$combination_formula, "\n")
  }

  if (!is.null(x$pooled_bias_sp)) {
    cat("Pooled bias surface available as a terra::SpatRaster (pooled_bias_sp).\n")
  }

  if (!is.null(x$directional_bias_stack)) {
    cat("Directional (standardized, direction-adjusted) bias stack is also stored.\n")
  }

  invisible(x)
}
