#' Apply sampling bias to prediction surfaces
#'
#' @description
#' Applies a prepared composite sampling bias surface to a prediction raster
#' by multiplication. The bias surface is aligned to the prediction grid when
#' needed and the result is cropped and masked to the prediction domain. The
#' output is a product of the prediction and the bias and is therefore no
#' longer interpretable as a probability.
#'
#' @usage apply_bias(prepared_bias, prediction, prediction_layer = NULL,
#'                   effect_direction = "direct", verbose = TRUE)
#'
#' @param prepared_bias A single-layer \code{SpatRaster} composite bias
#'   surface, or the list output from \code{\link{prepare_bias}} containing
#'   a \code{composite_surface} element.
#' @param prediction A \code{SpatRaster} containing one or more prediction
#'   layers, for example suitability or Mahalanobis distance.
#' @param prediction_layer Character. Name of the layer to extract from
#'   \code{prediction} when it contains multiple layers. If \code{NULL}
#'   (default) and \code{prediction} has a single layer, that layer is used.
#' @param effect_direction Character. How the prediction layer enters the
#'   product. \code{"direct"} (default) uses the prediction as is, so high
#'   prediction values increase sampling probability. \code{"inverse"}
#'   reflects the prediction around the midpoint of its own range,
#'   \eqn{(\max + \min) - x}, so high values decrease sampling probability.
#'   To reverse the bias surface instead, use \code{effect_direction} in
#'   \code{\link{prepare_bias}}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the composite bias surface from \code{prepared_bias}.
#'   \item Verifies the bias values are within \code{[0, 1]} and that the
#'   prediction values are non-negative.
#'   \item Aligns the bias surface to the prediction grid if geometries
#'   differ, using \code{terra::resample()} with nearest-neighbor
#'   interpolation.
#'   \item Reflects the prediction when \code{effect_direction = "inverse"}.
#'   \item Multiplies the (possibly reflected) prediction by the bias surface.
#'   \item Crops and masks the output to the prediction domain.
#' }
#'
#' The reflection is \eqn{(\max + \min) - x}, with \eqn{\min} and \eqn{\max}
#' taken from the layer being inverted. It preserves the units and the range
#' of the layer, does not rescale it, and reduces to \eqn{1 - x} when the
#' layer runs from 0 to 1. Because the bounds come from the raster supplied,
#' an inverted layer depends on the extent of that raster. Unbounded layers
#' such as Mahalanobis distance are the case where this matters most, since
#' the maximum is set by the most extreme cell in the study area.
#'
#' @return
#' A named list of class \code{"nicheR_biased_surface"} containing:
#' \itemize{
#'   \item One \code{SpatRaster} per input prediction layer, named
#'   \code{"<layer>_biased"}. The raster layer name includes the applied
#'   direction (e.g., \code{"suitability_biased_direct"}).
#'   \item \code{combination_formula}: a character string describing the
#'   operation applied (e.g., \code{"suitability * bias"} or
#'   \code{"(1 - suitability) * bias"}), including the reflection constant
#'   actually used.
#' }
#'
#' @seealso \code{\link{prepare_bias}} to build the composite bias surface,
#'   \code{\link{sample_biased_data}} to sample occurrences from the output.
#'
#' @importFrom terra compareGeom resample crop mask global nlyr
#'
#' @examples
#' pred_rast <- terra::rast(system.file("extdata/predictions_rast.tif",
#'                                      package = "nicheR"))
#'
#' bias_rast <- terra::rast(system.file("extdata/ma_biases.tif",
#'                                      package = "nicheR"))
#'
#' # 1. Prepare and standardize bias layers
#' bias <- prepare_bias(bias_surface = bias_rast[[1]],
#'                      effect_direction = "direct")
#'
#' # 2. Apply bias to the suitability layer
#' biased_pred <- apply_bias(prepared_bias = bias,
#'                           prediction = pred_rast,
#'                           prediction_layer = "suitability")
#'
#' terra::plot(biased_pred$suitability_biased)
#'
#' # 3. Sample away from the niche instead
#' biased_inv <- apply_bias(prepared_bias = bias,
#'                          prediction = pred_rast,
#'                          prediction_layer = "suitability",
#'                          effect_direction = "inverse")
#'
#' biased_inv$combination_formula
#'
#' @export
apply_bias <- function(prepared_bias,
                       prediction,
                       prediction_layer = NULL,
                       effect_direction = "direct",
                       verbose = TRUE){


  gc()
  verbose_message(verbose, "Starting: apply_bias()\n")

  # Basic Input checks --------------------------------------------------------

  if(missing(prepared_bias) || is.null(prepared_bias)){
    stop("'prepared_bias' must be provided.")
  }

  if(missing(prediction) || is.null(prediction)){
    stop("'prediction' must be provided.")
  }

  if(!inherits(prediction, "SpatRaster")){
    stop("'prediction' must be a terra::SpatRaster.")
  }

  prediction <- resolve_prediction(prediction, prediction_layer)$rast

  effect_direction <- match.arg(effect_direction,
                                choices = c("direct", "inverse"))


  # 1. Extract composite bias surface ----------------------------------------

  if(inherits(prepared_bias, "SpatRaster")){
    bias_rast <- prepared_bias

  }else if(is.list(prepared_bias)){

    if(!is.null(prepared_bias$composite_surface) &&
       inherits(prepared_bias$composite_surface, "SpatRaster")){
      bias_rast <- prepared_bias$composite_surface

    }else{
      # fallback: allow a named SpatRaster layer if user passed it directly
      stop("If 'prepared_bias' is a list, it must contain a SpatRaster named 'composite_surface'.")
    }

  }else{
    stop("'prepared_bias' must be a SpatRaster or the list output from prepare_bias().")
  }

  if(terra::nlyr(bias_rast) != 1){
    stop("'prepared_bias' must contain exactly 1 layer (a composite bias surface).")
  }

  # 2. Check bias range [0, 1] ----------------------------------------------

  bias_rng <- terra::global(bias_rast,
                            fun = c("min", "max"),
                            na.rm = TRUE)

  bias_min <- as.numeric(bias_rng[1, "min"])
  bias_max <- as.numeric(bias_rng[1, "max"])

  if(is.finite(bias_min) && bias_min < 0){
    stop("'prepared_bias' has values < 0. Bias must be standardized to [0, 1].")
  }

  if(is.finite(bias_max) && bias_max > 1){
    stop("'prepared_bias' has values > 1. Bias must be standardized to [0, 1].")
  }

  # 3. Prediction range ------------------------------------------------------

  # global() returns one row per layer, with min and max as columns. No upper
  # bound is imposed, the prediction can be suitability in [0, 1] or an
  # unbounded layer such as Mahalanobis distance. Negatives are rejected
  # because the product is used directly as a sampling weight.
  env_rng <- terra::global(prediction,
                           fun = c("min", "max"),
                           na.rm = TRUE)

  env_min <- suppressWarnings(min(env_rng[, "min"], na.rm = TRUE))

  if(is.finite(env_min) && env_min < 0){
    stop("'prediction' has values < 0. Expected non-negative prediction surfaces.")
  }

  # 4. Align bias to prediction grid -----------------------------------------

  same_grid <- terra::compareGeom(bias_rast,
                                  prediction[[1]],
                                  stopOnError = FALSE)

  if(!isTRUE(same_grid)){
    verbose_message(verbose, "Step: resampling prepared bias to match prediction grid...\n")
    bias_rast <- terra::resample(bias_rast,
                                 prediction[[1]],
                                 method = "near")
  }

  verbose_message(verbose, "Step: applying bias with '",
                  effect_direction,
                  "' effect to \"", paste(names(prediction), collapse = "\", \""),
                  "\" layer(s)...\n")

  # 5. Apply bias to each prediction layer -----------------------------------

  out_list <- vector("list", terra::nlyr(prediction))
  formula_entries <- character(terra::nlyr(prediction))


  for(i in 1:terra::nlyr(prediction)){

    s <- prediction[[i]]

    # Safe bias name
    bias_name <- names(bias_rast)
    if(is.null(bias_name) || length(bias_name) == 0L || !nzchar(bias_name[1]) || bias_name[1] == "lyr.1"){
      bias_name <- "bias"
    }else{
      bias_name <- bias_name[1]
    }

    # Safe prediction name
    suit_name <- names(s)
    if(is.null(suit_name) || length(suit_name) == 0L || !nzchar(suit_name[1]) || suit_name[1] == "lyr.1"){
      suit_name <- paste0("prediction_", i)
    }else{
      suit_name <- suit_name[1]
    }

    # Direction applies to the prediction, not to the bias. Inverting the
    # bias surface is prepare_bias()'s job, doing it here as well meant the
    # same flip could be applied twice without warning.
    if(effect_direction == "inverse"){

      lims <- c(env_rng[i, "min"], env_rng[i, "max"])

      if(any(!is.finite(lims))){
        stop("Cannot invert layer '", suit_name, "', its range is not finite.")
      }

      if(lims[1] == lims[2]){
        stop("Cannot invert layer '", suit_name, "', it is constant.")
      }

      # Reflection around the midpoint of the layer's own range. Preserves
      # units and range, no rescaling, and equals 1 - x on a [0, 1] layer.
      reflect <- lims[2] + lims[1]

      suit_effect <- reflect - s
      dir_tag <- "inverse"

      # The constant is data dependent, so it is recorded rather than
      # described, otherwise the formula cannot be reproduced later.
      suit_label <- paste0("(", signif(reflect, 6), " - ", suit_name, ")")

    }else{

      suit_effect <- s
      dir_tag <- "direct"
      suit_label <- suit_name
    }

    # Formula entry
    formula_entries[i] <- paste0(suit_label, " * ", bias_name)

    # Compute output raster
    out_r <- suit_effect * bias_rast
    out_r <- terra::crop(out_r, prediction[[1]], mask = TRUE)

    # Name list element + raster layer safely
    list_name <- paste0(suit_name, "_biased")
    if(!nzchar(list_name)){
      list_name <- paste0("prediction_", i, "_biased")
    }
    names(out_r) <- paste0(list_name, "_", dir_tag)

    out_list[[i]] <- out_r
    names(out_list)[i] <- list_name
  }


  # 6. Attach message / metadata ---------------------------------------------
  out_list$combination_formula <- formula_entries

  class(out_list) <- "nicheR_biased_surface"

  verbose_message(verbose, "Done: apply_bias(). Note: values are no longer probabilities\n")

  gc()

  out_list
}
