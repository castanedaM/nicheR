#' Apply Sampling Bias to a Prediction Surface
#'
#' Modifies a prediction surface by applying a spatial bias layer.
#' The bias layer is dynamically cropped and masked to match the prediction
#' layer, adjusted for effect direction, and combined multiplicatively.
#'
#' @param prepared_bias A \code{nicheR_bias_surface} object (output from
#'   \code{prepare_bias}) or a \code{terra::SpatRaster}.
#' @param pred A \code{terra::SpatRaster} representing the prediction surface
#'   (e.g., habitat suitability).
#' @param effect_direction Numeric vector. \code{1} for a direct effect (\code{pred * bias}),
#'   or \code{-1} for an indirect/inverse effect (\code{pred * (1 - bias)}).
#'   Can contain multiple values (e.g., \code{c(1, -1)}). If multiple values are provided,
#'   the function returns layers for every combination of valid prediction layer
#'   and effect direction. Default is \code{1}.
#' @param layers Character vector. Names of the layers in \code{pred} to which
#'   the bias should be applied. Defaults to \code{c("suitability", "suitability_trunc")}.
#'
#' @return A \code{terra::SpatRaster} containing the bias-adjusted prediction layers.
#' @export
apply_bias <- function(prepared_bias,
                       pred,
                       effect_direction = 1,
                       layers = c("suitability", "suitability_trunc")) {

  # 1. Input Validation
  if (!inherits(pred, "SpatRaster")) {
    stop("'pred' must be a terra::SpatRaster.")
  }

  # Ensure unique effect directions to avoid redundant layers
  effect_direction <- unique(effect_direction)
  if (!all(effect_direction %in% c(1, -1))) {
    stop("'effect_direction' must only contain 1 (direct) and/or -1 (indirect).")
  }

  # 2. Subset prediction to target layers dynamically
  available_layers <- names(pred)
  missing_layers <- setdiff(layers, available_layers)

  if (length(missing_layers) > 0) {
    warning(sprintf("The following layers were not found in 'pred' and will be skipped: %s",
                    paste(missing_layers, collapse = ", ")))
  }

  valid_layers <- intersect(layers, available_layers)

  if (length(valid_layers) == 0) {
    stop("None of the target 'layers' were found in 'pred'.")
  }

  pred_sub <- pred[[valid_layers]]

  # 3. Extract the bias raster
  if (inherits(prepared_bias, "nicheR_bias_surface")) {
    if (!is.null(prepared_bias$pooled_bias)) {
      bias_ras <- prepared_bias$pooled_bias
    } else {
      stop("'prepared_bias' is missing the pooled layer. Did you use out_bias = 'standardized'?")
    }
  } else if (inherits(prepared_bias, "SpatRaster")) {
    bias_ras <- prepared_bias
  } else {
    stop("'prepared_bias' must be a 'nicheR_bias_surface' list or a 'terra::SpatRaster'.")
  }

  # 4. Apply Effects (Nested Iteration)
  result_list <- list()
  counter <- 1

  for (i in seq_along(valid_layers)) {
    layer_name <- valid_layers[i]
    current_pred <- pred_sub[[layer_name]]

    # Match the specific pooled_bias layer from prepare_bias, or fallback to index matching
    expected_bias_name <- paste0("pooled_bias_", layer_name)
    if (expected_bias_name %in% names(bias_ras)) {
      current_bias <- bias_ras[[expected_bias_name]]
    } else {
      bias_idx <- min(i, terra::nlyr(bias_ras))
      current_bias <- bias_ras[[bias_idx]]
    }

    # Crop and Mask
    b_cropped <- terra::crop(current_bias, current_pred, mask = TRUE)

    # Apply Effect Directions
    for (eff in effect_direction) {
      if (eff == 1) {
        applied <- current_pred * b_cropped
        names(applied) <- paste0(layer_name, "_biased_direct")
      } else {
        applied <- current_pred * (1 - b_cropped)
        names(applied) <- paste0(layer_name, "_biased_indirect")
      }

      result_list[[counter]] <- applied
      counter <- counter + 1
    }
  }

  # 5. Re-combine into a single SpatRaster
  applied_bias <- do.call(c, result_list)

  return(applied_bias)
}
