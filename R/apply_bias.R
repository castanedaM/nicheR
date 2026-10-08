#' Apply sampling bias to prediction surfaces
#'
#' @description
#' Combines a prepared composite sampling bias surface with a prediction
#' layer. The prediction is first turned into a sampling weight, following
#' the same \code{sampling} and \code{method} rules as
#' \code{\link{sample_data}}, and that weight is multiplied by the bias
#' surface. The bias surface is aligned to the prediction grid when needed
#' and the result is cropped and masked to the prediction domain. The output
#' is a relative sampling weight and is no longer interpretable as a
#' probability.
#'
#' @usage apply_bias(prepared_bias, prediction, prediction_layer = NULL,
#'                   sampling = "centroid", method = "suitability",
#'                   verbose = TRUE)
#'
#' @param prepared_bias A single-layer \code{SpatRaster} composite bias
#'   surface, or the list output from \code{\link{prepare_bias}} containing
#'   a \code{composite_surface} element.
#' @param prediction A \code{SpatRaster} containing one or more prediction
#'   layers, for example suitability or Mahalanobis distance.
#' @param prediction_layer Character. Name of the layer to use from
#'   \code{prediction}. Required when \code{prediction} contains multiple
#'   layers. If \code{NULL} (default) and \code{prediction} has a single
#'   layer, that layer is used.
#' @param sampling Character. Sampling strategy. One of \code{"centroid"}
#'   (default), \code{"edge"}, or \code{"uniform"}. Controls where within the
#'   niche the weight is highest.
#' @param method Character. Weighting method. One of \code{"suitability"}
#'   (default) or \code{"mahalanobis"}. Must match the type of values in
#'   \code{prediction_layer}: suitability values must be in \code{[0, 1]},
#'   Mahalanobis distances must be non-negative.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @details
#' The \code{sampling} and \code{method} arguments define the weight that is
#' multiplied by the bias surface, the same way they define the sampling
#' weights in \code{\link{sample_data}}:
#' \itemize{
#'   \item \code{sampling = "centroid"}, \code{method = "suitability"}:
#'   the suitability itself, higher near the niche center.
#'   \item \code{sampling = "edge"}, \code{method = "suitability"}:
#'   \eqn{1 - \text{suitability}}, higher near the niche boundary.
#'   \item \code{sampling = "centroid"}, \code{method = "mahalanobis"}:
#'   \eqn{1 / D^2}, with \eqn{D^2} the Mahalanobis distance, higher near the
#'   centroid.
#'   \item \code{sampling = "edge"}, \code{method = "mahalanobis"}:
#'   \eqn{D^2}, higher near the boundary.
#'   \item \code{sampling = "uniform"}: a weight of 1 in every cell that has a
#'   prediction, so the output is the bias surface limited to those cells.
#' }
#'
#' Cells where the prediction is zero or \code{NA} are set to \code{NA} in
#' the output, so they cannot be sampled. This is what \code{strict = TRUE}
#' does in \code{\link{sample_data}}, and here it is always applied. In a
#' truncated layer these are the cells outside the niche, where
#' \code{"edge"} and \code{"uniform"} would otherwise place their highest
#' weights.
#'
#' The function stops when \code{method} disagrees with the layer name, for
#' example \code{method = "mahalanobis"} on a layer called
#' \code{"suitability_trunc"}. The two kinds of layer run in opposite
#' directions, so the wrong pairing weights the wrong part of the niche. It
#' warns when \code{method = "mahalanobis"} is used on an unnamed layer whose
#' values all fall in \code{[0, 1]}.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the composite bias surface from \code{prepared_bias}.
#'   \item Verifies the bias values are within \code{[0, 1]} and that the
#'   prediction values agree with \code{method}.
#'   \item Aligns the bias surface to the prediction grid if geometries
#'   differ, using \code{terra::resample()} with nearest-neighbor
#'   interpolation.
#'   \item Removes the cells where the prediction is zero.
#'   \item Turns the prediction into a weight according to \code{sampling}
#'   and \code{method}.
#'   \item Multiplies the weight by the bias surface.
#'   \item Crops and masks the output to the prediction domain.
#' }
#'
#' @return
#' A named list of class \code{"nicheR_biased_surface"} containing:
#' \itemize{
#'   \item One single-layer \code{SpatRaster} named
#'   \code{"<layer>_<sampling>_biased"} (e.g.,
#'   \code{"suitability_centroid_biased"}). The raster layer carries the
#'   same name.
#'   \item \code{combination_formula}: a character string describing the
#'   operation applied: \code{"suitability * bias"},
#'   \code{"(1 - suitability) * bias"}, \code{"(1 / Mahalanobis) * bias"},
#'   \code{"Mahalanobis * bias"}, or \code{"bias"} for
#'   \code{sampling = "uniform"}.
#' }
#'
#' @seealso \code{\link{prepare_bias}} to build the composite bias surface,
#'   \code{\link{sample_biased_data}} to sample data points from the output,
#'   \code{\link{sample_data}} for the same weights without a bias surface.
#'
#' @importFrom terra compareGeom resample crop global nlyr ifel clamp
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
#' # 2. Apply bias to the suitability layer, weighting toward the centroid
#' biased_pred <- apply_bias(prepared_bias = bias,
#'                           prediction = pred_rast,
#'                           prediction_layer = "suitability")
#'
#' terra::plot(biased_pred$suitability_centroid_biased)
#'
#' # 3. Weight toward the niche edge instead, using Mahalanobis distance
#' biased_edge <- apply_bias(prepared_bias = bias,
#'                           prediction = pred_rast,
#'                           prediction_layer = "Mahalanobis_trunc",
#'                           sampling = "edge",
#'                           method = "mahalanobis")
#'
#' biased_edge$combination_formula
#'
#' @export
apply_bias <- function(prepared_bias,
                       prediction,
                       prediction_layer = NULL,
                       sampling = "centroid",
                       method = "suitability",
                       verbose = TRUE){

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

  sampling <- match.arg(tolower(sampling),
                        choices = c("centroid", "edge", "uniform"),
                        several.ok = FALSE)

  method <- match.arg(tolower(method),
                      choices = c("suitability", "mahalanobis"),
                      several.ok = FALSE)

  eps <- 1e-8
  tol <- 1e-8

  # Always a single layer after this
  s <- resolve_prediction(prediction, prediction_layer)$rast

  # Safe prediction name
  pred_name <- names(s)
  if(is.null(pred_name) || length(pred_name) == 0L || !nzchar(pred_name[1]) || pred_name[1] == "lyr.1"){
    pred_name <- "prediction"
  }else{
    pred_name <- pred_name[1]
  }

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

  # Safe bias name
  bias_name <- names(bias_rast)
  if(is.null(bias_name) || length(bias_name) == 0L || !nzchar(bias_name[1]) || bias_name[1] == "lyr.1"){
    bias_name <- "bias"
  }else{
    bias_name <- bias_name[1]
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

  # 3. Check the prediction against the method -------------------------------

  pred_rng <- terra::global(s,
                            fun = c("min", "max"),
                            na.rm = TRUE)

  pred_min <- as.numeric(pred_rng[1, "min"])
  pred_max <- as.numeric(pred_rng[1, "max"])

  if(!is.finite(pred_min) || !is.finite(pred_max)){
    stop("Layer '", pred_name, "' has no finite values to apply bias to.")
  }

  # The layer name and the method have to describe the same kind of layer
  check_method_layer(method, pred_name, max_value = pred_max)

  if(method == "suitability"){

    if(pred_min < (0 - tol) || pred_max > (1 + tol)){
      stop(
        "method = 'suitability' requires prediction values in [0, 1]. ",
        "Found range [", format(pred_min), ", ", format(pred_max), "]. ",
        "Either rescale your prediction to [0,1] or set method = 'mahalanobis'."
      )
    }

    # clamp tiny numeric drift
    s <- terra::clamp(s, lower = 0, upper = 1, values = TRUE)

  }else{

    if(pred_min < 0){
      stop("method = 'mahalanobis' requires non-negative prediction values. Found values < 0.")
    }
  }

  # 4. Align bias to prediction grid -----------------------------------------

  same_grid <- terra::compareGeom(bias_rast,
                                  s,
                                  stopOnError = FALSE)

  if(!isTRUE(same_grid)){
    verbose_message(verbose, "Step: resampling prepared bias to match prediction grid...\n")
    bias_rast <- terra::resample(bias_rast,
                                 s,
                                 method = "near")
  }

  # 5. Remove cells where the prediction is zero -----------------------------

  # The same as strict = TRUE in sample_data(). In a truncated layer these
  # are the cells outside the niche, and "edge" and "uniform" would
  # otherwise give them the highest weight.
  if(pred_min <= 0){
    verbose_message(verbose, "Step: removing cells where the prediction is zero...\n")
    s <- terra::ifel(s == 0, NA, s)
  }

  # 6. Turn the prediction into a sampling weight ----------------------------

  verbose_message(verbose, "Step: applying bias to \"", pred_name,
                  "\" with sampling = '", sampling,
                  "' and method = '", method, "'...\n")

  # Same weights as sample_data()
  if(sampling == "uniform"){

    # Every cell with a prediction gets the same weight
    weight <- (s * 0) + 1
    formula_entry <- bias_name

  }else if(method == "suitability"){

    if(sampling == "centroid"){
      weight <- s
      formula_entry <- paste0(pred_name, " * ", bias_name)
    }else{ # edge
      weight <- 1 - s
      formula_entry <- paste0("(1 - ", pred_name, ") * ", bias_name)
    }

  }else{ # method == "mahalanobis"

    if(sampling == "centroid"){
      weight <- 1 / (s + eps)
      formula_entry <- paste0("(1 / ", pred_name, ") * ", bias_name)
    }else{ # edge
      weight <- s
      formula_entry <- paste0(pred_name, " * ", bias_name)
    }
  }

  # 7. Multiply by the bias surface ------------------------------------------

  out_r <- weight * bias_rast
  out_r <- terra::crop(out_r, s, mask = TRUE)

  # The list element and the raster layer carry the same name
  out_name <- paste0(pred_name, "_", sampling, "_biased")
  names(out_r) <- out_name

  out_list <- list(out_r)
  names(out_list) <- out_name

  # 8. Attach message / metadata ---------------------------------------------
  out_list$combination_formula <- formula_entry

  class(out_list) <- "nicheR_biased_surface"

  verbose_message(verbose, "Done: apply_bias(). Note: values are no longer probabilities\n")

  out_list
}
