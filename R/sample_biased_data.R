#' Sample data points from a bias-weighted prediction surface
#'
#' @description
#' Samples \code{n} virtual data points using the bias-weighted
#' prediction values directly as sampling probabilities. Unlike
#' \code{sample_data()}, there is no sampling strategy argument. The
#' prediction layer values themselves define where points are drawn from,
#' making this function suited for simulating realistically biased
#' data points.
#'
#' @usage sample_biased_data(n, prediction, prediction_layer = NULL,
#'                           sampling_mask = NULL, seed = 1, verbose = TRUE,
#'                           strict = NULL)
#'
#' @param n Integer. Number of data points to sample.
#' @param prediction A \code{SpatRaster} or data frame containing the
#'   bias-weighted prediction surface to sample from.
#' @param prediction_layer Character. Name of the layer or column to use as
#'   sampling weights. Required when \code{prediction} contains multiple
#'   layers or columns.
#' @param sampling_mask A \code{SpatRaster} or \code{SpatVector} used to
#'   restrict sampling to a geographic area. Only supported when
#'   \code{prediction} is a \code{SpatRaster}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{1}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#' @param strict Logical or \code{NULL}. If \code{TRUE}, removes \code{NA} and
#'   zero-valued cells before sampling. If \code{NULL} (default),
#'   auto-detected from the layer name and the proportion of zeros and
#'   \code{NA}s in the prediction values.
#'
#' @details
#' Prediction values are used directly as sampling weights, so they must be
#' non-negative. Higher values correspond to higher sampling probability,
#' reflecting areas of greater bias (e.g., higher detectability or observer
#' effort). This is in contrast to \code{sample_data()}, which transforms
#' prediction values according to a \code{sampling} and \code{method}
#' argument.
#'
#' Auto-detection of \code{strict} follows the same logic as
#' \code{sample_data()}: it is set to \code{TRUE} if the layer name contains
#' \code{"trunc"} or if the proportion of zeros or \code{NA}s exceeds 25\%.
#'
#' @return
#' A data frame of sampled data points with the same columns as the
#' input \code{prediction} (minus the internal \code{pred} column). If
#' \code{prediction} is a \code{SpatRaster}, the output includes \code{x}
#' and \code{y} coordinate columns.
#'
#' @seealso \code{\link{sample_data}} for unbiased sampling with explicit
#'   strategy and method control, \code{\link{apply_bias}} for generating
#'   the bias-weighted prediction surface used as input here.
#'
#' @examples
#' biased_pred <- terra::rast(system.file("extdata/applied_bias_rast.tif",
#'                                      package = "nicheR"))
#'
#' # Sample points from bias surface (not probability surface)
#' pts_biased <- sample_biased_data(n = 100,
#'                                  prediction = biased_pred,
#'                                  prediction_layer = "suitability_centroid_biased")
#'
#' head(pts_biased)
#'
#' @export
sample_biased_data <- function(n,
                               prediction,
                               prediction_layer = NULL,
                               sampling_mask = NULL,
                               seed = 1,
                               verbose = TRUE,
                               strict = NULL){

  verbose_message(verbose, "Starting: sample_biased_data()\n")

  if(!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0){
    stop("'n' must be a single positive number.")
  }
  n <- as.integer(n)

  if(!is.numeric(seed) || length(seed) != 1L || is.na(seed)){
    stop("'seed' must be a single number.")
  }

  if(!is.logical(verbose) || length(verbose) != 1L){
    stop("'verbose' must be TRUE or FALSE.")
  }

  if(!is.null(strict) && (!is.logical(strict) || length(strict) != 1L || is.na(strict))){
    stop("'strict' must be TRUE, FALSE, or NULL.")
  }

  zero_prop_threshold <- 0.25
  na_prop_threshold <- 0.25

  resolved <- resolve_prediction(prediction, prediction_layer)

  if(!is.null(sampling_mask) && resolved$type != "SpatRaster"){
    stop("'sampling_mask' is only supported when 'prediction' resolves to a SpatRaster.")
  }

  if(resolved$type == "SpatRaster"){
    r <- resolved$rast

    if(!is.null(sampling_mask)){
      if(inherits(sampling_mask, "SpatRaster")){
        r <- terra::mask(r, sampling_mask)
      } else if(inherits(sampling_mask, "SpatVector")){
        r <- terra::mask(r, sampling_mask)
      } else {
        stop("'sampling_mask' must be a SpatRaster or SpatVector.")
      }
    }

    df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)

    if(nrow(df) == 0L){
      stop("No non-NA cells available for sampling.")
    }

    pred_cols <- setdiff(names(df), c("x", "y"))
    if(length(pred_cols) != 1L){
      stop("Unexpected: raster extraction did not produce exactly one prediction column.")
    }

    pred_name <- pred_cols[1]
    df$pred <- df[[pred_name]]

  } else {

    df <- resolved$df
    pred_name <- resolved$pred_name
    df$pred <- df[[pred_name]]

    if(nrow(df) == 0L){
      stop("No prediction values available for sampling.")
    }
  }

  if(is.null(strict)){

    layer_name_flag <- FALSE
    if(!is.null(prediction_layer)){
      layer_name_flag <- grepl("trunc", prediction_layer, ignore.case = TRUE)
    }
    if(!is.null(pred_name) && nzchar(pred_name)){
      layer_name_flag <- layer_name_flag || grepl("trunc", pred_name, ignore.case = TRUE)
    }

    na_prop <- mean(is.na(df$pred))
    zero_prop <- mean(!is.na(df$pred) & df$pred == 0)

    strict <- isTRUE(layer_name_flag) ||
      isTRUE(na_prop >= na_prop_threshold) ||
      isTRUE(zero_prop >= zero_prop_threshold)

    if(isTRUE(strict)){
      verbose_message(verbose,
                      paste0(
                        "Step: auto-detected a likely truncated prediction surface. ",
                        "Setting 'strict = TRUE' and removing NA and zero values. ",
                        "You can override this behavior with the 'strict' argument...\n"
                      ))
    } else {
      strict <- FALSE
    }
  }

  if(isTRUE(strict)){
    df <- df[!is.na(df$pred) & is.finite(df$pred) & df$pred != 0, , drop = FALSE]
  } else {
    df <- df[is.finite(df$pred), , drop = FALSE]
  }

  if(nrow(df) == 0L){
    stop("No valid prediction values available for sampling after filtering.")
  }

  if(n > nrow(df)){
    stop("'n' is larger than the number of available samples.")
  }

  w <- df$pred

  if(any(w < 0, na.rm = TRUE)){
    stop("Weighted sampling requires non-negative values. Found values < 0.")
  }

  w[!is.finite(w)] <- 0

  if(sum(w) <= 0){
    stop("Sampling weights are all zero. Check your inputs.")
  }

  set.seed(seed)
  idx <- sample.int(nrow(df), size = n, replace = FALSE, prob = w)

  df$pred <- NULL
  out <- df[idx, , drop = FALSE]

  verbose_message(verbose, "Done: sampled ", nrow(out), " points from biased prediction layer\n")

  out
}
