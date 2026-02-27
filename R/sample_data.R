#' Sample Occurrence Data from a Prediction Surface
#'
#' Samples \code{n_occ} occurrence locations from a prediction surface using
#' weighted random sampling. The prediction can be provided as a
#' \code{terra::SpatRaster} or a data.frame containing prediction values.
#'
#' The user must specify whether the prediction values represent
#' \code{"suitability"} (expected range [0, 1]) or \code{"mahalanobis"}
#' (distance or distance-squared values).
#'
#' @param n_occ Integer. Number of occurrences to sample.
#' @param prediction Prediction object to sample from. Typically a
#' \code{terra::SpatRaster} or a data.frame. This is resolved by
#' \code{resolve_prediction()}.
#' @param prediction_layer Optional. Layer name or index used by
#' \code{resolve_prediction()} when \code{prediction} contains multiple layers.
#' @param sampling Character. Sampling bias direction. One of:
#' \itemize{
#'   \item \code{"centroid"}: sample more from high suitability or low distance
#'   \item \code{"edge"}: sample more from low suitability or high distance
#'   \item \code{"random"}: sample uniformly
#' }
#' @param method Character. Interpretation of \code{prediction} values. One of:
#' \itemize{
#'   \item \code{"suitability"}: values must be in [0, 1]
#'   \item \code{"mahalanobis"}: values treated as a distance-like measure (>= 0)
#' }
#' @param sampling_mask Optional. Restriction mask applied only when
#' \code{prediction} resolves to a \code{terra::SpatRaster}. Must be a
#' \code{terra::SpatRaster} or \code{terra::SpatVector}.
#' @param seed Numeric. Random seed.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return A data.frame of sampled points.
#' If raster input is used, includes \code{x} and \code{y} columns plus any
#' extracted columns (excluding the internal prediction column).
#'
#' @details
#' Weighting rules:
#' \itemize{
#'   \item \code{sampling = "random"}: uniform weights.
#'   \item \code{method = "suitability"}:
#'     \code{"centroid"} uses \eqn{w = pred + eps}, \code{"edge"} uses
#'     \eqn{w = (1 - pred) + eps}.
#'   \item \code{method = "mahalanobis"}:
#'     \code{"centroid"} uses \eqn{w = 1 / (pred + eps)}, \code{"edge"} uses
#'     \eqn{w = pred + eps}.
#' }
#'
#' If \code{method = "suitability"} and any finite prediction values fall outside
#' [0, 1] (with a small tolerance), the function stops with an error advising the
#' user to switch \code{method} or rescale their prediction.
#'
#' @export
sample_data <- function(n_occ,
                        prediction,
                        prediction_layer = NULL,
                        sampling = "centroid",
                        method = "suitability",
                        sampling_mask = NULL,
                        seed = 1,
                        verbose = TRUE){

  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)
  gc()

  verbose_message("Starting: sample_data()\n")

  sampling <- match.arg(sampling,
                        choices = c("centroid", "edge", "random"),
                        several.ok = FALSE)

  method <- match.arg(method,
                      choices = c("suitability", "mahalanobis"),
                      several.ok = FALSE)

  if(!is.numeric(n_occ) || length(n_occ) != 1L || is.na(n_occ) || n_occ <= 0){
    stop("'n_occ' must be a single positive number.")
  }
  n_occ <- as.integer(n_occ)

  if(!is.numeric(seed) || length(seed) != 1L || is.na(seed)){
    stop("'seed' must be a single number.")
  }

  if(!is.logical(verbose) || length(verbose) != 1L){
    stop("'verbose' must be TRUE or FALSE.")
  }

  eps <- .Machine$double.eps
  tol <- 1e-8

  # Resolve prediction input --------------------------------------------------

  resolved <- resolve_prediction(prediction, prediction_layer)

  # Enforce sampling_mask rules ----------------------------------------------

  if(!is.null(sampling_mask) && resolved$type != "SpatRaster"){
    stop("'sampling_mask' is only supported when 'prediction' resolves to a SpatRaster.")
  }

  # Apply sampling mask (raster only) ----------------------------------------

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

    df$pred <- df[[pred_cols[1]]]

  } else {

    df <- resolved$df
    df$pred <- df[[resolved$pred_name]]
    df <- df[is.finite(df$pred), , drop = FALSE]

    if(nrow(df) == 0L){
      stop("No finite prediction values available for sampling.")
    }
  }

  df <- df[is.finite(df$pred), , drop = FALSE]
  if(nrow(df) == 0L){
    stop("No finite prediction values available for sampling.")
  }

  # Method-specific value checks ---------------------------------------------

  if(method == "suitability"){
    rng <- range(df$pred, na.rm = TRUE)
    if(rng[1] < (0 - tol) || rng[2] > (1 + tol)){
      stop(
        "method = 'suitability' requires prediction values in [0, 1]. ",
        "Found range [", format(rng[1]), ", ", format(rng[2]), "]. ",
        "Either rescale your prediction to [0,1] or set method = 'mahalanobis'."
      )
    }
    # clamp tiny numeric drift
    df$pred[df$pred < 0] <- 0
    df$pred[df$pred > 1] <- 1
  } else {
    if(any(df$pred < 0, na.rm = TRUE)){
      stop("method = 'mahalanobis' requires non-negative prediction values. Found values < 0.")
    }
  }

  # Ensure sample size possible ----------------------------------------------

  if(n_occ > nrow(df)){
    stop("'n_occ' is larger than the number of available samples.")
  }

  # Weights -------------------------------------------------------------------

  if(sampling == "random"){

    w <- rep(1, nrow(df))

  } else if(method == "suitability"){

    if(sampling == "centroid"){
      w <- df$pred + eps
    } else { # edge
      w <- (1 - df$pred) + eps
    }

  } else { # method == "mahalanobis"

    if(sampling == "centroid"){
      w <- 1 / (df$pred + eps)
    } else { # edge
      w <- df$pred + eps
    }
  }

  w[!is.finite(w)] <- 0

  if(sum(w) <= 0){
    stop("Sampling weights are all zero. Check your inputs.")
  }

  # Sample --------------------------------------------------------------------

  set.seed(seed)
  idx <- sample.int(nrow(df), size = n_occ, replace = FALSE, prob = w)

  # Return sampled points (drop internal pred column)
  df$pred <- NULL
  out <- df[idx, , drop = FALSE]

  verbose_message("Done: sampled ", nrow(out), " points.\n")
  gc()

  out
}
