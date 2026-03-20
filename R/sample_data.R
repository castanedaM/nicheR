#' Sample occurrence data from a prediction surface
#'
#' @description
#' Samples \code{n_occ} virtual occurrence points from a suitability or
#' Mahalanobis distance prediction surface. Supports centroid, edge, and
#' random sampling strategies, and accepts both raster (\code{SpatRaster})
#' and data frame inputs.
#'
#' @param n_occ Integer. Number of occurrence points to sample.
#' @param prediction A \code{SpatRaster} or data frame containing the
#'   prediction surface to sample from.
#' @param prediction_layer Character. Name of the layer or column to use as
#'   the prediction values. Required when \code{prediction} contains multiple
#'   layers or columns.
#' @param sampling Character. Sampling strategy. One of \code{"centroid"}
#'   (default), \code{"edge"}, or \code{"random"}. Controls where within the
#'   niche points are preferentially drawn from.
#' @param method Character. Weighting method. One of \code{"suitability"}
#'   (default) or \code{"mahalanobis"}. Must match the type of values in
#'   \code{prediction_layer}: suitability values must be in \code{[0, 1]},
#'   Mahalanobis values must be non-negative.
#' @param sampling_mask A \code{SpatRaster} or \code{SpatVector} used to
#'   restrict sampling to a geographic area. Only supported when
#'   \code{prediction} is a \code{SpatRaster}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{1}.
#' @param strict Logical or \code{NULL}. If \code{TRUE}, removes \code{NA} and
#'   zero-valued cells before sampling (recommended with truncated prediction
#'   layers). If \code{NULL} (default), auto-detected from the layer name and
#'   the proportion of zeros and \code{NA}s in the prediction values.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @details
#' The \code{sampling} and \code{method} arguments interact to define the
#' probability weights used when drawing points:
#' \itemize{
#'   \item \code{sampling = "centroid"}, \code{method = "suitability"}:
#'   weights proportional to suitability — higher near the niche center.
#'   \item \code{sampling = "edge"}, \code{method = "suitability"}:
#'   weights proportional to \eqn{1 - \text{suitability}} — higher near the
#'   niche boundary.
#'   \item \code{sampling = "centroid"}, \code{method = "mahalanobis"}:
#'   weights inversely proportional to Mahalanobis distance — higher near
#'   the centroid.
#'   \item \code{sampling = "edge"}, \code{method = "mahalanobis"}:
#'   weights proportional to Mahalanobis distance — higher near the boundary.
#'   \item \code{sampling = "random"}: equal weights regardless of method.
#' }
#'
#' When \code{strict = NULL}, the function auto-detects truncation by checking
#' whether the layer name contains \code{"trunc"} or whether the proportion of
#' zeros or \code{NA}s exceeds 25\%.
#'
#' @return
#' A data frame of sampled occurrence points with the same columns as the
#' input \code{prediction} (minus the internal \code{pred} column). If
#' \code{prediction} is a \code{SpatRaster}, the output includes \code{x}
#' and \code{y} coordinate columns.
#'
#' @export
sample_data <- function(n_occ,
                        prediction,
                        prediction_layer = NULL,
                        sampling = "centroid",
                        method = "suitability",
                        sampling_mask = NULL,
                        seed = 1,
                        strict = NULL,
                        verbose = TRUE){

  gc()
  verbose_message(verbose, "Starting: sample_data()\n")

  sampling <- match.arg(tolower(sampling),
                        choices = c("centroid", "edge", "random"),
                        several.ok = FALSE)

  method <- match.arg(tolower(method),
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

  if(!is.null(strict) && (!is.logical(strict) || length(strict) != 1L || is.na(strict))){
    stop("'strict' must be TRUE, FALSE, or NULL.")
  }

  eps <- 1e-8
  tol <- 1e-8

  # thresholds for auto-detection of truncation
  zero_prop_threshold <- 0.25
  na_prop_threshold <- 0.25

  # Resolve prediction input --------------------------------------------------

  resolved <- resolve_prediction(prediction, prediction_layer)

  # Enforce sampling_mask rules ----------------------------------------------

  if(!is.null(sampling_mask) && resolved$type != "SpatRaster"){
    stop("'sampling_mask' is only supported when 'prediction' resolves to a SpatRaster.")
  }

  # Apply sampling mask and extract data -------------------------------------

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

    df <- terra::as.data.frame(prediction, xy = TRUE)

    if(nrow(df) == 0L){
      stop("No cells available after extracting the prediction surface.")
    }

    pred_name <- names(r)
    df$pred <- df[[pred_name]]

  } else {

    df <- resolved$df
    pred_name <- resolved$pred_name
    df$pred <- df[[pred_name]]

    if(nrow(df) == 0L){
      stop("No prediction values available for sampling.")
    }
  }

  # Auto-detect strict if needed ---------------------------------------------

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

  # Filter prediction values --------------------------------------------------

  if(isTRUE(strict)){
    df <- df[!is.na(df$pred) & is.finite(df$pred) & df$pred != 0, , drop = FALSE]
  } else {
    df <- df[is.finite(df$pred), , drop = FALSE]
  }

  if(nrow(df) == 0L){
    stop("No valid prediction values available for sampling after filtering.")
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

  # Return sampled points
  df$pred <- NULL
  out <- df[idx, , drop = FALSE]

  verbose_message(verbose, "Done: sampled ", nrow(out), " points.\n")
  gc()

  out
}
