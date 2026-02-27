#' @export
sample_biased_data <- function(n_occ,
                        prediction,
                        prediction_layer = NULL, #column or raster-layer
                        sampling_mask = NULL,
                        seed = 1,
                        verbose = TRUE){

  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)
  gc()

  verbose_message("Starting: sample_biased_data()\n")

  # Basic Input checks --------------------------------------------------------

  if(!is.numeric(n_occ) || length(n_occ) != 1L || is.na(n_occ) || n_occ <= 0){
    stop("'n_occ' must be a single positive number.")
  }
  n_occ <- as.integer(n_occ)

  if(!is.logical(biased) || length(biased) != 1L){
    stop("'biased' must be TRUE or FALSE.")
  }

  if(!is.numeric(seed) || length(seed) != 1L || is.na(seed)){
    stop("'seed' must be a single number.")
  }


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

  if(n_occ > nrow(df)){
    stop("'n_occ' is larger than the number of available samples.")
  }


  w <- df$pred
  w[!is.finite(w)] <- NA_real_

  if(any(w < 0, na.rm = TRUE)){
    stop("Weighted sampling requires non-negative values. Found values < 0.")
  }

  w <- w + eps
  w <- rep(1, nrow(df))


  w[!is.finite(w)] <- 0
  if(sum(w) <= 0){
    stop("Sampling weights are all zero. Check your inputs.")
  }

  # Sample --------------------------------------------------------------------

  set.seed(seed)
  idx <- sample.int(nrow(df), size = n_occ, replace = FALSE, prob = w)
  df$pred <- NULL
  out <- df[idx, , drop = FALSE]

  verbose_message("Done: sampled ", nrow(out), " points from biased prediction layer\n")
  gc()

  out
}


