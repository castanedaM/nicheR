#' @export
sample_virtual_data <- function(n_occ,
                                object,
                                virtual_prediction = NULL,
                                prediction_layer = NULL,
                                sampling = "centroid",
                                method = "suitability",
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

  eps <- 1e-8
  tol <- 1e-8

  if(is.null(virtual_prediction)){
    virtual_prediction <- as.data.frame(unclass(predict(object = object,
                                  newdata = virtual_data(object = object))))
  }

  # Resolve prediction input --------------------------------------------------

  resolved <- resolve_prediction(virtual_prediction, prediction_layer)

  df <- resolved$df
  df$pred <- df[[resolved$pred_name]]
  df <- df[is.finite(df$pred), , drop = FALSE]

  if(nrow(df) == 0L){
    stop("No finite prediction values available for sampling.")
  }


  # Method-specific value checks ---------------------------------------------

  if(method == "suitability"){
    rng <- range(df$pred, na.rm = TRUE)
    if(rng[1] < (0 - tol) || rng[2] > (1 + tol)){
      stop("method = 'suitability' requires prediction values in [0, 1]. ",
           "Found range [", format(rng[1]), ", ", format(rng[2]), "]. ",
           "Either rescale your prediction to [0,1] or set method = 'mahalanobis'.")
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
