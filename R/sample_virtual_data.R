#' @export
sample_virtual_data <- function(n_occ,
                                object,
                                virtual_prediction = NULL,
                                prediction_layer = NULL,
                                sampling = "centroid",
                                method = "suitability",
                                seed = 1,
                                verbose = TRUE,
                                strict = NULL){

  gc()

  verbose_message(verbose, "Starting: sample_virtual_data()\n")

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

  if(!is.null(strict) && (!is.logical(strict) || length(strict) != 1L || is.na(strict))){
    stop("'strict' must be TRUE, FALSE, or NULL.")
  }

  eps <- 1e-8
  tol <- 1e-8

  # thresholds for auto-detection of truncation
  zero_prop_threshold <- 0.25
  na_prop_threshold <- 0.25

  # Resolve prediction input --------------------------------------------------

  resolved <- resolve_prediction(virtual_prediction, prediction_layer)

  df <- resolved$df
  pred_name <- resolved$pred_name
  df$pred <- df[[pred_name]]

  if(nrow(df) == 0L){
    stop("No prediction values available for sampling.")
  }

  # Auto-detect strict --------------------------------------------------------

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
    } else {
      w <- (1 - df$pred) + eps
    }

  } else {

    if(sampling == "centroid"){
      w <- 1 / (df$pred + eps)
    } else {
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

  df$pred <- NULL
  out <- df[idx, , drop = FALSE]

  verbose_message(verbose, "Done: sampled ", nrow(out), " points.\n")
  gc()

  out
}
