#' Print message if verbose is TRUE
#'
#' Internal helper used to conditionally print progress messages.
#'
#' @param verbose Logical; whether to print the message.
#' @param ... Objects passed to \code{message()}.
#'
#' @return Invisibly \code{NULL}.
#' @noRd
#' @keywords internal
verbose_message <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

#' Compute ellipsoid hypervolume
#'
#' Computes the geometric volume (area in 2D, volume in 3D, hypervolume in higher
#' dimensions) of a \eqn{p}-dimensional ellipsoid defined by its semi-axis lengths.
#'
#' For semi-axes \eqn{a_1, \dots, a_p}, the volume is:
#' \deqn{
#' V_p = \frac{\pi^{p/2}}{\Gamma(p/2 + 1)} \prod_{i=1}^{p} a_i
#' }
#' where \eqn{\pi^{p/2} / \Gamma(p/2 + 1)} is the volume of the unit
#' \eqn{p}-dimensional ball.
#'
#' In probabilistic niche models, semi-axis lengths are typically derived from
#' covariance eigenvalues and a chi-square cutoff.
#'
#' @param n_dimensions Integer. Number of dimensions \eqn{p}.
#' @param semi_axes_lengths Numeric vector of length \code{n_dimensions}
#'   containing the ellipsoid semi-axis lengths.
#'
#' @return Numeric. Geometric volume (or hypervolume) of the ellipsoid.
#'
#' @seealso \code{\link{build_ellipsoid}},
#'   \code{\link{ellipsoid_calculator}}
#'
#' @export
ellipsoid_volume <- function(n_dimensions, semi_axes_lengths) {

  # Checks
  if(missing(n_dimensions))
    stop("Argument 'n_dimensions' must be defined.")

  if(missing(semi_axes_lengths))
    stop("Argument 'semi_axes_lengths' must be defined.")

  if(!is.numeric(n_dimensions) || length(n_dimensions) != 1)
    stop("'n_dimensions' must be a single numeric value.")

  if(!is.numeric(semi_axes_lengths))
    stop("'semi_axes_lengths' must be numeric.")

  if(length(semi_axes_lengths) != n_dimensions)
    stop("'semi_axes_lengths' must have length equal to 'n_dimensions'.")

  if(any(semi_axes_lengths <= 0))
    stop("All semi-axis lengths must be positive.")

  # dimensions-dimensional n volume
  unit_ball_volume <- pi^(n_dimensions / 2) /
    gamma(n_dimensions / 2 + 1)

  # ellipsoid volume
  volume <- unit_ball_volume * prod(semi_axes_lengths)

  return(volume)
}


#' Update covariances in a nicheR ellipsoid and recompute metrics
#'
#' @param x A `nicheR_ellipsoid` object from `build_ellipsoid()`.
#' @param covariance Either a single numeric (applied to all off-diagonals) or a
#'   named numeric vector with names like `"var1-var2"`.
#' @param tol Small positive number used by covariance limit helpers.
#' @param verbose Logical; print brief progress messages.
#'
#' @return A new `nicheR_ellipsoid` object with updated covariance matrix and
#'   recomputed ellipsoid metrics. Also includes remaining safe limits for any
#'   still-zero covariances.
#' @export
update_ellipsoid_covariance <- function(object,
                                        covariance,
                                        tol = 1e-6,
                                        verbose = TRUE){

  stopifnot(inherits(object, "nicheR_ellipsoid"))

  verbose_message(verbose, "Starting: updating covariance values...\n")

  up <- update_covariance(object$cov_matrix, covariance = covariance, tol = tol)

  object <- ellipsoid_calculator(cov_matrix = up$updated_matrix,
                                 centroid = object$centroid,
                                 cl = object$cl,
                                 verbose = verbose)

  # new fields
  object$cov_limits_remaining <- up$remaining_limits

  object
}


#' Resolves user input when to many columns or layers
#'
#'
#' @export
resolve_prediction <- function(prediction, prediction_layer){

  # ---- Case 1: prediction is a data.frame -----------------------------------
  if(is.data.frame(prediction)){

    if(is.null(prediction_layer) || !nzchar(prediction_layer)){
      stop("If 'prediction' is a data.frame, 'prediction_layer' must be provided (column name).")
    }

    if(!all(c("x", "y") %in% names(prediction))){
      warning("'prediction' is a data.frame, and it is missing 'x' and 'y', results wont show geographical connections.")
    }

    if(!(prediction_layer %in% names(prediction))){
      stop("Column '", prediction_layer, "' not found in data.frame prediction.")
    }

    return(list(type = "data.frame",
                df = prediction,
                pred_name = prediction_layer))
  }

  # ---- Case 2: prediction is a SpatRaster -----------------------------------
  if(inherits(prediction, "SpatRaster")){

    if(terra::nlyr(prediction) == 1L){

      if(!is.null(prediction_layer) && nzchar(prediction_layer)){
        if(is.null(names(prediction)) || !(prediction_layer %in% names(prediction))){
          stop("Layer '", prediction_layer, "' not found in SpatRaster prediction.")
        }
        prediction <- prediction[[prediction_layer]]
      }

      return(list(type = "SpatRaster", rast = prediction, pred_name = names(prediction)[1]))
    }

    # multi-layer raster: require layer name
    if(is.null(prediction_layer) || !nzchar(prediction_layer)){
      stop("If 'prediction' is a multi-layer SpatRaster, 'prediction_layer' must be provided (layer name).")
    }

    if(is.null(names(prediction)) || !(prediction_layer %in% names(prediction))){
      stop("Layer '", prediction_layer, "' not found in SpatRaster prediction. Available: ",
           paste(names(prediction), collapse = ", "))
    }

    r <- prediction[[prediction_layer]]
    return(list(type = "SpatRaster", rast = r, pred_name = prediction_layer))
  }

  # ---- Case 3: prediction is a list -----------------------------------------
  if(is.list(prediction)){

    if(length(prediction) == 0L){
      stop("'prediction' is an empty list.")
    }

    # If list has one element: repeat rules above
    if(length(prediction) == 1L){
      return(resolve_prediction(prediction[[1]], prediction_layer))
    }

    # list has multiple elements: require prediction_layer
    if(is.null(prediction_layer) || !nzchar(prediction_layer)){
      stop("If 'prediction' is a list with multiple elements, 'prediction_layer' must be provided.")
    }

    # (A) First: check if prediction_layer matches list names
    if(!is.null(names(prediction)) && prediction_layer %in% names(prediction)){
      return(resolve_prediction(prediction[[prediction_layer]], prediction_layer))
    }

    # (B) Otherwise: go element-by-element to find prediction_layer
    for(i in seq_along(prediction)){
      obj <- prediction[[i]]
      if(is.data.frame(obj)){
        if(prediction_layer %in% names(obj)){
          return(list(type = "data.frame",
                      df = obj,
                      pred_name = prediction_layer))
        }
      }

      if(inherits(obj, "SpatRaster")){
        if(!is.null(names(obj)) && prediction_layer %in% names(obj)){
          return(list(type = "SpatRaster",
                      rast = obj[[prediction_layer]],
                      pred_name = prediction_layer))
        }
      }
    }

    stop("Could not find 'prediction_layer' in list names or within any list element.")
  }

  stop("'prediction' must be a SpatRaster, a data.frame, or a list.")
}




