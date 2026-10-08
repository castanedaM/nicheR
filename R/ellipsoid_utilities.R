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

#' Calculate n-dimensional ellipsoid metrics
#'
#' Computes geometric and probabilistic metrics for an n-dimensional ellipsoid
#' defined by a centroid and covariance matrix, including semi-axis lengths,
#' axis vertices, and hypervolume for a chi-square confidence contour.
#'
#' @usage ellipsoid_calculator(cov_matrix, centroid,
#'                            cl, verbose = TRUE)
#'
#' @details
#' The ellipsoid boundary is defined by the constant Mahalanobis distance contour:
#' \deqn{(x - \mu)^\top \Sigma^{-1} (x - \mu) = c^2,}
#' where \eqn{\mu} is the centroid, \eqn{\Sigma} is the covariance matrix, and
#' \eqn{c^2 = \chi^2_{n}(\mathrm{cl})} is the chi-square cutoff with \eqn{n}
#' degrees of freedom.
#'
#' The covariance matrix must be symmetric positive definite (SPD). The inverse
#' covariance is computed via the Cholesky factorization. Semi-axis lengths are
#' computed from covariance eigenvalues \eqn{\lambda_i} as:
#' \deqn{a_i = \sqrt{\lambda_i c^2}.}
#'
#' Axis vertices are computed along each eigenvector direction as
#' \eqn{\mu \pm a_i v_i}. Hypervolume is computed with
#' \code{\link{ellipsoid_volume}}, and covariance-derived limits with
#' \code{\link{covariance_limits}}.
#'
#' @param cov_matrix A square, numeric covariance matrix \eqn{\Sigma}. Must be SPD.
#'   Row/column names (if provided) are used as variable names in the output.
#' @param centroid Numeric vector giving the centroid \eqn{\mu}. Must have length
#'   equal to \code{ncol(cov_matrix)}.
#' @param cl Numeric confidence level in (0, 1). Used to compute the chi-square
#'   cutoff defining the ellipsoid contour.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'
#' @return
#' An object of class \code{"nicheR_ellipsoid"} created by
#' \code{\link{new_nicheR_ellipsoid}}, containing ellipsoid geometry and
#' associated quantities (e.g., centroid, covariance matrix, chi-square cutoff,
#' semi-axis lengths, axis vertex coordinates, volume, and covariance limits).
#'
#' @examples
#' cm <- matrix(c(11.11, 0,
#'                0, 17777.78),
#'              nrow = 2,
#'              byrow = TRUE)
#' colnames(cm) <- rownames(cm) <- c("var1", "var2")
#' ctr <- c(20, 600)
#' ell <- ellipsoid_calculator(cov_matrix = cm, centroid = ctr, cl = 0.95, verbose = FALSE)
#'
#' @importFrom stats qchisq
#' @export
ellipsoid_calculator <- function(cov_matrix,
                                 centroid,
                                 cl,
                                 verbose = TRUE) {

  if (missing(cov_matrix) || missing(centroid) || missing(cl)) {
    stop("All arguments ('cov_matrix', 'centroid', 'cl') must be provided.")
  }
  if (!is.matrix(cov_matrix) || !is.numeric(cov_matrix)) {
    stop("'cov_matrix' must be a numeric matrix.")
  }
  if (nrow(cov_matrix) != ncol(cov_matrix)) {
    stop("'cov_matrix' must be a square matrix.")
  }
  if (!is.numeric(centroid) || length(centroid) != ncol(cov_matrix)) {
    stop("'centroid' must be a numeric vector,
    with length equal to the number of columns in 'cov_matrix'.")
  }
  if (!is.numeric(cl) || cl <= 0 || cl >= 1) {
    stop("'cl' must be a numeric value between 0 and 1.")
  }

  verbose_message(verbose, "Step: computing ellipsoid metrics...\n")

  # SPD + inverse
  chol_Sigma <- tryCatch(chol(cov_matrix), error = function(e) NULL)
  if (is.null(chol_Sigma)) stop("Updated covariance matrix is not SPD.")
  Sigma_inv <- chol2inv(chol_Sigma)

  # Cutoff
  chi2_cutoff <- stats::qchisq(cl, df = ncol(cov_matrix))

  # Eigen + axes
  eig <- eigen(cov_matrix, symmetric = TRUE)

  semi_axes_lengths <- sqrt(pmax(eig$values, 0) * chi2_cutoff)
  names(semi_axes_lengths) <- paste0("axis_", seq_along(semi_axes_lengths))

  # Axis coordinates
  axes_coordinates <- lapply(seq_len(ncol(cov_matrix)), function(i) {
    vec <- eig$vectors[, i] * semi_axes_lengths[i]

    ## Create the matrix of vertex coordinates for the current axis
    rbind(vertex_a = centroid - vec, vertex_b = centroid + vec)
  })
  names(axes_coordinates) <- paste0("axis_", seq_along(axes_coordinates))

  volume <- ellipsoid_volume(n_dimensions = ncol(cov_matrix),
                             semi_axes_lengths = semi_axes_lengths)

  cov_limits <- covariance_limits(cov_matrix)

  verbose_message(verbose, "Done: updated ellipsoidal niche metrics")

  new_nicheR_ellipsoid(
    ranges = data.frame(),
    dimensions = ncol(cov_matrix),
    var_names = colnames(cov_matrix),
    centroid = centroid,
    cov_matrix = cov_matrix,
    Sigma_inv = Sigma_inv,
    chol_Sigma = chol_Sigma,
    eigen = list(vectors = eig$vectors, values = eig$values),
    cl = cl,
    chi2_cutoff = chi2_cutoff,
    semi_axes_lengths = semi_axes_lengths,
    axes_coordinates = axes_coordinates,
    volume = volume,
    cov_limits = cov_limits,
    cov_limits_remaining = data.frame()
  )

}




#' Compute ellipsoid hypervolume
#'
#' @description
#'  Computes the geometric volume (area in 2D, volume in 3D, hypervolume in higher
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
#' @usage ellipsoid_volume(n_dimensions, semi_axes_lengths)
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
#' @examples
#' range_df <- data.frame(bio_1 = c(22, 28),
#'                        bio_12 = c(1000, 3500))
#' ell <- nicheR::build_ellipsoid(ranges = range_df)
#'
#' ell$volume
#'
#' # or recalculate
#' nicheR::ellipsoid_volume(n_dimensions = ell$dimensions, semi_axes_lengths = ell$semi_axes_lengths)
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
#' @description
#' Updates one or more off-diagonal covariance values in a
#' \code{nicheR_ellipsoid} object and recomputes all ellipsoid metrics
#' (centroid, semi-axes, volume, etc.) from the new covariance matrix. This
#' allows iterative niche shaping by adjusting the rotation and correlation
#' structure of the ellipsoid without rebuilding it from scratch.
#'
#' @usage update_ellipsoid_covariance(object, covariance,
#'                                    tol = 1e-6,
#'                                    verbose = TRUE)
#'
#' @param object A \code{nicheR_ellipsoid} object, typically created with
#'   \code{\link{build_ellipsoid}}.
#' @param covariance Either a single numeric value applied to all off-diagonal
#'   elements, or a named numeric vector where names identify the variable pair
#'   in the format \code{"var1-var2"} (e.g., \code{c("bio_1-bio_12" = 0.3)}).
#' @param tol Small positive number used as tolerance when computing safe
#'   covariance limits for positive definiteness. Default is \code{1e-6}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @details
#' Covariance values control the orientation and correlation structure of the
#' ellipsoid in environmental space. Setting a positive covariance between two
#' variables tilts the ellipsoid so that high values of one variable tend to
#' co-occur with high values of the other. Negative covariance tilts it in the
#' opposite direction.
#'
#' The updated covariance matrix must remain positive definite. If the
#' requested value would violate this, use \code{\link{covariance_limits}} to
#' find the safe range before calling this function.
#'
#' @return A \code{nicheR_ellipsoid} object with the updated covariance matrix
#'   and recomputed ellipsoid metrics. An additional element
#'   \code{cov_limits_remaining} is attached, giving the remaining safe
#'   covariance limits for any variable pairs that still have zero covariance.
#'
#' @seealso \code{\link{build_ellipsoid}} to create the initial ellipsoid,
#'   \code{\link{covariance_limits}} to find valid covariance ranges before
#'   updating.
#'
#' @examples
#' range_df <- data.frame(bio_1  = c(22, 28),
#'                        bio_12 = c(1000, 3500))
#' ell <- build_ellipsoid(ranges = range_df)
#'
#' # Check covariance allowed covariance range
#' ell$cov_limits
#'
#' # Introduce a negative correlation between bio_1 and bio_12
#' cov_limits <- c("bio_1-bio_12" = -100)
#' ell_corr <- update_ellipsoid_covariance(object = ell,
#'                                           covariance =  cov_limits)
#' ell_corr
#'
#' @export
update_ellipsoid_covariance <- function(object,
                                        covariance,
                                        tol = 1e-6,
                                        verbose = TRUE){

  stopifnot(inherits(object, "nicheR_ellipsoid"))

  verbose_message(verbose, "Starting: updating covariance values...\n")

  up <- update_covariance(object$cov_matrix, covariance = covariance, tol = tol)

  r <- object$ranges
  object <- ellipsoid_calculator(cov_matrix = up$updated_matrix,
                                 centroid = object$centroid,
                                 cl = object$cl,
                                 verbose = verbose)

  # new fields
  object$ranges <- r
  object$cov_limits_remaining <- up$remaining_limits

  object
}


#' Move the centroid of an ellipsoidal niche
#'
#' Translates a `nicheR_ellipsoid` through environmental space, leaving its
#' shape, size, and orientation unchanged.
#'
#' The covariance matrix is not touched, so the semi-axis lengths, the niche
#' volume, and the chi-square cutoff are identical before and after. Only the
#' position changes. The ranges are shifted by the same amount as the centroid,
#' so they continue to describe the extent of the ellipsoid rather than the
#' extent it was originally built from.
#'
#' Use this to ask what a species with the same tolerance breadth would look
#' like centered on different conditions, for example a congener occupying a
#' cooler part of the same environmental space. To change the breadth or the
#' orientation instead, rebuild with [build_ellipsoid()] or rotate with
#' [update_ellipsoid_covariance()].
#'
#' Nothing prevents moving a centroid outside the range of the available
#' environmental data. The resulting ellipsoid is valid, but predicting it onto
#' that data may return little or no suitable area.
#'
#' @usage update_ellipsoid_centroid(ell, new_centroid, verbose = FALSE)
#'
#' @param ell A `nicheR_ellipsoid` object, as returned by [build_ellipsoid()].
#' @param new_centroid A named numeric vector giving the new centroid position.
#' Names must include every entry of `ell$var_names`; extra names are ignored
#' and the order does not matter.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default is
#' `FALSE`.
#'
#' @returns A `nicheR_ellipsoid` object with `centroid` and `ranges` shifted by
#' the difference between `new_centroid` and the original centroid. All other
#' components are unchanged.
#'
#' @seealso [build_ellipsoid()], [update_ellipsoid_covariance()]
#'
#' @examples
#' range_df <- data.frame(bio_1 = c(22, 28),
#'                        bio_12 = c(1000, 3500))
#' ell <- build_ellipsoid(ranges = range_df)
#' ell$centroid
#' ell$ranges
#' rownames(ell$ranges) <- c("min", "max")
#'
#' # Same niche breadth, shifted to cooler and wetter conditions
#' ell_shifted <- update_ellipsoid_centroid(
#'   ell,
#'   new_centroid = c(bio_1 = 19, bio_12 = 2800)
#' )
#'
#' ell_shifted$centroid
#' ell_shifted$ranges
#'
#' # The shape is untouched
#' identical(ell$covariance_matrix, ell_shifted$covariance_matrix)
#'
#' # Both niches in the same environmental space
#' ma_bios <- terra::rast(
#'   system.file("extdata/ma_bios.tif", package = "nicheR"))
#' back_df <- as.data.frame(ma_bios[[c("bio_1", "bio_12")]])
#'
#' plot_ellipsoid(object = ell, background = back_df)
#' add_ellipsoid(object = ell_shifted, col = "red")
#'
#' @export
update_ellipsoid_centroid <- function(ell, new_centroid, verbose = FALSE){

  if(!inherits(ell, "nicheR_ellipsoid")){
    stop("ell must be a nicheR_ellipsoid object.")
  }

  vars <- ell$var_names

  if(is.null(names(new_centroid)) || !all(vars %in% names(new_centroid))){
    stop("new_centroid must be a named numeric vector with names matching ",
         "ell$var_names: ", paste(vars, collapse = ", "))
  }

  new_centroid <- new_centroid[vars]

  if(any(!is.finite(new_centroid))){
    stop("new_centroid contains non-finite values.")
  }

  if(verbose) message("Starting: moving ellipsoid centroid...")

  old_centroid <- ell$centroid
  delta <- new_centroid - old_centroid

  # Shift ranges and axis lengths by the same delta
  new_ranges <- ell$ranges

  for(v in vars){
    new_ranges[, v] <- sort(new_ranges[, v] + delta[[v]])
  }

  rownames(new_ranges) <- c("min", "max")

  ell$axes_coordinates <- lapply(ell$axes_coordinates,
                                 function(a){
                                   sweep(a, 2L, delta, "+")
                                 })

  ell$centroid <- new_centroid
  ell$ranges <- new_ranges

  if(verbose) message("Done: centroid moved.")

  ell
}



#' Resolve prediction input to a canonical form
#'
#' @description
#' Internal helper that normalizes the \code{prediction} argument accepted by
#' sampling functions into a consistent list with a type tag, a resolved data
#' frame or raster, and the name of the prediction column or layer. Handles
#' \code{data.frame}, single- and multi-layer \code{SpatRaster}, and list
#' inputs.
#'
#' @param prediction A \code{SpatRaster}, \code{data.frame}, or list containing
#'   the prediction surface.
#' @param prediction_layer Character. Name of the column or layer to extract.
#'   Required when \code{prediction} contains multiple layers or columns.
#'
#' @return A named list with three elements:
#' \itemize{
#'   \item \code{type}: Character. One of \code{"data.frame"} or
#'   \code{"SpatRaster"}.
#'   \item \code{df} or \code{rast}: The resolved data frame or single-layer
#'   \code{SpatRaster}.
#'   \item \code{pred_name}: Character. The name of the resolved prediction
#'   column or layer.
#' }
#'
#' @keywords internal
#' @noRd
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


#' Check that a weighting method agrees with its prediction layer
#'
#' @description
#' Internal helper shared by \code{sample_data()} and \code{apply_bias()}.
#' Suitability is highest at the niche centroid and Mahalanobis distance is
#' lowest there, so pairing a layer with the wrong \code{method} weights the
#' wrong part of the niche without raising any other error. This stops when
#' the layer name says one kind of layer and \code{method} says the other,
#' and warns when \code{method = "mahalanobis"} is used on a layer whose name
#' gives no hint and whose values all fall in \code{[0, 1]}.
#'
#' @param method Character. \code{"suitability"} or \code{"mahalanobis"}.
#' @param layer_name Character. Name of the prediction layer or column.
#' @param max_value Numeric. Largest value in the layer, or \code{NULL} to
#'   skip the value check.
#'
#' @return Invisibly \code{TRUE}. Called for its checks.
#'
#' @keywords internal
#' @noRd
check_method_layer <- function(method, layer_name, max_value = NULL){

  shown <- ""
  if(!is.null(layer_name) && length(layer_name) > 0L && !is.na(layer_name[1])){
    shown <- as.character(layer_name[1])
  }

  nm <- tolower(shown)
  looks_suit <- grepl("suit", nm, fixed = TRUE)
  looks_mahal <- grepl("mahal", nm, fixed = TRUE)

  if(method == "mahalanobis" && looks_suit && !looks_mahal){
    stop(
      "Layer '", shown, "' looks like a suitability layer, but ",
      "method = 'mahalanobis'. Use method = 'suitability', or rename the ",
      "layer if it holds Mahalanobis distances."
    )
  }

  if(method == "suitability" && looks_mahal && !looks_suit){
    stop(
      "Layer '", shown, "' looks like a Mahalanobis distance layer, but ",
      "method = 'suitability'. Use method = 'mahalanobis', or rename the ",
      "layer if it holds suitability values."
    )
  }

  if(method == "mahalanobis" && !looks_mahal && !is.null(max_value) &&
     is.finite(max_value) && max_value <= 1){
    warning(
      "method = 'mahalanobis', but every value in layer '", shown,
      "' is in [0, 1], which is what a suitability layer looks like. ",
      "Check that this layer holds Mahalanobis distances."
    )
  }

  invisible(TRUE)
}




