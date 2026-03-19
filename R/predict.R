#' Predict suitability and Mahalanobis distance from a nicheR ellipsoid
#'
#' @description
#' Computes Mahalanobis distance and optional suitability values from a
#' \code{nicheR_ellipsoid} for either (1) environmental samples provided as a
#' \code{data.frame} or \code{matrix}, (2) a \code{SpatRaster} stack of
#' predictors, or (3) virtual samples drawn in environmental space when
#' \code{newdata = NULL}.
#'
#' @param object A \code{nicheR_ellipsoid} object produced by
#'   \code{\link{build_ellipsoid}}.
#' @param newdata Environmental predictors. One of:
#'   \itemize{
#'     \item A \code{SpatRaster} (or legacy \code{raster} classes, coerced
#'     automatically).
#'     \item A \code{data.frame}, \code{tibble}, or \code{matrix} with columns
#'     named to match \code{object$var_names}.
#'   }
#' @param adjust_truncation_level Optional numeric confidence level in
#'   \code{(0, 1)} to override \code{object$cl} when computing truncated
#'   outputs. Default is \code{NULL} (uses the level stored in \code{object}).
#' @param include_suitability Logical. If \code{TRUE} (default), returns
#'   suitability values (\eqn{\exp(-0.5 D^2)}).
#' @param suitability_truncated Logical. If \code{TRUE}, returns a truncated
#'   suitability layer where values outside the chi-square contour are set to
#'   \code{0}. Default is \code{FALSE}.
#' @param include_mahalanobis Logical. If \code{TRUE} (default), returns
#'   Mahalanobis distance (\eqn{D^2}).
#' @param mahalanobis_truncated Logical. If \code{TRUE}, returns a truncated
#'   Mahalanobis layer where values outside the chi-square contour are set to
#'   \code{NA}. Default is \code{FALSE}.
#' @param keep_data Logical or \code{NULL}. If \code{TRUE}, includes the
#'   original predictors in the output. Default is \code{NULL}: \code{FALSE}
#'   for \code{SpatRaster} input, \code{TRUE} for tabular input.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @details
#' Suitability is computed as \eqn{\exp(-0.5 D^2)}, where \eqn{D^2} is the
#' squared Mahalanobis distance from the niche centroid. Truncated outputs use
#' a chi-square cutoff based on the ellipsoid confidence level (\code{cl}).
#'
#' For tabular inputs, coordinate columns (e.g., \code{x}, \code{y},
#' \code{lon}, \code{lat}) are detected and retained when
#' \code{keep_data = TRUE}. Extra non-predictor columns are ignored.
#'
#' @return
#' If \code{newdata} is a \code{SpatRaster}, returns a \code{SpatRaster} with
#' the requested prediction layers (and optionally the original predictor
#' layers if \code{keep_data = TRUE}).
#'
#' If \code{newdata} is tabular, returns a \code{data.frame} of class
#' \code{"nicheR_prediction"} with the requested prediction columns (and
#' optionally the original predictor columns if \code{keep_data = TRUE}).
#'
#' @method predict nicheR_ellipsoid
#' @importFrom stats qchisq complete.cases
#' @importFrom terra rast app
#'
#' @export
predict.nicheR_ellipsoid <- function(object,
                                     newdata,
                                     adjust_truncation_level = NULL,
                                     include_suitability = TRUE,
                                     suitability_truncated = FALSE,
                                     include_mahalanobis = TRUE,
                                     mahalanobis_truncated = FALSE,
                                     keep_data = NULL,
                                     verbose = TRUE){

  # Basic object checks -----------------------------------------------------

  if(!inherits(object, "nicheR_ellipsoid")){
    stop("'object' must be a nicheR_ellipsoid produced by build_ellipsoid().")
  }

  required_fields <- c("dimensions", "centroid", "cov_matrix", "Sigma_inv",
                       "cl", "var_names")
  missing_fields <- required_fields[!required_fields %in% names(object)]
  if(length(missing_fields) > 0){
    stop("object is missing required fields: ", paste(missing_fields, collapse = ", "))
  }

  if(!is.numeric(object$centroid) || length(object$centroid) != object$dimensions){
    stop("object$centroid must be a numeric vector of length object$dimensions.")
  }
  if(!is.matrix(object$cov_matrix) ||
     any(dim(object$cov_matrix) != c(object$dimensions, object$dimensions))){
    stop("object$cov_matrix must be a square matrix with dimensions x dimensions.")
  }
  if(!is.matrix(object$Sigma_inv) ||
     any(dim(object$Sigma_inv) != c(object$dimensions, object$dimensions))){
    stop("object$Sigma_inv must be a square matrix with dimensions x dimensions.")
  }
  if(!is.character(object$var_names) || length(object$var_names) != object$dimensions){
    stop("object$var_names must be a character vector of length object$dimensions.")
  }

  # Sanity checks for inclusions -------------------------------------------

  if(!is.logical(include_suitability) || length(include_suitability) != 1L){
    stop("include_suitability must be TRUE/FALSE.")
  }
  if(!is.logical(include_mahalanobis) || length(include_mahalanobis) != 1L){
    stop("include_mahalanobis must be TRUE/FALSE.")
  }
  if(!is.logical(suitability_truncated) || length(suitability_truncated) != 1L){
    stop("suitability_truncated must be TRUE/FALSE.")
  }
  if(!is.logical(mahalanobis_truncated) || length(mahalanobis_truncated) != 1L){
    stop("mahalanobis_truncated must be TRUE/FALSE.")
  }

  # Core object attributes --------------------------------------------------

  dimensions <- object$dimensions
  mu <- object$centroid
  Sigma_inv <- object$Sigma_inv
  var_names <- object$var_names
  truncation_level <- object$cl

  verbose_message(verbose,
                  paste0("Starting: suitability prediction using newdata of class: ",
                         paste(class(newdata), collapse = ", "),
                         "...\n"))

  # Cutoff handling ---------------------------------------------------------

  if(is.null(adjust_truncation_level)){
    truncation_threshold <- stats::qchisq(truncation_level, df = dimensions)
  }else{
    if(!is.numeric(adjust_truncation_level) || length(adjust_truncation_level) != 1L ||
       !is.finite(adjust_truncation_level) || adjust_truncation_level <= 0 ||
       adjust_truncation_level >= 1){
      stop("'adjust_truncation_level' must be a single finite number strictly between 0 and 1.")
    }
    truncation_threshold <- stats::qchisq(adjust_truncation_level, df = dimensions)
    truncation_level <- adjust_truncation_level
  }

  # Coerce newdata ----------------------------------------------------------

  if(inherits(newdata, c("RasterLayer", "RasterStack", "RasterBrick"))){
    newdata <- terra::rast(newdata)
  }else if(inherits(newdata, "SpatRaster")){
    # keep as-is
  }else if(inherits(newdata, c("data.frame", "matrix", "tbl_df"))){
    newdata <- as.data.frame(newdata)
  }else{
    stop("'newdata' must be a SpatRaster, Raster*, data.frame, tibble, or matrix.")
  }

  # Variable matching -------------------------------------------------------

  if(inherits(newdata, "SpatRaster")){
    nd_names <- names(newdata)
    if(is.null(nd_names)) stop("SpatRaster newdata must have named layers.")
  }else{
    nd_names <- names(newdata)
    if(is.null(nd_names)) stop("newdata must have column names matching object$var_names.")
  }

  used_vars <- intersect(var_names, nd_names)
  missing_vars <- setdiff(var_names, nd_names)
  extra_vars <- setdiff(nd_names, var_names)

  if(length(used_vars) == 0){
    stop("No matching predictor variables found.\n",
         "Variables expected: ", paste(var_names, collapse = ", "))
  }
  if(length(missing_vars) > 0){
    stop("newdata is missing required predictor variables: ",
         paste(missing_vars, collapse = ", "))
  }

  spatial_names <- c("x", "y", "lon", "lat", "longitude", "latitude")
  spatial_cols <- character(0)

  if(!inherits(newdata, "SpatRaster")){
    coords_lower <- tolower(colnames(newdata))
    spatial_cols <- colnames(newdata)[coords_lower %in% spatial_names]

    if(length(spatial_cols) > 0){
      verbose_message(verbose, "Step: Identified spatial columns: ",
                      paste(spatial_cols, collapse = ", "), "\n")
      extra_vars <- setdiff(extra_vars, spatial_cols)
    }
  }

  if(length(extra_vars) > 0){
    verbose_message(verbose, "Step: Ignoring extra predictor columns: ",
                    paste(extra_vars, collapse = ", "), "\n")
  }

  verbose_message(verbose, "Step: Using ", length(var_names),
                  " predictor variables: ",
                  paste(var_names, collapse = ", "), "\n")

  # Subset and reorder ------------------------------------------------------

  if(inherits(newdata, "SpatRaster")){
    newdata <- newdata[[var_names]]
  }else{
    newdata <- newdata[, c(spatial_cols, var_names), drop = FALSE]
  }

  # keep_data defaulting ----------------------------------------------------

  if(is.null(keep_data)){
    keep_data <- !inherits(newdata, "SpatRaster")
  }
  if(!is.logical(keep_data) || length(keep_data) != 1L){
    stop("keep_data must be TRUE or FALSE.")
  }

  # Predict: SpatRaster -----------------------------------------------------

  if(inherits(newdata, "SpatRaster")){
    D2 <- terra::app(newdata, fun = function(v){
      if(anyNA(v) || any(!is.finite(v))) return(NA_real_)
      d <- v - mu
      as.numeric(t(d) %*% Sigma_inv %*% d)
    })
    names(D2) <- "Mahalanobis"

    out_rast <- list()

    if(isTRUE(include_mahalanobis)) out_rast$Mahalanobis <- D2

    if(isTRUE(include_suitability)){
      S <- exp(-0.5 * D2)
      names(S) <- "suitability"
      out_rast$suitability <- S
    }

    if(isTRUE(mahalanobis_truncated)){
      Mt <- D2
      Mt[!is.na(Mt) & Mt > truncation_threshold] <- NA_real_
      names(Mt) <- "Mahalanobis_trunc"
      out_rast$Mahalanobis_trunc <- Mt
    }

    if(isTRUE(suitability_truncated)){
      St <- exp(-0.5 * D2)
      St[!is.na(St) & D2 > truncation_threshold] <- 0
      names(St) <- "suitability_trunc"
      out_rast$suitability_trunc <- St
    }

    out_rast <- if(isTRUE(keep_data)){
      c(newdata, terra::rast(out_rast))
    }else{
      terra::rast(out_rast)
    }

    verbose_message(verbose,
                    "Done: Prediction completed successfully. Returned raster layers: ",
                    paste(names(out_rast), collapse = ", "), "\n")
    return(out_rast)
  }

  # Predict: data.frame -----------------------------------------------------

  if(length(spatial_cols) > 0){
    out_df <- if(isTRUE(keep_data)){
      newdata[, c(spatial_cols, var_names), drop = FALSE]
    }else{
      newdata[, spatial_cols, drop = FALSE]
    }
  }else{
    out_df <- if(isTRUE(keep_data)){
      newdata[, var_names, drop = FALSE]
    }else{
      data.frame(row_id = seq_len(nrow(newdata)))
    }
  }

  cc  <- stats::complete.cases(newdata[, var_names, drop = FALSE])
  if(!any(cc)) stop("All rows contain NA in predictor columns.")

  pts   <- as.matrix(newdata[cc, var_names, drop = FALSE])
  diffs <- sweep(pts, 2, mu, "-")
  D2v   <- rowSums((diffs %*% Sigma_inv) * diffs)

  if(isTRUE(include_mahalanobis)){
    out_df$Mahalanobis        <- NA_real_
    out_df$Mahalanobis[cc]    <- D2v
  }

  if(isTRUE(include_suitability)){
    out_df$suitability        <- NA_real_
    out_df$suitability[cc]    <- exp(-0.5 * D2v)
  }

  if(isTRUE(mahalanobis_truncated)){
    out_df$Mahalanobis_trunc     <- NA_real_
    out_df$Mahalanobis_trunc[cc] <- ifelse(D2v <= truncation_threshold,
                                           D2v, NA_real_)
  }

  if(isTRUE(suitability_truncated)){
    s <- exp(-0.5 * D2v)
    out_df$suitability_trunc     <- NA_real_
    out_df$suitability_trunc[cc] <- ifelse(D2v <= truncation_threshold, s, 0)
  }

  if("row_id" %in% names(out_df)) out_df$row_id <- NULL

  verbose_message(verbose,
                  "Done: Prediction completed successfully. Returned columns: ",
                  paste(colnames(out_df), collapse = ", "), "\n")

  class(out_df) <- c("nicheR_prediction", class(out_df))
  out_df
}
