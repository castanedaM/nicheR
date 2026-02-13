#' Predict Mahalanobis distance and suitability from a nicheR ellipsoid
#'
#' Computes squared Mahalanobis distance \eqn{D^2} and a monotonic suitability
#' transform for new environmental data under a probabilistic ellipsoid niche.
#' The ellipsoid is defined by a multivariate normal (MVN) distribution with
#' mean vector \eqn{\mu} and inverse covariance matrix \eqn{\Sigma^{-1}} stored
#' in a \code{nicheR_ellipsoid} object.
#'
#' The squared Mahalanobis distance is:
#' \deqn{
#' D^2(x) = (x - \mu)^T \Sigma^{-1} (x - \mu)
#' }
#'
#' Suitability is computed as a convenient monotonic transform:
#' \deqn{
#' s(x) = \exp\left(-\frac{1}{2} D^2(x)\right)
#' }
#'
#' A point is considered inside the ellipsoid at probability \code{cl} when:
#' \deqn{
#' D^2(x) \le \chi^2_p(\text{cl})
#' }
#' where \eqn{p} is the number of predictors (dimensions).
#'
#' When \code{newdata} is a \code{terra::SpatRaster}, computations are performed
#' per cell and the function returns a multi-layer \code{SpatRaster}. When
#' \code{newdata} is a \code{data.frame} or \code{matrix}, the function returns a
#' \code{data.frame} with added output columns.
#'
#' @param object A \code{nicheR_ellipsoid} object produced by
#'   \code{\link{build_ellipsoid}}. Must contain \code{dimensions},
#'   \code{centroid}, \code{cov_matrix}, \code{Sigma_inv}, \code{chi2_cutoff}, and
#'   \code{var_names}.
#' @param newdata Environmental data to predict to. Supported inputs are:
#'   \code{terra::SpatRaster}, \code{raster::RasterLayer},
#'   \code{raster::RasterStack}, \code{raster::RasterBrick},
#'   \code{data.frame}, \code{tibble}, or \code{matrix}.
#' @param adjust_level Numeric scalar in (0, 1). Optional probability level used
#'   to compute a new chi-square cutoff for truncated outputs.
#' @param include_suitability Logical; if \code{TRUE}, include suitability
#'   \eqn{s(x)} in the output.
#' @param suitability_truncated Logical; if \code{TRUE}, include a truncated
#'   suitability layer/column where values outside the chi-square cutoff are set
#'   to 0 (and NA values remain NA).
#' @param include_mahalanobis Logical; if \code{TRUE}, include squared Mahalanobis
#'   distance \eqn{D^2} in the output.
#' @param mahalanobis_truncated Logical; if \code{TRUE}, include a binary inside
#'   mask for Mahalanobis distance at \code{adjust_level} (or object level if
#'   \code{adjust_level} is \code{NULL}). Values are 1 for
#'   \eqn{D^2 \le \chi^2_p(level)}, 0 otherwise (NA values remain NA).
#' @param verbose Logical; if \code{TRUE}, print brief progress messages.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' If \code{newdata} is a \code{terra::SpatRaster} (or a \code{raster::Raster*}
#' converted to \code{SpatRaster}), returns a multi-layer \code{SpatRaster}
#' containing the requested outputs.
#'
#' If \code{newdata} is a \code{data.frame} or \code{matrix}, returns a
#' \code{data.frame} with the same number of rows as \code{newdata} and added
#' columns matching the requested outputs. Rows with missing predictor values
#' receive NA for numeric outputs.

#' @method predict nicheR_ellipsoid
#' @export
predict.nicheR_ellipsoid <- function(object,
                                     newdata,
                                     adjust_level = NULL,
                                     include_suitability = TRUE,
                                     suitability_truncated = FALSE,
                                     include_mahalanobis = TRUE,
                                     mahalanobis_truncated = FALSE,
                                     verbose = TRUE,
                                     ...){

  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)

  # Basic object checks -----------------------------------------------------

  if(!inherits(object, "nicheR_ellipsoid")){
    stop("'object' must be a nicheR_ellipsoid produced by build_ellipsoid().")
  }

  required_fields <- c("dimensions", "centroid", "cov_matrix", "Sigma_inv", "chi2_cutoff", "var_names")
  missing_fields <- required_fields[!required_fields %in% names(object)]
  if(length(missing_fields) > 0){
    stop("object is missing required fields: ", paste(missing_fields, collapse = ", "))
  }

  if(!is.numeric(object$centroid) || length(object$centroid) != object$dimensions){
    stop("object$centroid must be a numeric vector of length object$dimensions.")
  }
  if(!is.matrix(object$cov_matrix) || any(dim(object$cov_matrix) != c(object$dimensions, object$dimensions))){
    stop("object$cov_matrix must be a square matrix with dimensions x dimensions.")
  }
  if(!is.matrix(object$Sigma_inv) || any(dim(object$Sigma_inv) != c(object$dimensions, object$dimensions))){
    stop("object$Sigma_inv must be a square matrix with dimensions x dimensions.")
  }
  if(!is.character(object$var_names) || length(object$var_names) != object$dimensions){
    stop("object$var_names must be a character vector of length object$dimensions.")
  }

  verbose_message("Starting: predict suitab...\n")


  # Sanity checks for inclusions
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

  # Core objects attributes
  dimensions <- object$dimensions
  mu <- object$centroid
  Sigma_inv <- object$Sigma_inv
  var_names <- object$var_names


  # Cutoff handling
  if(is.null(adjust_level)){
    chi2_cutoff <- object$chi2_cutoff
  }else{
    if(!is.numeric(adjust_level) || length(adjust_level) != 1L || !is.finite(adjust_level) || adjust_level <= 0 || adjust_level >= 1){
      stop("'adjust_level' must be a single finite number strictly between 0 and 1.")
    }
    chi2_cutoff <- stats::qchisq(adjust_level, df = dimensions)
  }


  # Coerce newdata
  if(inherits(newdata, c("RasterLayer", "RasterStack", "RasterBrick"))){
    newdata <- terra::rast(newdata)
  }else if(inherits(newdata, "SpatRaster")){
    # keep
  }else if(inherits(newdata, c("data.frame", "matrix", "tbl_df"))){
    newdata <- as.data.frame(newdata)
  }else{
    stop("'newdata' must be a SpatRaster, Raster*, data.frame, tibble, or matrix.")
  }

  # Variable matching
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
         "Variables expected: ",
         paste(var_names, collapse = ", "))
  }

  if(length(missing_vars) > 0){
    stop("newdata is missing required predictor variables: ",
         paste(missing_vars, collapse = ", "))
  }

  if(length(extra_vars) > 0){
    verbose_message("Step: Ignoring extra predictor columns: ",
                    paste(extra_vars, collapse = ", "),
                    "\n")
  }

  verbose_message("Step: Using ", length(var_names), " predictor variables: ",
                  paste(var_names, collapse = ", "),
                  "\n")

  # Subset and reorder
  if(inherits(newdata, "SpatRaster")){
    newdata <- newdata[[var_names]]
  }else{
    newdata <- newdata[, var_names, drop = FALSE]
  }

  # Predict: SpatRaster -----------------------------------------------------
#   To do: use ifelse instead of app, smae thing
  if(inherits(newdata, "SpatRaster")){

    D2 <- terra::app(newdata, fun = function(v){
      if(anyNA(v) || any(!is.finite(v))) return(NA_real_)
      d <- v - mu
      as.numeric(t(d) %*% Sigma_inv %*% d)
    })
    names(D2) <- "Mahalanobis"

    out_rast <- list( )

    if(isTRUE(include_mahalanobis)){
      out_rast$Mahalanobis <- D2
    }

    if(isTRUE(include_suitability) || isTRUE(suitability_truncated)){
      S <- terra::app(D2, fun = function(z){
        ifelse(is.na(z), NA_real_, exp(-0.5 * z))
      })
      names(S) <- "suitability"
      if(isTRUE(include_suitability)) out_rast$suitability <- S
    }

    if(isTRUE(mahalanobis_truncated)){
      Mt <- terra::app(D2, fun = function(z){
        ifelse(is.na(z), NA_real_, as.numeric(z <= chi2_cutoff))
      })
      names(Mt) <- "Mahalanobis_trunc"
      out_rast$Mahalanobis_trunc <- Mt
    }

    if(isTRUE(suitability_truncated)){
      St <- terra::app(D2, fun = function(z){
        ifelse(is.na(z), NA_real_,
               ifelse(z <= chi2_cutoff, exp(-0.5 * z), 0))
      })
      names(St) <- "suitability_trunc"
      out_rast$suitability_trunc <- St
    }

    verbose_message("Done: prediction completed.\n")

    return(terra::rast(out_rast))
  }

  # Predict: data.frame -----------------------------------------------------
  out_df <- data.frame(row_id = seq_len(nrow(newdata)))
  cc <- stats::complete.cases(newdata)

  if(!any(cc)) stop("All rows contain NA in predictor columns.")

  pts <- as.matrix(newdata[cc, var_names, drop = FALSE])
  diffs <- sweep(pts, 2, mu, "-")
  D2v <- rowSums((diffs %*% Sigma_inv) * diffs)

  if(isTRUE(include_mahalanobis)){
    out_df$Mahalanobis <- NA_real_
    out_df$Mahalanobis[cc] <- D2v
  }

  if(isTRUE(include_suitability)){
    out_df$suitability <- NA_real_
    out_df$suitability[cc] <- exp(-0.5 * D2v)
  }

  if(isTRUE(mahalanobis_truncated)){
    out_df$Mahalanobis_trunc <- NA_real_
    out_df$Mahalanobis_trunc[cc] <- as.numeric(D2v <= chi2_cutoff)
  }

  if(isTRUE(suitability_truncated)){
    out_df$suitability_trunc <- NA_real_
    s <- exp(-0.5 * D2v)
    out_df$suitability_trunc[cc] <- ifelse(D2v <= chi2_cutoff, s, 0)
  }

  out_df$row_id <- NULL

  verbose_message("Done: prediction completed.\n")

  out_df
}
