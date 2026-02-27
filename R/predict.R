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
#' Suitability is computed as a monotonic transform:
#' \deqn{
#' s(x) = \exp\left(-\frac{1}{2} D^2(x)\right)
#' }
#'
#' A point is considered inside the ellipsoid at probability level \code{cl} when:
#' \deqn{
#' D^2(x) \le \chi^2_p(\mathrm{cl})
#' }
#' where \eqn{p} is the number of predictors.
#'
#' Depending on user options, the following outputs may be returned:
#'
#' \itemize{
#'   \item \code{"Mahalanobis"}: squared Mahalanobis distance \eqn{D^2}
#'   \item \code{"suitability"}: \eqn{\exp(-\frac{1}{2} D^2)}
#'   \item \code{"Mahalanobis_trunc"}: Mahalanobis distance truncated at the
#'         chi square cutoff (values outside set to NA)
#'   \item \code{"suitability_trunc"}: suitability truncated at the chi square
#'         cutoff (values outside set to 0)
#' }
#'
#' When \code{newdata} is a \code{terra::SpatRaster} (or a \code{raster::Raster*}
#' object coerced to \code{SpatRaster}), predictions are computed per cell and
#' returned as a multi layer \code{SpatRaster}. When \code{newdata} is a
#' \code{data.frame} or \code{matrix}, the function returns a \code{data.frame}
#' with added output columns using the names listed above.
#'
#' @param object A \code{nicheR_ellipsoid} object produced by
#'   \code{\link{build_ellipsoid}}. Must contain \code{dimensions},
#'   \code{centroid}, \code{cov_matrix}, \code{Sigma_inv}, \code{truncation_level},
#'   and \code{var_names}.
#' @param newdata Environmental data to predict to. Supported inputs are:
#'   \code{terra::SpatRaster}, \code{raster::RasterLayer},
#'   \code{raster::RasterStack}, \code{raster::RasterBrick},
#'   \code{data.frame}, \code{tibble}, or \code{matrix}. Predictor names must
#'   match \code{object$var_names}.
#' @param adjust_truncation_level Optional numeric scalar in (0, 1). If provided,
#'   overrides \code{object$truncation_level} for truncated outputs.
#' @param include_suitability Logical; if \code{TRUE}, include the
#'   \code{"suitability"} output.
#' @param suitability_truncated Logical; if \code{TRUE}, include the
#'   \code{"suitability_trunc"} output.
#' @param include_mahalanobis Logical; if \code{TRUE}, include the
#'   \code{"Mahalanobis"} output.
#' @param mahalanobis_truncated Logical; if \code{TRUE}, include the
#'   \code{"Mahalanobis_trunc"} output.
#' @param keep_data Logical or \code{NULL}. Controls whether predictor data are
#'   retained in the returned object.
#'
#'   If \code{NULL} (default), behavior depends on the class of \code{newdata}:
#'   \itemize{
#'     \item For \code{data.frame} or \code{matrix} inputs, \code{keep_data} defaults to \code{TRUE}
#'           (predictor columns are retained).
#'     \item For \code{terra::SpatRaster} inputs, \code{keep_data} defaults to \code{FALSE}
#'           (only prediction layers are returned).
#'   }
#'
#'   If explicitly set to \code{TRUE} or \code{FALSE}, that value overrides the default behavior.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' If \code{newdata} is a \code{SpatRaster}, returns a multi layer
#' \code{SpatRaster} with layer names matching the selected outputs
#' (and optionally predictor layers when \code{keep_data = TRUE}).
#'
#' If \code{newdata} is a \code{data.frame} or \code{matrix}, returns a
#' \code{data.frame} with added columns named according to the selected outputs
#' (and optionally predictor columns when \code{keep_data = TRUE}).
#'
#' @method predict nicheR_ellipsoid
#' @export
predict.nicheR_ellipsoid <- function(object,
                                     newdata,
                                     adjust_truncation_level = NULL,
                                     include_suitability = TRUE,
                                     suitability_truncated = FALSE,
                                     include_mahalanobis = TRUE,
                                     mahalanobis_truncated = FALSE,
                                     keep_data = NULL,
                                     verbose = TRUE,
                                     ...){

  verbose_message <- function(...) if(isTRUE(verbose)) cat(...)

  # Basic object checks -----------------------------------------------------

  if(!inherits(object, "nicheR_ellipsoid")){
    stop("'object' must be a nicheR_ellipsoid produced by build_ellipsoid().")
  }

  required_fields <- c("dimensions", "centroid", "cov_matrix", "Sigma_inv", "truncation_level", "var_names")
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

  verbose_message("Starting: suitbaility prediction using newdata of class: ",
    paste(class(newdata), collapse = ", "),
    "...\n")

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
  truncation_level <- object$truncation_level

  # Cutoff handling
  if(is.null(adjust_truncation_level)){
    truncation_threshold <- stats::qchisq(truncation_level, df = dimensions)
  }else{
    if(!is.numeric(adjust_truncation_level) || length(adjust_truncation_level) != 1L || !is.finite(adjust_truncation_level) || adjust_truncation_level <= 0 || adjust_truncation_level >= 1){
      stop("'adjust_truncation_level' must be a single finite number strictly between 0 and 1.")
    }
    truncation_threshold <- stats::qchisq(adjust_truncation_level, df = dimensions)
    truncation_level <- adjust_truncation_level
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

  spatial_names <- c("x", "y", "lon", "lat", "longitude", "latitude")

  coords_lower  <- tolower(colnames(newdata))
  spatial_cols  <- colnames(newdata)[coords_lower %in% spatial_names]


  if(length(spatial_cols) > 0){

    verbose_message("Step: Identified spatial columns: ",
                    paste(spatial_cols, collapse = ", "),
                    "\n")

    # Remove spatial columns from extra vars
    extra_vars <- setdiff(extra_vars, spatial_cols)

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
    newdata <- newdata[, c(spatial_cols, var_names), drop = FALSE]
  }

  if (is.null(keep_data)){
    if (inherits(newdata, "SpatRaster")){
      keep_data <- FALSE
    }else if (inherits(newdata, c("data.frame", "matrix"))) {
      keep_data <- TRUE
    }else {
      keep_data <- FALSE
    }
  }

  if (!is.logical(keep_data) || length(keep_data) != 1L) {
    stop("keep_data must be TRUE or FALSE")
  }


  # Predict: SpatRaster -----------------------------------------------------
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
      S <- D2
      S[!is.na(S)] <- exp(-0.5 * S)

      names(S) <- "suitability"
      if(isTRUE(include_suitability)) out_rast$suitability <- S
    }

    if(isTRUE(mahalanobis_truncated)){
      Mt <- D2
      Mt[!is.na(Mt) & Mt > truncation_threshold] <- NA_real_

      names(Mt) <- "Mahalanobis_trunc_"
      out_rast$Mahalanobis_trunc <- Mt
    }

    if(isTRUE(suitability_truncated)){
      St  <- D2
      St [!is.na(St)] <- exp(-0.5 * St)

      St[!is.na(St) & St > exp(-0.5 * truncation_threshold)] <- 0

      names(St) <- "suitability_trunc"
      out_rast$suitability_trunc <- St
    }


    if(keep_data){
      out_rast <- c(newdata, terra::rast(out_rast))
    }else{
      out_rast <- terra::rast(out_rast)
    }

    verbose_message("Done: Prediction completed successfully. ",
                    "Returned raster layers: ",
                    paste(names(out_rast), collapse = ", "),
                    "\n")

    return(out_rast)
  }

  # Predict: data.frame -----------------------------------------------------
  # detect common spatial column names

  if(length(spatial_cols) > 0){
    if(keep_data){
      out_df <- newdata[ , c(spatial_cols, var_names), drop = FALSE]
    }else{
      out_df <- newdata[, spatial_cols, drop = FALSE]
    }
  }else{
    if(keep_data){
      out_df <- newdata[ , var_names, drop = FALSE]
    }else{
      out_df <- data.frame(row_id = seq_len(nrow(newdata)))
    }
  }

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
    out_df$Mahalanobis_trunc[cc] <- ifelse(D2v <= truncation_threshold, D2v, NA_real_)
  }

  if(isTRUE(suitability_truncated)){
    out_df$suitability_trunc <- NA_real_
    s <- exp(-0.5 * D2v)
    out_df$suitability_trunc[cc] <- ifelse(D2v <= truncation_threshold, s, 0)
  }

  out_df$row_id <- NULL

  verbose_message("Done: Prediction completed successfully. Returned columns: ",
                  paste(colnames(out_df), collapse = ", "),
                  "\n")

  out_df <- do.call(cbind, out_df)
  class(out_df) <- "nicheR_prediction"

  out_df
}
