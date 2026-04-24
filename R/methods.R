#' Print method for nicheR objects
#'
#' Provides a concise summary of \code{nicheR} objects.
#'
#' @name print
#' @aliases print,nicheR_nicheR_ellipsoid-method
#' @aliases print,nicheR_nicheR_community-method
#' 
#' @rdname print
#'
#' @param x An object of the classes \code{"nicheR_ellipsoid"} or
#'   \code{"nicheR_community"}.
#' @param digits Integer. Number of decimal places used when printing numeric
#'   values. Default is 3.
#' @param ... Additional arguments.
#'
#' @details
#' The function formats and rounds key quantities for readability but does
#' not modify the underlying object.
#'
#' @return
#' The input object \code{x}, returned invisibly.
#'
#' @seealso \code{\link{build_ellipsoid}} and \code{\link{generate_community}}
#'
#' @method print nicheR_ellipsoid
#' @export

print.nicheR_ellipsoid <- function(x, digits = 3, ...) {

  cat("nicheR Ellipsoid Object\n")
  cat("-----------------------\n")

  cat("Dimensions:        ", x$dimensions, "D\n", sep = "")
  cat("Chi-square cutoff: ", round(x$chi2_cutoff, digits), "\n", sep = "")

  cat("Centroid (mu):     ",
      paste(round(x$centroid, digits), collapse = ", "),
      "\n", sep = "")

  cat("\nCovariance matrix:\n")
  print(round(x$cov_matrix, digits))

  cat("\nCovariance Limits:\n")
  cov_lims <- x$cov_limits
  rownames(cov_lims) <- apply(
    combn(x$var_names, 2), 2,
    function(pair) paste(pair, collapse = "-")
  )
  print(round(cov_lims, digits))

  cat("\nEllipsoid semi-axis lengths:\n  ",
      paste(round(x$semi_axes_lengths, digits), collapse = ", "),
      "\n", sep = "")

  cat("\nEllipsoid axis endpoints:\n")

  for(i in seq_len(x$dimensions)){
    cat(ifelse(i == 1, "", "\n"), " Axis ", i, ":\n", sep = "")
    print(round(x$axes_coordinates[[i]], digits))
  }

  cat("\nEllipsoid volume:  ", round(x$volume, digits), "\n", sep = "")

  cat("\n")
  invisible(x)
}



#' @rdname print
#' @method print nicheR_community
#' @importFrom stats sd
#' @export

print.nicheR_community <- function(x, digits = 3, ...) {

  cat("nicheR Community Object\n")
  cat("-----------------------\n")

  # Generation details
  cat("Generation Metadata:\n")
  cat("  Pattern:            ", x$details$pattern, "\n", sep = "")
  cat("  Number of ellipses: ", x$details$n, "\n", sep = "")
  cat("  Smallest prop.:     ", round(x$details$smallest_proportion, digits),
      "\n", sep = "")

  ## Only print these if they aren't NA
  if (!is.na(x$details$largest_proportion)) {
    cat("  Largest prop.:      ", round(x$details$largest_proportion, digits),
        "\n", sep = "")
  }
  if (!is.na(x$details$bias)) {
    cat("  Bias exponent:      ", round(x$details$bias, digits), "\n", sep = "")
  }
  if (!is.na(x$details$thin_background)) {
    cat("  Thin background:    ", x$details$thin_background, "\n", sep = "")
  }
  if (!is.na(x$details$resolution)) {
    cat("  Resolution:         ", x$details$resolution, "\n", sep = "")
  }
  if (!is.na(x$details$seed)) {
    cat("  Random seed:        ", x$details$seed, "\n", sep = "")
  }

  # Reference ellipsoid summary
  cat("\nReference ellipsoid summary:\n")
  cat("  Dimensions:        ", x$reference$dimensions, "D\n", sep = "")
  cat("  Variables:         ", paste(x$reference$var_names, collapse = ", "),
      "\n", sep = "")
  cat("  Centroid (mu):     ",
      paste(round(x$reference$centroid, digits), collapse = ", "),
      "\n", sep = "")
  cat("  Ellipsoid volume:  ", round(x$reference$volume, digits),
      "\n", sep = "")

  # A few community summary statistics
  cat("\nCommunity summary (n =", x$details$n, "):\n")

  ## Extract centroids and volumes from the list of ellipsoids
  all_centroids <- do.call(rbind, lapply(x$ellipse_community, function(e) {
    e$centroid
  }))
  all_volumes <- vapply(x$ellipse_community, function(e) e$volume, numeric(1))

  ## Calculate descriptive stats
  mean_cent <- colMeans(all_centroids)
  sd_cent   <- apply(all_centroids, 2, sd)
  mean_vol  <- mean(all_volumes)
  sd_vol    <- sd(all_volumes)

  ## Print centroids
  cat("  Centroid positions | mean (+/-SD):\n")
  for (i in seq_along(mean_cent)) {
    cat("   ", names(mean_cent)[i], ": ",
        round(mean_cent[i], digits),
        " (+/-", round(sd_cent[i], digits), ")\n", sep = "")
  }

  ## Print volumes
  cat("\n  Ellipsoid volumes:\n")
  cat("   Mean: ", round(mean_vol, digits), "\n", sep = "")
  cat("   SD:   ", round(sd_vol, digits), "\n", sep = "")

  cat("\n")
  invisible(x)
}



# Predict methods ------------------------------------------------------

#' Predict suitability and Mahalanobis distance from a nicheR ellipsoid
#'
#' @name predict
#' @aliases predict,nicheR_nicheR_ellipsoid-method
#' @aliases predict,nicheR_nicheR_community-method
#'
#' @rdname predict
#'
#' @description
#' Computes Mahalanobis distance and suitability values deriving from a
#' \code{nicheR_ellipsoid} or \code{nicheR_community} object,
#' for \code{newdata} provided as a \code{data.frame}, \code{matrix}, or
#' \code{SpatRaster}.
#'
#' @param object An object of the classes \code{"nicheR_ellipsoid"} or
#'   \code{"nicheR_community"}.
#' @param newdata Environmental predictors. One of:
#'   \itemize{
#'     \item A \code{SpatRaster} (or legacy \code{raster} classes, coerced
#'     automatically).
#'     \item A \code{data.frame} or \code{matrix} with columns
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
#' @param ... Additional arguments, not currently used.
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
#' For \code{nicheR_ellipsoid} objects, if \code{newdata} is a
#' \code{SpatRaster}, returns a \code{SpatRaster} with the requested
#' prediction as layers (and optionally the original predictor
#' layers if \code{keep_data = TRUE}). If \code{newdata} is tabular,
#' returns a \code{data.frame} with the requested predictions as columns
#' (by default returns the original predictors as columns).
#'
#' For \code{nicheR_community} objects, if \code{newdata} is a
#' \code{SpatRaster}, returns a \code{SpatRaster} where each layer represents
#' predictions for each ellipse. If \code{newdata} is a \code{data.frame},
#' returns a \code{data.frame} with the original data plus one prediction
#' column per ellipse.
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
                                     verbose = TRUE,
                                     ...){
  
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
  
  cc <- stats::complete.cases(newdata[, var_names, drop = FALSE])
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
    out_df$Mahalanobis_trunc[cc] <- ifelse(D2v <= truncation_threshold,
                                           D2v, NA_real_)
  }
  
  if(isTRUE(suitability_truncated)){
    s <- exp(-0.5 * D2v)
    out_df$suitability_trunc <- NA_real_
    out_df$suitability_trunc[cc] <- ifelse(D2v <= truncation_threshold, s, 0)
  }
  
  if("row_id" %in% names(out_df)) out_df$row_id <- NULL
  
  verbose_message(verbose,
                  "Done: Prediction completed successfully. Returned columns: ",
                  paste(colnames(out_df), collapse = ", "), "\n")
  
  class(out_df) <- c("nicheR_prediction", class(out_df))
  out_df
}




#' @param prediction Character. The type of prediction to return. One of:
#'   \code{"Mahalanobis"} (default), \code{"suitability"},
#'   \code{"Mahalanobis_trunc"}, or \code{"suitability_trunc"}.
#'
#' @rdname predict
#' @method predict nicheR_community
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom terra rast
#' @export

predict.nicheR_community <- function(object,
                                     newdata,
                                     prediction = "Mahalanobis",
                                     verbose = TRUE,
                                     ...) {

  # Validation and Setup
  if (missing(object) || missing(newdata)) {
    stop("Arguments 'object' and 'newdata' must be specified.")
  }
  if (!inherits(object, "nicheR_community")) {
    stop("'object' must be of class 'nicheR_community'.")
  }
  if (!inherits(newdata, c("SpatRaster", "data.frame", "matrix"))) {
    stop("'newdata' must be a 'SpatRaster', 'data.frame', or 'matrix'.")
  }
  valid_preds <- c("Mahalanobis", "suitability",
                   "Mahalanobis_trunc", "suitability_trunc")
  if (!(prediction %in% valid_preds)) {
    stop("'prediction' must be one of: ", paste(valid_preds, collapse = ", "))
  }

  is_raster <- inherits(newdata, "SpatRaster")

  # Variable matching check
  var_names <- object$reference$var_names
  nd_names <- if (is_raster) names(newdata) else colnames(newdata)

  missing_vars <- setdiff(var_names, nd_names)

  if (length(missing_vars) > 0) {
    stop("'newdata' is missing required variables: ",
         paste(var_names, collapse = ", "))
  }

  verbose_message(
    verbose,
    paste0("Starting: using newdata of class: ", class(newdata)[1], "...\n")
  )

  # Map the 'prediction' argument to the internal flags
  inc_suit <- prediction == "suitability"
  inc_mahal <- prediction == "Mahalanobis"
  trunc_s <- prediction == "suitability_trunc"
  trunc_m <- prediction == "Mahalanobis_trunc"

  # Number of ellipses in the community
  n_ell <- length(object$ellipse_community)
  ell_names <- paste0("ell_", seq_len(n_ell))

  verbose_message(
    verbose,
    paste0("Predictions for a ", object$details$pattern, " community of ",
           n_ell, " ellipses...\n")
  )

  # Iterate Predictions
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n_ell, style = 3)
  }

  results_list <- lapply(seq_len(n_ell), function(i) {
    p <- predict.nicheR_ellipsoid(object = object$ellipse_community[[i]],
                                  newdata = newdata,
                                  include_suitability = inc_suit,
                                  suitability_truncated = trunc_s,
                                  include_mahalanobis = inc_mahal,
                                  mahalanobis_truncated = trunc_m,
                                  verbose = FALSE)

    if (verbose) utils::setTxtProgressBar(pb, i)

    res <- p[[prediction]]

    if (is_raster) {
      names(res) <- ell_names[i]
    }

    return(res)
  })

  if (verbose) close(pb)

  verbose_message(verbose, "\nFinalizing results...\n")

  # Format and Return Output
  if (is_raster) {
    return(terra::rast(results_list))
  } else {
    results_df <- do.call(cbind, results_list)
    colnames(results_df) <- ell_names
    return(cbind(as.data.frame(newdata), results_df))
  }
}