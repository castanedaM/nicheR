#' Print a nicheR Ellipsoid Object
#'
#' Provides a concise summary of a \code{nicheR_ellipsoid} object created by
#' \code{\link{build_ellipsoid}}. The printed output includes dimensionality,
#' chi-square cutoff, centroid, covariance matrix, principal semi-axis lengths,
#' axis endpoints, and ellipsoid volume.
#'
#' @param x A \code{nicheR_ellipsoid} object.
#' @param digits Integer. Number of decimal places used when printing numeric
#'   values. Default is 3.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This is an S3 method for objects of class \code{"nicheR_ellipsoid"}.
#' The function formats and rounds key quantities for readability but does
#' not modify the underlying object.
#'
#' @return
#' The input object \code{x}, returned invisibly.
#'
#' @seealso \code{\link{build_ellipsoid}}
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



#' Print a nicheR Community Object
#'
#' Provides a structured summary of a \code{nicheR_community} object.
#' Includes the generation metadata, a summary of the reference ellipsoid,
#' and descriptive statistics (mean and standard deviation) for the
#' centroids and volumes of the generated community.
#'
#' @param x A \code{nicheR_community} object.
#' @param digits Integer. Number of decimal places used when printing
#'   numeric values. Default is 3.
#' @param ... Additional arguments passed to the \code{nicheR_ellipsoid}
#'   print method.
#'
#' @method print nicheR_community
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
  cat("  Centroid positions | mean (±SD):\n")
  for (i in seq_along(mean_cent)) {
    cat("   ", names(mean_cent)[i], ": ",
        round(mean_cent[i], digits),
        " (±", round(sd_cent[i], digits), ")\n", sep = "")
  }

  ## Print volumes
  cat("\n  Ellipsoid volumes:\n")
  cat("   Mean: ", round(mean_vol, digits), "\n", sep = "")
  cat("   SD:   ", round(sd_vol, digits), "\n", sep = "")

  cat("\n")
  invisible(x)
}



# Predcit methods ------------------------------------------------------
#' Predict method for a nicheR Community
#'
#' @description
#' Iterates through all ellipses in a \code{nicheR_community} and calculates
#' predictions (suitability or Mahalanobis distance) using
#' \code{\link{predict.nicheR_ellipsoid}}.
#'
#' @param object A \code{nicheR_community} object.
#' @param newdata A \code{SpatRaster}, \code{data.frame}, or \code{matrix}
#'   containing at least the environmental variables used to create the
#'   reference ellipsoid in the community.
#' @param prediction Character. The type of prediction to return. One of:
#'   \code{"Mahalanobis"} (default), \code{"suitability"},
#'   \code{"Mahalanobis_trunc"}, or \code{"suitability_trunc"}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'   Default = \code{TRUE}.
#'
#' @return
#' If \code{newdata} is a \code{SpatRaster}, returns a \code{SpatRaster}
#' where each layer represents one ellipse. If \code{newdata} is a
#' \code{data.frame}, returns a \code{data.frame} with the original data
#' plus one column per ellipse.
#'
#' @method predict nicheR_community
#' @export

predict.nicheR_community <- function(object,
                                     newdata,
                                     prediction = "Mahalanobis",
                                     verbose = TRUE) {

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
  inc_suit  <- prediction == "suitability"
  inc_mahal <- prediction == "Mahalanobis"
  trunc_s   <- prediction == "suitability_trunc"
  trunc_m   <- prediction == "Mahalanobis_trunc"

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
    # Call the existing method for a single ellipsoid
    p <- predict.nicheR_ellipsoid(object = object$ellipse_community[[i]],
                                  newdata = newdata,
                                  include_suitability = inc_suit,
                                  suitability_truncated = trunc_s,
                                  include_mahalanobis = inc_mahal,
                                  mahalanobis_truncated = trunc_m,
                                  verbose = FALSE)

    if (verbose) utils::setTxtProgressBar(pb, i)

    # Extract the specific column/layer requested
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
    # Combine results into a data frame
    results_list <- do.call(cbind, results_list)
    colnames(results_list) <- ell_names
    
    # Return original data (with spatial cols if present) + new predictions
    return(cbind(as.data.frame(newdata), results_list))
  }
}