#' Extract Suitable Environmental Area from a Niche Ellipsoid
#'
#' This function identifies and extracts environmental grid cells or data
#' points that fall within a defined ellipsoid niche based on Mahalanobis distance.
#'
#' @param niche An object of class \code{ellipsoid} created by
#'   \code{build_ellipsoid()}.
#' @param env_bg A \code{terra::SpatRaster}, \code{data.frame}, or \code{matrix}
#'   of environmental predictor variables. It must contain the variables
#'   referenced by the \code{niche} object. If a \code{data.frame} is used and
#'   spatial output is requested, it should also contain \code{x} and \code{y}
#'   columns for spatial referencing.
#' @param out A character string specifying the desired output format. One of
#'   \code{"data.frame"}, \code{"spatial"}, or \code{"both"}:
#'   \itemize{
#'     \item \code{"data.frame"}: Returns a data frame of all suitable points.
#'     \item \code{"spatial"}: Returns a \code{terra::SpatRaster} where suitable
#'       cells are marked with 1 and unsuitable with 0. Requires \code{env_bg}
#'       to be a \code{terra::SpatRaster}.
#'     \item \code{"both"}: Returns a list containing both the spatial raster and
#'       the data frame of suitable points.
#'   }
#' @param distances Logical; if \code{TRUE}, an additional column named
#'   \code{dist_sq} is added to the output data frame containing the squared
#'   Mahalanobis distance for each suitable point.
#' @param verbose Logical; if \code{TRUE}, prints basic progress messages.
#'
#' @return The suitable environmental area in the specified format: a
#'   \code{data.frame}, a \code{terra::SpatRaster}, or a list containing both.
#'
#' @details
#' The function converts \code{env_bg} (if needed) to a data frame, calculates
#' the squared Mahalanobis distance for each point, and then filters for points
#' where this distance is less than or equal to 1. This identifies all locations
#' that fall within the ellipsoid niche boundary.
#'
#' @seealso \code{\link{build_ellipsoid}}, \code{\link{get_sample_occ}}
#'
#' @export
get_suitable_env <- function(niche,
                             env_bg,
                             out = c("data.frame", "spatial", "both"),
                             distances = FALSE,
                             verbose = TRUE) {

  gc()
  out <- tolower(match.arg(out))

  # --- 1) Check niche structure ----------------------------------------------

  if (!inherits(niche, "ellipsoid")) {
    stop("'niche' must be an object of class 'ellipsoid' produced by build_ellipsoid().")
  }

  need <- c("center", "Sigma_inv", "dimen")
  miss <- setdiff(need, names(niche))
  if (length(miss)) {
    stop(sprintf("'niche' is missing required fields: %s", paste(miss, collapse = ", ")))
  }

  if (!(is.numeric(niche$center) && length(niche$center) == niche$dimen)) {
    stop("'niche$center' must be numeric with length equal to 'niche$dimen'.")
  }

  if (!is.matrix(niche$Sigma_inv) || any(!is.finite(niche$Sigma_inv))) {
    stop("'niche$Sigma_inv' must be a finite numeric matrix.")
  }

  if (!is.logical(distances) || length(distances) != 1) {
    stop("'distances' must be a single logical value.")
  }

  if (missing(env_bg) || is.null(env_bg)) {
    stop("'env_bg' is required and cannot be NULL.")
  }

  # --- 2) Coerce env_bg and build env_bg_df ----------------------------------

  # tibble -> data.frame
  if (inherits(env_bg, "tbl_df")) {
    env_bg <- as.data.frame(env_bg)
  }

  # raster::Raster* -> SpatRaster
  if (inherits(env_bg, "Raster")) {
    env_bg <- terra::rast(env_bg)
  }

  # Track raster version for spatial output
  env_bg_rast <- NULL
  env_is_raster <- inherits(env_bg, "SpatRaster")

  if (!inherits(env_bg, c("SpatRaster", "data.frame", "matrix"))) {
    stop("'env_bg' must be a terra::SpatRaster, data.frame, or matrix.")
  }

  if (env_is_raster) {
    env_bg_rast <- env_bg
    # Use size-aware helper for large rasters
    env_bg_df <- as.data.frame.nicheR(env_bg_rast, verbose = verbose, use_cache = TRUE)
  } else {
    env_bg_df <- as.data.frame(env_bg)
  }

  # --- 3) Determine predictor columns ----------------------------------------

  # If x,y exist, predictors = all other columns; otherwise use all columns
  if (all(c("x", "y") %in% names(env_bg_df))) {
    niche_vars <- setdiff(names(env_bg_df), c("x", "y"))
  } else {
    niche_vars <- names(env_bg_df)
  }

  if (length(niche_vars) < niche$dimen) {
    stop("The number of predictor columns in 'env_bg' is less than 'niche$dimen'.")
  }

  # --- 4) Compute Mahalanobis distance and filter inside ---------------------

  cc <- stats::complete.cases(env_bg_df[, niche_vars, drop = FALSE])
  if (!any(cc)) {
    stop("All candidate rows contain NA in predictor columns. Provide complete predictors or impute values.")
  }

  env_bg_df_cc <- env_bg_df[cc, , drop = FALSE]

  pts   <- as.matrix(env_bg_df_cc[, niche_vars, drop = FALSE])
  diffs <- sweep(pts, 2, niche$center, "-")
  m_sq  <- rowSums((diffs %*% niche$Sigma_inv) * diffs)

  is_inside_cc <- is.finite(m_sq) & (m_sq <= 1)

  if (!any(is_inside_cc)) {
    stop("No points fall inside the ellipsoid after removing rows with NA predictors.")
  }

  inside_rows <- which(cc)[is_inside_cc]

  return_df <- env_bg_df[inside_rows, , drop = FALSE]
  if (isTRUE(distances)) {
    return_df$dist_sq <- m_sq[is_inside_cc]
  }

  # --- 5) Spatial output (binary raster) -------------------------------------

  return_ras <- NULL

  if (out %in% c("spatial", "both")) {

    if (!env_is_raster) {
      stop("For 'spatial' or 'both' output, 'env_bg' must be a terra::SpatRaster.")
    }

    if (!all(c("x", "y") %in% names(return_df))) {
      stop("For spatial output, 'env_bg' or its data.frame representation must contain 'x' and 'y' columns.")
    }

    # Start from first layer, keep geometry
    return_ras <- env_bg_rast[[1]]
    vals <- terra::values(return_ras)
    vals[!is.na(vals)] <- 0
    terra::values(return_ras) <- vals

    xy <- stats::na.omit(return_df[, c("x", "y"), drop = FALSE])
    if (nrow(xy) == 0) {
      stop("No valid XY coordinates found for suitable environments.")
    }

    inside_cells <- terra::cellFromXY(env_bg_rast, as.matrix(xy))
    inside_cells <- inside_cells[!is.na(inside_cells)]

    if (length(inside_cells)) {
      return_ras[as.integer(inside_cells)] <- 1
    }

    names(return_ras) <- "suitable"
  }

  # --- 6) Pack and return -----------------------------------------------------

  res <- switch(
    out,
    "data.frame" = return_df,
    "spatial"    = return_ras,
    "both"       = {
      tmp <- list(
        suitable_env_sp = return_ras,
        suitable_env_df = return_df
      )
      class(tmp) <- c("suitable_env", class(tmp))
      tmp
    }
  )

  gc()
  return(res)
}

#' @export
print.suitable_env <- function(x, ...) {
  if (is.list(x) && all(c("suitable_env_sp", "suitable_env_df") %in% names(x))) {
    cat("Suitable environment object:\n")
    cat("Spatial layer (SpatRaster):\n")
    print(x$suitable_env_sp)
    cat("\nData frame (showing first 6 rows):\n")
    print(utils::head(x$suitable_env_df))
  } else if (inherits(x, "SpatRaster")) {
    cat("Suitable environment raster:\n")
    print(x)
  } else if (is.data.frame(x)) {
    cat("Suitable environment data frame (showing first 6 rows):\n")
    print(utils::head(x))
  } else {
    NextMethod()
  }
  invisible(x)
}
