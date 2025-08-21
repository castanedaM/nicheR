# Extract from the enviroemntal raster layer the suitble area based  on the
# ellipsoid

get_suitable_env <- function(niche,
                             env_bg,
                             out = c("data.frame", "spatial", "both"),
                             distances = FALSE) {

  out <- tolower(match.arg(out))

  # --- 1) Input validation and coercion ---

  # 1a) Check niche structure early
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

  # 1b) Accept tibble -> data.frame (keeps names and types)
  if (inherits(env_bg, "tbl_df")) {
    env_bg <- as.data.frame(env_bg)
  }

  # 1c) Accept raster::Raster* by converting to terra::SpatRaster
  if (inherits(env_bg, "Raster")) {
    env_bg <- terra::rast(env_bg)
  }

  # 1d) Validate env_bg type now
  if (!inherits(env_bg, c("SpatRaster", "data.frame", "matrix"))) {
    stop("'env_bg' must be a terra::SpatRaster, data.frame, or matrix.")
  }

  # 1e) Branch on desired output type
  if (out %in% c("spatial", "both")) {
    # For spatial outputs we REQUIRE a SpatRaster so geometry is preserved
    if (!inherits(env_bg, "SpatRaster")) {
      stop("For spatial output, 'env_bg' must be a terra::SpatRaster.")
    }

    # Build a data.frame with XY and layer values for lookups later
    env_bg_df <- terra::as.data.frame(env_bg, xy = TRUE, na.rm = FALSE)

    # Ensure XY names exist
    if (!all(c("x", "y") %in% names(env_bg_df))) {
      stop("Could not find 'x' and 'y' columns after converting raster to data.frame.")
    }

    # Predictor columns will be all raster layers (exclude x,y)
    niche_vars <- setdiff(names(env_bg_df), c("x", "y"))

    if (length(niche_vars) < niche$dimen) {
      stop("Raster has fewer predictor layers than 'niche$dimen'.")
    }

  } else {
    # data.frame or matrix path (non-spatial outputs)
    if (!(is.matrix(env_bg) || is.data.frame(env_bg))) {
      stop("For 'data.frame' output, 'env_bg' must be a matrix or data.frame.")
    }

    env_bg_df <- as.data.frame(env_bg)
    niche_vars <- setdiff(names(env_bg_df), c("x", "y"))

    if (ncol(env_bg_df) < niche$dimen) {
      stop("The number of predictor columns in 'env_bg' is less than 'niche$dimen'.")
    }
  }

  # 1f) Optional distances flag type check
  if (!is.logical(distances) || length(distances) != 1) {
    stop("'distances' must be a single logical value.")
  }


  # --- 2. Calculate Mahalanobis Distance and Filter ---

  # Clean predictors first ---
  # niche_vars should be the predictor columns you selected earlier
  cc <- stats::complete.cases(env_bg_df[, niche_vars, drop = FALSE])

  if (!any(cc)) {
    stop("All candidate rows contain NA in predictor columns. Provide complete predictors or impute values.")
  }

  env_bg_df_cc <- env_bg_df[cc, , drop = FALSE]

  ## --- Distances on clean data ---
  pts <- as.matrix(env_bg_df_cc[, niche_vars, drop = FALSE])
  diffs <- sweep(pts, 2, niche$center, "-")
  m_sq_dist <- rowSums((diffs %*% niche$Sigma_inv) * diffs)

  # Treat any non-finite distance as outside
  is_inside_cc <- is.finite(m_sq_dist) & (m_sq_dist <= 1)

  if (!any(is_inside_cc)) {
    stop("No points fall inside the ellipsoid after removing rows with NA predictors.")
  }

  # Map back to original row indices if needed later (e.g., raster cell lookup)
  inside_rows <- which(cc)[is_inside_cc]

  # Data-frame return: keep original columns, not just predictors
  return_df <- env_bg_df[inside_rows, , drop = FALSE]

  if (isTRUE(distances)) {
    return_df$dist_sq <- m_sq_dist[is_inside_cc]
  }

  if (out %in% c("spatial", "both")) {

    return_ras <- env_bg[[1]]                    # keep geometry
    return_ras[!is.na(terra::values(return_ras))] <- 0

    # get XY for cells to set = 1
    if (all(c("x", "y") %in% names(return_df))) {
      # Already have coordinates, no join needed
      xy <- return_df[, c("x", "y"), drop = FALSE]

    } else {
      stop("Something went wrong with the spatial calcualtions")
    }

    # Clean XY and map to cells
    xy <- stats::na.omit(xy)

    if (nrow(xy) == 0) {
      stop("No valid XY coordinates found for suitable environments.")
    }

    inside_cells <- terra::cellFromXY(env_bg, as.matrix(xy))
    inside_cells <- inside_cells[!is.na(inside_cells)]

    if (length(inside_cells)) {
      return_ras[as.integer(inside_cells)] <- 1
    }

    names(return_ras) <- "suitable"
  }



  # --- 3) Return results ---
  res <- switch(out,
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

  return(res)

}

# ---- Pretty print helper ----
#' @export
print.suitable_env <- function(x, ...) {
  if (is.list(x) && all(c("suitable_env_sp", "suitable_env_df") %in% names(x))) {
    cat("Suitable environment object:\n")
    cat("  • Spatial layer (SpatRaster):\n")
    print(x$suitable_env_sp)
    cat("\n  • Data frame (showing first 6 rows):\n")
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
