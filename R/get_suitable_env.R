# Extract from the enviroemntal raster layer the suitble area based  on the
# ellipsoid

get_suitable_environment <- function(niche,
                                     env_bg,
                                     out = c("data.frame", "spatial", "both"),
                                     distances = FALSE) {
  out <- tolower(match.arg(out))

  # EDIT: add the option to accept raster, but transform to SpatRaster. Also for
  # tibble

  # --- 1. Input Validation ---
  if(!class(env_bg) %in% c("data.frame", "SpatRaster", "matrix")){
    stop("'env_bg' must be either a matrix, data.frame or a terra::SpatRaster.")
  }

  if(out %in% c("spatial", "both")){
    if(class(env_bg) != "SpatRaster"){
      stop("If desire output contains a spatial layer 'env_bg' must be a terra::SpatRaster")
    }
    env_bg_df <- as.data.frame(env_bg, xy = TRUE)
  } else {
    if(!is.matrix(env_bg) && !is.data.frame(env_bg)) {
      stop("'env_bg' must be a matrix or data frame.")
    }
    env_bg_df <- env_bg

    if (ncol(env_bg_df) != niche$dimen) {
      stop("The number of columns in 'env_bg' does not match the ellipsoid's dimensions. Specify columns, coordinate coolumns x and y are not necessary.")
    }
  }

  if (!inherits(niche, "ellipsoid")) {
    stop("'niche' must be an object of class 'ellipsoid' from build_ellipsoid().")
  }


  # --- 2. Calculate Mahalanobis Distance and Filter ---
  pts <- as.matrix(env_bg_df[, names(env_bg)])

  diffs <- sweep(pts, 2, niche$center, "-")
  m_sq_dist <- rowSums((diffs %*% niche$Sigma_inv) * diffs)

  is_inside <- m_sq_dist <= 1

  return_df <- as.data.frame(pts[is_inside, , drop = FALSE])


  if(isTRUE(distances)){
    return_df$dist_sq <- m_sq_dist[m_sq_dist <= 1]
  }

  if(out %in% c("spatial", "both")){

    return_ras <- env_bg[[1]]   # keep geometry
    return_ras[!is.na(values(return_ras))] <- 0  # set all cells to 0

    xy_match <- dplyr::left_join(return_df, env_bg_df,
                                 by = names(env_bg),
                                 relationship = "many-to-many")

    # mark suitable cells
    # return_ras[as.numeric(rownames(return_df))] <- 1
    inside_cells <- terra::cellFromXY(env_bg, xy_match[, c("x","y")])
    return_ras[as.numeric(inside_cells)] <- 1

    names(return_ras) <- "suitable"

  }


  if(out == "data.frame"){
    return(return_df)
  }else if(out == "spatial"){
    return(return_ras)
  }else{
    return(list(suitable_env_sp = return_ras,
                suitable_env_df = return_df))

    # EDIT: to only print the top 6 of the df, create print helper function
  }
}
