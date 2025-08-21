# Sample Occurrence Points from a Suitable Environment

get_sample_occurrences <- function(n_occ,
                                   env_bg, # accepts Spat.Rasters and data.frames
                                   niche,
                                   method = c("random", "center", "edge"),
                                   seed = NULL) {

  method <- tolower(match.arg(method))


  # EDIT: add the option to accept raster, but transform to SpatRaster. Also for
  # tibble

  # Make env_bg a data.frame if other
  if(class(env_bg) != "data.frame"){
    env_bg <- as.data.frame(env_bg, xy = TRUE, na.rm = TRUE)
    niche_vars <- names(env_bg)[which(!names(env_bg) %in% c("x", "y"))]

  }else{

    if (ncol(env_bg) != niche$dimen + 2) {
      stop("'evn_bg' must include spatial columns, x and y are necessary. Update data.frame or use Spat.Raster")
    }

    niche_vars <- names(env_bg)[which(!names(env_bg) %in% c("x", "y"))]
  }

  # Calculate mahalanobis distance
  suitable_pool <- get_suitable_environment(niche,
                                            env_bg[ , c(niche_vars)],
                                            out = "data.frame",
                                            distances = TRUE)

  # If there is no suitable area send message
  if (nrow(suitable_pool) < 1)
    stop("No suitable environments were found in the provided niche space.")

  # If suitable are exists, continue

  # --- 3. Sampling weights by method ---
  # d in [0, 1] (distance from center to boundary in Mahalanobis metric)
  # Making sure that there are numbers not larger than 1 or less than 0?
  d <- sqrt(pmax(0, suitable_pool$dist_sq))
  d <- pmin(d, 1)

  w <- switch(
    method,
    "random" = rep(1, nrow(suitable_pool)),
    "center" = 1 - d,   # higher near center
    "edge"   = d        # higher near boundary
  )

  w[!is.finite(w)] <- 0

  if (sum(w) == 0) {
    warning("All weights are zero; falling back to uniform sampling.",
            call. = FALSE)
    w <- rep(1, length(w))
  }

  # --- 4. Sample indices (ensure we return exactly n rows) ---
  if (!is.null(seed)) set.seed(seed)
  replace_flag <- n_occ > nrow(suitable_pool)

  idx <- sample.int(n = nrow(suitable_pool),
                    size = n_occ,
                    replace = replace_flag,
                    prob = w)

  # --- 5. Return sampled rows (drop temp column) ---
  occ <- suitable_pool[idx, , drop = FALSE]

  # Remove distances
  occ$dist_sq <- NULL
  rownames(occ) <- NULL

  occ_xy <- left_join(occ, env_bg, by = niche_vars)
  # table(duplicated(suitable_pool[,1:3]))

  message("Done sampling occurrences")
  return(occ_xy)
}
