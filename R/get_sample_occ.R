#' Sample occurrence points from a suitable environment pool
#'
#' Sample occurrence points from a precomputed pool of suitable environments.
#' Sampling can be uniform, biased toward the niche center, or biased toward
#' the niche edge. Weights can be computed from Mahalanobis distance
#' (\code{dist_sq}) or from a multivariate normal (MVN) density.
#'
#' @param n_occ Integer; number of occurrence points to sample.
#' @param suitable_env A pool of suitable environments. Can be:
#' \itemize{
#'   \item a \code{terra::SpatRaster} or \code{raster::Raster*} (typically a
#'         distance raster such as \code{dist_sq}),
#'   \item a \code{data.frame} or \code{matrix} containing \code{x} and \code{y}
#'         (and for non-random methods, \code{dist_sq}),
#'   \item a \code{suitable_env} object returned by
#'         \code{get_suitable_env(..., out.suit = "both")}.
#' }
#' @param method Character; one of \code{"random"}, \code{"center"}, or
#'   \code{"edge"}.
#' \itemize{
#'   \item \code{"random"} samples uniformly from the pool.
#'   \item \code{"center"} gives higher weight to points nearer the niche center.
#'   \item \code{"edge"} gives higher weight to points nearer the niche boundary.
#' }
#' @param sampling Character; how to compute weights. One of:
#' \itemize{
#'   \item \code{"mahalanobis"}: use \code{dist_sq} as the weight base.
#'   \item \code{"mvn"}: use MVN density via \code{mvtnorm::dmvnorm}. Requires
#'         \code{niche}.
#' }
#' @param niche Optional \code{ellipsoid} object. Required if
#'   \code{sampling = "mvn"}. Must include \code{center} and \code{Sigma}.
#' @param bias_surface Optional 1-layer \code{terra::SpatRaster} (or
#'   \code{raster::Raster*}) scaled to [0, 1], used to modulate sampling
#'   probability. Bias is extracted at \code{x,y} and multiplied into the
#'   sampling weights. \code{NA} bias values are treated as 0.
#' @param seed Optional integer seed for reproducible sampling.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'
#' @return A \code{data.frame} with \code{n_occ} sampled rows from the suitable
#'   pool. The returned object includes \code{x}, \code{y}, and any predictor
#'   columns. The \code{dist_sq} column is removed before returning.
#'
#' @details
#' This function assumes suitability has already been computed (e.g. with
#' \code{get_suitable_env()}). It only performs the sampling step:
#' \enumerate{
#'   \item Convert \code{suitable_env} to a data.frame with \code{x}, \code{y},
#'         and optionally \code{dist_sq}.
#'   \item Compute weights based on \code{sampling} and \code{method}.
#'   \item Optionally multiply weights by \code{bias_surface}.
#'   \item Sample \code{n_occ} rows with probability proportional to weights.
#' }
#'
#' @seealso \code{\link{get_suitable_env}}, \code{\link{set_bias_surface}}
#'
#' @export
get_sample_occ <- function(n_occ,
                           suitable_env,
                           method = c("random", "center", "edge"),
                           sampling = c("mahalanobis", "mvn"),
                           niche = NULL,
                           bias_surface = NULL,
                           seed = NULL,
                           verbose = TRUE) {

  gc()
  method   <- tolower(match.arg(method))
  sampling <- tolower(match.arg(sampling))

  if (isTRUE(verbose)) message("Starting get_sample_occ()...")

  # --- 0) Basic checks -------------------------------------------------------

  if (missing(n_occ) || length(n_occ) != 1 || !is.numeric(n_occ)) {
    stop("'n_occ' must be a single numeric value.")
  }
  n_occ <- as.integer(n_occ)
  if (!is.finite(n_occ) || n_occ <= 0) stop("'n_occ' must be a positive integer.")

  if (is.null(seed) && isTRUE(verbose)) {
    warning("Sampling seed not set; results will differ each time this function is run.",
            call. = FALSE, immediate. = TRUE)
  }

  if (missing(suitable_env) || is.null(suitable_env)) {
    stop("'suitable_env' must be provided.")
  }

  # --- 1) Coerce suitable_env to data.frame ---------------------------------

  suitable_pool <- NULL
  suitable_rast <- NULL

  if (isTRUE(verbose)) {
    message("Coercing 'suitable_env' to a data.frame of candidate environments...")
  }

  # Accept suitable_env outputs from get_suitable_env(..., out.suit="both")
  if (inherits(suitable_env, "suitable_env") ||
      (is.list(suitable_env) &&
       any(c("suitable_env_df", "suitable_env_sp") %in% names(suitable_env)))) {

    if ("suitable_env_df" %in% names(suitable_env) &&
        is.data.frame(suitable_env$suitable_env_df)) {

      suitable_pool <- suitable_env$suitable_env_df

    } else if ("suitable_env_sp" %in% names(suitable_env)) {

      sp <- suitable_env$suitable_env_sp
      if (inherits(sp, "Raster")) sp <- terra::rast(sp)

      # sp can be a SpatRaster or a list of SpatRaster objects
      if (inherits(sp, "SpatRaster")) {
        suitable_rast <- sp
      } else if (is.list(sp) &&
                 length(sp) > 0 &&
                 all(vapply(sp, inherits, logical(1), "SpatRaster"))) {

        # For center/edge we need dist_sq; for random any layer is ok
        if (method == "random") {
          suitable_rast <- if ("suitable" %in% names(sp)) sp[["suitable"]] else sp[[1]]
        } else {
          if ("dist_sq" %in% names(sp)) {
            suitable_rast <- sp[["dist_sq"]]
          } else {
            stop("For method = '", method,
                 "', 'suitable_env' does not contain a 'dist_sq' raster.\n",
                 "Run get_suitable_env(..., distances = TRUE) or pass a data.frame with 'dist_sq'.")
          }
        }

      } else {
        stop("'suitable_env$suitable_env_sp' must be a SpatRaster or a list of SpatRaster objects.")
      }

      suitable_pool <- as.data.frame.nicheR(suitable_rast)

    } else {
      stop("'suitable_env' object must contain either 'suitable_env_df' or 'suitable_env_sp'.")
    }

  } else {

    if (inherits(suitable_env, "Raster")) suitable_env <- terra::rast(suitable_env)

    if (inherits(suitable_env, "SpatRaster")) {
      suitable_rast <- suitable_env
      suitable_pool <- as.data.frame.nicheR(suitable_rast)
    } else if (is.data.frame(suitable_env) || is.matrix(suitable_env)) {
      suitable_pool <- as.data.frame(suitable_env)
    } else {
      stop("'suitable_env' must be a 'suitable_env' object, a terra::SpatRaster, raster::Raster*, data.frame, or matrix.")
    }
  }

  # Need x, y always
  required_xy <- c("x", "y")
  missing_xy  <- setdiff(required_xy, names(suitable_pool))
  if (length(missing_xy) > 0) {
    stop("'suitable_env' (after coercion) must contain columns: ",
         paste(required_xy, collapse = ", "),
         ". Missing: ", paste(missing_xy, collapse = ", "), ".")
  }

  if (nrow(suitable_pool) < 1) stop("No suitable environments were found in 'suitable_env'.")

  # For method center/edge, need dist_sq for mahalanobis sampling
  if (sampling == "mahalanobis" && method %in% c("center", "edge") && !"dist_sq" %in% names(suitable_pool)) {
    stop("For sampling = 'mahalanobis' and method = '", method,
         "', 'suitable_env' must include a 'dist_sq' column.\n",
         "Run get_suitable_env(..., distances = TRUE).")
  }

  if (isTRUE(verbose)) {
    message("Suitable pool contains ", nrow(suitable_pool), " candidate environments.")
  }

  # --- 2) Optional bias surface ---------------------------------------------

  bias_values <- NULL

  if (!is.null(bias_surface)) {

    if (inherits(bias_surface, "Raster")) bias_surface <- terra::rast(bias_surface)

    if (!inherits(bias_surface, "SpatRaster")) {
      stop("'bias_surface' must be a 1-layer terra::SpatRaster or raster::Raster*.")
    }
    if (terra::nlyr(bias_surface) != 1) {
      stop("'bias_surface' must have exactly one layer (0–1 bias values).")
    }

    if (isTRUE(verbose)) message("Extracting bias values at suitable locations...")

    vals <- terra::values(bias_surface)
    rng  <- range(vals, na.rm = TRUE)
    if (rng[1] < 0 - 1e-6 || rng[2] > 1 + 1e-6) {
      stop("'bias_surface' must be scaled to [0, 1]. Observed range: [",
           signif(rng[1], 3), ", ", signif(rng[2], 3), "].")
    }

    coords <- cbind(suitable_pool$x, suitable_pool$y)
    ext_df <- terra::extract(bias_surface, coords)
    bias_values <- ext_df[, 1]
    bias_values[is.na(bias_values)] <- 0

  } else if (isTRUE(verbose)) {
    message("No 'bias_surface' supplied; sampling based only on 'method'.")
  }

  # --- 3) Compute weights ----------------------------------------------------

  if (isTRUE(verbose)) {
    message("Computing sampling weights with sampling = '", sampling,
            "' and method = '", method, "'...")
  }

  n <- nrow(suitable_pool)
  w <- rep(1, n)

  if (method == "random") {
    w <- rep(1, n)

  } else if (sampling == "mahalanobis") {

    # dist_sq is Mahalanobis squared distance; smaller = center, larger = edge
    r <- suitable_pool$dist_sq
    r[!is.finite(r)] <- NA_real_

    r_max <- max(r, na.rm = TRUE)
    if (!is.finite(r_max) || r_max <= 0) {
      stop("Invalid 'dist_sq' values: cannot compute weights.")
    }

    r01 <- r / r_max
    r01 <- pmax(0, pmin(r01, 1))

    w <- switch(
      method,
      "center" = 1 - r01,  # highest at center (r01 ~ 0)
      "edge"   = r01       # highest at edge   (r01 ~ 1)
    )

    w[!is.finite(w)] <- 0

  } else if (sampling == "mvn") {

    # Require niche
    if (is.null(niche) || !inherits(niche, "ellipsoid")) {
      stop("For sampling = 'mvn', please provide 'niche' as an 'ellipsoid' object.")
    }
    need_n <- c("center", "Sigma", "dimen")
    miss_n <- setdiff(need_n, names(niche))
    if (length(miss_n)) {
      stop("Provided 'niche' is missing required fields: ", paste(miss_n, collapse = ", "))
    }
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      stop("Package 'mvtnorm' is required for sampling = 'mvn'. Please install it.")
    }

    # Need predictor columns to evaluate MVN pdf
    exclude <- c("x", "y", "dist_sq")
    pred_cols <- setdiff(names(suitable_pool), exclude)

    if (length(pred_cols) < niche$dimen) {
      stop("For sampling = 'mvn', 'suitable_env' must include the predictor columns used by the niche.")
    }

    # Prefer matching by names if possible
    if (!is.null(names(niche$center)) && all(names(niche$center) %in% pred_cols)) {
      use_cols <- names(niche$center)
    } else {
      use_cols <- pred_cols[seq_len(niche$dimen)]
    }

    X <- as.matrix(suitable_pool[, use_cols, drop = FALSE])
    if (!is.numeric(X)) storage.mode(X) <- "double"
    if (any(!is.finite(X))) stop("Predictor columns contain non-finite values; remove/clean before sampling.")

    pdf <- mvtnorm::dmvnorm(X, mean = niche$center, sigma = niche$Sigma)
    pdf[!is.finite(pdf)] <- 0

    if (method == "center") {
      # higher probability near center
      w <- pdf
    } else if (method == "edge") {
      # emphasize edge: inverse of pdf (stable with epsilon)
      eps <- max(.Machine$double.eps, max(pdf, na.rm = TRUE) * 1e-12)
      w <- 1 / (pdf + eps)
    }

    w[!is.finite(w)] <- 0
  }

  # Apply bias
  if (!is.null(bias_values)) {
    if (isTRUE(verbose)) message("Applying bias surface to sampling weights...")
    w <- w * bias_values
  }

  if (sum(w, na.rm = TRUE) == 0) {
    warning("All sampling weights are zero; falling back to uniform sampling.", call. = FALSE)
    w <- rep(1, n)
  }

  # --- 4) Draw samples -------------------------------------------------------

  if (!is.null(seed)) {
    if (isTRUE(verbose)) message("Setting seed to ", seed, " and drawing ", n_occ, " samples...")
    set.seed(seed)
  } else if (isTRUE(verbose)) {
    message("Drawing ", n_occ, " samples without a fixed seed...")
  }

  replace_flag <- n_occ > nrow(suitable_pool)

  idx <- sample.int(
    n       = nrow(suitable_pool),
    size    = n_occ,
    replace = replace_flag,
    prob    = w
  )

  occ <- suitable_pool[idx, , drop = FALSE]
  if ("dist_sq" %in% names(occ)) occ$dist_sq <- NULL
  rownames(occ) <- NULL

  if (isTRUE(verbose)) message("Finished get_sample_occ(): sampled ", n_occ, " occurrences.")
  gc()
  return(occ)
}

