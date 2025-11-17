#' Sample Occurrence Points from a Suitable Environment
#'
#' This function samples occurrence points from a background environmental
#' dataset based on a defined ellipsoid niche. Sampling can be biased towards
#' the center, the edge, or be purely random. Optionally, a sampling bias
#' surface can be supplied to modulate the probability of selecting each
#' suitable location.
#'
#' @param n_occ Integer; the number of occurrence points to sample.
#' @param env_bg Optional. A \code{terra::SpatRaster}, \code{data.frame}, or
#'   \code{matrix} of environmental predictor variables. It must contain the
#'   variables referenced by the \code{niche} object. If a \code{data.frame}
#'   is used, it should also contain \code{x} and \code{y} columns for spatial
#'   referencing. Required if \code{suitable_env} is not provided, or if
#'   \code{bias_surface} is used.
#' @param niche An object of class \code{ellipsoid} created by
#'   \code{build_ellipsoid()}.
#' @param method A character string specifying the sampling strategy. One of
#'   \code{"random"}, \code{"center"}, or \code{"edge"}:
#'   \itemize{
#'     \item \code{"random"}: Samples are drawn with uniform probability from
#'       the suitable environmental space.
#'     \item \code{"center"}: Sampling probability is inversely proportional to
#'       the Mahalanobis distance, favoring points closer to the niche center.
#'     \item \code{"edge"}: Sampling probability is proportional to the
#'       Mahalanobis distance, favoring points closer to the niche boundary.
#'   }
#' @param bias_surface Optional. One or more bias layers describing sampling
#'   bias. Passed to \code{set_bias_surface()} using either \code{env_bg} or
#'   \code{suitable_env} as the raster template. Can be:
#'   \itemize{
#'     \item a \code{terra::SpatRaster} (single or multi layer),
#'     \item a \code{terra::SpatVector},
#'     \item a character path to a raster or vector file,
#'     \item a \code{list} of any of the above.
#'   }
#'   The resulting pooled bias surface is used to modulate sampling weights
#'   over the suitable pool of environments. Only supported when at least one
#'   of \code{env_bg} or \code{suitable_env} is a raster object.
#' @param suitable_env Optional. A precomputed suitability object describing
#'   the suitable environment pool. May be:
#'   \itemize{
#'     \item a \code{terra::SpatRaster} or \code{raster::Raster*} already
#'           containing suitability and Mahalanobis distances, or
#'     \item a \code{data.frame} / \code{matrix} with columns \code{x},
#'           \code{y}, and \code{dist_sq}, typically produced by
#'           \code{get_suitable_env(..., out = "data.frame", distances = TRUE)}.
#'   }
#'   When supplied, \code{get_suitable_env()} is not called internally.
#' @param seed An integer used to set the random number generator seed for
#'   reproducible results.
#'
#' @return A \code{data.frame} containing \code{n_occ} sampled rows from the
#'   suitable environment pool, including the \code{x} and \code{y} coordinates
#'   and the environmental predictor values.
#'
#' @details
#' The function first obtains a pool of suitable environments either by:
#' \enumerate{
#'   \item calling \code{get_suitable_env()} on \code{env_bg}, or
#'   \item using a user-supplied \code{suitable_env} object that already
#'         contains suitability and Mahalanobis distances.
#' }
#'
#' It then applies a weighting scheme based on \code{method} and, if provided,
#' a pooled bias surface from \code{set_bias_surface()} to probabilistically
#' select \code{n_occ} points from this suitable environment pool.
#'
#' @seealso \code{\link{build_ellipsoid}}, \code{\link{get_suitable_env}},
#'   \code{\link{set_bias_surface}}
#'
#' @export
get_sample_occ <- function(n_occ,
                           env_bg = NULL,
                           niche,
                           method = c("random", "center", "edge"),
                           bias_surface = NULL,
                           suitable_env = NULL,
                           seed = NULL) {

  method <- tolower(match.arg(method))

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

  # --- 2) Ensure we have at least env_bg or suitable_env ---------------------

  if (is.null(env_bg) && is.null(suitable_env)) {
    stop("Either 'env_bg' or 'suitable_env' must be provided.")
  }

  # --- 3) Coerce env_bg if present -------------------------------------------

  env_bg_rast <- NULL
  if (!is.null(env_bg)) {

    if (inherits(env_bg, "tbl_df")) {
      env_bg <- as.data.frame(env_bg)
    }

    if (inherits(env_bg, "Raster")) {
      env_bg <- terra::rast(env_bg)
    }

    if (!inherits(env_bg, c("SpatRaster", "data.frame", "matrix"))) {
      stop("'env_bg' must be a terra::SpatRaster, data.frame, or matrix.")
    }

    if (inherits(env_bg, "SpatRaster")) {
      env_bg_rast <- env_bg

      if (terra::nlyr(env_bg_rast) < niche$dimen) {
        stop("Raster has fewer predictor layers than 'niche$dimen'.")
      }

    } else {
      env_bg <- as.data.frame(env_bg)
      if (!all(c("x", "y") %in% names(env_bg))) {
        stop("'env_bg' must include spatial columns 'x' and 'y' when supplied as a data.frame or matrix.")
      }
    }
  }

  # --- 4) Coerce suitable_env / build suitable_pool --------------------------

  suitable_env_rast <- NULL
  if (is.null(suitable_env)) {

    if (is.null(env_bg)) {
      stop("When 'suitable_env' is NULL, 'env_bg' must be provided.")
    }

    suitable_pool <- get_suitable_env(
      niche     = niche,
      env_bg    = env_bg,
      out       = "data.frame",
      distances = TRUE
    )

  } else {

    # keep raster version (if any) for potential bias template
    if (inherits(suitable_env, "Raster")) {
      suitable_env <- terra::rast(suitable_env)
    }
    if (inherits(suitable_env, "SpatRaster")) {
      suitable_env_rast <- suitable_env
      suitable_pool <- as.data.frame.nicheR(suitable_env_rast)
    } else if (is.data.frame(suitable_env) || is.matrix(suitable_env)) {
      suitable_pool <- as.data.frame(suitable_env)
    } else {
      stop(
        "'suitable_env' must be a terra::SpatRaster, raster::Raster*, ",
        "data.frame, or matrix when provided directly."
      )
    }

    required_cols <- c("x", "y", "dist_sq")
    missing_cols  <- setdiff(required_cols, names(suitable_pool))
    if (length(missing_cols) > 0) {
      stop(
        "'suitable_env' (after coercion) must contain columns: ",
        paste(required_cols, collapse = ", "),
        ". Missing: ", paste(missing_cols, collapse = ", "),
        ". If needed, run get_suitable_env(..., distances = TRUE) first."
      )
    }
  }

  # --- 5) Check suitable pool size -------------------------------------------

  if (nrow(suitable_pool) < 1) {
    stop("No suitable environments were found in the provided niche space.")
  }

  # --- 6) Optional sampling bias ---------------------------------------------

  bias_values <- NULL

  if (!is.null(bias_surface)) {

    # Choose template: prefer env_bg_rast, otherwise suitable_env_rast
    template_rast <- if (!is.null(env_bg_rast)) {
      env_bg_rast
    } else {
      suitable_env_rast
    }

    if (is.null(template_rast)) {
      stop(
        "A 'bias_surface' was supplied, but neither 'env_bg' nor 'suitable_env' ",
        "is a spatial raster object.\n",
        "Bias layers can only be applied when 'env_bg' or 'suitable_env' is a SpatRaster, ",
        "so they can be aligned and extracted at occurrence locations."
      )
    }

    bias_res <- set_bias_surface(
      bias_surface = bias_surface,
      template     = template_rast,
      out          = "default",
      verbose      = FALSE
    )

    pooled_bias_sp <- bias_res$pooled_bias_sp

    coords  <- cbind(suitable_pool$x, suitable_pool$y)
    ext_df  <- terra::extract(pooled_bias_sp, coords)
    bias_values <- ext_df[, 2]
    bias_values[is.na(bias_values)] <- 0
  }

  # --- 7) Sampling weights by method -----------------------------------------

  d <- sqrt(pmax(0, suitable_pool$dist_sq))
  d <- pmin(d, 1)

  w <- switch(
    method,
    "random" = rep(1, nrow(suitable_pool)),
    "center" = 1 - d,   # higher near center
    "edge"   = d        # higher near boundary
  )

  w[!is.finite(w)] <- 0

  if (!is.null(bias_values)) {
    w <- w * bias_values
  }

  if (sum(w) == 0) {
    warning("All sampling weights are zero; falling back to uniform sampling.",
            call. = FALSE)
    w <- rep(1, length(w))
  }

  # --- 8) Sample indices ------------------------------------------------------

  if (!is.null(seed)) set.seed(seed)
  replace_flag <- n_occ > nrow(suitable_pool)

  idx <- sample.int(
    n       = nrow(suitable_pool),
    size    = n_occ,
    replace = replace_flag,
    prob    = w
  )

  # --- 9) Return sampled rows -------------------------------------------------

  occ <- suitable_pool[idx, , drop = FALSE]
  occ$dist_sq <- NULL
  rownames(occ) <- NULL

  message(sprintf("Done sampling %d occurrences", n_occ))
  return(occ)
}
