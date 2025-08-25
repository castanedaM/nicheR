#' Validate arguments for plot_e_space
#'
#' Performs input checks for \code{plot_e_space} and emits informative errors
#' or warnings. Ensures required columns exist and are numeric, labels are correctly
#' specified, sampling size is valid, and optional objects have the needed structure.
#' This function is designed for internal use within the `NicheR` package
#' to ensure robust input handling for plotting functions.
#'
#' @param env_bg A data.frame with at least 3 numeric predictor columns.
#' @param x,y,z Column names or integer indices in \code{env_bg} for the three predictors,
#'  in the same order used to define the niche/ellipsoid.
#' @param labels Character vector of length 3 for axis labels.
#' @param n_bg Positive number of background rows to plot at most.
#' @param niche Optional list with elements \code{center} (numeric length 3)
#'  and \code{axes} (numeric length 3). If `niche` is provided, it must also
#'  contain an `angles` element (numeric length 3) for 2D projections to work
#'  correctly in `plot_e_space`.
#' @param show.pts.in Logical. If TRUE, expects \code{niche} to be provided.
#' @param occ_pts Optional data.frame with the same predictor columns present.
#' @param show.occ.density Logical. If TRUE, expects \code{occ_pts}.
#' @param plot.3d Logical. Indicates whether the main plot will be 3D. While
#'   this validator doesn't change behavior based on `plot.3d`, it's included
#'   for completeness of the arguments being validated from `plot_e_space`.
#'
#' @return An invisible list with \code{col_names}, the resolved predictor names
#'   in the order `c(x, y, z)`.
#' @keywords internal
#' @seealso \code{\link{plot_e_space}}
#' @export
validate_plot_e_space_args <- function(env_bg, x, y, z,
                                       labels, n_bg,
                                       niche, show.pts.in,
                                       occ_pts, show.occ.density,
                                       plot.3d) { # Added plot.3d to the validated arguments
  # --- env_bg checks ---
  if (!is.data.frame(env_bg)) stop("'env_bg' must be a data.frame.")
  if (ncol(env_bg) < 3) stop("'env_bg' must have at least 3 columns.")
  if (nrow(env_bg) < 1) stop("'env_bg' has no rows to plot.")

  # Helper: resolve a single spec (name or index) to a valid column name
  resolve_one <- function(spec) {
    if (is.character(spec) && length(spec) == 1) {
      if (!spec %in% names(env_bg)) stop(sprintf("Column '%s' not found in 'env_bg'.", spec))
      spec
    } else if (is.numeric(spec) && length(spec) == 1) {
      idx <- as.integer(spec)
      if (is.na(idx) || idx < 1 || idx > ncol(env_bg)) stop("'x', 'y', 'z' index out of range for 'env_bg'.")
      names(env_bg)[idx]
    } else {
      stop("'x', 'y', and 'z' must be column names or single integer indices.")
    }
  }

  col_x <- resolve_one(x)
  col_y <- resolve_one(y)
  col_z <- resolve_one(z)
  col_names <- c(col_x, col_y, col_z)

  # Distinctness
  if (any(duplicated(col_names))) {
    stop("'x', 'y', and 'z' must refer to three distinct columns in 'env_bg'.")
  }

  # All numeric
  non_numeric <- col_names[!vapply(env_bg[, col_names, drop = FALSE], is.numeric, logical(1))]
  if (length(non_numeric)) {
    stop(sprintf("These columns must be numeric in 'env_bg': %s", paste(non_numeric, collapse = ", ")))
  }

  # --- labels checks ---
  if (!(is.character(labels) && length(labels) == 3 && all(!is.na(labels)))) {
    stop("'labels' must be a non-NA character vector of length 3.")
  }

  # --- n_bg checks ---
  if (!(length(n_bg) == 1 && is.finite(n_bg) && n_bg > 0)) {
    stop("'n_bg' must be a single positive number.")
  }
  if (n_bg > 100000) {
    warning("'n_bg' is very large. Plotting may be slow.", call. = FALSE, immediate. = TRUE)
  } else if (n_bg > 10000) {
    warning("Selected number of background points is large and may slow plotting.",
            call. = FALSE, immediate. = TRUE)
  }
  if (is.numeric(n_bg) && n_bg %% 1 != 0) {
    warning("'n_bg' is not an integer; it will be truncated when sampling.", call. = FALSE, immediate. = TRUE)
  }
  if (nrow(env_bg) < n_bg) {
    warning(sprintf("'n_bg' (%d) exceeds available rows in 'env_bg' (%d). All rows will be used.",
                    as.integer(n_bg), nrow(env_bg)), call. = FALSE, immediate. = TRUE)
  }

  # --- niche object checks ---
  if (!is.null(niche)) {
    need <- c("center", "axes", "angles") # angles is explicitly required for 2D projections
    miss <- setdiff(need, names(niche))
    if (length(miss)) {
      stop(sprintf("'niche' is missing required fields: %s", paste(miss, collapse = ", ")))
    }
    if (!(is.numeric(niche$center) && length(niche$center) == 3)) {
      stop("'niche$center' must be numeric of length 3.")
    }
    if (!(is.numeric(niche$axes) && length(niche$axes) == 3)) {
      stop("'niche$axes' must be numeric of length 3.")
    }
    if (any(!is.finite(niche$axes)) || any(niche$axes <= 0)) {
      stop("'niche$axes' must contain finite, positive values.")
    }
    if (!(is.numeric(niche$angles) && length(niche$angles) == 3)) {
      stop("'niche$angles' must be numeric of length 3.")
    }
    if (any(!is.finite(niche$angles))) {
      stop("'niche$angles' must contain finite values.")
    }
  }

  # --- show.pts.in requires niche ---
  if (isTRUE(show.pts.in) && is.null(niche)) {
    warning("'show.pts.in' is TRUE but 'niche' is NULL. Points inside cannot be computed.",
            call. = FALSE, immediate. = TRUE)
  }

  # --- occ_pts and density flags checks ---
  if (!is.null(occ_pts)) {
    if (!is.data.frame(occ_pts)) stop("'occ_pts' must be a data.frame.")
    missing_occ <- setdiff(col_names, names(occ_pts))
    if (length(missing_occ)) {
      stop(sprintf("'occ_pts' is missing required columns: %s", paste(missing_occ, collapse = ", ")))
    }
    non_num_occ <- col_names[!vapply(occ_pts[, col_names, drop = FALSE], is.numeric, logical(1))]
    if (length(non_num_occ)) {
      stop(sprintf("These 'occ_pts' columns must be numeric: %s", paste(non_num_occ, collapse = ", ")))
    }
    if (nrow(occ_pts) < 1) warning("'occ_pts' has zero rows; nothing to plot.", call. = FALSE, immediate. = TRUE)
  } else if (isTRUE(show.occ.density)) {
    warning("'show.occ.density' is TRUE but 'occ_pts' is NULL. Density panels will be skipped.",
            call. = FALSE, immediate. = TRUE)
  }

  invisible(list(col_names = col_names))
}
