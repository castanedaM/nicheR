#' Validate arguments for plot_e_space
#'
#' Performs input checks for \code{plot_e_space} and emits informative errors
#' or warnings. Ensures required columns exist and are numeric, labels are correctly
#' specified, sampling size is valid, and optional objects have the needed structure.
#' This function is designed for internal use within the \code{NicheR} package
#' to ensure robust input handling for plotting functions.
#'
#' @param env_bg A data.frame with at least 3 numeric predictor columns.
#' @param x,y,z Column names or integer indices in \code{env_bg} for the three predictors,
#'   in the same order used to define the niche/ellipsoid.
#' @param labels Character vector of length 3 for axis labels.
#' @param n_bg Positive number of background rows to plot at most.
#' @param niche Optional list/object with elements \code{center} (numeric length 3)
#'   and \code{axes} (numeric length 3). If \code{niche} is provided and
#'   2D projections will be drawn in \code{plot_e_space}, it is recommended that
#'   \code{niche} also contain an \code{angles} element (numeric length 3) so
#'   projected ellipses are correctly oriented.
#' @param occ_pts Optional data.frame with the same predictor columns present.
#' @param show.occ.density Logical. If \code{TRUE}, expects \code{occ_pts} to
#'   be provided; otherwise density panels are skipped with a warning.
#'
#' @return An invisible list with \code{col_names}, the resolved predictor names
#'   in the order \code{c(x, y, z)}.
#' @keywords internal
#' @seealso \code{\link{plot_e_space}}
#' @export
validate_plot_e_space_args <- function(env_bg, x, y, z,
                                       labels, n_bg,
                                       niche,
                                       occ_pts, show.occ.density) {

  if (!is.data.frame(env_bg)) {
    stop("'env_bg' must be a data.frame.")
  }

  # Helper: resolve a single spec (name or index) to a valid column name
  resolve_one <- function(spec) {
    if (is.character(spec) && length(spec) == 1) {
      if (!spec %in% names(env_bg)) {
        stop(sprintf("Column '%s' not found in 'env_bg'.", spec))
      }
      spec
    } else if (is.numeric(spec) && length(spec) == 1) {
      idx <- as.integer(spec)
      if (is.na(idx) || idx < 1 || idx > ncol(env_bg)) {
        stop("'x', 'y', and 'z' indices must be between 1 and ncol(env_bg).")
      }
      names(env_bg)[idx]
    } else {
      stop("'x', 'y', and 'z' must be column names or single integer indices.")
    }
  }

  col_x <- resolve_one(x)
  col_y <- resolve_one(y)
  col_z <- resolve_one(z)
  col_names <- c(col_x, col_y, col_z)

  # All numeric
  non_numeric <- col_names[!vapply(env_bg[, col_names, drop = FALSE], is.numeric, logical(1))]
  if (length(non_numeric)) {
    stop(
      sprintf(
        "These columns must be numeric in 'env_bg': %s",
        paste(non_numeric, collapse = ", ")
      )
    )
  }

  # labels
  if (!(is.character(labels) && length(labels) == 3)) {
    stop("'labels' must be a character vector of length 3.")
  }

  # n_bg
  if (!(length(n_bg) == 1 && is.numeric(n_bg) && is.finite(n_bg) && n_bg > 0)) {
    stop("'n_bg' must be a single positive number.")
  }
  if (n_bg > 100000) {
    warning("'n_bg' is very large; plotting may be slow.",
            call. = FALSE, immediate. = TRUE)
  } else if (n_bg > 10000) {
    warning(
      "Selected number of background points is relatively large; ",
      "plots may take longer to render.",
      call. = FALSE, immediate. = TRUE
    )
  }

  # niche object if provided
  if (!is.null(niche)) {
    need <- c("center", "axes")
    miss <- setdiff(need, names(niche))
    if (length(miss)) {
      stop(
        sprintf(
          "'niche' is missing required fields: %s",
          paste(miss, collapse = ", ")
        )
      )
    }

    if (!(is.numeric(niche$center) && length(niche$center) == 3)) {
      stop("'niche$center' must be numeric of length 3.")
    }
    if (!(is.numeric(niche$axes) && length(niche$axes) == 3)) {
      stop("'niche$axes' must be numeric of length 3.")
    }
    if (any(niche$axes <= 0)) {
      stop("'niche$axes' must contain positive values.")
    }

    # Angles are strongly recommended for 2D projections
    if (!"angles" %in% names(niche)) {
      warning(
        "'niche' does not contain an 'angles' element. ",
        "2D ellipsoid projections in plot_e_space() will assume no rotation.",
        call. = FALSE, immediate. = TRUE
      )
    } else if (!(is.numeric(niche$angles) && length(niche$angles) == 3)) {
      stop("'niche$angles' must be numeric of length 3 when supplied.")
    }
  }

  # occ_pts and density flags
  if (!is.null(occ_pts)) {
    if (!is.data.frame(occ_pts)) {
      stop("'occ_pts' must be a data.frame.")
    }

    missing_occ <- setdiff(col_names, names(occ_pts))
    if (length(missing_occ)) {
      stop(
        sprintf(
          "'occ_pts' is missing required columns: %s",
          paste(missing_occ, collapse = ", ")
        )
      )
    }

    non_num_occ <- col_names[!vapply(occ_pts[, col_names, drop = FALSE], is.numeric, logical(1))]
    if (length(non_num_occ)) {
      stop(
        sprintf(
          "These 'occ_pts' columns must be numeric: %s",
          paste(non_num_occ, collapse = ", ")
        )
      )
    }
  } else if (isTRUE(show.occ.density)) {
    warning(
      "'show.occ.density' is TRUE but 'occ_pts' is NULL. ",
      "Density panels will be skipped.",
      call. = FALSE, immediate. = TRUE
    )
  }

  invisible(list(col_names = col_names))
}
