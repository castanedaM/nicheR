#' Plot a nicheR ellipsoid in environmental space
#'
#' @description
#' Plots the 2D boundary of a \code{nicheR_ellipsoid} object in environmental
#' space for a chosen pair of dimensions. Optionally overlays background points
#' or a prediction surface colored by a continuous variable (e.g., suitability).
#' Use \code{\link{add_data}} and \code{\link{add_ellipsoid}} to layer additional
#' data onto the plot after calling this function.
#'
#' @param object A \code{nicheR_ellipsoid} object containing at least
#'   \code{centroid}, \code{cov_matrix}, \code{chi2_cutoff}, and
#'   \code{var_names}.
#' @param background Optional data frame or matrix of background points to plot
#'   behind the ellipsoid. Rows are observations, columns are environmental
#'   variables. If provided, \code{prediction} is ignored.
#' @param prediction Optional data frame or matrix of prediction values to plot.
#'   Used when \code{background} is \code{NULL}. Can be colored by a continuous
#'   variable using \code{col_layer}.
#' @param dim Integer vector of length 2. Indices of the two dimensions to plot.
#'   Default is \code{c(1, 2)}.
#' @param col_layer Character or \code{NULL}. Name of a column in
#'   \code{prediction} to use for coloring points by a continuous variable.
#'   If \code{NULL} (default), all prediction points are drawn with
#'   \code{col_bg}.
#' @param pal A color palette function or character vector used when
#'   \code{col_layer} is provided. Default is \code{hcl.colors(100, palette = "Viridis")}.
#' @param rev_pal Logical. If \code{TRUE}, reverses the color palette. Default
#'   is \code{FALSE}.
#' @param bg_sample Integer or \code{NULL}. If provided and the number of
#'   background or prediction rows exceeds this value, a random subsample of
#'   this size is drawn before plotting. Useful for large data frames. Default
#'   is \code{NULL} (plot all points).
#' @param lty Integer. Line type for the ellipsoid boundary. Default is
#'   \code{1} (solid).
#' @param lwd Numeric. Line width for the ellipsoid boundary. Default is
#'   \code{1}.
#' @param col_ell Character. Color of the ellipsoid boundary line. Default is
#'   \code{"#000000"} (black).
#' @param col_bg Character. Color of background or prediction points when
#'   \code{col_layer} is \code{NULL}. Default is \code{"#8A8A8A"} (grey).
#' @param pch Integer or character. Point symbol for background or prediction
#'   points. Default is \code{1}.
#' @param alpha_bg Numeric in \code{[0, 1]}. Transparency of background or
#'   prediction points. Default is \code{1} (fully opaque).
#' @param alpha_ell Numeric in \code{[0, 1]}. Transparency of the ellipsoid
#'   boundary line. Default is \code{1} (fully opaque).
#' @param cex_ell Numeric. Size scaling for the ellipsoid boundary line.
#'   Default is \code{1}.
#' @param cex_bg Numeric. Size scaling for background or prediction points.
#'   Default is \code{1}.
#' @param fixed_lims A named list with elements \code{xlim} and \code{ylim},
#'   each a numeric vector of length 2. When provided, overrides the limits
#'   computed by \code{safe_lims()}. Intended for use by
#'   \code{\link{plot_ellipsoid_pairs}} to enforce consistent axis scales across
#'   panels, but can also be set manually by the user. Default is \code{NULL}
#'   (limits computed from data).
#' @param ... Additional graphical parameters passed to
#'   \code{\link[graphics]{plot}}.
#'
#' @details
#' The function has three display modes depending on what is provided:
#' \enumerate{
#'   \item \strong{Background only} (\code{background} is not \code{NULL}):
#'   plots background points in \code{col_bg} with the ellipsoid boundary
#'   overlaid.
#'   \item \strong{Prediction surface} (\code{background} is \code{NULL},
#'   \code{prediction} is not \code{NULL}): plots prediction points, optionally
#'   colored by \code{col_layer} using values mapped onto \code{pal}.
#'   When \code{col_layer} is provided, points outside the ellipsoid (zero or
#'   \code{NA} in \code{col_layer}, as produced by truncated prediction types)
#'   are drawn in \code{col_bg} behind the colored interior points. Axis limits
#'   are computed from the full \code{prediction} extent so the view is never
#'   collapsed to the ellipsoid interior.
#'   \item \strong{Ellipsoid only} (both \code{NULL}): plots the ellipsoid
#'   boundary alone with no background.
#' }
#'
#' @return Called for its side effect of creating a plot. Returns \code{NULL}
#'   invisibly.
#'
#' @seealso \code{\link{add_data}} to overlay occurrence points,
#'   \code{\link{add_ellipsoid}} to overlay additional ellipsoid boundaries,
#'   \code{\link{plot_ellipsoid_pairs}} for pairwise plots of all dimensions.
#'
#' @importFrom graphics plot lines
#' @importFrom grDevices adjustcolor
#'
#' @export
plot_ellipsoid <- function(object,
                           background = NULL,
                           prediction = NULL,
                           dim = c(1, 2),
                           col_layer = NULL,
                           pal = hcl.colors(100, palette = "Viridis"),
                           rev_pal = FALSE,
                           bg_sample = NULL,
                           lty = 1,
                           lwd = 1,
                           col_ell = "#000000",
                           col_bg = "#8A8A8A",
                           pch = 1,
                           alpha_bg = 1,
                           alpha_ell = 1,
                           cex_ell = 1,
                           cex_bg = 1,
                           fixed_lims = NULL,
                           ...) {

  if (missing(object) || !inherits(object, "nicheR_ellipsoid")) {
    stop("Please provide a valid 'nicheR_ellipsoid' object.")
  }

  if (!is.null(background) && !is.data.frame(background) &&
      !is.matrix(background)) {
    stop("Background must be a 'data.frame' or 'matrix'.")
  }

  if (!is.null(prediction) && !is.data.frame(prediction) &&
      !is.matrix(prediction)) {
    stop("Prediction must be a 'data.frame' or 'matrix'.")
  }

  if (!is.null(fixed_lims)) {
    if (!is.list(fixed_lims) || !all(c("xlim", "ylim") %in% names(fixed_lims))) {
      stop("'fixed_lims' must be a named list with elements 'xlim' and 'ylim'.")
    }
  }

  # Calculate ellipse boundaries
  ell_points <- ellipsoid_boundary_2d(object = object,
                                      n_segments = 100,
                                      dim = dim)

  # Background sampling
  if (!is.null(background)) {

    if (!is.null(bg_sample) && nrow(background) > bg_sample) {
      pts_indx <- sample(seq_len(nrow(background)), bg_sample)
    } else {
      pts_indx <- seq_len(nrow(background))
    }

    pts_xy <- background[pts_indx, c(object$var_names[dim]), drop = FALSE]
    lims   <- if (!is.null(fixed_lims)) fixed_lims else safe_lims(pts_xy, ell_points)

    plot(pts_xy[[1L]], pts_xy[[2L]],
         col = adjustcolor(col_bg, alpha.f = alpha_bg),
         pch = pch,
         cex = cex_bg,
         xlim = lims$xlim,
         ylim = lims$ylim, ...)

    lines(ell_points,
          lty = lty,
          lwd = lwd,
          col = adjustcolor(col_ell, alpha.f = alpha_ell),
          cex = cex_ell)

  } else if (!is.null(prediction)){

    if(is.null(col_layer)){

      if (!is.null(bg_sample) && nrow(prediction) > bg_sample) {
        pts_indx <- sample(seq_len(nrow(prediction)), bg_sample)
      } else {
        pts_indx <- seq_len(nrow(prediction))
      }

      pts_xy <- prediction[pts_indx, c(object$var_names[dim]), drop = FALSE]
      lims   <- if (!is.null(fixed_lims)) fixed_lims else safe_lims(pts_xy, ell_points)

      plot(pts_xy[[1L]], pts_xy[[2L]],
           col = adjustcolor(col_bg, alpha.f = alpha_bg),
           pch = pch,
           cex = cex_bg,
           xlim = lims$xlim,
           ylim = lims$ylim, ...)

      lines(ell_points,
            lty = lty,
            lwd = lwd,
            col = adjustcolor(col_ell, alpha.f = alpha_ell),
            cex = cex_ell)


    } else {

      if(is.function(pal)){
        pal <- pal(100)
      }

      if(isTRUE(rev_pal)){
        pal <- rev(pal)
      }

      # Split on col_layer: inside = non-NA and non-zero, outside = zero or NA.
      # Zeros arise from suitability_trunc outside the ellipsoid;
      # NAs arise from Mahalanobis_trunc outside the ellipsoid.
      # Both cases get col_bg. Limits use the full prediction so the view
      # never collapses to the ellipsoid interior when data are truncated.
      col_vals   <- prediction[ , col_layer]
      is_outside <- is.na(col_vals) | col_vals == 0
      is_inside  <- !is_outside

      pred_outside <- prediction[is_outside, c(object$var_names[dim]), drop = FALSE]
      pred_inside  <- prediction[is_inside,  c(object$var_names[dim]), drop = FALSE]

      # Subsample outside points if requested
      if (!is.null(bg_sample) && nrow(pred_outside) > bg_sample) {
        pts_indx_out <- sample(seq_len(nrow(pred_outside)), bg_sample)
      } else {
        pts_indx_out <- seq_len(nrow(pred_outside))
      }

      # Subsample inside points if requested
      if (!is.null(bg_sample) && nrow(pred_inside) > bg_sample) {
        pts_indx_in <- sample(seq_len(nrow(pred_inside)), bg_sample)
      } else {
        pts_indx_in <- seq_len(nrow(pred_inside))
      }


      col_indx <- map_to_pal(col_vals[is_inside][pts_indx_in], length(pal))

      # Limits: fixed_lims takes priority, otherwise span full prediction extent
      lims <- if (!is.null(fixed_lims)) fixed_lims else
        safe_lims(prediction[ , c(object$var_names[dim]), drop = FALSE], ell_points)

      # Open the device with plot.default using explicit x/y vectors so
      # xlim and ylim are always respected. plot() on a data frame dispatches
      # to plot.data.frame which ignores xlim/ylim and sizes axes to the data
      # it receives, collapsing the view to just the outside or inside region.
      plot(lims$xlim, lims$ylim,
           type = "n",
           xlim = lims$xlim,
           ylim = lims$ylim, ...)

      # Outside points first (grey, behind), then inside (colored)
      if (nrow(pred_outside) > 0L) {
        points(pred_outside[pts_indx_out, 1L],
               pred_outside[pts_indx_out, 2L],
               col = adjustcolor(col_bg, alpha.f = alpha_bg),
               pch = pch,
               cex = cex_bg)
      }

      if (nrow(pred_inside) > 0L) {
        points(pred_inside[pts_indx_in, 1L],
               pred_inside[pts_indx_in, 2L],
               col = pal[col_indx],
               pch = pch,
               cex = cex_bg)
      }

      lines(ell_points,
            lty = lty,
            lwd = lwd,
            col = adjustcolor(col_ell, alpha.f = alpha_ell),
            cex = cex_ell)
    }
  }else{
    # Ellipsoid boundary only, no background
    plot(ell_points, type = "l",
         lty = lty, lwd = lwd,
         col = adjustcolor(col_ell, alpha.f = alpha_ell),
         cex = cex_ell, ...)
  }

}



#' Add occurrence points or other data to an existing E-space plot
#'
#' @description
#' Adds points to an existing environmental space plot created with
#' \code{plot_ellipsoid()}. Points can be plotted with a single color or
#' colored by a continuous variable (e.g., suitability) using a color palette.
#'
#' @param data A data frame containing the points to plot. Must include
#'   columns matching \code{x} and \code{y}, and \code{col_layer} if provided.
#' @param x Character. Name of the column to use as the x-axis variable.
#' @param y Character. Name of the column to use as the y-axis variable.
#' @param pts_col Character. Color for all points when \code{col_layer} is
#'   \code{NULL}. Default is \code{"#000000"} (black).
#' @param pts_alpha Numeric in \code{[0, 1]}. Transparency of points when
#'   \code{col_layer} is \code{NULL}. Default is \code{1} (fully opaque).
#' @param col_layer Character or \code{NULL}. Name of a column in \code{data}
#'   to use for coloring points by a continuous variable. If \code{NULL}
#'   (default), all points are drawn with \code{pts_col}.
#' @param pal A color palette function or character vector of colors used when
#'   \code{col_layer} is provided. Default is \code{hcl.colors(100, palette = "Viridis")}.
#' @param rev_pal Logical. If \code{TRUE}, reverses the color palette before
#'   applying it. Default is \code{FALSE}.
#' @param pch Integer or character. Point symbol. Default is \code{1}.
#' @param cex Numeric. Size scaling for points. Default is \code{1}.
#' @param bg_sample Integer or \code{NULL}. If provided and \code{nrow(data)}
#'   exceeds this value, a random subsample of this size is drawn before
#'   plotting. Useful for large data frames. Default is \code{NULL} (plot all).
#' @param ... Additional arguments passed to \code{\link[graphics]{points}}.
#'
#' @details
#' When \code{col_layer} is provided, points are colored by the values of that
#' column mapped onto the palette. \code{NA}s in \code{col_layer} are removed
#' before plotting; zeros are retained as valid values (e.g., truncated
#' suitability predictions outside the ellipsoid boundary).
#'
#' @return Called for its side effect of adding points to the current plot.
#'   Returns \code{NULL} invisibly.
#'
#' @seealso \code{\link{plot_ellipsoid}}, \code{\link{add_ellipsoid}}
#'
#' @importFrom graphics points
#' @importFrom grDevices adjustcolor
#'
#' @export
add_data <- function(data, x, y,
                     pts_col = "#000000",
                     pts_alpha  = 1,
                     col_layer = NULL,
                     pal = hcl.colors(100, palette = "Viridis"),
                     rev_pal = FALSE,
                     pch = 1,
                     cex = 1,
                     bg_sample = NULL,
                     ...) {

  if(is.null(col_layer)){

    if (!is.null(bg_sample) && nrow(data) > bg_sample) {
      pts_indx <- sample(seq_len(nrow(data)), bg_sample)
    } else {
      pts_indx <- seq_len(nrow(data))
    }

    points(data[pts_indx, c(x, y)],
           col = adjustcolor(pts_col, alpha.f = pts_alpha),
           pch = pch, cex = cex, ...)

  } else {

    # Remove NAs only; zeros are valid (e.g., truncated suitability outside the ellipsoid)
    data_clean <- data[!is.na(data[ , col_layer]), ]

    if (!is.null(bg_sample) && nrow(data_clean) > bg_sample) {
      pts_indx <- sample(seq_len(nrow(data_clean)), bg_sample)
    } else {
      pts_indx <- seq_len(nrow(data_clean))
    }

    if(is.function(pal)){
      pal <- pal(100)
    }

    if(isTRUE(rev_pal)){
      pal <- rev(pal)
    }

    col_indx <- map_to_pal(data_clean[pts_indx, col_layer], length(pal))

    points(data_clean[pts_indx, c(x, y)],
           col = pal[col_indx],
           pch = pch, cex = cex, ...)

  }
}


#' Add an ellipsoid boundary to an existing E-space plot
#'
#' @description
#' Draws the 2D boundary of a \code{nicheR_ellipsoid} object onto an existing
#' environmental space plot created with \code{plot_ellipsoid()}. The boundary
#' is computed as a cross-section of the ellipsoid at the chosen pair of
#' dimensions.
#'
#' @param object A \code{nicheR_ellipsoid} object.
#' @param dim Integer vector of length 2. Indices of the two dimensions to
#'   plot. Default is \code{c(1, 2)}.
#' @param lty Integer. Line type. Default is \code{1} (solid).
#' @param lwd Numeric. Line width. Default is \code{1}.
#' @param col_ell Character. Color of the ellipsoid boundary line. Default is
#'   \code{"#000000"} (black).
#' @param alpha_ell Numeric in \code{[0, 1]}. Transparency of the ellipsoid
#'   boundary line. Default is \code{1} (fully opaque).
#' @param cex_ell Numeric. Size scaling for the ellipsoid boundary. Default
#'   is \code{1}.
#' @param ... Additional arguments passed to \code{\link[graphics]{lines}}.
#'
#' @return Called for its side effect of adding lines to the current plot.
#'   Returns \code{NULL} invisibly.
#'
#' @seealso \code{\link{plot_ellipsoid}}, \code{\link{add_data}}
#'
#' @importFrom graphics lines
#' @importFrom grDevices adjustcolor
#'
#' @export
add_ellipsoid <- function(object,
                          dim = c(1, 2),
                          lty = 1,
                          lwd = 1,
                          col_ell = "#000000",
                          alpha_ell = 1,
                          cex_ell = 1, ...){

  if (!inherits(object, "nicheR_ellipsoid"))
    stop("'object' must be a nicheR_ellipsoid.")

  ell_points <- ellipsoid_boundary_2d(object = object,
                                      n_segments = 100,
                                      dim = dim)

  lines(ell_points,
        lty = lty,
        lwd = lwd,
        col = adjustcolor(col_ell, alpha.f = alpha_ell),
        cex = cex_ell, ...)

}

#' Generate 2D ellipsoid boundary points for plotting
#'
#' Computes ordered boundary points for a two-dimensional slice of a
#' \code{nicheR_ellipsoid}, suitable for plotting with \code{lines()} or
#' \code{plot(type = "l")}. The boundary is derived from the selected
#' covariance submatrix and the stored chi-square cutoff.
#'
#' @param object A \code{nicheR_ellipsoid} object.
#' @param n_segments Integer. Number of boundary points to generate
#'   (must be >= 4).
#' @param dim Integer vector of length 2 indicating which dimensions
#'   (indices of the original variables) to use for the 2D slice.
#'
#' @return A \code{data.frame} with \code{n_segments} ordered boundary points
#'   in the selected dimensions.
#'
#' @keywords internal
ellipsoid_boundary_2d <- function(object,
                                  n_segments = 50,
                                  dim = c(1, 2)) {

  if (!inherits(object, "nicheR_ellipsoid")) stop("'object' must be a nicheR_ellipsoid.")
  if (!is.numeric(n_segments) || length(n_segments) != 1L || !is.finite(n_segments) || n_segments < 4) {
    stop("'n_segments' must be a single number >= 4.")
  }

  n_segments <- as.integer(n_segments)

  mu2 <- object$centroid[dim]
  Sigma2 <- object$cov_matrix[dim, dim, drop = FALSE]
  chi2_cutoff <- object$chi2_cutoff

  eig2 <- eigen(Sigma2, symmetric = TRUE)
  vals2 <- pmax(eig2$values, 0)
  a2 <- sqrt(vals2 * chi2_cutoff)
  rot2 <- eig2$vectors

  t <- seq(0, 2 * pi, length.out = n_segments)
  unit <- rbind(a2[1] * cos(t), a2[2] * sin(t))

  pts <- t(rot2 %*% unit)
  pts <- sweep(pts, 2, mu2, "+")
  pts <- as.data.frame(pts)

  # name columns by variable names
  if (!is.null(object$var_names) && length(object$var_names) >= max(dim)) {
    colnames(pts) <- object$var_names[dim]
  } else {
    colnames(pts) <- paste0("dim", dim)
  }

  pts
}

#' Plot all pairwise 2D ellipsoid projections
#'
#' Plots all pairwise two-dimensional slices of a \code{nicheR_ellipsoid}
#' in a multi-panel layout using \code{plot_ellipsoid()}. When \code{background}
#' or \code{prediction} is supplied, axis limits are computed once from the
#' global range of all variables and shared across every panel, so projections
#' are directly comparable without distortion from per-panel rescaling.
#'
#' @param object A \code{nicheR_ellipsoid} object.
#' @param background Optional data frame or matrix of background points passed
#'   to each \code{plot_ellipsoid()} call. When provided, global axis limits
#'   are computed from the range of all variables in \code{background} combined
#'   with all pairwise ellipsoid boundaries.
#' @param prediction Optional data frame or matrix of prediction values passed
#'   to each \code{plot_ellipsoid()} call. Used when \code{background} is
#'   \code{NULL}. Global limits are computed from the range of all variables
#'   in \code{prediction}.
#' @param ... Additional graphical arguments passed to \code{plot_ellipsoid()}.
#'
#' @details
#' Global limits are computed per variable across the full data and all
#' ellipsoid boundary projections, then passed to each panel via the
#' \code{fixed_lims} argument of \code{plot_ellipsoid()}. This prevents
#' individual panels from rescaling to their own data extent, which would make
#' niche widths appear identical across dimensions even when they differ. If
#' neither \code{background} nor \code{prediction} is provided, each panel
#' shows only the ellipsoid boundary and limits come from that boundary alone,
#' which is the intended behavior for a boundary-only view.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @importFrom graphics par
#' @importFrom utils combn
#'
#' @export
plot_ellipsoid_pairs <- function(object,
                                 background = NULL,
                                 prediction = NULL,
                                 ...) {

  if (!inherits(object, "nicheR_ellipsoid"))
    stop("'object' must be a nicheR_ellipsoid.")

  pairs_idx <- t(combn(seq_len(object$dimensions), 2))
  n_pairs   <- nrow(pairs_idx)
  n_cols    <- ceiling(sqrt(n_pairs))
  n_rows    <- ceiling(n_pairs / n_cols)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(n_rows, n_cols))

  # Pre-compute per-variable ranges across data and all ellipsoid boundaries,
  # so every panel shares the same axis scale for each variable.
  data_ref <- if (!is.null(background)) background else prediction

  if (!is.null(data_ref)) {
    var_ranges <- lapply(object$var_names, function(v) {
      # Collect ellipsoid boundary points across all pairs that include this variable
      ell_vals <- unlist(lapply(seq_len(n_pairs), function(i) {
        if (v %in% object$var_names[pairs_idx[i, ]]) {
          ell <- ellipsoid_boundary_2d(object, n_segments = 100, dim = pairs_idx[i, ])
          ell[[v]]
        }
      }))
      data_vals <- if (v %in% colnames(data_ref)) data_ref[[v]] else numeric(0)
      range(c(data_vals, ell_vals), na.rm = TRUE)
    })
    names(var_ranges) <- object$var_names
  } else {
    var_ranges <- NULL
  }

  for (i in seq_len(n_pairs)) {
    v1 <- object$var_names[pairs_idx[i, 1]]
    v2 <- object$var_names[pairs_idx[i, 2]]

    fixed_lims <- if (!is.null(var_ranges)) {
      list(xlim = var_ranges[[v1]], ylim = var_ranges[[v2]])
    } else {
      NULL
    }

    plot_ellipsoid(object     = object,
                   background = background,
                   prediction = prediction,
                   dim        = pairs_idx[i, ],
                   fixed_lims = fixed_lims,
                   main       = paste0(v1, " vs. ", v2),
                   ...)
  }

  invisible(NULL)
}




# Helper: map a numeric vector to palette indices via linear min-max scaling
#
#' Map numeric values to palette indices
#'
#' Rescales a numeric vector linearly from its observed range onto integer
#' indices in \code{[1, pal_len]} for use in palette lookups. When all values
#' are identical, returns the middle index for every element.
#'
#' @param vals Numeric vector of values to map.
#' @param pal_len Integer. Length of the target palette.
#'
#' @return Integer vector of the same length as \code{vals}, with values in
#'   \code{[1, pal_len]}.
#'
#' @keywords internal
map_to_pal <- function(vals, pal_len) {
  vmin <- min(vals, na.rm = TRUE)
  vmax <- max(vals, na.rm = TRUE)
  if (vmax == vmin) {
    return(rep(ceiling(pal_len / 2L), length(vals)))
  }
  as.integer(ceiling((vals - vmin) / (vmax - vmin) * (pal_len - 1L))) + 1L
}



# Helper: compute xlim/ylim that covers both a data matrix and ell_points
#
#' Compute safe axis limits covering data and ellipsoid boundary
#'
#' Returns x and y ranges that span both a matrix of data points and a
#' data frame of ellipsoid boundary points, so the ellipsoid boundary is
#' never clipped by the data extent when passed to \code{plot()}.
#'
#' @param pts_xy A data frame or matrix with at least two columns. The first
#'   column is used for x, the second for y.
#' @param ell_xy A data frame or matrix of ellipsoid boundary points with the
#'   same column structure as \code{pts_xy}.
#'
#' @return A named list with elements \code{xlim} and \code{ylim}, each a
#'   numeric vector of length 2.
#'
#' @keywords internal
safe_lims <- function(pts_xy, ell_xy) {
  all_x <- c(pts_xy[, 1], ell_xy[, 1])
  all_y <- c(pts_xy[, 2], ell_xy[, 2])
  list(
    xlim = range(all_x, na.rm = TRUE),
    ylim = range(all_y, na.rm = TRUE)
  )
}
