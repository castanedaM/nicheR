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
#' in a multi-panel layout using \code{plot_ellipsoid()}.
#'
#' @param object A \code{nicheR_ellipsoid} object.
#' @param ... Additional graphical arguments passed to \code{plot_ellipsoid()}.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @importFrom graphics par
#' @importFrom utils combn
#'
#' @export
plot_ellipsoid_pairs <- function(object, ...) {

  if (!inherits(object, "nicheR_ellipsoid"))
    stop("'object' must be a nicheR_ellipsoid.")

  pairs_idx <- t(combn(seq_len(object$dimensions), 2))
  n_pairs <- nrow(pairs_idx)

  n_cols <- ceiling(sqrt(n_pairs))
  n_rows <- ceiling(n_pairs / n_cols)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_rows, n_cols))

  for (i in seq_len(n_pairs)) {

    plot_ellipsoid(object = object,
                   dim = pairs_idx[i, ],
                   main = paste0(object$var_names[pairs_idx[i, 1]],
                                 "vs.", object$var_names[pairs_idx[i, 2]]),
                   ...)
  }

  invisible(NULL)
}





