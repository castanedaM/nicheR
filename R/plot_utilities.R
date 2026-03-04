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
                   ...)
  }

  invisible(NULL)
}




#' Plot 2d base R just for testing right now
#'
#' @export
plot_nicheR <- function(ell_list = NULL,
                        dims = c(1, 2),
                        n_points = 200,
                        labels = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        main = "Ellipsoid",
                        add = FALSE,
                        legend_pos = "topright",
                        occ_list = NULL,
                        occ_pch = 1,
                        occ_col = "green4",
                        occ_cex = 0.6,
                        ...){

  if(is.null(ell_list) || length(ell_list) < 1){
    stop("'ell_list' must be a non-empty list of ellipsoid objects.")
  }
  if(!is.numeric(dims) || length(dims) != 2L){
    stop("'dims' must be an integer vector of length 2.")
  }
  if(dims[1] == dims[2]){
    stop("'dims' must contain two different dimensions.")
  }
  if(!is.numeric(n_points) || length(n_points) != 1L || n_points < 50){
    stop("'n_points' must be a single number >= 50.")
  }

  # check var_names exists
  if(is.null(ell_list[[1]]$var_names)){
    stop("ell_list[[1]] must contain 'var_names' to plot environmental points in E-space.")
  }

  var_names <- ell_list[[1]]$var_names

  if(max(dims) > length(var_names)){
    stop("Requested 'dims' exceed the number of variables in ellipsoid var_names.")
  }

  xvar <- var_names[dims[1]]
  yvar <- var_names[dims[2]]

  # collect boundary points + centroids
  boundary_points <- vector("list", length(ell_list))
  mu_mat <- matrix(NA_real_, nrow = length(ell_list), ncol = 2)

  for(i in seq_along(ell_list)){
    ell <- ell_list[[i]]

    if(is.null(ell$centroid) || is.null(ell$cov_matrix) || is.null(ell$chi2_cutoff)){
      stop("Each element of 'ell_list' must contain centroid, cov_matrix, and chi2_cutoff.")
    }
    if(is.null(ell$var_names)){
      stop("Each element of 'ell_list' must contain var_names.")
    }
    if(!identical(ell$var_names, var_names)){
      stop("All ellipsoids in 'ell_list' must share the same var_names in the same order.")
    }
    if(length(ell$centroid) < max(dims)){
      stop("An ellipsoid has fewer dimensions than requested in 'dims'.")
    }

    mu2 <- ell$centroid[dims]
    Sig2 <- ell$cov_matrix[dims, dims, drop = FALSE]

    boundary_points[[i]] <- ellipsoid_boundary_2d(object = object,
                                                  dim= dim)

    mu_mat[i, ] <- mu2
  }

  # axis labels
  if(is.null(xlab)){
    xlab <- xvar
  }

  if(is.null(ylab)){
    ylab <- yvar
  }

  # plot ranges from ellipsoids
  xs <- unlist(lapply(boundary_points, function(p) p[, 1]))
  ys <- unlist(lapply(boundary_points, function(p) p[, 2]))

  rx <- range(c(xs, mu_mat[, 1]))
  ry <- range(c(ys, mu_mat[, 2]))

  # add padding (10% of span)
  padx <- 0.10 * diff(rx)
  pady <- 0.10 * diff(ry)

  rx <- c(rx[1] - padx, rx[2] + padx)
  ry <- c(ry[1] - pady, ry[2] + pady)

  # if occurrences exist, expand ranges to include them
  if(!is.null(occ_list)){

    occ_tmp <- if(is.data.frame(occ_list)) list(occ_list) else occ_list

    if(!is.list(occ_tmp)){
      stop("'occ_list' must be a data.frame or a list of data.frames.")
    }

    occ_x_all <- c()
    occ_y_all <- c()

    for(j in seq_along(occ_tmp)){

      occ <- occ_tmp[[j]]

      if(!is.data.frame(occ)){
        stop("Each element of 'occ_list' must be a data.frame.")
      }

      if(!(xvar %in% names(occ)) || !(yvar %in% names(occ))){
        stop("occ_list data.frames must include the environmental columns being plotted: ",
             xvar, " and ", yvar, ".")
      }

      occ_x_all <- c(occ_x_all, occ[[xvar]])
      occ_y_all <- c(occ_y_all, occ[[yvar]])
    }

    occ_x_all <- occ_x_all[is.finite(occ_x_all)]
    occ_y_all <- occ_y_all[is.finite(occ_y_all)]

    if(length(occ_x_all) > 0){
      rx <- range(c(rx, occ_x_all))
    }
    if(length(occ_y_all) > 0){
      ry <- range(c(ry, occ_y_all))
    }
  }

  if(!isTRUE(add)){
    plot(NA,
         xlim = rx,
         ylim = ry,
         xlab = xlab,
         ylab = ylab,
         main = main,
         ...)
  }

  # draw ellipsoids + centroids
  ltys <- rep_len(c(1, 2, 3, 4, 5, 6), length(ell_list))
  lwds <- rep_len(2, length(ell_list))
  pchs <- rep_len(c(16, 4, 17, 3, 15, 8), length(ell_list))

  for(i in seq_along(boundary_points)){
    p <- boundary_points[[i]]
    lines(p[, 1], p[, 2], lwd = lwds[i], lty = ltys[i])
    points(mu_mat[i, 1], mu_mat[i, 2], pch = pchs[i], cex = 1.1, lwd = 2, col = "red")
  }

  # overlay environmental occurrences (if provided) ------------------------

  if(!is.null(occ_list)){

    if(is.data.frame(occ_list)){
      occ_list <- list(occ_list)
    }

    if(!is.list(occ_list)){
      stop("'occ_list' must be a data.frame or a list of data.frames.")
    }

    for(j in seq_along(occ_list)){

      occ <- occ_list[[j]]

      if(!is.data.frame(occ)){
        stop("Each element of 'occ_list' must be a data.frame.")
      }

      if(!(xvar %in% names(occ)) || !(yvar %in% names(occ))){
        stop("occ_list data.frames must include the environmental columns being plotted: ",
             xvar, " and ", yvar, ".")
      }

      points(occ[[xvar]], occ[[yvar]], pch = occ_pch, cex = occ_cex, col = occ_col)
    }
  }

  # legend
  if(is.null(labels)){
    labels <- paste0("ellipsoid ", seq_along(ell_list))
  }
  if(length(labels) != length(ell_list)){
    labels <- rep_len(labels, length(ell_list))
  }

  legend(legend_pos,
         legend = labels,
         lty = ltys,
         lwd = lwds,
         bty = "n")

  invisible(list(boundary_points = boundary_points,
                 centroids = mu_mat,
                 dims = dims,
                 xvar = xvar,
                 yvar = yvar))
}



