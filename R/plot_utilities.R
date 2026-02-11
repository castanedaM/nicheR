
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

    boundary_points[[i]] <- ellipsoid_surface_points(mu_vec = mu2,
                                                     cov_matrix = Sig2,
                                                     chi2_cutoff = ell$chi2_cutoff,
                                                     n_points = n_points)

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




#' Generate Ellipsoid Surface Coordinates
#'
#' Computes boundary coordinates of a 2D or 3D ellipsoid defined by a covariance
#' matrix and chi-square probability contour.
#'
#' @param mu_vec Numeric vector. Ellipsoid centroid.
#' @param cov_matrix Covariance matrix.
#' @param chi2_cutoff Chi-square cutoff defining the ellipsoid contour.
#' @param n_points Number of points used to approximate the surface.
#'
#' @return A data.frame containing surface coordinates:
#' \itemize{
#'   \item \code{x1, x2} for 2D
#'   \item \code{x1, x2, x3} for 3D
#' }
#'
#' @details
#' Surface points are generated by rotating and scaling unit circles or spheres
#' using the eigen decomposition of the covariance matrix.
#'
#' @keywords internal
ellipsoid_surface_points <- function(mu_vec,
                                     cov_matrix,
                                     chi2_cutoff,
                                     n_points) {

  # Idea: For easier plotting, save the 2D surfaces pairs
  dimensions <- length(mu_vec)

  # eigen decomposition gives rotation + axis lengths
  eig <- eigen(cov_matrix, symmetric = TRUE)

  # semi-axis lengths at chi-square level
  semi_axes_lengths <- sqrt(eig$values * chi2_cutoff)

  # rotation matrix
  cor_matrix <- eig$vectors

  if (dimensions == 2) {

    t <- seq(0, 2 * pi, length.out = n_points)

    unit_ellipse <- rbind(
      semi_axes_lengths[1] * cos(t),
      semi_axes_lengths[2] * sin(t)
    )

    rotated_points <- t(cor_matrix %*% unit_ellipse)

    surface <- sweep(rotated_points, 2, mu_vec, "+")

    surface <- as.data.frame(surface)
    colnames(surface) <- c("x1", "x2")

    return(surface)
  }

  if (dimensions == 3) {

    grid <- expand.grid(
      u = seq(0, 2 * pi, length.out = n_points),
      v = seq(0, pi, length.out = n_points)
    )

    unit_ellipsoid <- rbind(
      semi_axes_lengths[1] * sin(grid$v) * cos(grid$u),
      semi_axes_lengths[2] * sin(grid$v) * sin(grid$u),
      semi_axes_lengths[3] * cos(grid$v)
    )

    rotated_points <- t(cor_matrix %*% unit_ellipsoid)

    surface <- sweep(rotated_points, 2, mu_vec, "+")

    surface <- as.data.frame(surface)
    colnames(surface) <- c("x1", "x2", "x3")

    return(surface)
  }

  stop("Surface points only implemented for 2D or 3D.")
}
