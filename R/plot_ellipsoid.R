#' Plot a nicheR Ellipsoid
#'
#' Plots the boundary of a probabilistic ellipsoid in environmental space,
#' optionally overlaying background points.
#'
#' @param object A nicheR ellipsoid object containing \code{centroid},
#'   \code{cov_matrix}, and \code{chi2_cutoff}.
#' @param background Optional data frame of background points (rows = observations,
#'   columns = environmental variables).
#' @param sample Integer. Number of background points to sample for plotting.
#' @param lty Line type for ellipsoid boundary.
#' @param lwd Line width for ellipsoid boundary.
#' @param col_ell Color for ellipsoid boundary.
#' @param col_bg Color for background points.
#' @param pch Plotting symbol for background points.
#' @param alpha_bg Transparency for background points.
#' @param alpha_ell Transparency for ellipsoid boundary.
#' @param cex_ell Expansion factor for ellipsoid line.
#' @param cex_bg Expansion factor for background points.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @return Invisibly returns \code{NULL}.
#' @export
plot_ellipsoid <- function(object,
                           background = NULL,
                           sample = 1000,
                           lty = 1,
                           lwd = 1,
                           col_ell = "#000000",
                           col_bg = "#8A8A8A",
                           pch = 1,
                           alpha_bg = 0.5,
                           alpha_ell = 0.7,
                           cex_ell = 1,
                           cex_bg = 1, ...){

#   Check for data frame
  ell_points <- ellipsoid_surface_points(mu_vec = object$centroid,
                                         cov_matrix = object$cov_matrix,
                                         chi2_cutoff = object$chi2_cutoff,
                                         n_point = 100)
  # to do: make sure name of vars does not disappear

  if(!is.null(background)){

    if(nrow(background) <= sample){
      sample <- nrow(background)
    }

    sample_bg <- sample(1:nrow(background), sample)

    plot(background[sample_bg, ],
         col = adjustcolor(col_bg, alpha.f = alpha_bg),
         pch = pch,
         cex = cex_bg, ...)

    lines(ell_points,
          lty = lty,
          lwd = lwd,
          col = adjustcolor(col_ell, alpha.f = alpha_ell),
          cex = cex_ell)

  }else{
    # Basic line for elliposid
    plot(ell_points, type = "l",
         lty = lty, lwd = lwd,
         col = adjustcolor(col_ell, alpha.f = alpha_ell),
         cex = cex_ell, ...)
  }

}

