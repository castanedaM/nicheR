#' Plot a nicheR Ellipsoid
#'
#' Plots the boundary of a probabilistic ellipsoid in environmental space,
#' optionally overlaying background points.
#'
#' @param object A nicheR ellipsoid object containing \code{centroid},
#'   \code{cov_matrix}, and \code{chi2_cutoff}.
#' @param background Optional data frame of background points (rows = observations,
#'   columns = environmental variables).
#' @param bg_sample Integer. Number of background points to bg_sample for plotting.
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
                           dim = c(1, 2),
                           bg_sample = 1000,
                           lty = 1,
                           lwd = 1,
                           col_ell = "#000000",
                           col_bg = "#8A8A8A",
                           pch = 1,
                           alpha_bg = 0.5,
                           alpha_ell = 0.7,
                           cex_ell = 1,
                           cex_bg = 1, ...){

  # Check for data frame
  ell_points <- ellipsoid_boundary_2d(object = object,
                                      n_segments = 50,
                                      dim = dim)

  # to do: make sure name of vars does not disappear

  if(!is.null(background)){

    if(nrow(background) <= bg_sample){
      bg_sample <- nrow(background)
    }

    bg_sample_bg <- sample(1:nrow(background), bg_sample)

    plot(background[bg_sample_bg, ],
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



#' @export
add_data <- function(x, y,
                     pts_col = "#000000",
                     pts_alpha  = 1,
                     col_layer = NULL,
                     pal = heat.colors(100),
                     rev_col = FALSE,
                     pch = 1,
                     pts_sample = 1000,
                     ...) {

  if(length(x) > pts_sample || length(y) > pts_sample ||
     length(col_layer) > pts_sample){

    pts_idx <- sample(1:length(x), pts_sample)

    x <- x[pts_idx]
    y <- y[pts_idx]

  }else{
    pts_idx <- 1:length(x)
  }


  # Now plot

  if(is.null(col_layer)){

    points(x = x,
           y = y,
           pch = pch,
           col = adjustcolor(pts_col, alpha.f = pts_alpha),
           ...)


  } else {

    col_layer <- col_layer[pts_idx]

    if(is.function(pal)){
      pal <- pal(100)
    }

    if(rev_col){
      pal <- rev(pal)
    }

    # Remove zeros and NAs
    x <- x[col_layer > 0 & !is.na(col_layer)]
    y <- y[col_layer > 0 & !is.na(col_layer)]

    col_layer <- col_layer[col_layer > 0 & !is.na(col_layer)]

    col_indx <- as.numeric(cut(log(col_layer),
                               breaks = length(pal),
                               include.lowest = TRUE))
    points(x = x,
           y = y,
           col = pal[col_indx],
           pch = pch,
           ...)
  }
}


#' @export
add_ellipsoid <- function(object,
                          dim = c(1, 2),
                          lty = 1,
                          lwd = 1,
                          col_ell = "#000000",
                          pch = 1,
                          alpha_ell = 1,
                          cex_ell = 1, ...){

  # Check for data frame
  ell_points <- ellipsoid_boundary_2d(object = object,
                                      n_segments = 50,
                                      dim = dim)

  # Basic line for elliposid
  lines(ell_points,
        lty = lty,
        lwd = lwd,
        col = adjustcolor(col_ell, alpha.f = alpha_ell),
        cex = cex_ell, ...)

}



