#' Plot a nicheR Ellipsoid
#'
#' Plots the boundary of a probabilistic ellipsoid in environmental space,
#' optionally overlaying background points.
#'
#' @param object A nicheR ellipsoid object containing \code{centroid},
#'   \code{cov_matrix}, and \code{chi2_cutoff}.
#' @param background Optional data frame of background points (rows = observations,
#'   columns = environmental variables).
#' @param bg_sample Integer. Number of background points to sample for plotting.
#'   The default, NULL, plots all background points.
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
                           prediction = NULL,
                           dim = c(1, 2),
                           col_layer = NULL,
                           pal = heat.colors(100),
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
                           cex_bg = 1, ...) {

  if (missing(object) || !inherits(object, "nicheR_ellipsoid")) {
    stop("Please provide a valid 'nicheR_ellipsoid' object.")
  }
  if (!is.null(background) && !is.data.frame(background) &&
        !is.matrix(background)) {
    stop("Background must be a 'data.frame' or 'matrix'.")
  }

  # Calculate ellipse boundaries
  ell_points <- ellipsoid_boundary_2d(object = object, n_segments = 100,
                                      dim = dim)

  # Background sampling
  if (!is.null(background)) {
    if (!is.null(bg_sample) && nrow(background) > bg_sample) {
      bg_sample_bg <- sample(seq_len(nrow(background)), bg_sample)
    } else {
      bg_sample_bg <- seq_len(nrow(background))
    }

    plot(background[bg_sample_bg, c(object$var_names[dim])],
         col = adjustcolor(col_bg, alpha.f = alpha_bg),
         pch = pch,
         cex = cex_bg, ...)

    lines(ell_points,
          lty = lty,
          lwd = lwd,
          col = adjustcolor(col_ell, alpha.f = alpha_ell),
          cex = cex_ell)

  } else if (!is.null(prediction)){

    if(is.null(col_layer)){

      if(nrow(prediction) <= bg_sample){
        pred_sample <- nrow(prediction)
      }

      pred_sample_indx <- sample(1:nrow(prediction), pred_sample)


      plot(prediction[pred_sample_indx, c(object$var_names[dim])],
           col = adjustcolor(col_bg, alpha.f = alpha_bg),
           pch = pch,
           cex = cex_bg, ...)

      lines(ell_points,
            lty = lty,
            lwd = lwd,
            col = adjustcolor(col_ell, alpha.f = alpha_ell),
            cex = cex_ell)


    } else {

      # Remove zeros and NAs
      col_layer_clean <- prediction[ , col_layer]
      col_layer_clean <- col_layer_clean[col_layer_clean > 0 & !is.na(col_layer_clean)]

      if(length(col_layer_clean) <= bg_sample){
        pred_sample <- length(col_layer_clean)
      }

      pred_sample_indx <- sample(1:length(col_layer_clean), pred_sample)

      if(is.function(pal)){
        pal <- pal(100)
      }

      if(isTRUE(rev_pal)){
        pal <- rev(pal)
      }


      col_indx <- as.numeric(cut(log(col_layer_clean[pred_sample_indx]),
                                 breaks = length(pal),
                                 include.lowest = TRUE))

      plot(prediction[pred_sample_indx, c(object$var_names[dim])],
           col = pal[col_indx],
           pch = pch,
           cex = cex_bg, ...)

      lines(ell_points,
            lty = lty,
            lwd = lwd,
            col = adjustcolor(col_ell, alpha.f = alpha_ell),
            cex = cex_ell)
    }
  }else{
    # Basic line for elliposid
    plot(ell_points, type = "l",
         lty = lty, lwd = lwd,
         col = adjustcolor(col_ell, alpha.f = alpha_ell),
         cex = cex_ell, ...)
  }

}



#' @export
add_data <- function(data, x, y,
                     pts_col = "#000000",
                     pts_alpha  = 1,
                     col_layer = NULL,
                     pal = heat.colors(100),
                     rev_pal = FALSE,
                     pch = 1,
                     pts_sample = 1000,
                     ...) {

  if(is.null(col_layer)){

    if(nrow(data) > pts_sample){
      pts_idx <- sample(1:nrow(data), pts_sample)
    }else{
      pts_idx <- 1:nrow(data)
    }


    points(data[pts_idx, c(x, y)],
           col = adjustcolor(pts_col, alpha.f = pts_alpha),
           pch = pch, ...)

  } else {

    # Remove zeros and NAs
    col_layer_clean <- data[  , col_layer]
    col_layer_clean <- col_layer_clean[col_layer_clean > 0 & !is.na(col_layer_clean)]

    if(length(col_layer_clean) <= pts_sample){
      pts_sample <- length(col_layer_clean)
    }

    pts_idx <- sample(1:length(col_layer_clean), pts_sample)

    if(is.function(pal)){
      pal <- pal(100)
    }

    if(isTRUE(rev_pal)){
      pal <- rev(pal)
    }


    col_indx <- as.numeric(cut(log(col_layer_clean[pts_idx]),
                               breaks = length(pal),
                               include.lowest = TRUE))

    points(data[pts_idx, c(x, y)],
           col = pal[col_indx],
           pch = pch, ...)

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



