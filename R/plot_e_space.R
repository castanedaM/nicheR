#' Helper: Retrieve environmental background for plot_e_space()
#'
#' Priority:
#'  1. User-supplied env_bg
#'  2. From vs via nr_get_env()
#'  3. If neither exists → message and use suitable_env or niche extents
#'
#' @keywords internal
.pe_get_env_bg <- function(env_bg, vs, suitable_env, niche) {

  # 1. User-supplied
  if (!is.null(env_bg)) {
    message("Using user-supplied env_bg.")
    return(env_bg)
  }

  # 2. Try to get from vs
  if (!is.null(vs)) {
    env_vs <- nr_get_env(vs)
    if (!is.null(env_vs)) {
      message("Using env_bg from NicheR_species object via nr_get_env().")
      return(env_vs)
    }
  }

  # 3. No direct env_bg found
  message(
    "No env_bg provided and none found in vs. ",
    "Will build environment using suitable_env or niche extents."
  )

  # 4. If suitable_env exists, convert to df
  if (!is.null(suitable_env)) {
    s_df <- nr_get_suitable_df(suitable_env)
    if (!is.null(s_df)) {
      message("Using suitable_env values as environmental background.")
      return(s_df)
    }
  }

  # 5. If niche exists, simulate bounding box grid
  if (!is.null(niche)) {
    message("Using ellipsoid bounding box as env_bg (no raster or df available).")
    bb <- .pe_make_bbox_df(niche)
    return(bb)
  }

  # 6. Nothing available → stop
  stop(
    "Could not determine env_bg. Supply env_bg explicitly ",
    "or provide a vs/suitable_env/niche object."
  )
}

#' Helper: Create bounding-box environmental grid from ellipsoid
#'
#' @keywords internal
.pe_make_bbox_df <- function(niche, n = 5000) {
  cx <- niche$center
  ax <- niche$axes

  # min/max ranges
  mins <- cx - ax
  maxs <- cx + ax

  grid <- as.data.frame(matrix(
    runif(n * length(cx), mins, maxs),
    ncol = length(cx)
  ))
  names(grid) <- names(cx)

  return(grid)
}

#' Helper: Get suitable-env points as data.frame for E-space plotting
#'
#' For E-space we only use data.frame-based suitability. If no
#' data.frame is available (e.g., only rasters are stored), suitable
#' points are not plotted to avoid unexpected large conversions.
#'
#' @keywords internal
.pe_get_suitable_pts <- function(suitable_env) {

  if (is.null(suitable_env)) return(NULL)

  # 1) Try via nr_get_suitable_df() (works on vs, suitable_env objects, etc.)
  df <- suppressWarnings(nr_get_suitable_df(suitable_env))
  if (!is.null(df)) return(df)

  # 2) If user directly passed a df / matrix, use that
  if (is.data.frame(suitable_env) || is.matrix(suitable_env)) {
    return(as.data.frame(suitable_env))
  }

  # 3) Otherwise, do NOT try to create a df from rasters here
  #    (avoid big implicit conversions in E-space).
  message(
    "No suitable_env data.frame found; suitable environments will not be ",
    "plotted in E-space. If you want them, re-run get_suitable_env() with ",
    "out.suit = 'both' (or 'data.frame') and pass that object."
  )

  return(NULL)
}


#' Helper: Extract occurrence points (df)
#'
#' @keywords internal
.pe_get_occ_pts <- function(occ_pts, vs) {

  # user-supplied
  if (!is.null(occ_pts)) return(as.data.frame(occ_pts))

  # from vs
  if (!is.null(vs)) {
    op <- suppressWarnings(nr_get_occ(vs))
    if (!is.null(op)) return(as.data.frame(op))
  }

  return(NULL)
}

#' Helper: Resolve x,y,z columns
#'
#' If user did not supply all, infer from env_bg
#'
#' @keywords internal
.pe_resolve_xyz <- function(env_bg, x, y, z) {

  auto_inf <- FALSE

  if (missing(x) || missing(y) || missing(z)) auto_inf <- TRUE

  # Extract predictor names (remove x/y if present)
  if (all(c("x", "y") %in% names(env_bg))) {
    preds <- setdiff(names(env_bg), c("x", "y"))
  } else preds <- names(env_bg)

  if (length(preds) < 3)
    stop("Cannot infer x,y,z from env_bg (need ≥3 predictors). Provide x,y,z explicitly.")

  # Auto inference
  if (auto_inf) {
    message("Auto-inferring x,y,z from first 3 predictor columns: ",
            paste(preds[1:3], collapse=", "))
    return(preds[1:3])
  }

  return(c(x, y, z))
}

#' Helper: build 2D ellipsoid objects for each pair
#'
#' @keywords internal
.pe_build_ell2d <- function(niche) {
  list(
    y_x = build_ellps(center = niche$center[c(2,1)],
                      axes   = niche$axes[c(2,1)],
                      angles = niche$angles[c(2,1)]),
    z_x = build_ellps(center = niche$center[c(3,1)],
                      axes   = niche$axes[c(3,1)],
                      angles = niche$angles[c(3,1)]),
    z_y = build_ellps(center = niche$center[c(3,2)],
                      axes   = niche$axes[c(3,2)],
                      angles = niche$angles[c(3,2)])
  )
}



#' Plot niche and environments in E-space (2D or 3D)
#'
#' Automatically retrieves:
#'  - env_bg (environmental background),
#'  - suitable_env (inside ellipsoid),
#'  - occurrences,
#'  - x,y,z predictors
#' from: user args → vs → nr_get().
#'
#' @export
plot_e_space <- function(env_bg = NULL,
                         x, y, z,
                         labels = c("ENV 1", "ENV 2", "ENV 3"),
                         n_bg = 10000,
                         niche = NULL,
                         suitable_env = NULL,
                         occ_pts = NULL,
                         rand_seed = 1234,
                         show.occ.density = FALSE,
                         plot.3d = FALSE,
                         colors = NULL,
                         palette = "default",
                         vs = NULL) {

  ## -------------------------------------------------------------------------
  ## 0. Pull information from vs if provided
  ## -------------------------------------------------------------------------

  if (!is.null(vs)) {

    if (!inherits(vs, "NicheR_species"))
      stop("'vs' must be a NicheR_species object produced by create_virtual_species().")

    # niche
    if (is.null(niche))
      niche <- nr_get_niche(vs)

    # suitable_env
    if (is.null(suitable_env))
      suitable_env <- nr_get_suitable_df(vs)

    # occurrences
    if (is.null(occ_pts))
      occ_pts <- nr_get_occ(vs)
  }

  ## -------------------------------------------------------------------------
  ## 1. Resolve environmental background correctly
  ## -------------------------------------------------------------------------
  env_bg <- .pe_get_env_bg(env_bg, vs, suitable_env, niche)

  # Coerce tibble
  if (inherits(env_bg, "tbl_df")) env_bg <- as.data.frame(env_bg)

  # Coerce Raster → SpatRaster → df (with memory check)
  if (inherits(env_bg, "Raster")) env_bg <- terra::rast(env_bg)

  if (inherits(env_bg, "SpatRaster")) {

    ncell <- terra::ncell(env_bg)
    nlyr  <- terra::nlyr(env_bg)
    est_mb <- ncell * nlyr * 8 / 1024^2

    if (est_mb > 5000) {
      stop("env_bg raster too large to convert safely (~", round(est_mb,1),
           " MB). Convert manually:\n  as.data.frame.nicheR(env_bg)\nand pass the df.")
    }

    env_bg <- as.data.frame.nicheR(env_bg)
  }

  if (!is.data.frame(env_bg))
    stop("'env_bg' must be coercible to a data.frame.")

  ## -------------------------------------------------------------------------
  ## 2. Resolve predictor columns x,y,z
  ## -------------------------------------------------------------------------
  xyz <- .pe_resolve_xyz(env_bg, x, y, z)
  col_x <- xyz[1]; col_y <- xyz[2]; col_z <- xyz[3]

  ## -------------------------------------------------------------------------
  ## 3. Retrieve inside points and occurrences
  ## -------------------------------------------------------------------------
  pts_in <- .pe_get_suitable_pts(suitable_env)
  occ_pts <- .pe_get_occ_pts(occ_pts, vs)

  ## -------------------------------------------------------------------------
  ## 4. Colors and palettes
  ## -------------------------------------------------------------------------
  palettes <- list(
    default = list(
      bg           = "#9093A2FF",
      ellipsoid    = "#2A363BFF",
      centroid     = "#D72000FF",
      tolerance    = "#EE6100FF",
      suitable_env = "#FED789FF",
      occ          = "#B4BF3AFF"
    )
  )

  if (!palette %in% names(palettes))
    stop("Unknown palette '", palette, "'.")

  base_colors <- palettes[[palette]]

  if (is.null(colors)) {
    colors <- base_colors
  } else {
    if (!is.list(colors)) colors <- as.list(colors)
    colors <- utils::modifyList(base_colors, colors)
  }

  ## -------------------------------------------------------------------------
  ## 5. Legend metadata
  ## -------------------------------------------------------------------------

  # background label changes depending on source
  bg_label <- "Background environments"
  if (missing(env_bg) && !is.null(vs)) {
    bg_label <- "Background (from vs)"
  }

  legend_items <- data.frame(
    id    = c("bg","ell","cen","tol","suit","occ"),
    type  = c("point","line","point","line","point","point"),
    label = c(
      bg_label,
      "Niche boundary",
      "Centroid",
      "Tolerance range",
      "Suitable environments",
      "Occurrences"
    ),
    color = c(
      colors$bg,
      colors$ellipsoid,
      colors$centroid,
      colors$tolerance,
      colors$suitable_env,
      colors$occ
    ),
    linetype = c(NA,1,NA,2,NA,NA),
    stringsAsFactors = FALSE
  )

  # deactivate unavailable parts
  active <- c(
    TRUE,
    !is.null(niche),
    !is.null(niche),
    !is.null(niche),
    !is.null(pts_in),
    !is.null(occ_pts)
  )
  legend_items <- legend_items[active, , drop=FALSE]

  # legend plot
  if (nrow(legend_items)==0) {
    legend_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  } else {

    top_y <- 2
    spacing <- 0.25

    legend_items$y <- top_y - (seq_len(nrow(legend_items))-1)*spacing

    legend_plot <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim=c(-0.5,0.5), ylim=c(0,2.5), clip="off") +
      ggplot2::theme_void()

    # points
    pt <- legend_items[legend_items$type=="point",]
    if (nrow(pt)>0) {
      legend_plot <- legend_plot +
        ggplot2::geom_point(
          data=pt,
          aes(x=0, y=y, colour=color),
          size=2
        )
    }

    # lines
    ln <- legend_items[legend_items$type=="line",]
    if (nrow(ln)>0) {
      legend_plot <- legend_plot +
        ggplot2::geom_segment(
          data=ln,
          aes(x=-0.2, xend=0, y=y, yend=y,
              colour=color, linetype=linetype),
          linewidth=0.5
        )
    }

    legend_plot <- legend_plot +
      ggplot2::geom_text(
        data=legend_items,
        aes(x=0.1, y=y, label=label),
        hjust=0
      ) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_linetype_identity()
  }

  ## -------------------------------------------------------------------------
  ## 6. Downsample environmental background
  ## -------------------------------------------------------------------------
  if (nrow(env_bg) > n_bg) {
    set.seed(rand_seed)
    env_bg <- env_bg[sample.int(nrow(env_bg), n_bg), ]
  }

  ## -------------------------------------------------------------------------
  ## 7. 3D PLOTTING BRANCH
  ## -------------------------------------------------------------------------
  if (plot.3d) {

    p3 <- plotly::plot_ly(
      data = env_bg,
      x = env_bg[[col_x]],
      y = env_bg[[col_y]],
      z = env_bg[[col_z]],
      type = "scatter3d",
      mode = "markers",
      marker = list(color=colors$bg, size=2),
      name="Background"
    ) %>%
      plotly::layout(
        scene=list(
          xaxis=list(title=labels[1]),
          yaxis=list(title=labels[2]),
          zaxis=list(title=labels[3])
        )
      )

    # Add ellipsoid surface
    if (!is.null(niche) && !is.null(niche$surface)) {
      surf <- as.data.frame(niche$surface)
      if (ncol(surf)>=3) {
        p3 <- p3 %>%
          plotly::add_trace(
            data=surf, x=surf[[1]], y=surf[[2]], z=surf[[3]],
            type="scatter3d", mode="lines",
            line=list(color=colors$ellipsoid),
            name="Boundary"
          )
      }

      # centroid
      if (!is.null(niche$center)) {
        p3 <- p3 %>%
          plotly::add_markers(
            x=niche$center[1], y=niche$center[2], z=niche$center[3],
            marker=list(size=5, color=colors$centroid),
            name="Centroid"
          )
      }

      # suitable
      if (!is.null(pts_in)) {
        p3 <- p3 %>%
          plotly::add_markers(
            data=pts_in,
            x=pts_in[[col_x]], y=pts_in[[col_y]], z=pts_in[[col_z]],
            marker=list(size=3, color=colors$suitable_env),
            name="Suitable"
          )
      }
    }

    # occurrences
    if (!is.null(occ_pts)) {
      p3 <- p3 %>%
        plotly::add_markers(
          data=occ_pts,
          x=occ_pts[[col_x]], y=occ_pts[[col_y]], z=occ_pts[[col_z]],
          marker=list(size=4, color=colors$occ),
          name="Occurrences"
        )
    }

    return(p3)
  }

  ## -------------------------------------------------------------------------
  ## 8. BUILD GRID OF MAIN 2D PANELS
  ## -------------------------------------------------------------------------

  make_scatter <- function(df, xcol, ycol, color) {
    ggplot2::ggplot(df, aes(x=.data[[xcol]], y=.data[[ycol]])) +
      ggplot2::geom_point(alpha=0.5, color=color, pch=".") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title=element_blank())
  }

  p_yx <- make_scatter(env_bg, col_y, col_x, colors$bg)
  p_zx <- make_scatter(env_bg, col_z, col_x, colors$bg)
  p_zy <- make_scatter(env_bg, col_z, col_y, colors$bg)

  # axis labels
  x_name <- ggplot2::ggplot() + ggplot2::theme_void() +
    ggplot2::geom_text(aes(0,0,label=labels[1]), fontface="bold")
  y_name <- ggplot2::ggplot() + ggplot2::theme_void() +
    ggplot2::geom_text(aes(0,0,label=labels[2]), fontface="bold")
  z_name <- ggplot2::ggplot() + ggplot2::theme_void() +
    ggplot2::geom_text(aes(0,0,label=labels[3]), fontface="bold")

  ## -------------------------------------------------------------------------
  ## 9. Add ellipsoid + suitable + occurrences
  ## -------------------------------------------------------------------------

  if (!is.null(niche)) {

    ell2d <- .pe_build_ell2d(niche)

    # suitable
    if (!is.null(pts_in)) {
      p_yx <- p_yx + geom_point(data=pts_in, aes(x=.data[[col_y]], y=.data[[col_x]]),
                                color=colors$suitable_env, pch=".")
      p_zx <- p_zx + geom_point(data=pts_in, aes(x=.data[[col_z]], y=.data[[col_x]]),
                                color=colors$suitable_env, pch=".")
      p_zy <- p_zy + geom_point(data=pts_in, aes(x=.data[[col_z]], y=.data[[col_y]]),
                                color=colors$suitable_env, pch=".")
    }

    # occurrences
    if (!is.null(occ_pts)) {
      p_yx <- p_yx + geom_point(data=occ_pts, aes(x=.data[[col_y]], y=.data[[col_x]]),
                                color=colors$occ, pch=".")
      p_zx <- p_zx + geom_point(data=occ_pts, aes(x=.data[[col_z]], y=.data[[col_x]]),
                                color=colors$occ, pch=".")
      p_zy <- p_zy + geom_point(data=occ_pts, aes(x=.data[[col_z]], y=.data[[col_y]]),
                                color=colors$occ, pch=".")
    }

    # ellipsoid boundary lines
    add_ell <- function(p, ell, col1=colors$ellipsoid, col2=colors$tolerance, center_col=colors$centroid) {

      p +
        geom_path(data=ell$surface, aes(x,y), color=col1, linewidth=0.5) +
        annotate("segment",
                 x=ell$center[1]-ell$axes[1], xend=ell$center[1]+ell$axes[1],
                 y=ell$center[2], yend=ell$center[2],
                 color=col2, linetype="dashed") +
        annotate("segment",
                 y=ell$center[2]-ell$axes[2], yend=ell$center[2]+ell$axes[2],
                 x=ell$center[1], xend=ell$center[1],
                 color=col2, linetype="dashed") +
        annotate("point", x=ell$center[1], y=ell$center[2],
                 color=center_col, size=2)
    }

    p_yx <- add_ell(p_yx, ell2d$y_x)
    p_zx <- add_ell(p_zx, ell2d$z_x)
    p_zy <- add_ell(p_zy, ell2d$z_y)
  }

  ## -------------------------------------------------------------------------
  ## 10. Assemble the main 3x3 grid
  ## -------------------------------------------------------------------------
  main_plot <- ggpubr::ggarrange(
    x_name,     p_yx,       p_zx,
    legend_plot,y_name,     p_zy,
    NULL,       NULL,       z_name,
    ncol=3, nrow=3,
    widths = c(0.12,0.44,0.44),
    heights = c(0.44,0.44,0.12)
  )

  ## -------------------------------------------------------------------------
  ## OPTIONAL: DENSITY PANELS (Option A)
  ## -------------------------------------------------------------------------
  if (show.occ.density && !is.null(occ_pts)) {

    rngx <- range(env_bg[[col_x]], na.rm=TRUE)
    rngy <- range(env_bg[[col_y]], na.rm=TRUE)
    rngz <- range(env_bg[[col_z]], na.rm=TRUE)

    dens_z <- ggplot2::ggplot(occ_pts, aes(x=.data[[col_z]])) +
      geom_density(fill=colors$occ, alpha=0.6) +
      scale_x_continuous(limits=rngz) +
      scale_y_continuous(n.breaks=3) +
      theme_minimal() +
      theme(axis.text.x=element_blank(),
            axis.title=element_blank())

    dens_y <- ggplot2::ggplot(occ_pts, aes(x=.data[[col_y]])) +
      geom_density(fill=colors$occ, alpha=0.6) +
      scale_x_continuous(limits=rngy) +
      scale_y_continuous(n.breaks=3) +
      theme_minimal() +
      theme(axis.text.x=element_blank(),
            axis.title=element_blank())

    dens_x_right <- ggplot2::ggplot(occ_pts, aes(x=.data[[col_x]])) +
      geom_density(fill=colors$occ, alpha=0.6) +
      scale_x_continuous(limits=rngx) +
      coord_flip() +
      scale_y_continuous(n.breaks=3) +
      theme_minimal() +
      theme(axis.text.y=element_blank(),
            axis.title=element_blank())

    dens_y_right <- ggplot2::ggplot(occ_pts, aes(x=.data[[col_y]])) +
      geom_density(fill=colors$occ, alpha=0.6) +
      scale_x_continuous(limits=rngy) +
      coord_flip() +
      scale_y_continuous(n.breaks=3) +
      theme_minimal() +
      theme(axis.text.y=element_blank(),
            axis.title=element_blank())

    # layout 4x4
    full_plot <- ggpubr::ggarrange(
      NULL,       dens_y,    dens_z,    NULL,
      x_name,     p_yx,      p_zx,      dens_x_right,
      legend_plot,y_name,    p_zy,      dens_y_right,
      NULL,       NULL,      z_name,    NULL,
      ncol=4, nrow=4,
      widths  = c(0.1,0.4,0.4,0.1),
      heights = c(0.1,0.4,0.4,0.1)
    )

    return(full_plot)
  }

  return(main_plot)
}
