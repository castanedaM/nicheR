# Title: Helper functions

# Description: Functions shared across the server scripts. Grouped by the
# stage of the workflow they belong to. Several read `input`, `session`, or
# `session_data` from the enclosing environment and only work because this
# file is sourced with local = TRUE from server.R; those are marked below.

# Date last updated: 10/06/2026


#' Pairwise covariances as a flat named vector
#'
#' Values are ordered to match the covariance sliders.
#'
#' @param ell A nicheR_ellipsoid object, or NULL.
#'
#' @returns A named numeric vector of upper-triangle covariances, or NULL.
#'
#' @noRd
cov_upper <- function(ell){
  if(is.null(ell)) return(NULL)
  pairs <- t(combn(ell$var_names, 2))
  setNames(sapply(seq_len(nrow(pairs)),
                  function(i) ell$cov_matrix[pairs[i, 1], pairs[i, 2]]),
           apply(pairs, 1, paste, collapse = "-"))
}


# FILE READERS ------------------------------------------------------------

#' Read a raster file and return a SpatRaster
#'
#' @param path File path to read from.
#' @param ext Lowercase file extension, one of tif, tiff, or rds.
#'
#' @returns A SpatRaster.
#'
#' @noRd
load_raster_file <- function(path, ext){

  obj <- switch(ext,
                "tif" = terra::rast(path),
                "tiff" = terra::rast(path),
                "rds" = readRDS(path),
                stop("Unsupported file type: ", ext))

  if(!inherits(obj, c("SpatRaster", "RasterStack", "RasterBrick", "RasterLayer"))){
    stop("File does not contain a raster object. Found class: ",
         paste(class(obj), collapse = ", "))
  }

  # Convert raster package objects to terra if needed
  if(inherits(obj, c("RasterStack", "RasterBrick", "RasterLayer"))){
    obj <- terra::rast(obj)
  }

  obj
}


#' Read a tabular file and return a data frame
#'
#' @param path File path to read from.
#' @param ext Lowercase file extension, one of csv or rds.
#'
#' @returns A data frame.
#'
#' @noRd
load_df_file <- function(path, ext){

  obj <- switch(ext,
                "csv" = read.csv(path),
                "rds" = readRDS(path),
                stop("Unsupported file type: ", ext))

  if(!inherits(obj, c("data.frame", "tbl_df", "tbl"))){
    stop("File does not contain a data frame. Found class: ",
         paste(class(obj), collapse = ", "))
  }

  as.data.frame(obj)
}


# ELLIPSOID IDENTITY ------------------------------------------------------

#' Session-level ellipsoid counter
#'
#' Incremented every time an ellipsoid is tagged. Saved and restored with
#' the session so ids never collide after a reload.
#'
#' @noRd
ell_id_counter <- reactiveVal(0L)

#' Generate the next ellipsoid id
#'
#' Increments the session counter as a side effect.
#'
#' @returns A character id of the form "E3_060826".
#'
#' @noRd
make_ell_id <- function(){
  n <- ell_id_counter() + 1L
  ell_id_counter(n)
  paste0("E", n, "_", format(Sys.time(), "%d%m%y"))
}


#' Give an ellipsoid a fresh id and name
#'
#' The default name is read after make_ell_id() has already incremented the
#' counter, so the number in the name matches the number in the id.
#'
#' @param ell A nicheR_ellipsoid object.
#' @param name Optional name. Defaults to "ellipsoid_<n>".
#'
#' @returns The ellipsoid with ell_id and ell_name set.
#'
#' @noRd
tag_ellipsoid <- function(ell, name = NULL){
  ell$ell_id <- make_ell_id()
  ell$ell_name <- if(is.null(name)) paste0("ellipsoid_", ell_id_counter()) else name
  ell
}


# ELLIPSOID LINEAGE -------------------------------------------------------

#' Resolve an ellipsoid's parent
#'
#' Reads session_data from the server environment. Note that list[[NULL]]
#' errors, so parent_id has to be checked before indexing.
#'
#' @param ell A nicheR_ellipsoid object, or NULL.
#'
#' @returns The parent ellipsoid, or NULL if it has none or the parent has
#'   been deleted.
#'
#' @noRd
ell_parent <- function(ell){

  if(is.null(ell)) return(NULL)

  pid <- ell$parent_id
  if(is.null(pid) || !nzchar(pid)) return(NULL)
  if(!pid %in% names(session_data$ellipsoid_list)) return(NULL)

  session_data$ellipsoid_list[[pid]]
}


#' Find the state a reset returns to
#'
#' A copy resets to its parent; an ellipsoid with no parent resets to what
#' its own ranges produce, which is zero covariance and the centroid at the
#' range midpoint. Reads session_data from the server environment.
#'
#' @param ell A nicheR_ellipsoid object.
#'
#' @returns A nicheR_ellipsoid object, or NULL if neither source is available.
#'
#' @noRd
ell_reset_target <- function(ell){

  parent <- ell_parent(ell)
  if(!is.null(parent)) return(parent)

  inputs <- ell$range_inputs
  if(is.null(inputs)) return(NULL)

  range_df <- as.data.frame(rbind(unlist(inputs$min), unlist(inputs$max)),
                            row.names = c("min", "max"))

  tryCatch(build_ellipsoid(range = range_df, cl = ell$cl, verbose = FALSE),
           error = function(e) NULL)
}


#' Label for the centroid reset link
#'
#' @param ell A nicheR_ellipsoid object.
#'
#' @returns A character label naming the parent, or a generic fallback.
#'
#' @noRd
ell_reset_label <- function(ell){
  parent <- ell_parent(ell)
  if(!is.null(parent)) paste0("Reset to ", parent$ell_name) else "Reset to original"
}


#' Lineage label shown in the ellipsoid library
#'
#' @param ell A nicheR_ellipsoid object.
#'
#' @returns "from <parent name>", or "root" when it has no parent.
#'
#' @noRd
ell_lineage_label <- function(ell){
  parent <- ell_parent(ell)
  if(!is.null(parent)) paste0("from ", parent$ell_name) else "root"
}


# BUILD TAB ---------------------------------------------------------------

#' Default values for the range panel
#'
#' Prefers the stored range inputs of the working ellipsoid, then quartiles
#' of the background data, then a plain 0 to 1 scale for virtual mode. Reads
#' session_data from the server environment.
#'
#' @param vars Character vector of variable names.
#' @param has_bg_df Logical, whether background data is available.
#'
#' @returns A list with min, max, mean, sd, expand_min, and expand_max, each
#'   a named list over `vars`.
#'
#' @noRd
build_range_defaults <- function(vars, has_bg_df){

  if(!is.null(session_data$session_range)){
    return(session_data$session_range)
  }

  zeros <- setNames(as.list(rep(0, length(vars))), vars)
  ones <- setNames(as.list(rep(1, length(vars))), vars)

  if(has_bg_df){
    q <- round(apply(session_data$bg_df[, vars, drop = FALSE], 2,
                     quantile, na.rm = TRUE), 2)

    return(list(min = as.list(q[2, ]),
                max = as.list(q[4, ]),
                mean = as.list((q[2, ] + q[4, ]) / 2),
                sd = as.list((q[4, ] - q[2, ]) / 4),
                expand_min = zeros,
                expand_max = zeros))
  }

  list(min = zeros,
       max = ones,
       mean = zeros,
       sd = ones,
       expand_min = zeros,
       expand_max = zeros)
}


# BIAS TAB ----------------------------------------------------------------

#' Apply bias to one layer and merge the result into the biased list
#'
#' If ell_id already exists, the new layer is appended to the existing
#' SpatRaster stack. If the layer is already in the stack, it is skipped.
#'
#' @param biased_list Named list of SpatRasters keyed by ellipsoid id.
#' @param ell_id Ellipsoid id to merge into.
#' @param pred_rast The unbiased prediction SpatRaster.
#' @param layer Name of the prediction layer to bias.
#' @param prepared_bias Output of prepare_bias().
#' @param direction Either "direct" or "inverse".
#'
#' @returns The biased list with the new layer merged in, unchanged on failure.
#'
#' @noRd
apply_bias_to_list <- function(biased_list, ell_id, pred_rast,
                               layer, prepared_bias, direction = "direct"){

  result <- tryCatch(
    apply_bias(prepared_bias = prepared_bias,
               prediction = pred_rast,
               prediction_layer = layer,
               effect_direction = direction,
               verbose = FALSE),
    error = function(e) NULL
  )

  if(is.null(result)) return(biased_list)

  # Extract the SpatRaster from the result, dropping combination_formula
  new_rast <- result[[paste0(layer, "_biased")]]
  if(is.null(new_rast) || !inherits(new_rast, "SpatRaster")) return(biased_list)

  if(!ell_id %in% names(biased_list)){
    biased_list[[ell_id]] <- new_rast
  } else {
    existing <- biased_list[[ell_id]]
    new_lyr_nm <- names(new_rast)

    if(new_lyr_nm %in% names(existing)){
      message("Layer '", new_lyr_nm, "' already exists for ", ell_id, ". Skipping.")
    } else {
      biased_list[[ell_id]] <- c(existing, new_rast)
    }
  }

  biased_list
}


# GENERATE TAB ------------------------------------------------------------


#' Parameters that define an occurrence set
#'
#' Two sets are the same set only if all of these match.
#'
#' @noRd
occ_set_fields <- c("mode", "layer", "n_occ", "seed", "sampling",
                    "strict", "mask", "truncate", "effect")

#' Read one metadata attribute from an occurrence set
#'
#' @param df An occurrence set data frame.
#' @param field Attribute name, normally one of occ_set_fields.
#' @param default Value returned when the attribute is absent.
#'
#' @returns The attribute value, or `default`.
#'
#' @noRd
occ_meta <- function(df, field, default = NA){
  val <- attr(df, field, exact = TRUE)
  if(is.null(val)) default else val
}


#' Sampling label for sets drawn from a biased surface
#'
#' The surface values are the weights, so no centroid or edge strategy
#' applies. The direction used when the bias was applied is already part of
#' the layer name, so it is not repeated here.
#'
#' @noRd
occ_sampling_biased <- "bias weighted"

#' Sampling labels for the virtual_data() effects
#'
#' Virtual sets keep the effect passed to virtual_data(). Their sampling
#' value uses the same words as the raster sets, so both modes read alike in
#' the summary and in the downloads.
#'
#' @noRd
virtual_sampling_labels <- c(direct = "centroid",
                             inverse = "edge",
                             uniform = "random")

#' Sampling label for an occurrence set
#'
#' Returns the stored value when there is one. Sets generated before biased
#' and virtual sets recorded a sampling value fall back to the "biased"
#' attribute or the virtual effect, so older sessions still read correctly.
#'
#' @param df An occurrence set data frame.
#'
#' @returns A single character label, or NA when nothing is recorded.
#'
#' @noRd
occ_sampling <- function(df){

  val <- occ_meta(df, "sampling")
  if(!is.na(val)) return(as.character(val))

  if(isTRUE(attr(df, "biased", exact = TRUE))) return(occ_sampling_biased)

  effect <- as.character(occ_meta(df, "effect"))
  if(effect %in% names(virtual_sampling_labels)){
    return(unname(virtual_sampling_labels[effect]))
  }

  NA_character_
}


#' Canonical key for an occurrence set
#'
#' Never shown to the user, only used for identity. Any parameter change
#' produces a different key, so sets accumulate rather than overwrite.
#'
#' @param attrs Named list covering every entry in occ_set_fields.
#'
#' @returns A single character signature.
#'
#' @noRd
occ_set_signature <- function(attrs){
  paste(vapply(occ_set_fields,
               function(f) paste0(f, "=", as.character(attrs[[f]])),
               character(1)),
        collapse = "|")
}


#' Short display labels for a list of occurrence sets
#'
#' Always names the layer, then adds only the fields that vary across the
#' sets present, so a run differing only by seed reads
#' "suitability_trunc, seed 123" rather than repeating shared parameters.
#'
#' Virtual sets also always show sampling and truncation. Their layer is
#' just "virtual", so unlike a raster layer name it says nothing about how
#' the points were drawn.
#'
#' @param occ_list Named list of occurrence set data frames.
#'
#' @returns A named character vector: names are labels, values are keys.
#'
#' @noRd
occ_set_labels <- function(occ_list){

  if(length(occ_list) == 0) return(character(0))

  meta <- lapply(occ_list, function(df){
    m <- setNames(lapply(occ_set_fields,
                         function(f) as.character(occ_meta(df, f, "NA"))),
                  occ_set_fields)

    # Same label the downloads use, including the fallback for older sets
    samp <- occ_sampling(df)
    m$sampling <- if(is.na(samp)) "NA" else samp
    m
  })

  varies <- vapply(occ_set_fields, function(f){
    length(unique(vapply(meta, function(m) m[[f]], character(1)))) > 1
  }, logical(1))

  # Effect is left out because sampling already says the same thing
  show <- setdiff(unique(c("layer", occ_set_fields[varies])), "effect")

  labs <- vapply(meta, function(m){

    show_set <- if(identical(m$mode, "virtual")){
      unique(c("layer", "sampling", "truncate", show))
    } else {
      show
    }

    parts <- vapply(show_set, function(f){
      v <- m[[f]]
      switch(f,
             "layer" = v,
             "n_occ" = paste0("n = ", v),
             "seed" = paste0("seed ", v),
             "sampling" = v,
             "strict" = if(identical(v, "TRUE")) "strict" else "not strict",
             "mask" = if(identical(v, "none")) "no mask" else paste0("mask: ", v),
             "truncate" = if(identical(v, "TRUE")){
               "truncated"
             } else if(identical(v, "FALSE")){
               "not truncated"
             } else {
               v
             },
             v)
    }, character(1))
    paste(parts, collapse = ", ")
  }, character(1))


  # Names are what the user sees, values are the keys
  setNames(names(occ_list), labs)
}


#' Sample occurrences for one ellipsoid across several layers
#'
#' Each layer is pulled from either the prediction or the biased prediction.
#' Biased layers go to sample_biased_data(), which already carries the
#' sampling weights in the surface itself and so takes no sampling strategy
#' or method. Unbiased layers go to sample_data().
#'
#' @param ell_id Ellipsoid id.
#' @param pred_list Named list of unbiased prediction SpatRasters.
#' @param biased_list Named list of biased prediction SpatRasters.
#' @param layers Character vector of layer names to sample from.
#' @param n_occ Number of occurrences to draw per layer.
#' @param sampling One of "centroid", "edge", or "random". Unbiased only.
#' @param strict Logical, or NULL to let the function auto-detect.
#' @param sampling_mask Optional SpatRaster restricting where points fall.
#' @param seed Integer random seed.
#'
#' @returns A named list of data frames, one per layer that succeeded. Each
#'   carries a "biased" attribute recording which function produced it.
#'
#' @noRd
generate_occ_for_ell <- function(ell_id, pred_list, biased_list,
                                 layers, n_occ, sampling,
                                 strict, sampling_mask, seed = 123L){

  results <- list()

  for(layer in layers){

    pred <- pred_list[[ell_id]]
    bias <- biased_list[[ell_id]]

    # Biased list checked first, since a biased layer name never appears
    # in the unbiased stack
    is_biased <- !is.null(bias) && inherits(bias, "SpatRaster") &&
      layer %in% names(bias)

    source_rast <- if(is_biased){
      bias
    } else if(!is.null(pred) && inherits(pred, "SpatRaster") &&
              layer %in% names(pred)){
      pred
    } else {
      NULL
    }

    if(is.null(source_rast)){
      message("Layer '", layer, "' not found for ", ell_id, ". Skipping.")
      next
    }

    occ <- if(is_biased){

      tryCatch(
        sample_biased_data(n_occ = n_occ,
                           prediction = source_rast,
                           prediction_layer = layer,
                           sampling_mask = sampling_mask,
                           seed = seed,
                           strict = strict,
                           verbose = FALSE),
        error = function(e){
          message("Biased generate failed for ", ell_id, " (", layer, "): ", e$message)
          NULL
        }
      )

    } else {

      method <- if(grepl("mahalanobis", layer, ignore.case = TRUE)){
        "mahalanobis"
      } else {
        "suitability"
      }

      tryCatch(
        sample_data(n_occ = n_occ,
                    prediction = source_rast,
                    prediction_layer = layer,
                    sampling = sampling,
                    method = method,
                    sampling_mask = sampling_mask,
                    seed = seed,
                    strict = strict,
                    verbose = FALSE),
        error = function(e){
          message("Generate failed for ", ell_id, " (", layer, "): ", e$message)
          NULL
        }
      )
    }

    if(!is.null(occ)){
      out <- occ[, c("x", "y"), drop = FALSE]
      attr(out, "biased") <- is_biased
      results[[layer]] <- out
    }
  }

  results
}


#' Sample virtual occurrences directly from an ellipsoid
#'
#' Virtual mode has no prediction surface and no geography, so points come
#' from virtual_data(), which draws from the multivariate normal defined by
#' the ellipsoid's centroid and covariance. The result is environmental
#' values only, with no coordinates.
#'
#' @param ell A nicheR_ellipsoid object.
#' @param n Number of points to generate.
#' @param truncate Logical, constrain points within the confidence limit.
#' @param effect One of "direct", "inverse", or "uniform". The latter two
#'   require truncate = TRUE.
#' @param seed Integer random seed.
#'
#' @returns A data frame of environmental values, or NULL on failure.
#'
#' @noRd
generate_virtual_for_ell <- function(ell, n, truncate, effect, seed = 123L){

  # virtual_data enforces this too, but failing here gives a clearer message
  if(effect %in% c("inverse", "uniform") && !isTRUE(truncate)){
    message("Effect '", effect, "' requires truncation.")
    return(NULL)
  }

  res <- tryCatch(
    virtual_data(object = ell,
                 n = n,
                 truncate = truncate,
                 effect = effect,
                 seed = seed),
    error = function(e){
      message("Virtual generate failed: ", e$message)
      NULL
    }
  )

  if(is.null(res)) return(NULL)

  as.data.frame(res)
}

# PLOTTING ----------------------------------------------------------------

#' Read a plot input with a NULL-safe default
#'
#' Reads `input` from the server environment.
#'
#' @param id Input id.
#' @param default Value returned when the input has not reported yet.
#'
#' @returns The input value, or `default`.
#'
#' @noRd
get_input <- function(id, default){
  if(!is.null(input[[id]])){
    input[[id]]
  } else {
    default
  }
}


#' Keep two axis selectors mutually exclusive
#'
#' Stops the same variable being chosen on both axes of a 2D plot. Reads
#' `input` and `session` from the server environment.
#'
#' @param x_id Input id of the x-axis selector.
#' @param y_id Input id of the y-axis selector.
#' @param vars Character vector of available variables.
#'
#' @returns Invisibly NULL. Called for the update side effect.
#'
#' @noRd
update_axis_selectors <- function(x_id, y_id, vars){

  x_sel <- input[[x_id]]
  y_sel <- input[[y_id]]

  if(is.null(x_sel) || !x_sel %in% vars) x_sel <- vars[1]

  y_choices <- setdiff(vars, x_sel)
  if(is.null(y_sel) || !y_sel %in% y_choices) y_sel <- y_choices[1]

  x_choices <- setdiff(vars, y_sel)

  updateSelectInput(session, x_id, choices = x_choices, selected = x_sel)
  updateSelectInput(session, y_id, choices = y_choices, selected = y_sel)
}


#' Continuous palette for prediction layers
#'
#' Uses base R's hcl.colors, so no extra dependency is needed. An unknown
#' palette name falls back to viridis rather than erroring, which matters
#' when loading a session saved under a different R version.
#'
#' @param name Palette name, as listed by hcl.pals().
#' @param reverse Logical, reverse the ramp.
#' @param n Number of colors.
#'
#' @returns A character vector of colors.
#'
#' @noRd
pred_palette <- function(name = "viridis", reverse = FALSE, n = 100){
  if(!name %in% hcl.pals()) name <- "viridis"
  cols <- grDevices::hcl.colors(n, palette = name)
  if(isTRUE(reverse)) rev(cols) else cols
}

#' Map numeric values onto palette colors
#'
#' @param vals Numeric vector.
#' @param pal Character vector of colors.
#' @param rng Optional range to scale against, so several panels can share
#'   one scale. Defaults to the range of `vals`.
#'
#' @returns A character vector of colors the same length as `vals`.
#'
#' @noRd
pred_colors <- function(vals, pal, rng = NULL){
  if(is.null(rng)) rng <- range(vals, na.rm = TRUE)
  if(!all(is.finite(rng)) || diff(rng) == 0) return(rep(pal[1], length(vals)))
  idx <- cut(vals, breaks = seq(rng[1], rng[2], length.out = length(pal) + 1),
             labels = FALSE, include.lowest = TRUE)
  pal[idx]
}

#' Value range for an overlay layer
#'
#' @param vals_df Data frame of extracted prediction values.
#' @param layer Layer name.
#'
#' @returns A length-two numeric range, or NULL when unavailable.
#'
#' @noRd
pred_layer_range <- function(vals_df, layer){
  if(is.null(vals_df) || is.null(layer)) return(NULL)
  if(!layer %in% names(vals_df)) return(NULL)
  rng <- range(vals_df[[layer]], na.rm = TRUE)
  if(!all(is.finite(rng))) return(NULL)
  rng
}

#' Draw a shared color bar in its own panel
#'
#' Called after a plot grid so every panel shares one legend.
#'
#' @param rng Length-two numeric range the bar spans.
#' @param pal Character vector of colors.
#' @param label Optional title above the bar.
#'
#' @returns Invisibly NULL. Called for the plotting side effect.
#'
#' @noRd
pred_legend_panel <- function(rng, pal, label = ""){

  if(is.null(rng) || !all(is.finite(rng))) return(invisible(NULL))

  old_mar <- par("mar")
  on.exit(par(mar = old_mar))
  par(mar = c(2.5, 4, 0.5, 4))

  plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "")

  xs <- seq(0.15, 0.85, length.out = length(pal) + 1)
  rect(xs[-length(xs)], 0.45, xs[-1], 0.8, col = pal, border = NA)
  rect(0.15, 0.45, 0.85, 0.8, border = "#888", lwd = 0.5)

  ticks <- seq(0.15, 0.85, length.out = 5)
  vals <- seq(rng[1], rng[2], length.out = 5)
  segments(ticks, 0.45, ticks, 0.38, col = "#888", lwd = 0.5)
  text(ticks, 0.32, format(round(vals, 2)), cex = 0.65, adj = c(0.5, 1))

  if(nzchar(label)) text(0.5, 0.9, label, cex = 0.75, adj = c(0.5, 0))
}

#' Fill a 2x2 block with panel ids
#'
#' Spans cells when there are fewer than four panels: one fills the block,
#' two give two rows, three or four give a 2x2. Used to lay out each half of
#' a combined plot.
#'
#' @param ids Integer vector of up to four layout panel numbers.
#'
#' @returns A 2x2 integer matrix for layout().
#'
#' @noRd
plot_panel_block <- function(ids){

  n <- length(ids)

  if(n <= 0) return(matrix(0L, 2, 2))
  if(n == 1) return(matrix(ids[1], 2, 2))
  if(n == 2) return(matrix(c(ids[1], ids[1],
                             ids[2], ids[2]), 2, 2, byrow = TRUE))
  if(n == 3) return(matrix(c(ids[1], ids[2],
                             ids[3], 0L), 2, 2, byrow = TRUE))

  matrix(ids[1:4], 2, 2, byrow = TRUE)
}


# EXPORT ------------------------------------------------------------------

#' Open a graphics device for a figure export
#'
#' Shared by every tab's export handler. The caller is responsible for
#' closing the device, normally with on.exit(dev.off()).
#'
#' @param file Output path supplied by downloadHandler.
#' @param ext One of "png", "pdf", or "svg".
#' @param w_val Width in `unit`.
#' @param h_val Height in `unit`.
#' @param unit One of "mm", "in", or "px".
#' @param res Resolution in dpi, used for png and for px conversion.
#' @param cex_val Base text scaling applied to every cex parameter.
#'
#' @returns Invisibly NULL. Called for the device side effect.
#'
#' @noRd
plot_open_device <- function(file, ext, w_val, h_val, unit, res, cex_val){

  if(is.null(unit)) unit <- "mm"
  if(is.null(w_val)) w_val <- 166
  if(is.null(h_val)) h_val <- 166
  if(is.null(res)) res <- 300
  if(is.null(cex_val)) cex_val <- 1

  to_inches <- function(val){
    switch(unit, "mm" = val / 25.4, "in" = val, "px" = val / res)
  }
  to_px <- function(val){
    switch(unit,
           "mm" = round(val / 25.4 * res),
           "in" = round(val * res),
           "px" = round(val))
  }

  if(ext == "png"){
    png(file, width = to_px(w_val), height = to_px(h_val), res = res)
  } else if(ext == "pdf"){
    pdf(file, width = to_inches(w_val), height = to_inches(h_val))
  } else if(ext == "svg"){
    svg(file, width = to_inches(w_val), height = to_inches(h_val))
  }

  # Set every cex parameter so individual plot calls do not override them
  par(cex = cex_val,
      cex.axis = cex_val,
      cex.lab = cex_val,
      cex.main = cex_val * 1.1,
      cex.sub = cex_val * 0.9)
}
