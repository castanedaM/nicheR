#' Create a virtual species end to end
#'
#' Orchestrates the full NicheR workflow by routing `...` to
#' [build_ellps()], [get_suitable_env()], [set_bias_surface()], and
#' [get_sample_occ()]. Returns an S3 object of class **NicheR_species**
#' containing the niche object, suitability, (optional) bias surface,
#' and sampled occurrences.
#'
#' @section Workflow:
#' 1. Calls [build_ellps()] to define the ellipsoid in E-space.
#' 2. Calls [get_suitable_env()] to compute suitability in G-space.
#' 3. Optionally calls [set_bias_surface()] to build a pooled bias raster
#'    when bias arguments are supplied.
#' 4. Calls [get_sample_occ()] to sample occurrences from the suitable
#'    (and optionally biased) environment.
#'
#' @param ... Named arguments passed to the component functions.
#'   Arguments are routed by name to the matching formal parameters of
#'   [build_ellps()], [get_suitable_env()], [set_bias_surface()], and
#'   [get_sample_occ()]. Unknown names are ignored with a warning.
#' @param out.file Logical. If `TRUE`, save the returned object to an `.rds`
#'   file in the working directory.
#' @param out.file.name Optional base name (without extension) for the saved
#'   file. If `NULL` or empty, a timestamped name is used.
#' @param verbose Logical. If `TRUE`, print progress messages.
#'
#' @details
#' If any component function requires `env_bg` and it is supplied only at the
#' top level, it is forwarded to that function automatically. If a required
#' argument (such as `env_bg`) is missing for a component, the function stops
#' with an informative error.
#'
#' If bias-related arguments (e.g. `bias_surface`, `bias_dir`, `bias_template`,
#' `bias_res`, `bias_ext`) are provided, a pooled bias surface is built via
#' [set_bias_surface()] and stored in the output. This precomputed bias object
#' is also passed into [get_sample_occ()] so that bias is not recomputed.
#'
#' @return
#' A list of class **NicheR_species** with elements:
#' \itemize{
#'   \item \code{niche}: the ellipsoid object returned by [build_ellps()].
#'   \item \code{suitability}: suitability surface from [get_suitable_env()].
#'   \item \code{bias_surface}: optional \code{nicheR_bias_surface} object
#'         returned by [set_bias_surface()] (or \code{NULL} if no bias was used).
#'   \item \code{occurrences}: sampled occurrences from [get_sample_occ()].
#'   \item \code{call_args}: the original `...` captured as a named list.
#'   \item \code{routed_args}: a list showing which args went to each function.
#'   \item \code{save_path}: file path if \code{out.file = TRUE}, otherwise `NULL`.
#' }
#'
#' @seealso [build_ellps()], [get_suitable_env()], [set_bias_surface()],
#'   [get_sample_occ()]
#'
#' @export
create_virtual_species <- function(...,
                                   out.file = FALSE,
                                   out.file.name = NULL,
                                   verbose = TRUE) {

  # capture all args once
  args <- list(...)

  # ensure required functions exist
  needed_funs <- c("build_ellps", "get_suitable_env",
                   "set_bias_surface", "get_sample_occ")
  missing_funs <- needed_funs[!vapply(needed_funs, exists, logical(1), mode = "function")]

  if (length(missing_funs)) {
    stop("These functions are not available in the current environment: ",
         paste(missing_funs, collapse = ", "))
  }

  # pull formals for routing
  f_build <- names(formals(build_ellps))
  f_suit  <- names(formals(get_suitable_env))
  f_bias  <- names(formals(set_bias_surface))
  f_occ   <- names(formals(get_sample_occ))

  # split args by target function
  args_build <- args[names(args) %in% f_build]
  args_suit  <- args[names(args) %in% f_suit]
  # bias is special – we mostly care about its arguments, but
  # we won't call set_bias_surface() purely from these; instead
  # we use them together with env_bg.
  args_bias  <- args[names(args) %in% f_bias]
  args_occ   <- args[names(args) %in% f_occ]

  # pass env_bg through if user supplied it at top level but it was not routed
  if ("env_bg" %in% names(args)) {
    if (!("env_bg" %in% names(args_suit)) && ("env_bg" %in% f_suit)) {
      args_suit$env_bg <- args$env_bg
    }
    if (!("env_bg" %in% names(args_occ)) && ("env_bg" %in% f_occ)) {
      args_occ$env_bg <- args$env_bg
    }
  }

  # check that env_bg will be available when needed for suitability
  if (("env_bg" %in% f_suit) && !("env_bg" %in% names(args_suit))) {
    stop("env_bg is required by get_suitable_env but was not supplied.")
  }

  # build ellipsoid niche
  niche_obj <- tryCatch(
    do.call(build_ellps, args_build),
    error = function(e) stop("build_ellps failed: ", e$message)
  )

  if (verbose) message("Built niche object.")

  # suitability in G space
  suit_env <- tryCatch(
    do.call(get_suitable_env, c(list(niche = niche_obj), args_suit)),
    error = function(e) stop("get_suitable_env failed: ", e$message)
  )

  if (verbose) message("Computed suitable environments.")

  # Optionally build bias surface once, if bias args were supplied -----------
  bias_obj <- NULL

  # We treat presence of 'bias_surface' in either args_bias or args_occ as the signal
  # that the user wants a bias surface.
  has_bias_surface <- "bias_surface" %in% c(names(args_bias), names(args_occ))

  if (has_bias_surface) {

    # Reconstruct a bias call using the shared arguments:
    # priority: explicit args_bias first, then fall back to args_occ
    combined_bias_args <- list()
    for (nm in f_bias) {
      if (nm %in% names(args_bias)) {
        combined_bias_args[[nm]] <- args_bias[[nm]]
      } else if (nm %in% names(args_occ)) {
        # map get_sample_occ names to set_bias_surface names where appropriate
        if (nm == "bias_surface") {
          combined_bias_args[[nm]] <- args_occ[["bias_surface"]]
        } else if (nm == "bias_dir" && "bias_dir" %in% names(args_occ)) {
          combined_bias_args[[nm]] <- args_occ[["bias_dir"]]
        } else if (nm == "template" && "bias_template" %in% names(args_occ)) {
          combined_bias_args[[nm]] <- args_occ[["bias_template"]]
        } else if (nm == "res" && "bias_res" %in% names(args_occ)) {
          combined_bias_args[[nm]] <- args_occ[["bias_res"]]
        } else if (nm == "ext" && "bias_ext" %in% names(args_occ)) {
          combined_bias_args[[nm]] <- args_occ[["bias_ext"]]
        }
      }
    }

    # If no explicit template/res/ext are supplied, try to infer from env_bg
    if (!("template" %in% names(combined_bias_args)) ||
        is.null(combined_bias_args$template)) {

      if ("env_bg" %in% names(args_suit) &&
          inherits(args_suit$env_bg, c("SpatRaster", "Raster"))) {

        combined_bias_args$template <- if (inherits(args_suit$env_bg, "Raster")) {
          terra::rast(args_suit$env_bg)
        } else {
          args_suit$env_bg
        }

      }
    }

    if (!"bias_surface" %in% names(combined_bias_args) ||
        is.null(combined_bias_args$bias_surface)) {
      warning("Bias arguments were detected but 'bias_surface' is missing; ",
              "skipping bias surface construction.", call. = FALSE)
    } else {

      bias_obj <- tryCatch(
        do.call(set_bias_surface, combined_bias_args),
        error = function(e) stop("set_bias_surface failed: ", e$message)
      )

      if (verbose) message("Constructed pooled bias surface.")

      # Now, ensure get_sample_occ uses this precomputed bias object
      args_occ$bias_surface <- bias_obj
    }
  }

  # sample occurrences --------------------------------------------------------
  occ <- tryCatch(
    do.call(get_sample_occ, c(list(niche = niche_obj), args_occ)),
    error = function(e) stop("get_sample_occ failed: ", e$message)
  )

  if (verbose) message("Sampled occurrences.")

  # warn on unused args
  used_names <- union(
    names(args_build),
    union(names(args_suit), union(names(args_occ), names(args_bias)))
  )
  unused <- setdiff(names(args), used_names)
  if (length(unused) && verbose) {
    warning("These arguments did not match any target function and were ignored: ",
            paste(unused, collapse = ", "))
  }

  # assemble S3 object
  out <- list(
    niche        = niche_obj,
    suitability  = suit_env,
    bias_surface = bias_obj,
    occurrences  = occ,
    call_args    = args,
    routed_args  = list(
      build_ellps      = args_build,
      get_suitable_env = args_suit,
      set_bias_surface = args_bias,
      get_sample_occ   = args_occ
    ),
    save_path    = NULL
  )

  class(out) <- c("NicheR_species", class(out))

  # ---- Saving logic ----
  if (isTRUE(out.file)) {

    dir_path <- getwd()

    if (is.null(out.file.name) || out.file.name == "") {
      timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      out.file.name <- paste0("NicheR_vs_", timestamp)
    }

    save_path <- file.path(dir_path, paste0(out.file.name, ".rds"))
    saveRDS(out, save_path)

    if (verbose) message("Virtual species saved to: ", normalizePath(save_path))

    out$save_path <- save_path

  } else if (verbose) {
    message("`out.file = FALSE`: object not saved to disk.")
  }

  return(out)
}

#' Print a NicheR virtual species
#'
#' @param x A \code{NicheR_species} object.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @method print NicheR_species
#' @export
print.NicheR_species <- function(x, ...) {
  cat("NicheR virtual species components:\n")
  cat("  niche:        ", paste(class(x$niche), collapse = "/"), "\n", sep = "")
  cat("  suitability:  ", paste(class(x$suitability), collapse = "/"), "\n", sep = "")
  cat("  bias_surface: ", if (is.null(x$bias_surface)) "NULL" else
    paste(class(x$bias_surface), collapse = "/"), "\n", sep = "")
  cat("  occurrences:  ", paste(class(x$occurrences), collapse = "/"), "\n", sep = "")
  invisible(x)
}
