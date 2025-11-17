#' Create a virtual species end to end
#'
#' Orchestrates the full NicheR workflow by routing \code{...} to
#' \code{\link{build_ellps}}, \code{\link{get_suitable_env}}, and
#' \code{\link{get_sample_occ}}. Returns an S3 object of class
#' \code{"NicheR_species"} containing the niche object, suitability, and
#' sampled occurrences.
#'
#' @section Workflow:
#' \enumerate{
#'   \item Calls \code{build_ellps()} to define the ellipsoid in E-space.
#'   \item Calls \code{get_suitable_env()} to compute the suitable environment
#'         in G-space (always requesting distances and a data.frame output).
#'   \item Calls \code{get_sample_occ()} to sample occurrences from the
#'         suitable environment pool, passing \code{suitable_env} automatically
#'         unless the user has supplied it explicitly.
#' }
#'
#' @param ... Named arguments passed to the component functions.
#'   Arguments are routed by name to the matching formal parameters of
#'   \code{build_ellps()}, \code{get_suitable_env()}, and
#'   \code{get_sample_occ()}. Unknown names are ignored with a warning.
#' @param out.file Logical. If \code{TRUE}, save the returned object to an
#'   \code{.rds} file in the working directory.
#' @param out.file.name Optional base name (without extension) for the saved
#'   file. If \code{NULL} or empty, a timestamped name is used.
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @details
#' Each low-level function (\code{build_ellps()}, \code{get_suitable_env()},
#' \code{get_sample_occ()}) can be used independently. This wrapper function
#' simply wires them together and standardizes a few details:
#' \itemize{
#'   \item \code{get_suitable_env()} is always called with
#'         \code{distances = TRUE}.
#'   \item If \code{out} is not supplied for \code{get_suitable_env()}, it
#'         defaults to \code{"both"} so that both a raster and a data.frame
#'         are available.
#'   \item The resulting suitable environment data.frame is passed to
#'         \code{get_sample_occ()} as \code{suitable_env} if the user did not
#'         already supply one.
#'   \item If \code{env_bg} is supplied only at the top level, it is forwarded
#'         to \code{get_suitable_env()} (and to \code{get_sample_occ()}
#'         only if the user routed it explicitly).
#' }
#'
#' @return
#' A list of class \code{"NicheR_species"} with elements:
#' \itemize{
#'   \item \code{niche}: the ellipsoid object returned by \code{build_ellps()}.
#'   \item \code{suitability}: the object returned by \code{get_suitable_env()}
#'         (data.frame, SpatRaster, or \code{"suitable_env"} list depending
#'         on \code{out}).
#'   \item \code{occurrences}: sampled occurrences returned by
#'         \code{get_sample_occ()}.
#'   \item \code{call_args}: the original \code{...} captured as a named list.
#'   \item \code{routed_args}: a list showing which args went to each function.
#'   \item \code{save_path}: file path if \code{out.file = TRUE}, otherwise
#'         \code{NULL}.
#' }
#'
#' @seealso \code{\link{build_ellps}}, \code{\link{get_suitable_env}},
#'   \code{\link{get_sample_occ}}
#'
#' @examples
#' \dontrun{
#' vs <- create_virtual_species(
#'   # build_ellps() args
#'   center = c(temp = 20, precip = 140, seas = 25),
#'   axes   = c(8, 45, 6.5),
#'
#'   # get_suitable_env() args
#'   env_bg = my_env_spatraster,   # SpatRaster of predictors
#'
#'   # get_sample_occ() args
#'   n_occ        = 500,
#'   method       = "center",
#'   bias_surface = my_bias_raster  # optional
#' )
#' print(vs)
#' }
#'
#' @export
create_virtual_species <- function(...,
                                   out.file = FALSE,
                                   out.file.name = NULL,
                                   verbose = TRUE) {

  # capture all args once
  args <- list(...)

  # ensure required functions exist
  needed_funs <- c("build_ellps", "get_suitable_env", "get_sample_occ")
  missing_funs <- needed_funs[!vapply(needed_funs, exists, logical(1), mode = "function")]

  if (length(missing_funs)) {
    stop(
      "These functions are not available in the current environment: ",
      paste(missing_funs, collapse = ", ")
    )
  }

  # pull formals for routing
  f_build <- names(formals(build_ellps))
  f_suit  <- names(formals(get_suitable_env))
  f_occ   <- names(formals(get_sample_occ))

  # split args by target function
  args_build <- args[names(args) %in% f_build]
  args_suit  <- args[names(args) %in% f_suit]
  args_occ   <- args[names(args) %in% f_occ]

  # pass env_bg through to get_suitable_env if supplied top-level
  if ("env_bg" %in% names(args) &&
      !("env_bg" %in% names(args_suit)) &&
      ("env_bg" %in% f_suit)) {

    args_suit$env_bg <- args$env_bg
  }

  # check that env_bg will be available for get_suitable_env
  if (("env_bg" %in% f_suit) && !("env_bg" %in% names(args_suit))) {
    stop("env_bg is required by get_suitable_env but was not supplied.")
  }

  # --- normalize get_suitable_env control args (out, distances) --------------

  # force distances = TRUE (warn if user tried FALSE)
  if ("distances" %in% names(args_suit)) {
    if (isFALSE(args_suit$distances)) {
      warning(
        "create_virtual_species() requires distances = TRUE for get_suitable_env(); ",
        "overriding user-supplied distances = FALSE."
      )
    }
  }
  args_suit$distances <- TRUE

  # default out = "both" if user did not specify
  if (!"out" %in% names(args_suit)) {
    args_suit$out <- "both"
  }

  # --- build ellipsoid niche --------------------------------------------------

  niche_obj <- tryCatch(
    do.call(build_ellps, args_build),
    error = function(e) stop("build_ellps failed: ", e$message)
  )

  if (verbose) message("Built niche object.")

  # --- compute suitable environment ------------------------------------------

  suit_env <- tryCatch(
    do.call(get_suitable_env, c(list(niche = niche_obj), args_suit)),
    error = function(e) stop("get_suitable_env failed: ", e$message)
  )

  if (verbose) message("Computed suitable environments.")

  # Extract a suitable_env data.frame for get_sample_occ if needed
  suitable_df <- NULL

  if (inherits(suit_env, "suitable_env") &&
      is.list(suit_env) &&
      "suitable_env_df" %in% names(suit_env)) {

    suitable_df <- suit_env$suitable_env_df

  } else if (is.data.frame(suit_env)) {

    suitable_df <- suit_env

  } else if (inherits(suit_env, "SpatRaster")) {

    # last resort: convert raster to df; must contain dist_sq already
    suitable_df <- as.data.frame.nicheR(suit_env)
  }

  if (!is.null(suitable_df)) {
    if (!all(c("x", "y", "dist_sq") %in% names(suitable_df))) {
      stop(
        "Suitable environment returned by get_suitable_env() does not contain ",
        "the required columns 'x', 'y', and 'dist_sq'. ",
        "Ensure get_suitable_env(..., distances = TRUE, out = 'data.frame' or 'both')."
      )
    }
  }

  # --- wire suitable_env into get_sample_occ if user did not supply it -------

  if (!"suitable_env" %in% names(args_occ) && !is.null(suitable_df)) {
    args_occ$suitable_env <- suitable_df
  }

  # always force the same niche into get_sample_occ, even if user passed another
  occ <- tryCatch(
    do.call(get_sample_occ, c(list(niche = niche_obj), args_occ)),
    error = function(e) stop("get_sample_occ failed: ", e$message)
  )

  if (verbose) message("Sampled occurrences.")

  # warn on unused args
  used_names <- union(
    names(args_build),
    union(names(args_suit), names(args_occ))
  )
  unused <- setdiff(names(args), used_names)
  if (length(unused) && verbose) {
    warning(
      "These arguments did not match any target function and were ignored: ",
      paste(unused, collapse = ", ")
    )
  }

  # assemble S3 object
  out <- list(
    niche       = niche_obj,
    suitability = suit_env,
    occurrences = occ,
    call_args   = args,
    routed_args = list(
      build_ellps      = args_build,
      get_suitable_env = args_suit,
      get_sample_occ   = args_occ
    ),
    save_path   = NULL
  )

  class(out) <- c("NicheR_species", class(out))

  # ---- Saving logic ----------------------------------------------------------

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
    message("out.file = FALSE: object not saved to disk.")
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
  cat("  niche:       ", paste(class(x$niche), collapse = "/"), "\n", sep = "")
  cat("  suitability: ", paste(class(x$suitability), collapse = "/"), "\n", sep = "")
  cat("  occurrences: ", paste(class(x$occurrences), collapse = "/"), "\n", sep = "")
  invisible(x)
}
