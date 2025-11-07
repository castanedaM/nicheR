#' Create a virtual species end to end
#'
#' Orchestrates the full NicheR workflow by routing `...` to
#' [build_ellps()], [get_suitable_env()], and [get_sample_occ()].
#' Returns an S3 object of class **NicheR_species** containing the
#' niche object, suitability raster, and sampled occurrences.
#'
#' @section Workflow:
#' 1. Calls [build_ellps()] to define the ellipsoid in E-space.
#' 2. Calls [get_suitable_env()] to compute suitability in G-space.
#' 3. Calls [get_sample_occ()] to sample occurrences from the constrained surface.
#'
#' @param ... Named arguments passed to the component functions.
#'   Arguments are routed by name to the matching formal parameters of
#'   [build_ellps()], [get_suitable_env()], and [get_sample_occ()].
#'   Unknown names are ignored with a warning.
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
#' @return
#' A list of class **NicheR_species** with elements:
#' \itemize{
#'   \item \code{niche}: the ellipsoid object returned by [build_ellps()].
#'   \item \code{suitability}: suitability surface from [get_suitable_env()].
#'   \item \code{occurrences}: sampled occurrences from [get_sample_occ()].
#'   \item \code{call_args}: the original `...` captured as a named list.
#'   \item \code{routed_args}: a list showing which args went to each function.
#'   \item \code{save_path}: file path if \code{out.file = TRUE}, otherwise `NULL`.
#' }
#'
#' @seealso [build_ellps()], [get_suitable_env()], [get_sample_occ()]
#'
#' @examples
#' \dontrun{
#' vs <- create_virtual_species(
#'   # build_ellps() args
#'   centroid = c(temp = 20, precip = 140, seas = 25),
#'   axes     = c(8, 45, 6.5),
#'   # get_suitable_env() args
#'   env_bg   = my_env_spatraster,
#'   bias     = my_bias_raster,
#'   # get_sample_occ() args
#'   n_points = 500
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
    stop("These functions are not available in the current environment: ",
         paste(missing_funs, collapse = ", "))
  }

  # pull formals for routing
  f_build <- names(formals(build_ellps))
  f_suit  <- names(formals(get_suitable_env))
  f_occ   <- names(formals(get_sample_occ))

  # split args by target function
  args_build <- args[names(args) %in% f_build]
  args_suit  <- args[names(args) %in% f_suit]
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

  # check that env_bg will be available when needed
  if (("env_bg" %in% f_suit) && !("env_bg" %in% names(args_suit))) {
    stop("env_bg is required by get_suitable_env but was not supplied.")
  }
  if (("env_bg" %in% f_occ) && !("env_bg" %in% names(args_occ))) {
    stop("env_bg is required by get_sample_occ but was not supplied.")
  }

  # build ellipsoid niche
  niche_obj <- tryCatch(
    do.call(build_ellps, args_build),
    error = function(e) stop("build_ellps failed: ", e$message)
  )

  if(verbose) message("Built niche object.")

  # suitability in G space
  # force niche into the call even if user passed one in args_suit
  suit_env <- tryCatch(
    do.call(get_suitable_env, c(list(niche = niche_obj), args_suit)),
    error = function(e) stop("get_suitable_env failed: ", e$message)
  )

  if(verbose) message("Computed suitable environments.")

  # sample occurrences
  occ <- tryCatch(
    do.call(get_sample_occ, c(list(niche = niche_obj), args_occ)),
    error = function(e) stop("get_sample_occ failed: ", e$message)
  )
  if (verbose) message("Sampled occurrences.")

  # warn on unused args
  used_names <- union(names(args_build), union(names(args_suit), names(args_occ)))
  unused <- setdiff(names(args), used_names)
  if (length(unused) && verbose) {
    warning("These arguments did not match any target function and were ignored: ",
            paste(unused, collapse = ", "))
  }

  # assemble S3 object
  out <- list(
    niche       = niche_obj,
    suitability = suit_env,
    occurrences = occ,
    call_args   = args,
    routed_args = list(
      build_ellps = args_build,
      get_suitable_env = args_suit,
      get_sample_occ = args_occ
    )
  )

  class(out) <- c("NicheR_species", class(out))


  # ---- Saving logic ----
  if (isTRUE(out.file)) {
    # directory = current working directory
    dir_path <- getwd()

    # auto-generate file name if not provided
    if (is.null(out.file.name) || out.file.name == "") {
      timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      out.file.name <- paste0("NicheR_vs_", timestamp)
    }

    # full save path
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
  cat("  niche:       ", paste(class(x$niche), collapse = "/"), "\n", sep = "")
  cat("  suitability: ", paste(class(x$suitability), collapse = "/"), "\n", sep = "")
  cat("  occurrences: ", paste(class(x$occurrences), collapse = "/"), "\n", sep = "")
  invisible(x)
}
