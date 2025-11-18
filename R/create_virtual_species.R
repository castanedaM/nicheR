#' Create a virtual species end to end
#'
#' Orchestrates the full NicheR workflow by routing `...` to
#' [build_ellps()], [get_suitable_env()], optionally [set_bias_surface()],
#' and [get_sample_occ()]. Returns an S3 object of class **NicheR_species**
#' containing the niche object, suitability, (optional) bias surface, and
#' sampled occurrences.
#'
#' @section Workflow:
#' 1. Calls [build_ellps()] to define the ellipsoid in E-space.
#' 2. Calls [get_suitable_env()] to compute suitability in G-space.
#' 3. Optionally calls [set_bias_surface()] to build a pooled bias raster
#'    when bias-construction arguments are supplied (e.g., `bias_surface`,
#'    `bias_dir`, `suitable_env`, `out.bias`).
#' 4. Calls [get_sample_occ()] to sample occurrences from the suitable
#'    (and optionally biased) environment.
#'
#' @param ... Named arguments passed to the component functions.
#'   Arguments are routed by name to the matching formal parameters of
#'   [build_ellps()], [get_suitable_env()], [set_bias_surface()], and
#'   [get_sample_occ()]. Unknown names are ignored with a warning.
#'   \itemize{
#'     \item \strong{build_ellps()}: e.g. `center`, `axes`, `angles`, etc.
#'     \item \strong{get_suitable_env()}: e.g. `env_bg`, `out.suit`,
#'           `distances`, `var_names`.
#'     \item \strong{set_bias_surface()}: e.g. `bias_surface`, `bias_dir`,
#'           `suitable_env`, `out.bias`.
#'     \item \strong{get_sample_occ()}: e.g. `n_occ`, `method`,
#'           `bias_surface`, `suitable_env`, `seed`, `verbose`.
#'   }
#' @param out.file Logical. If `TRUE`, save the returned object to an `.rds`
#'   file in the working directory.
#' @param out.file.name Optional base name (without extension) for the saved
#'   file. If `NULL` or empty, a timestamped name is used.
#' @param verbose Logical. If `TRUE`, print progress messages.
#'
#' @details
#' If `env_bg` is supplied only at the top level, it is forwarded to both
#' [get_suitable_env()] and [get_sample_occ()] (if those functions accept an
#' `env_bg` argument). If a required argument (such as `env_bg`) is missing
#' for a component, the function stops with an informative error.
#'
#' When available in [get_suitable_env()]'s formals and not explicitly set
#' by the user, `create_virtual_species()` sets:
#' \itemize{
#'   \item `out.suit = "both"` so that both raster and data.frame outputs are
#'         available for plotting and downstream use.
#'   \item `distances = TRUE` so that distance-to-centroid information is
#'         computed and can be used in geographic plots.
#' }
#'
#' If [get_sample_occ()] has a `suitable_env` argument and the user does not
#' supply it, the suitability object returned by [get_suitable_env()] is
#' automatically forwarded as `suitable_env`.
#'
#' \strong{Bias handling:}
#' \itemize{
#'   \item If the user passes only `bias_surface` (and no additional
#'         [set_bias_surface()] arguments), it is assumed to be a
#'         \emph{precomputed single-layer raster} scaled to \eqn{[0, 1]} and
#'         is passed directly to [get_sample_occ()] (which expects a 0–1 layer).
#'   \item If the user passes `bias_surface` \emph{plus} any
#'         [set_bias_surface()] arguments (e.g. `bias_dir`, `suitable_env`,
#'         `out.bias`), then [set_bias_surface()] is called once, and its
#'         pooled bias raster is passed to [get_sample_occ()] as
#'         `bias_surface`. The full bias object is stored in the output
#'         as `bias_surface`.
#'   \item When [set_bias_surface()] is invoked and `suitable_env` is not
#'         provided to it explicitly, the suitability object from
#'         [get_suitable_env()] is used as the mask/template.
#' }
#'
#' @return
#' A list of class **NicheR_species** with elements:
#' \itemize{
#'   \item \code{niche}: the ellipsoid object returned by [build_ellps()].
#'   \item \code{suitability}: suitability object from [get_suitable_env()].
#'   \item \code{bias_surface}: either `NULL`, a \code{nicheR_bias_surface}
#'         object (if [set_bias_surface()] was used), or a user-supplied
#'         single-layer bias raster.
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

  if (isTRUE(verbose)) {
    message("[create_virtual_species] Starting virtual species workflow.")
  }

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

  # ensure env_bg will be available when needed for suitability
  if (("env_bg" %in% f_suit) && !("env_bg" %in% names(args_suit))) {
    stop("env_bg is required by get_suitable_env but was not supplied.")
  }

  # set sensible defaults for get_suitable_env behavior
  if ("out.suit" %in% f_suit && !("out.suit" %in% names(args_suit))) {
    args_suit$out.suit <- "both"
    if (isTRUE(verbose)) {
      message("[create_virtual_species] get_suitable_env(): 'out.suit' not supplied; using 'both'.")
    }
  }
  if ("distances" %in% f_suit && !("distances" %in% names(args_suit))) {
    args_suit$distances <- TRUE
    if (isTRUE(verbose)) {
      message("[create_virtual_species] get_suitable_env(): 'distances' not supplied; computing distances = TRUE.")
    }
  }

  # force inner functions to be quiet if they have a 'verbose' argument
  if ("verbose" %in% f_build) args_build$verbose <- FALSE
  if ("verbose" %in% f_suit)  args_suit$verbose  <- FALSE
  if ("verbose" %in% f_bias)  args_bias$verbose  <- FALSE
  if ("verbose" %in% f_occ)   args_occ$verbose   <- FALSE

  # ---------------------------------------------------------------------------
  # 1) build ellipsoid niche
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    message("[create_virtual_species] Running build_ellps()...")
  }
  niche_obj <- tryCatch(
    do.call(build_ellps, args_build),
    error = function(e) stop("build_ellps failed: ", e$message)
  )
  if (isTRUE(verbose)) {
    message("[create_virtual_species] Finished build_ellps().")
  }

  # ---------------------------------------------------------------------------
  # 2) suitability in G space
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    message("[create_virtual_species] Running get_suitable_env()...")
  }
  suit_env <- tryCatch(
    do.call(get_suitable_env, c(list(niche = niche_obj), args_suit)),
    error = function(e) stop("get_suitable_env failed: ", e$message)
  )
  if (isTRUE(verbose)) {
    message("[create_virtual_species] Finished get_suitable_env().")
  }

  # If get_sample_occ has a 'suitable_env' arg and user didn't set it,
  # wire through the suitability object from get_suitable_env().
  if ("suitable_env" %in% f_occ && !("suitable_env" %in% names(args_occ))) {
    args_occ$suitable_env <- suit_env
    if (isTRUE(verbose)) {
      message("[create_virtual_species] get_sample_occ(): using suitability from get_suitable_env() as 'suitable_env'.")
    }
  }

  # ---------------------------------------------------------------------------
  # 3) Optional bias construction
  # ---------------------------------------------------------------------------
  bias_obj <- NULL

  # "Bias in play" if user provided bias_surface at all (in top-level args).
  has_any_bias_arg <- "bias_surface" %in% names(args)

  if (has_any_bias_arg) {

    # Decide if we should CALL set_bias_surface() or just forward bias_surface.
    # We call set_bias_surface() if any set_bias_surface-specific args
    # beyond 'bias_surface' and 'verbose' are present (e.g., bias_dir,
    # suitable_env, out.bias).
    bias_construction_args <- setdiff(names(args_bias),
                                      c("bias_surface", "verbose"))

    wants_bias_build <- length(bias_construction_args) > 0

    if (wants_bias_build) {
      # We will build a bias surface with set_bias_surface()
      if (isTRUE(verbose)) {
        message("[create_virtual_species] Running set_bias_surface()...")
      }

      bias_call_args <- args_bias

      # If no suitable_env given to set_bias_surface, use the suitability
      # object from get_suitable_env() as mask/template.
      if (!("suitable_env" %in% names(bias_call_args))) {
        bias_call_args$suitable_env <- suit_env
      }

      # make sure inner verbose is FALSE
      if ("verbose" %in% f_bias) {
        bias_call_args$verbose <- FALSE
      }

      bias_obj <- tryCatch(
        do.call(set_bias_surface, bias_call_args),
        error = function(e) stop("set_bias_surface failed: ", e$message)
      )

      if (isTRUE(verbose)) {
        message("[create_virtual_species] Finished set_bias_surface().")
      }

      # get_sample_occ() expects a single-layer bias raster;
      # pass the pooled layer if present.
      if (!is.null(bias_obj$pooled_bias_sp)) {
        args_occ$bias_surface <- bias_obj$pooled_bias_sp
      } else {
        warning("set_bias_surface() returned no 'pooled_bias_sp'; ",
                "bias will not be applied in get_sample_occ().",
                call. = FALSE)
      }

    } else {
      # No extra set_bias_surface args: assume the user passed a precomputed
      # 0–1 bias raster as 'bias_surface' that should go straight into
      # get_sample_occ().
      if ("bias_surface" %in% names(args)) {
        args_occ$bias_surface <- args$bias_surface
        bias_obj <- args$bias_surface  # store in output for transparency
        if (isTRUE(verbose)) {
          message("[create_virtual_species] Using user-supplied 'bias_surface' directly in get_sample_occ().")
        }
      }
    }
  } else if (isTRUE(verbose)) {
    message("[create_virtual_species] No 'bias_surface' supplied; skipping bias construction.")
  }

  # ---------------------------------------------------------------------------
  # 4) sample occurrences
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    message("[create_virtual_species] Running get_sample_occ()...")
  }
  occ <- tryCatch(
    do.call(get_sample_occ, args_occ),
    error = function(e) stop("get_sample_occ failed: ", e$message)
  )
  if (isTRUE(verbose)) {
    message("[create_virtual_species] Finished get_sample_occ().")
  }

  # ---------------------------------------------------------------------------
  # 5) warn on unused args
  # ---------------------------------------------------------------------------
  used_names <- union(
    names(args_build),
    union(names(args_suit), union(names(args_occ), names(args_bias)))
  )
  unused <- setdiff(names(args), used_names)
  if (length(unused) && verbose) {
    warning("[create_virtual_species] These arguments did not match any target function and were ignored: ",
            paste(unused, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # 6) assemble S3 object
  # ---------------------------------------------------------------------------
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

    if (verbose) message("[create_virtual_species] Virtual species saved to: ",
                         normalizePath(save_path))

    out$save_path <- save_path

  } else if (verbose) {
    message("[create_virtual_species] out.file = FALSE: object not saved to disk.")
  }

  if (isTRUE(verbose)) {
    message("[create_virtual_species] Workflow completed.")
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
  cat("  bias_surface: ",
      if (is.null(x$bias_surface)) "NULL" else paste(class(x$bias_surface), collapse = "/"),
      "\n", sep = "")
  cat("  occurrences:  ", paste(class(x$occurrences), collapse = "/"), "\n", sep = "")
  invisible(x)
}
