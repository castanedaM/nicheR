#' Create a virtual species end to end
#' Routes ... to build_ellps(), get_suitable_env(), and get_sample_occ()
#' Returns an S3 object of class "NicheR_species"

create_virtual_species <- function(..., verbose = TRUE) {
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

  return(out)
}

#' Print method for NicheR_species
print.NicheR_species <- function(x, ...) {
  cat("NicheR virtual species components:\n")
  cat("  niche:       ", paste(class(x$niche), collapse = "/"), "\n", sep = "")
  cat("  suitability: ", paste(class(x$suitability), collapse = "/"), "\n", sep = "")
  cat("  occurrences: ", paste(class(x$occurrences), collapse = "/"), "\n", sep = "")
  invisible(x)
}
