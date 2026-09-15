# Title: Session report builder

# Description: Turns the current session into an Rmd holding a methods
# narrative and the nicheR code that reproduces it. The generated document
# uses real code chunks, so one file is both the runnable script the user
# takes away and the source knitted to produce the HTML report.

# Layout: shared primitives first, then one block per workflow step. Each step
# keeps its inspection helpers, code emitters, prose, and section assembly
# together, so adding Bias and Generate means appending two more blocks in the
# same shape rather than threading through what is here.

# Ellipsoids are emitted self contained, built from their own range_inputs
# rather than chained off a parent. The library records current state, not
# click history, so a chained script stops reproducing the session as soon as a
# parent is edited after being copied. Lineage appears in the prose instead.

# Date last updated: 08/13/2026


# FORMATTING --------------------------------------------------------------

#' Format a number for emitted code
#'
#' Seven significant digits keeps typed values exact and computed values short.
#' trim = TRUE suppresses padding to a common width.
#'
#' @param x Numeric vector.
#'
#' @returns A character vector the same length as `x`.
#'
#' @noRd
report_num <- function(x){
  format(x, trim = TRUE, digits = 7, scientific = FALSE)
}


#' Quote a name that is not a syntactic R identifier
#'
#' Variable names come from raster layers and covariance pair labels contain a
#' hyphen, so backticks are needed more often than not.
#'
#' @param nm Character vector of names.
#'
#' @returns `nm`, with non-syntactic entries wrapped in backticks.
#'
#' @noRd
report_name <- function(nm){
  ifelse(nm == make.names(nm), nm, paste0("`", nm, "`"))
}


#' Character vector as an R literal
#'
#' @param x Character vector.
#'
#' @returns A single string of the form "c(\"a\", \"b\")".
#'
#' @noRd
report_chr <- function(x){
  paste0("c(", paste0("\"", x, "\"", collapse = ", "), ")")
}


#' Named numeric vector as an R literal
#'
#' @param x Named numeric vector.
#'
#' @returns A single string of the form "c(a = 1, b = 2)".
#'
#' @noRd
report_vec <- function(x){
  paste0("c(",
         paste0(report_name(names(x)), " = ", report_num(x),
                collapse = ", "),
         ")")
}


#' Join a character vector into a readable list
#'
#' @param x Character vector.
#'
#' @returns A single string.
#'
#' @noRd
report_and <- function(x){
  if(length(x) < 2) return(paste(x, collapse = ""))
  paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])
}


#' Object names for the emitted script
#'
#' Ellipsoid names are user supplied, so they need sanitising, and two
#' ellipsoids may share a name.
#'
#' @param nms Character vector of display names.
#' @param prefix Fallback prefix when a name sanitises to nothing.
#'
#' @returns A character vector of unique syntactic names.
#'
#' @noRd
report_obj_names <- function(nms, prefix = "object"){
  out <- gsub("[^A-Za-z0-9]+", "_", nms)
  out <- gsub("^_+|_+$", "", out)
  out[!nzchar(out)] <- prefix
  out[grepl("^[0-9]", out)] <- paste0(prefix, "_", out[grepl("^[0-9]", out)])
  make.unique(out, sep = "_")
}


# CODE EMISSION -----------------------------------------------------------

#' Format a call as assignment code
#'
#' Every emitter builds calls, so comma and indent handling lives here once.
#' Values in args must already be formatted as code.
#'
#' @param obj Name to assign to.
#' @param fn Function name.
#' @param args Named list of formatted argument values.
#'
#' @returns A character vector of code lines.
#'
#' @noRd
report_call <- function(obj, fn, args){

  # An argument formatted as more than one string would recycle against the
  # names below and emit duplicated arguments, which reads as valid code until
  # it is run
  bad <- names(args)[lengths(args) != 1]
  if(length(bad) > 0){
    stop("report_call: arguments must be single strings, got several for ",
         paste(bad, collapse = ", "), call. = FALSE)
  }

  vals <- unlist(args, use.names = FALSE)
  c(paste0(obj, " <- ", fn, "("),
    paste0("  ", names(args), " = ", vals,
           c(rep(",", length(args) - 1), "")),
    ")")
}


#' Wrap emitted code as a knitr chunk
#'
#' @param lines Character vector of code lines.
#' @param label Chunk label, unique within the document.
#' @param eval Logical, whether the chunk runs when the report is knitted.
#'
#' @returns A character vector including the chunk delimiters.
#'
#' @noRd
report_chunk <- function(lines, label, eval = TRUE){
  c(paste0("```{r ", label, ", eval = ", eval, "}"), lines, "```", "")
}


# SESSION DATA ------------------------------------------------------------

#' Whether the emitted script can run at report time
#'
#' Sessions built on uploaded files refer to a folder only the user can point
#' at, so their chunks are emitted unevaluated. Example and virtual sessions
#' are self contained and are evaluated, which checks that the script runs.
#' Reads session_data from the enclosing server environment.
#'
#' @returns A single logical.
#'
#' @noRd
report_can_eval <- function(){
  isTRUE(session_data$input_mode %in% c("example", "virtual"))
}


#' Name the emitted script gives the background data
#'
#' @returns A single string.
#'
#' @noRd
report_data_object <- function(){
  if(identical(session_data$file_type, "df") &&
     is.null(session_data$bg_raster)) "bg" else "bios"
}


#' Object names for every ellipsoid in the library
#'
#' Computed over the whole library rather than a subset, so the name an
#' ellipsoid is given in the Build section is the same one later sections use
#' to refer to it. Reads session_data from the enclosing server environment.
#'
#' @returns A character vector of object names, named by ell_id.
#'
#' @noRd
report_ell_objects <- function(){
  ells <- session_data$ellipsoid_list
  if(length(ells) == 0) return(character(0))
  setNames(report_obj_names(vapply(ells, function(e) e$ell_name, character(1)),
                            prefix = "ellipsoid"),
           vapply(ells, function(e) e$ell_id, character(1)))
}


#' Code that reloads the environmental data
#'
#' Uploaded files cannot be located from the app. The browser sends only a file
#' name and Shiny's copy of the upload is deleted when the session ends, so the
#' emitted script declares the folder as a variable for the user to set rather
#' than writing a path that will not exist. Reads session_data from the
#' enclosing server environment.
#'
#' @returns A character vector of code lines.
#'
#' @noRd
report_data_code <- function(){

  vars <- session_data$vars

  if(identical(session_data$input_mode, "virtual")){
    return(c("# Virtual mode: the niche is defined in environmental space",
             "# only, so there are no environmental layers to load.",
             paste0("vars <- ", report_chr(vars))))
  }

  if(identical(session_data$input_mode, "example")){
    return(c("# Bundled example data",
             "bios <- terra::rast(",
             "  system.file(\"extdata/ma_bios.tif\", package = \"nicheR\"))",
             paste0("bios <- terra::subset(bios, ", report_chr(vars), ")")))
  }

  src <- session_data$data_source
  if(is.null(src) || length(src) == 0){
    src <- if(identical(session_data$file_type, "df")){
      "your_data.csv"
    } else {
      "your_layers.tif"
    }
  }

  head <- c("# Set this to the folder holding the files you uploaded to the app",
            "data_dir <- \".\"",
            "")

  if(identical(session_data$file_type, "df")){
    return(c(head,
             paste0("bg <- read.csv(file.path(data_dir, \"", src[1], "\"))"),
             paste0("bg <- bg[, ", report_chr(vars), "]")))
  }

  load <- if(length(src) == 1){
    paste0("bios <- terra::rast(file.path(data_dir, \"", src, "\"))")
  } else {
    c("bios <- terra::rast(file.path(data_dir,",
      paste0("                            ", report_chr(src), "))"))
  }

  c(head, load,
    paste0("bios <- terra::subset(bios, ", report_chr(vars), ")"))
}


# ELLIPSOID INSPECTION ----------------------------------------------------

#' Off diagonal covariance entries
#'
#' Returned as a table rather than a named vector so that callers read the two
#' variables from their own columns. Splitting a "v1-v2" label back apart
#' breaks as soon as a raster layer name contains a hyphen.
#'
#' @param ell A nicheR ellipsoid.
#'
#' @returns A data frame with columns v1, v2, value. Zero rows when every pair
#' is zero.
#'
#' @noRd
report_cov_pairs <- function(ell){

  v <- ell$var_names
  if(length(v) < 2) return(data.frame(v1 = character(0), v2 = character(0),
                                      value = numeric(0)))

  idx <- which(upper.tri(ell$cov_matrix), arr.ind = TRUE)
  out <- data.frame(v1 = v[idx[, 1]], v2 = v[idx[, 2]],
                    value = ell$cov_matrix[idx],
                    stringsAsFactors = FALSE)

  out[out$value != 0, , drop = FALSE]
}


#' Covariance entries formatted for update_ellipsoid_covariance()
#'
#' @param cov A table from report_cov_pairs().
#'
#' @returns A named numeric vector.
#'
#' @noRd
report_cov_vec <- function(cov){
  setNames(cov$value, paste0(cov$v1, "-", cov$v2))
}


#' Centroid an ellipsoid had when it was built from its ranges
#'
#' Derived from range_inputs rather than stored, since build_ellipsoid() places
#' the centroid at the midpoint of each range and range_inputs is frozen at
#' build time.
#'
#' @param ell A nicheR ellipsoid carrying the app fields.
#'
#' @returns A named numeric vector.
#'
#' @noRd
report_centroid_origin <- function(ell){
  vars <- ell$var_names
  ri <- ell$range_inputs
  (unlist(ri$min[vars]) + unlist(ri$max[vars])) / 2
}


#' Whether the centroid differs meaningfully from a reference position
#'
#' The centroid sliders round to two decimals, so no move a user can make is
#' smaller than 0.01. A difference below that came from rounding.
#'
#' @param x Named numeric vector, the current centroid.
#' @param ref Named numeric vector, the position to compare against.
#' @param ell A nicheR ellipsoid.
#'
#' @returns A logical vector, one entry per variable.
#'
#' @noRd
report_centroid_diff <- function(x, ref, ell){
  vars <- ell$var_names
  abs(x[vars] - ref[vars]) >= 0.01
}


#' Whether the centroid was moved away from where the ranges put it
#'
#' Compared against the range midpoint, not against a parent, because this
#' decides whether the emitted script needs an update call and the script
#' rebuilds each ellipsoid from its own range_inputs.
#'
#' @param ell A nicheR ellipsoid.
#'
#' @returns A single logical.
#'
#' @noRd
report_centroid_moved <- function(ell){
  any(report_centroid_diff(ell$centroid, report_centroid_origin(ell), ell))
}


# BUILD, CODE -------------------------------------------------------------

#' Code that rebuilds one range specification
#'
#' The three branches mirror the three ways the Build tab sets ranges. Manual
#' and data derived ranges become a data frame literal, stats ranges become a
#' call, so the script documents the method where it can.
#'
#' @param ell A nicheR ellipsoid carrying the app fields.
#' @param obj Name to give the range object.
#'
#' @returns A list with `pre`, lines to emit first, and `code`, the assignment.
#'
#' @noRd
report_range_code <- function(ell, obj){

  ri <- ell$range_inputs
  vars <- ell$var_names
  emin <- unlist(ri$expand_min[vars])
  emax <- unlist(ri$expand_max[vars])

  if(identical(ri$method, "stats")){
    return(list(pre = character(0),
                code = report_call(obj, "ranges_from_stats",
                                   list(mean = report_vec(unlist(ri$mean[vars])),
                                        sd = report_vec(unlist(ri$sd[vars])),
                                        cl = report_num(ell$cl),
                                        expand_min = report_vec(emin),
                                        expand_max = report_vec(emax)))))
  }

  mins <- unlist(ri$min[vars])
  maxs <- unlist(ri$max[vars])
  args <- as.list(setNames(paste0("c(", report_num(mins), ", ",
                                  report_num(maxs), ")"),
                           report_name(vars)))

  # The occurrence table behind a data derived range is not carried out of the
  # app, so the resolved bounds are given directly and the method is recorded
  # in a comment instead
  pre <- if(identical(ri$method, "df")){
    c("# Ranges derived in the app from an uploaded occurrence table with",
      paste0("# ranges_from_data(), expanding by ",
             report_and(paste0(vars, " ", report_num(emin), "/",
                               report_num(emax), "%")), "."),
      "# The resolved bounds are given directly so this runs without that table.")
  } else {
    character(0)
  }

  list(pre = pre, code = report_call(obj, "data.frame", args))
}


#' Code that rebuilds one ellipsoid
#'
#' Emitted in the order the app applies them: build, rotate, translate.
#' update_ellipsoid_covariance() rebuilds the object from its components, so a
#' centroid move written before it would be discarded.
#'
#' @param ell A nicheR ellipsoid.
#' @param obj Name to give the ellipsoid object.
#' @param range_obj Name of the range object it is built from.
#' @param parent_obj Name of the parent's object, or NULL for a root.
#'
#' @returns A character vector of code lines.
#'
#' @noRd
report_ell_code <- function(ell, obj, range_obj, parent_obj = NULL){

  out <- c(if(is.null(parent_obj)){
    paste0("# ", ell$ell_name, ", built from ranges")
  } else {
    paste0("# ", ell$ell_name, ", copied from ", parent_obj)
  },
  report_call(obj, "build_ellipsoid",
              list(range = range_obj, cl = report_num(ell$cl))))

  cov <- report_cov_pairs(ell)
  if(nrow(cov) > 0){
    out <- c(out, "",
             report_call(obj, "update_ellipsoid_covariance",
                         list(object = obj,
                              covariance = report_vec(report_cov_vec(cov)))))
  }

  if(report_centroid_moved(ell)){
    out <- c(out, "",
             report_call(obj, "update_ellipsoid_centroid",
                         list(ell = obj,
                              new_centroid = report_vec(round(ell$centroid[ell$var_names], 2)))))
  }

  out
}


# BUILD, PROSE ------------------------------------------------------------

#' Sentence describing where the ranges came from
#'
#' @param ell A nicheR ellipsoid carrying the app fields.
#'
#' @returns A single string.
#'
#' @noRd
report_range_prose <- function(ell){

  ri <- ell$range_inputs
  vars <- ell$var_names
  mins <- unlist(ri$min[vars])
  maxs <- unlist(ri$max[vars])
  spans <- report_and(paste0(vars, " from ", report_num(mins), " to ",
                             report_num(maxs), " (a span of ",
                             report_num(maxs - mins), ")"))

  lead <- switch(ri$method,
                 "man" = "Its bounds were entered directly",
                 "stats" = paste0("Its bounds were derived from a mean and standard ",
                                  "deviation per variable, so they describe a tolerance ",
                                  "reported as a distribution rather than one observed ",
                                  "directly"),
                 "df" = paste0("Its bounds were derived from the observed spread of an ",
                               "occurrence table and then expanded, so they describe the ",
                               "recorded conditions plus a margin"),
                 "Its bounds were set")

  paste0(lead, ", giving ", spans, ".")
}


#' Sentence describing the confidence level
#'
#' @param ell A nicheR ellipsoid.
#'
#' @returns A single string.
#'
#' @noRd
report_cl_prose <- function(ell){
  paste0("The confidence level of ", report_num(ell$cl), " sets where the ",
         "boundary falls, so conditions beyond it are treated as outside the ",
         "niche.")
}


#' Sentence describing the covariance, if any was set
#'
#' Describes the sign and magnitude of each pair and what that implies about
#' the orientation of the niche. No claim is made about what the variables
#' represent, since the app knows only their names.
#'
#' @param ell A nicheR ellipsoid.
#'
#' @returns A single string, or NULL when every pair is zero.
#'
#' @noRd
report_cov_prose <- function(ell){

  cov <- report_cov_pairs(ell)
  if(nrow(cov) == 0) return(NULL)

  parts <- paste0("a ", ifelse(cov$value > 0, "positive", "negative"),
                  " covariance of ", report_num(cov$value),
                  " between ", cov$v1, " and ", cov$v2)

  paste0("Covariance was set for ", nrow(cov),
         if(nrow(cov) == 1) " variable pair: " else " variable pairs: ",
         report_and(parts), ". ",
         "A non-zero covariance tilts the ellipsoid away from the variable ",
         "axes, so the range of ", cov$v1[1], " the niche includes depends on ",
         "the value of ", cov$v2[1], " rather than being the same throughout. ",
         "A positive value means the two increase together across the niche; ",
         "a negative value means one increases as the other decreases. ",
         "Because the covariance matrix determines the volume, this also ",
         "changes how much environmental space the niche occupies.")
}


#' Sentence describing the centroid move, if any
#'
#' Compared against the parent where there is one, since for a copy the useful
#' number is how far it sits from what it was copied from. This differs from
#' report_centroid_moved(), which decides what the script has to emit and so
#' must measure from the range midpoint.
#'
#' @param ell A nicheR ellipsoid.
#'
#' @returns A single string, or NULL when the centroid was not moved.
#'
#' @noRd
report_centroid_prose <- function(ell){

  vars <- ell$var_names
  parent <- ell_parent(ell)

  ref <- if(!is.null(parent)){
    parent$centroid[vars]
  } else {
    report_centroid_origin(ell)
  }

  ref_txt <- if(!is.null(parent)){
    paste0("relative to ", parent$ell_name)
  } else {
    "relative to the midpoint of its ranges"
  }

  moved <- vars[report_centroid_diff(ell$centroid, ref, ell)]
  if(length(moved) == 0) return(NULL)

  dv <- ell$centroid[moved] - ref[moved]
  parts <- paste0(moved, " by ", ifelse(dv > 0, "+", ""),
                  report_num(dv), " (from ",
                  report_num(ref[moved]), " to ",
                  report_num(ell$centroid[moved]), ")")

  paste0("The centroid sits ", ref_txt, " shifted in ", report_and(parts),
         ". Translating the centroid changes only where the niche sits, so ",
         "the ellipsoid keeps whatever shape and orientation it had before ",
         "the move. It describes the same set of tolerances centred on ",
         "different values of ", report_and(moved), ".")
}


#' Paragraph introducing one ellipsoid
#'
#' @param ell A nicheR ellipsoid.
#' @param parent_obj Name of the parent's object, or NULL for a root.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_ell_prose <- function(ell, parent_obj = NULL){

  open <- if(is.null(parent_obj)){
    paste0("**", ell$ell_name, "** is defined directly from environmental ",
           "ranges, so it does not inherit anything from another niche. ")
  } else {
    paste0("**", ell$ell_name, "** was derived from `", parent_obj,
           "`, which means the two are alternative descriptions of a niche ",
           "rather than independent hypotheses. ")
  }

  vol <- if(!is.null(ell$volume)){
    paste0(" Its niche volume is ", report_num(signif(ell$volume, 4)),
           ", which is the measure to compare across niches: a larger volume ",
           "means a broader set of tolerable conditions.")
  } else {
    ""
  }

  cov_txt <- report_cov_prose(ell)
  cen_txt <- report_centroid_prose(ell)

  c(paste0(open, report_range_prose(ell), " ", report_cl_prose(ell), vol), "",
    if(!is.null(cov_txt)) c(cov_txt, ""),
    if(!is.null(cen_txt)) c(cen_txt, ""))
}


#' Paragraph describing the library as a whole
#'
#' Generalisations live here and not in the code, because a description that is
#' slightly loose is still true while a loop that is slightly wrong produces a
#' script that does not reproduce the session.
#'
#' @param ells The ellipsoid list, in creation order.
#' @param n_range Number of distinct range specifications.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_build_prose <- function(ells, n_range){

  n <- length(ells)
  vars <- ells[[1]]$var_names
  cls <- unique(vapply(ells, function(e) e$cl, numeric(1)))
  n_root <- sum(vapply(ells, function(e) is.null(e$parent_id), logical(1)))

  out <- paste0(
    "We defined ", n, " ellipsoidal niche", if(n == 1) "" else "s", " in ",
    length(vars), " environmental dimensions (", paste(vars, collapse = ", "),
    ") using ",
    if(length(cls) == 1){
      paste0("a confidence level of ", report_num(cls))
    } else {
      paste0("confidence levels between ", report_num(min(cls)), " and ",
             report_num(max(cls)))
    },
    ". ",
    if(n_range == 1){
      "All were built from a single set of environmental ranges. "
    } else {
      paste0("These were built from ", n_range,
             " distinct sets of environmental ranges. ")
    },
    n_root, " of the ", n, " were defined directly from ranges",
    if(n_root < n){
      paste0(", and the remaining ", n - n_root, " were derived from these.")
    } else {
      "."
    })

  cov_n <- vapply(ells, function(e) nrow(report_cov_pairs(e)) > 0, logical(1))
  if(any(cov_n)){
    rng <- range(unlist(lapply(ells[cov_n], function(e) report_cov_pairs(e)$value)))
    out <- c(out, "",
             paste0("Covariance between environmental variables was set in ",
                    sum(cov_n), " of the ", n,
                    " niches, with off-diagonal values from ",
                    report_num(rng[1]), " to ", report_num(rng[2]),
                    ", rotating the ellipsoid relative to the variable axes."))
  }

  moved <- vapply(ells, report_centroid_moved, logical(1))
  if(any(moved)){
    out <- c(out, "",
             paste0("The centroid was translated in ", sum(moved), " of the ",
                    n, " niches, shifting the position of the niche in ",
                    "environmental space."))
  }

  c(out, "")
}


# BUILD, SECTION ----------------------------------------------------------

#' Rmd source for the Build section
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A character vector of Rmd lines.
#'
#' @noRd
report_build_section <- function(){

  ells <- session_data$ellipsoid_list
  if(length(ells) == 0){
    return(c("## Building ellipsoidal niches", "",
             "No ellipsoids were saved in this session.", ""))
  }

  obj <- report_ell_objects()
  ev <- report_can_eval()

  # Ranges are deduplicated on the code they produce, so two ellipsoids share a
  # range object exactly when rebuilding them would emit the same lines
  drafts <- lapply(ells, report_range_code, obj = "range")
  keys <- vapply(drafts, function(d) paste(c(d$pre, d$code), collapse = "\n"),
                 character(1))
  uniq <- unique(keys)
  range_obj <- paste0("range_", seq_along(uniq))
  range_of <- range_obj[match(keys, uniq)]

  range_lines <- unlist(lapply(seq_along(uniq), function(i){
    d <- report_range_code(ells[[match(uniq[i], keys)]], range_obj[i])
    c(d$pre, d$code, "")
  }))

  # One chunk per ellipsoid, each introduced by its own paragraph
  blocks <- unlist(lapply(seq_along(ells), function(i){
    p <- ells[[i]]$parent_id
    p_obj <- if(!is.null(p) && p %in% names(obj)) obj[[p]] else NULL
    c(paste0("### ", ells[[i]]$ell_name), "",
      report_ell_prose(ells[[i]], p_obj),
      report_chunk(report_ell_code(ells[[i]], obj[[i]], range_of[i], p_obj),
                   label = paste0("build-", i), eval = ev))
  }))

  c("## Building ellipsoidal niches", "",
    report_build_prose(ells, length(uniq)),
    report_chunk(report_data_code(), label = "build-data", eval = ev),
    paste0("The ", length(uniq), " set", if(length(uniq) == 1) "" else "s",
           " of environmental ranges below ",
           if(length(uniq) == 1) "is" else "are",
           " shared by the niches that follow."),
    "",
    report_chunk(range_lines, label = "build-ranges", eval = ev),
    blocks)
}


# PREDICT, INSPECTION -----------------------------------------------------

#' Output layers present in a stored prediction
#'
#' Coordinates and predictor variables are dropped, leaving only what predict()
#' added. With keep_data = TRUE the predictors are stacked in front of the
#' outputs, so this filter is what separates the two.
#'
#' @param pred A stored prediction, SpatRaster or data frame.
#' @param ell The ellipsoid it belongs to.
#'
#' @returns A character vector.
#'
#' @noRd
report_pred_layer_names <- function(pred, ell){
  if(is.null(pred)) return(character(0))
  nms <- if(inherits(pred, "SpatRaster")) names(pred) else colnames(pred)
  setdiff(nms, c("x", "y", ell$var_names))
}


#' Layer arguments implied by a stored prediction
#'
#' Recovered from the layers that came back, which is exact: a truncated layer
#' appears only when it was asked for.
#'
#' @param pred A stored prediction.
#' @param ell The ellipsoid it belongs to.
#'
#' @returns A named list of logicals.
#'
#' @noRd
report_pred_layers <- function(pred, ell){
  nms <- report_pred_layer_names(pred, ell)
  list(include_mahalanobis = "Mahalanobis" %in% nms,
       include_suitability = "suitability" %in% nms,
       mahalanobis_truncated = "Mahalanobis_trunc" %in% nms,
       suitability_truncated = "suitability_trunc" %in% nms)
}


#' Predict arguments for one stored prediction
#'
#' The layer arguments are recovered from the output. The truncation level
#' leaves no trace there, so it is read from what the app recorded, and is
#' omitted for sessions saved before that was stored.
#'
#' @param id The ell_id the prediction is stored under.
#' @param pred The stored prediction.
#' @param ell The ellipsoid it belongs to.
#'
#' @returns A named list of formatted argument values.
#'
#' @noRd
report_pred_args <- function(id, pred, ell){

  out <- lapply(report_pred_layers(pred, ell),
                function(x) if(isTRUE(x)) "TRUE" else "FALSE")

  lvl <- session_data$prediction_settings[[id]]$adjust_truncation_level
  if(length(lvl) == 1 && is.finite(lvl)){
    out$adjust_truncation_level <- report_num(lvl)
  }

  out
}


#' Value range of each layer a prediction produced
#'
#' Read from the stored surface rather than recomputed, so the numbers are the
#' ones the app showed.
#'
#' @param pred A stored prediction.
#' @param ell The ellipsoid it belongs to.
#'
#' @returns A data frame with columns layer, min, max.
#'
#' @noRd
report_pred_layer_stats <- function(pred, ell){

  lyrs <- report_pred_layer_names(pred, ell)
  if(length(lyrs) == 0){
    return(data.frame(layer = character(0), min = numeric(0),
                      max = numeric(0), stringsAsFactors = FALSE))
  }

  rng <- vapply(lyrs, function(l){
    v <- if(inherits(pred, "SpatRaster")){
      as.numeric(terra::minmax(pred[[l]]))
    } else {
      range(pred[[l]], na.rm = TRUE)
    }
    if(length(v) != 2 || any(!is.finite(v))) c(NA_real_, NA_real_) else v
  }, numeric(2))

  data.frame(layer = lyrs, min = rng[1, ], max = rng[2, ],
             stringsAsFactors = FALSE)
}


#' Proportion of the study area a prediction marks as suitable
#'
#' Read from the truncated suitability layer where it exists, since that is the
#' layer with a hard niche boundary. Returns NA when no such layer was
#' produced, rather than inventing a threshold on the continuous surface.
#'
#' @param pred A stored prediction.
#'
#' @returns A single numeric between 0 and 1, or NA.
#'
#' @noRd
report_pred_suitable <- function(pred){

  if(!inherits(pred, "SpatRaster")) return(NA_real_)
  if(!"suitability_trunc" %in% names(pred)) return(NA_real_)

  r <- pred[["suitability_trunc"]]
  tot <- terra::global(!is.na(r), "sum", na.rm = TRUE)[1, 1]
  if(is.na(tot) || tot == 0) return(NA_real_)

  terra::global(r > 0, "sum", na.rm = TRUE)[1, 1] / tot
}


#' Predict groups and the object each one is emitted as
#'
#' The name of the list a prediction lands in depends on how many argument
#' groups the session produced, and later sections have to index that list by
#' name. Deciding it here rather than inside the section builder means Bias and
#' Generate read the same answer rather than recomputing the grouping and
#' drifting from it. A wrong answer subsets a list that exists, so the emitted
#' script would run and bias the wrong surface.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A data frame with columns id, group, list_obj and key, one row per
#' predicted ellipsoid, in the order predictions are stored. Zero rows when
#' nothing was predicted.
#'
#' @noRd
report_pred_groups <- function(){

  empty <- data.frame(id = character(0), group = integer(0),
                      list_obj = character(0), key = character(0),
                      stringsAsFactors = FALSE)

  preds <- session_data$ellipsoid_prediction_list

  if(identical(session_data$input_mode, "virtual") || length(preds) == 0){
    return(empty)
  }

  ids <- names(preds)
  keep <- !vapply(session_data$ellipsoid_list[ids], is.null, logical(1))
  ids <- ids[keep]

  if(length(ids) == 0) return(empty)

  # Grouped on names as well as values, since a group carrying the truncation
  # argument and one without differ in length and values alone could collide
  keys <- vapply(ids, function(id){
    a <- report_pred_args(id, preds[[id]], session_data$ellipsoid_list[[id]])
    paste(names(a), unlist(a), collapse = "|")
  }, character(1))

  uniq <- unique(keys)
  g <- match(keys, uniq)

  data.frame(id = ids,
             group = g,
             list_obj = if(length(uniq) > 1) paste0("predictions_", g) else "predictions",
             key = keys,
             row.names = NULL,
             stringsAsFactors = FALSE)
}


# PREDICT, CODE -----------------------------------------------------------

#' Code that predicts a group of ellipsoids
#'
#' Prediction applies the same call to each ellipsoid, so a loop here states
#' only what the stored arguments already show. Ellipsoids predicted
#' differently land in different groups and so in different calls.
#'
#' @param objs Object names of the ellipsoids in the group.
#' @param args The predict arguments they share, already formatted.
#' @param out Name to give the resulting list.
#'
#' @returns A character vector of code lines.
#'
#' @noRd
report_pred_code <- function(objs, args, out){
  c(paste0(out, "_ellipsoids <- list(",
           paste0(objs, " = ", objs, collapse = ", "), ")"),
    "",
    paste0(out, " <- lapply(", out, "_ellipsoids, function(e){"),
    "  predict(e,",
    paste0("          newdata = ", report_data_object(), ","),
    paste0("          ", names(args), " = ", unlist(args), ","),
    "          verbose = FALSE)",
    "})")
}


# PREDICT, PROSE ----------------------------------------------------------

#' Paragraphs introducing the predict step
#'
#' @param ells The predicted ellipsoids.
#' @param preds Their stored predictions.
#' @param n_group Number of distinct argument groups.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_predict_prose <- function(ells, preds, n_group){

  n <- length(ells)
  spatial <- inherits(preds[[1]], "SpatRaster")

  out <- paste0(
    "Each of the ", n, " niche", if(n == 1) "" else "s",
    " defined above was projected ",
    if(spatial) "across the study area" else "onto the background records",
    ". For every ", if(spatial) "cell" else "record",
    ", the environmental values are compared against the ellipsoid to give a ",
    "measure of how close those conditions sit to the centre of the niche.",
    if(n_group > 1){
      paste0(" Not every niche was projected the same way, so the calls below ",
             "are grouped by the arguments they share.")
    } else {
      ""
    })

  layers <- unique(unlist(lapply(seq_along(ells), function(i){
    report_pred_layer_names(preds[[i]], ells[[i]])
  })))

  meanings <- c(
    "Mahalanobis" = paste0("**Mahalanobis** is the distance from the centre ",
                           "of the niche, measured in units that account for ",
                           "the spread and correlation of the variables. ",
                           "Larger values are further from the conditions the ",
                           "niche describes."),
    "suitability" = paste0("**suitability** rescales that distance to run ",
                           "from 1 at the centre toward 0 further out. It is ",
                           "continuous everywhere, so conditions outside the ",
                           "niche boundary still receive a small non-zero ",
                           "value."),
    "Mahalanobis_trunc" = paste0("**Mahalanobis_trunc** is the same distance ",
                                 "with everything beyond the niche boundary ",
                                 "set to `NA`, so only conditions inside the ",
                                 "niche carry a value."),
    "suitability_trunc" = paste0("**suitability_trunc** is the suitability ",
                                 "surface with everything beyond the boundary ",
                                 "set to zero. This is the layer with a hard ",
                                 "edge, and the one to use when the question ",
                                 "is whether conditions fall inside the niche ",
                                 "rather than how close to its centre they sit."))

  c(out, "",
    "The outputs produced were:", "",
    paste0("- ", meanings[names(meanings) %in% layers]), "")
}


#' Paragraph describing what one prediction produced
#'
#' Describes the surface that came back, not the call that produced it. The
#' call is shared across a group and is shown once; what differs between the
#' niches is the result.
#'
#' @param ell The ellipsoid.
#' @param pred Its stored prediction.
#' @param id The ell_id, used to read the recorded predict settings.
#' @param obj Its object name in the emitted script.
#' @param areas Named numeric vector of suitable area proportions for every
#' predicted ellipsoid, used for the comparison sentence.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_pred_ell_prose <- function(ell, pred, id, obj, areas){

  spatial <- inherits(pred, "SpatRaster")
  stats <- report_pred_layer_stats(pred, ell)
  p <- report_pred_suitable(pred)

  n_units <- if(spatial){
    tot <- terra::global(!is.na(pred[[stats$layer[1]]]), "sum",
                         na.rm = TRUE)[1, 1]
    paste0(format(tot, big.mark = ","), " cells")
  } else {
    paste0(format(nrow(pred), big.mark = ","), " records")
  }

  open <- paste0("Projecting **", ell$ell_name, "** onto ", n_units,
                 " produced ", nrow(stats), " layer",
                 if(nrow(stats) == 1) "" else "s", ".")

  # The truncation override leaves no trace in the returned layers, so without
  # this sentence a reader has no way to know the boundary was moved
  lvl <- session_data$prediction_settings[[id]]$adjust_truncation_level
  trunc_txt <- if(length(lvl) == 1 && is.finite(lvl)){
    paste0(" The truncation level was set to ", report_num(lvl),
           " rather than the confidence level of ", report_num(ell$cl),
           " this niche was built with, so the boundary of the prediction sits ",
           if(lvl < ell$cl) "tighter" else "wider", " than the ellipsoid.")
  } else {
    ""
  }

  ranges <- paste0("- `", stats$layer, "` ranges from ",
                   report_num(signif(stats$min, 4)), " to ",
                   report_num(signif(stats$max, 4)))

  area <- if(is.na(p)){
    paste0("No truncated suitability layer was produced for this niche, so ",
           "the extent of the suitable area is not reported. Suitability ",
           "falls off continuously and never reaches zero, which means there ",
           "is no boundary to measure without one.")
  } else {

    txt <- paste0(report_num(round(p * 100, 1)), "% of the ",
                  if(spatial) "study area" else "background records",
                  " falls inside the niche boundary.")

    rest <- areas[names(areas) != id]
    rest <- rest[!is.na(rest)]

    if(length(rest) >= 2){
      txt <- paste0(txt, if(p >= max(rest)){
        " This is the broadest of the niches projected here."
      } else if(p <= min(rest)){
        " This is the narrowest of the niches projected here."
      } else {
        paste0(" The others projected here range from ",
               report_num(round(min(rest) * 100, 1)), "% to ",
               report_num(round(max(rest) * 100, 1)), "%.")
      })
    } else if(length(rest) == 1){
      other <- session_data$ellipsoid_list[[names(rest)]]$ell_name
      txt <- paste0(txt, " ", other, " covers ",
                    report_num(round(rest * 100, 1)), "%.")
    }

    txt
  }

  c(paste0(open, trunc_txt), "",
    ranges, "",
    area, "",
    paste0("These surfaces are in `predictions[[\"", obj, "\"]]`."), "")
}


# PREDICT, SECTION --------------------------------------------------------

#' Rmd source for the Predict section
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A character vector of Rmd lines.
#'
#' @noRd
report_predict_section <- function(){

  preds <- session_data$ellipsoid_prediction_list

  # Virtual sessions never reach this step, and a session can be saved before
  # anything has been predicted
  if(identical(session_data$input_mode, "virtual") || length(preds) == 0){
    return(character(0))
  }

  ids <- names(preds)
  ells <- session_data$ellipsoid_list[ids]
  keep <- !vapply(ells, is.null, logical(1))
  ids <- ids[keep]
  ells <- ells[keep]
  preds <- preds[ids]

  if(length(ells) == 0) return(character(0))

  objs <- report_ell_objects()[ids]

  args <- lapply(seq_along(ells), function(i){
    report_pred_args(ids[i], preds[[i]], ells[[i]])
  })

  # The grouping and the name each group's list is given are decided by
  # report_pred_groups(), so the Bias section indexes the same object this
  # section creates rather than recomputing the rule and drifting from it
  grp <- report_pred_groups()
  if(nrow(grp) == 0) return(character(0))

  keys <- grp$key[match(ids, grp$id)]
  uniq <- unique(keys)

  # Computed once so each paragraph can place its niche against the others
  areas <- setNames(vapply(preds, report_pred_suitable, numeric(1)), ids)

  # Each group is followed by the niches it covers, so a reader moves from a
  # call to its results rather than back through every call first
  multi <- length(uniq) > 1

  blocks <- unlist(lapply(seq_along(uniq), function(g){

    sel <- which(keys == uniq[g])
    out <- grp$list_obj[match(ids[sel[1]], grp$id)]

    # A rule and a heading with several niches under it, so the code reads as
    # belonging to what follows rather than to the summaries above. Niche
    # headings drop a level to sit inside the group.
    open <- if(multi){
      c("***", "",
        paste0("### Prediction ", g), "",
        paste0("The following ", length(sel), " niche",
               if(length(sel) == 1) " was" else "s were",
               " projected together, since they were run with the same ",
               "settings."), "")
    } else {
      character(0)
    }

    # Shown once. Repeating it under every group would say the same thing
    # about a differently named list each time.
    how <- if(g == 1){
      lyr <- report_pred_layer_names(preds[[sel[1]]], ells[[sel[1]]])[1]
      c(paste0("Each element of `", out, "` is named after the niche it came ",
               "from, so a single surface can be pulled out with, for ",
               "example, `", out, "[[\"", objs[sel[1]], "\"]][[\"", lyr,
               "\"]]`."), "")
    } else {
      character(0)
    }

    detail <- unlist(lapply(sel, function(i){
      c(paste0(if(multi) "#### " else "### ", ells[[i]]$ell_name), "",
        report_pred_ell_prose(ells[[i]], preds[[i]], ids[i], objs[i], areas))
    }))

    c(open,
      report_chunk(report_pred_code(objs[sel], args[[sel[1]]], out),
                   label = paste0("predict-", g), eval = report_can_eval()),
      how,
      detail)
  }))

  c("## Projecting the niches", "",
    report_predict_prose(ells, preds, length(uniq)),
    blocks)
}


# BIAS, INSPECTION --------------------------------------------------------

#' Whether the Bias chunks can run at report time
#'
#' Evaluability has to be monotonic down the document. A chunk emitted with
#' eval = FALSE leaves its objects undefined, so any later chunk that runs and
#' references them fails the render outright rather than being skipped. Bias is
#' therefore evaluable only when Build and Predict were, and only when the bias
#' raster came from the bundled example rather than an upload.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A single logical.
#'
#' @noRd
report_bias_can_eval <- function(){
  report_can_eval() && identical(session_data$bias_source, "example")
}


#' Layers and effect directions behind a prepared bias surface
#'
#' The directions the user selected are read from what the app recorded.
#' prepare_bias() also encodes them in combination_formula as "(1-name)", but
#' recovering them that way means splitting a single string on " * " and a layer
#' named with a leading "(1-" or containing " * " would be misread, so the parse
#' is a fallback for sessions saved before bias_settings existed.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @param prepared The stored prepare_bias() output.
#'
#' @returns A data frame with columns layer and direction. Zero rows when
#' neither source is available.
#'
#' @noRd
report_bias_parts <- function(prepared){

  empty <- data.frame(layer = character(0), direction = character(0),
                      stringsAsFactors = FALSE)

  stored <- session_data$bias_settings$effect_direction
  if(length(stored) > 0 && !is.null(names(stored))){
    return(data.frame(layer = names(stored), direction = stored,
                      stringsAsFactors = FALSE))
  }

  f <- prepared$combination_formula
  if(is.null(f) || length(f) != 1L || !nzchar(f)) return(empty)

  parts <- trimws(strsplit(f, " * ", fixed = TRUE)[[1]])
  inv <- grepl("^\\(1-.*\\)$", parts)

  data.frame(layer = ifelse(inv, sub("^\\(1-(.*)\\)$", "\\1", parts), parts),
             direction = ifelse(inv, "inverse", "direct"),
             stringsAsFactors = FALSE)
}


#' Which optional outputs a prepared bias surface carries
#'
#' Read from the result rather than from what the app currently passes, so the
#' emitted call keeps matching the stored object if the tab changes.
#'
#' @param prepared The stored prepare_bias() output.
#'
#' @returns A named list of logicals.
#'
#' @noRd
report_bias_outputs <- function(prepared){
  list(include_composite = inherits(prepared$composite_surface, "SpatRaster"),
       include_processed_layers = inherits(prepared$processed_layers,
                                           "SpatRaster"))
}


#' Value summary of the composite bias surface
#'
#' The share of cells left as NA is reported because it is the only visible
#' consequence of the union or intersection choice, which is otherwise a
#' setting the reader has to take on trust.
#'
#' @param prepared The stored prepare_bias() output.
#'
#' @returns A named numeric vector with min, max, mean and na_share, or NULL
#' when no composite was produced.
#'
#' @noRd
report_bias_stats <- function(prepared){

  r <- prepared$composite_surface
  if(!inherits(r, "SpatRaster")) return(NULL)

  n_all <- terra::ncell(r)
  n_valid <- terra::global(!is.na(r), "sum", na.rm = TRUE)[1, 1]
  g <- terra::global(r, c("min", "max", "mean"), na.rm = TRUE)

  c(min = as.numeric(g[1, "min"]),
    max = as.numeric(g[1, "max"]),
    mean = as.numeric(g[1, "mean"]),
    na_share = if(n_all > 0) 1 - (n_valid / n_all) else NA_real_)
}


#' Prediction layer and direction behind each stored biased layer
#'
#' apply_bias() names every output "<layer>_biased_<direction>", and the app
#' stacks those layers per ellipsoid across repeated applications, so the stored
#' raster records exactly which calls were made. Nothing here is inferred, which
#' is why the Apply step needs no recorded settings the way Predict did.
#'
#' @param biased One element of ellipsoid_prediction_list_biased, a SpatRaster.
#'
#' @returns A data frame with columns name, layer and direction, one row per
#' stored layer. layer and direction are NA for a name that does not parse.
#'
#' @noRd
report_bias_applied <- function(biased){

  if(!inherits(biased, "SpatRaster")){
    return(data.frame(name = character(0), layer = character(0),
                      direction = character(0), stringsAsFactors = FALSE))
  }

  nms <- names(biased)
  parsed <- regmatches(nms, regexec("^(.+)_biased_(direct|inverse)$", nms))

  part <- function(i){
    vapply(parsed, function(x) if(length(x) == 3L) x[i] else NA_character_,
           character(1))
  }

  data.frame(name = nms, layer = part(2), direction = part(3),
             stringsAsFactors = FALSE)
}


#' Ellipsoids grouped by the bias applications they received
#'
#' Two ellipsoids share an emitted call only when their stored layers name the
#' same prediction layers applied in the same directions, and only when their
#' predictions live in the same emitted list, since the call subsets that list
#' by name. Ellipsoids with no biased layers, or whose layer names did not
#' parse, are dropped rather than guessed at.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A data frame with columns id, pred_obj, key and group, one row per
#' ellipsoid carrying a biased prediction.
#'
#' @noRd
report_bias_groups <- function(){

  empty <- data.frame(id = character(0), pred_obj = character(0),
                      key = character(0), group = integer(0),
                      list_obj = character(0), stringsAsFactors = FALSE)

  biased <- session_data$ellipsoid_prediction_list_biased
  if(length(biased) == 0) return(empty)

  pg <- report_pred_groups()
  if(nrow(pg) == 0) return(empty)

  ids <- names(biased)

  keys <- vapply(ids, function(id){
    ap <- report_bias_applied(biased[[id]])
    if(nrow(ap) == 0 || any(is.na(ap$layer))) return(NA_character_)
    paste(ap$layer, ap$direction, collapse = "|")
  }, character(1))

  pred_obj <- pg$list_obj[match(ids, pg$id)]

  keep <- !is.na(keys) & !is.na(pred_obj) &
    !vapply(ids, function(id) is.null(session_data$ellipsoid_list[[id]]),
            logical(1))

  ids <- ids[keep]
  if(length(ids) == 0) return(empty)

  keys <- keys[keep]
  pred_obj <- pred_obj[keep]

  full <- paste(keys, pred_obj, sep = "||")
  g <- match(full, unique(full))

  data.frame(id = ids, pred_obj = pred_obj, key = keys, group = g,
             list_obj = if(length(unique(full)) > 1){
               paste0("biased_predictions_", g)
             } else {
               "biased_predictions"
             },
             stringsAsFactors = FALSE)
}


#' Value summary of each biased layer, against the layer it came from
#'
#' The unbiased mean is carried alongside so the prose can report how much the
#' composite reduced the surface. The product of a surface and a composite in
#' [0, 1] is never larger than the surface, so this is always a reduction.
#'
#' @param biased One element of ellipsoid_prediction_list_biased.
#' @param pred The unbiased prediction for the same ellipsoid.
#'
#' @returns The table from report_bias_applied() with min, max, mean and
#' source_mean columns added.
#'
#' @noRd
report_bias_layer_stats <- function(biased, pred){

  ap <- report_bias_applied(biased)

  if(nrow(ap) == 0){
    return(cbind(ap, data.frame(min = numeric(0), max = numeric(0),
                                mean = numeric(0), source_mean = numeric(0))))
  }

  vals <- t(vapply(seq_len(nrow(ap)), function(i){

    g <- terra::global(biased[[ap$name[i]]], c("min", "max", "mean"),
                       na.rm = TRUE)

    src <- if(inherits(pred, "SpatRaster") && !is.na(ap$layer[i]) &&
              ap$layer[i] %in% names(pred)){
      terra::global(pred[[ap$layer[i]]], "mean", na.rm = TRUE)[1, 1]
    } else {
      NA_real_
    }

    c(as.numeric(g[1, "min"]), as.numeric(g[1, "max"]),
      as.numeric(g[1, "mean"]), as.numeric(src))
  }, numeric(4)))

  cbind(ap, data.frame(min = vals[, 1], max = vals[, 2], mean = vals[, 3],
                       source_mean = vals[, 4]))
}


# BIAS, CODE --------------------------------------------------------------

#' Code that reloads the bias layers
#'
#' The bundled raster is reachable through system.file() and so runs anywhere.
#' Uploaded files cannot be located from the app for the same reason the
#' environmental layers cannot, so the folder is declared as a variable. It is
#' declared here only when the Build chunk did not already declare it, which
#' happens when an example session uploaded its own bias layers.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A character vector of code lines.
#'
#' @noRd
report_bias_data_code <- function(){

  src <- session_data$bias_source

  if(identical(src, "example")){
    return(c("# Bundled example bias layers",
             "bias_rast <- terra::rast(",
             "  system.file(\"extdata/ma_biases.tif\", package = \"nicheR\"))"))
  }

  if(is.null(src) || length(src) == 0){
    src <- "your_bias_layers.tif"
  }

  # Build declares data_dir for every mode that reads uploaded files, so
  # repeating it there would overwrite a folder the user had already set
  head <- if(session_data$input_mode %in% c("example", "virtual")){
    c("# Set this to the folder holding the bias files you uploaded to the app",
      "data_dir <- \".\"",
      "")
  } else {
    character(0)
  }

  load <- if(length(src) == 1){
    paste0("bias_rast <- terra::rast(file.path(data_dir, \"", src, "\"))")
  } else {
    c("bias_rast <- terra::rast(file.path(data_dir,",
      paste0("                                   ", report_chr(src), "))"))
  }

  c(head, load)
}


#' Code that rebuilds the composite bias surface
#'
#' effect_direction is always written out. prepare_bias() declares its default
#' as c("direct", "inverse") and match.arg(several.ok = TRUE) returns both
#' elements when the argument is left alone, so the documented default does not
#' resolve to a single direction and relying on it errors.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @param prepared The stored prepare_bias() output.
#'
#' @returns A character vector of code lines.
#'
#' @noRd
report_bias_prepare_code <- function(prepared){

  parts <- report_bias_parts(prepared)
  outs <- report_bias_outputs(prepared)

  dirs <- if(nrow(parts) > 0){
    paste0("c(",
           paste0(report_name(parts$layer), " = \"", parts$direction, "\"",
                  collapse = ", "),
           ")")
  } else {
    "\"direct\""
  }

  args <- list(bias_surface = "bias_rast",
               effect_direction = dirs,
               include_composite = if(outs$include_composite) "TRUE" else "FALSE",
               include_processed_layers = if(outs$include_processed_layers) "TRUE" else "FALSE")

  mask <- session_data$bias_settings$mask_na

  # Nothing in the returned object records which overlap rule was used, so a
  # session saved before it was stored leaves the argument at its default and
  # says so rather than guessing
  pre <- character(0)
  if(length(mask) == 1 && !is.na(mask)){
    args$mask_na <- if(isTRUE(mask)) "TRUE" else "FALSE"
  } else {
    pre <- c("# This session did not record the NA overlap setting, so mask_na",
             "# is left at its default of FALSE, the union of layer extents.")
  }

  args$verbose <- "FALSE"

  c(pre, report_call("prepared_bias", "prepare_bias", args))
}


#' Code that applies the composite to one group of predictions
#'
#' A loop is honest here because the grouping is read from the stored layer
#' names rather than detected: two ellipsoids share this call only when the same
#' prediction layers were biased in the same directions, and only when their
#' predictions live in the same emitted list.
#'
#' @param objs Object names of the ellipsoids in the group.
#' @param ap The applied table shared by the group, from report_bias_applied().
#' @param pred_obj Name of the prediction list these ellipsoids sit in.
#' @param out Name to give the resulting list.
#'
#' @returns A character vector of code lines.
#'
#' @noRd
report_bias_apply_code <- function(objs, ap, pred_obj, out){

  calls <- unlist(lapply(seq_len(nrow(ap)), function(i){
    c("    apply_bias(prepared_bias = prepared_bias,",
      paste0("               prediction = ", pred_obj, "[[nm]],"),
      paste0("               prediction_layer = \"", ap$layer[i], "\","),
      paste0("               effect_direction = \"", ap$direction[i], "\","),
      paste0("               verbose = FALSE)[[\"", ap$layer[i], "_biased\"]]",
             if(i < nrow(ap)) "," else ""))
  }))

  c(paste0(out, "_targets <- ", report_chr(objs)),
    "",
    paste0(out, " <- lapply(setNames(", out, "_targets, ", out,
           "_targets), function(nm){"),
    "",
    "  parts <- list(",
    calls,
    "  )",
    "",
    "  terra::rast(parts)",
    "})")
}


# BIAS, PROSE -------------------------------------------------------------

#' Paragraphs introducing the bias step
#'
#' @param n_ell Number of ellipsoids that received a biased surface.
#' @param n_group Number of distinct application groups.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_bias_prose <- function(n_ell, n_group){

  out <- paste0(
    "A prediction describes where conditions suit the species. It does not ",
    "describe where anyone looked. Records accumulate near roads, towns and ",
    "reserves, so the places a species is recorded are the places that are ",
    "both suitable and sampled. This step builds a surface representing that ",
    "second component and combines it with the predictions, so that ",
    "occurrences drawn later carry the same uneven effort a real dataset ",
    "would.")

  if(n_ell > 0){
    out <- c(out, "",
             paste0("The composite was applied to ", n_ell, " niche",
                    if(n_ell == 1) "" else "s",
                    if(n_group > 1){
                      paste0(", in ", n_group, " groups, since not every niche ",
                             "was biased the same way.")
                    } else {
                      "."
                    },
                    " The result is a product of two surfaces and is no longer ",
                    "a probability. It is a relative weight, useful for ",
                    "ranking cells against one another and for sampling, not ",
                    "for reading off a value."))
  }

  c(out, "")
}


#' Paragraphs describing the composite bias surface
#'
#' @param prepared The stored prepare_bias() output.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_bias_prepare_prose <- function(prepared){

  parts <- report_bias_parts(prepared)
  n <- nrow(parts)

  lead <- paste0(
    "The bias surface was built from ", n, " layer", if(n == 1) "" else "s",
    ". Each was rescaled to run from 0 to 1, so layers measured on different ",
    "units contribute on the same scale, and then given a direction of ",
    "effect. A direct layer is one where higher values mean more sampling ",
    "effort; an inverse layer is one where higher values mean less, and is ",
    "replaced by 1 minus itself before being combined.")

  detail <- if(n > 0){
    paste0("- `", parts$layer, "` was treated as ", parts$direction,
           ", so ", ifelse(parts$direction == "direct",
                           "higher values raise the chance of a record",
                           "higher values lower the chance of a record"))
  } else {
    character(0)
  }

  f <- prepared$combination_formula

  combine <- paste0(
    "The layers were multiplied together into one composite surface",
    if(!is.null(f) && length(f) == 1 && nzchar(f)){
      paste0(", following `", f, "`")
    } else {
      ""
    },
    ". Multiplying means a cell needs to score well on every layer to end up ",
    "with a high weight: one layer near zero pulls the product near zero ",
    "regardless of the others.")

  mask <- session_data$bias_settings$mask_na

  na_txt <- if(length(mask) == 1 && !is.na(mask)){
    if(isTRUE(mask)){
      paste0("Where the layers did not overlap, the intersection was kept, so ",
             "any cell missing a value in any layer was left as missing in the ",
             "composite. This is the conservative choice: the composite covers ",
             "only the area every layer describes.")
    } else {
      paste0("Where the layers did not overlap, the union was kept, so a cell ",
             "missing a value in one layer was still given a weight from the ",
             "layers that did cover it. This keeps the composite over the full ",
             "area at the cost of some cells being informed by fewer layers ",
             "than others.")
    }
  } else {
    ""
  }

  stats <- report_bias_stats(prepared)

  stats_txt <- if(!is.null(stats)){
    c(paste0("The composite runs from ", report_num(signif(stats[["min"]], 4)),
             " to ", report_num(signif(stats[["max"]], 4)), ", with a mean of ",
             report_num(signif(stats[["mean"]], 4)),
             ". The mean is the number to compare against 1: the further below ",
             "it sits, the more of the study area the composite treats as ",
             "poorly sampled.",
             if(!is.na(stats[["na_share"]]) && stats[["na_share"]] > 0){
               paste0(" ", report_num(round(stats[["na_share"]] * 100, 1)),
                      "% of cells carry no weight at all and are excluded from ",
                      "anything sampled later.")
             } else {
               ""
             }))
  } else {
    character(0)
  }

  c(lead, "",
    if(length(detail) > 0) c(detail, ""),
    paste0(combine, if(nzchar(na_txt)) paste0(" ", na_txt) else ""), "",
    if(length(stats_txt) > 0) c(stats_txt, ""))
}


#' Paragraph describing what bias did to one niche
#'
#' @param ell The ellipsoid.
#' @param biased Its stored biased raster.
#' @param pred Its unbiased prediction.
#' @param obj Its object name in the emitted script.
#' @param out Name of the list the biased raster is stored in.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_bias_ell_prose <- function(ell, biased, pred, obj, out){

  st <- report_bias_layer_stats(biased, pred)

  if(nrow(st) == 0){
    return(c(paste0("No biased surface was stored for **", ell$ell_name, "**."),
             ""))
  }

  dirs <- unique(st$direction)

  open <- paste0(
    "The composite was applied to ", nrow(st), " layer",
    if(nrow(st) == 1) "" else "s", " of the **", ell$ell_name,
    "** prediction, ", report_and(paste0("`", st$layer, "`")), ". ",
    if(length(dirs) == 1 && identical(dirs, "direct")){
      paste0("Each was multiplied by the composite directly, so cells the ",
             "composite scores highly keep most of their value and cells it ",
             "scores near zero lose almost all of theirs.")
    } else if(length(dirs) == 1 && identical(dirs, "inverse")){
      paste0("Each was multiplied by 1 minus the composite, so the ",
             "relationship is reversed and cells the composite scores highly ",
             "are the ones pushed toward zero.")
    } else {
      paste0("The direction differed between layers, so they are not directly ",
             "comparable with one another.")
    })

  lines <- paste0(
    "- `", st$name, "` ranges from ", report_num(signif(st$min, 4)), " to ",
    report_num(signif(st$max, 4)), ", mean ", report_num(signif(st$mean, 4)),
    ifelse(is.na(st$source_mean) | st$source_mean == 0, "",
           paste0(", against a mean of ",
                  report_num(signif(st$source_mean, 4)),
                  " before bias, a reduction of ",
                  report_num(round((1 - st$mean / st$source_mean) * 100, 1)),
                  "%")))

  reduction <- paste0(
    "The reduction is the share of the surface the sampling term removed. It ",
    "is not a statement about the niche, which is unchanged: the same ",
    "conditions are still suitable, they are simply less likely to be ",
    "recorded.")

  c(open, "",
    lines, "",
    reduction, "",
    paste0("These surfaces are in `", out, "[[\"", obj, "\"]]`."), "")
}


# BIAS, SECTION -----------------------------------------------------------

#' Rmd source for the Bias section
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A character vector of Rmd lines.
#'
#' @noRd
report_bias_section <- function(){

  prepared <- session_data$prepared_bias

  # The tab is optional, is unavailable in virtual mode, and a session can be
  # saved with layers uploaded but never prepared
  if(identical(session_data$input_mode, "virtual") || is.null(prepared)){
    return(character(0))
  }

  ev <- report_bias_can_eval()
  grp <- report_bias_groups()

  head <- c("## Introducing sampling bias", "",
            report_bias_prose(nrow(grp), length(unique(grp$group))),
            report_chunk(c(report_bias_data_code(), "",
                           report_bias_prepare_code(prepared)),
                         label = "bias-prepare", eval = ev),
            report_bias_prepare_prose(prepared))

  if(nrow(grp) == 0){
    return(c(head,
             paste0("The composite was prepared but not applied to any ",
                    "prediction, so the surfaces carried into the next step ",
                    "are the unbiased ones."),
             ""))
  }

  biased <- session_data$ellipsoid_prediction_list_biased
  objs <- report_ell_objects()
  multi <- length(unique(grp$group)) > 1

  blocks <- unlist(lapply(sort(unique(grp$group)), function(g){

    sel <- which(grp$group == g)
    ids <- grp$id[sel]
    out <- grp$list_obj[match(ids[1], grp$id)]
    ap <- report_bias_applied(biased[[ids[1]]])

    open <- if(multi){
      c("***", "",
        paste0("### Bias application ", g), "",
        paste0("The following ", length(sel), " niche",
               if(length(sel) == 1) " was" else "s were",
               " biased together, since the same prediction layers were ",
               "combined with the composite in the same direction."), "")
    } else {
      character(0)
    }

    detail <- unlist(lapply(seq_along(ids), function(i){
      id <- ids[i]
      c(paste0(if(multi) "#### " else "### ",
               session_data$ellipsoid_list[[id]]$ell_name), "",
        report_bias_ell_prose(session_data$ellipsoid_list[[id]],
                              biased[[id]],
                              session_data$ellipsoid_prediction_list[[id]],
                              objs[[id]], out))
    }))

    c(open,
      report_chunk(report_bias_apply_code(objs[ids], ap,
                                          grp$pred_obj[sel[1]], out),
                   label = paste0("bias-apply-", g), eval = ev),
      detail)
  }))

  c(head, blocks)
}
# GENERATE, INSPECTION ----------------------------------------------------

#' Weighting method implied by a layer name
#'
#' The app derives the method from the layer name rather than storing it, in
#' generate_occ_for_ell() and again in the tab's method message. This has to
#' give the same answer or the emitted script samples with different weights
#' than the session did.
#'
#' @param layer Character vector of layer names.
#'
#' @returns A character vector, "mahalanobis" or "suitability".
#'
#' @noRd
report_occ_method <- function(layer){
  ifelse(grepl("mahalanobis", layer, ignore.case = TRUE),
         "mahalanobis", "suitability")
}


#' Whether an occurrence set was drawn from a biased surface
#'
#' generate_occ_for_ell() records this on the set as a "biased" attribute, and
#' that is read first because it is what the app itself decided. Membership in
#' the stored biased raster is the fallback for sets written before that
#' attribute existed. The two samplers take different arguments, so this
#' decides both the emitted call and what the prose may say about strategy.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @param id The ell_id.
#' @param layer The layer the set was drawn from.
#' @param df The occurrence set itself.
#'
#' @returns A single logical.
#'
#' @noRd
report_occ_is_biased <- function(id, layer, df = NULL){

  flag <- if(!is.null(df)) occ_meta(df, "biased", NA) else NA
  if(length(flag) == 1 && !is.na(flag)) return(isTRUE(flag))

  if(length(layer) != 1 || is.na(layer)) return(FALSE)

  b <- session_data$ellipsoid_prediction_list_biased[[id]]
  inherits(b, "SpatRaster") && layer %in% names(b)
}


#' Name of the emitted list holding the surface a set was drawn from
#'
#' Read from the group tables rather than rebuilt, so this indexes the object
#' the earlier sections actually created. Virtual sets have no surface and
#' return NA.
#'
#' @param id The ell_id.
#' @param biased Whether the set came from a biased surface.
#' @param mode The recorded mode, "raster" or "virtual".
#'
#' @returns A single string, or NA.
#'
#' @noRd
report_occ_source_obj <- function(id, biased, mode){

  if(identical(mode, "virtual")) return(NA_character_)

  tab <- if(isTRUE(biased)) report_bias_groups() else report_pred_groups()
  if(nrow(tab) == 0) return(NA_character_)

  i <- match(id, tab$id)
  if(is.na(i)) return(NA_character_)

  tab$list_obj[i]
}


#' Whether the chunk sampling one set can run at report time
#'
#' A sampling mask is an uploaded file, so a set that used one cannot be
#' reproduced without the user pointing at it and its chunk is emitted
#' unevaluated even when everything upstream ran. Evaluability stays monotonic:
#' a set drawn from a biased surface inherits the Bias answer, since its chunk
#' references an object that chunk created.
#'
#' @param biased Whether the set came from a biased surface.
#' @param mask The recorded mask file name, or "none".
#' @param mode The recorded mode.
#'
#' @returns A single logical.
#'
#' @noRd
report_occ_can_eval <- function(biased, mask, mode){

  up <- if(identical(mode, "virtual")){
    report_can_eval()
  } else if(isTRUE(biased)){
    report_bias_can_eval()
  } else {
    report_can_eval()
  }

  isTRUE(up) && identical(as.character(mask), "none")
}


#' Every occurrence set in the session, flattened
#'
#' The parameters travel with each set as attributes, so this step has better
#' provenance than any other: nothing is inferred from the result. Sets whose
#' ellipsoid has been deleted are dropped rather than reported against a
#' missing niche.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A data frame with one row per set.
#'
#' @noRd
report_occ_sets <- function(){

  empty <- data.frame(id = character(0), set = character(0),
                      layer = character(0), mode = character(0),
                      n_occ = numeric(0), n = numeric(0), seed = numeric(0),
                      sampling = character(0), strict = logical(0),
                      mask = character(0), truncate = logical(0),
                      effect = character(0), biased = logical(0),
                      source_obj = character(0), stringsAsFactors = FALSE)

  occ <- session_data$ellipsoid_occurrence_list
  if(length(occ) == 0) return(empty)

  rows <- lapply(names(occ), function(id){

    if(is.null(session_data$ellipsoid_list[[id]])) return(NULL)

    sets <- occ[[id]]
    if(length(sets) == 0) return(NULL)

    parts <- lapply(names(sets), function(k){

      df <- sets[[k]]
      if(is.null(df) || nrow(df) == 0) return(NULL)

      lyr <- as.character(occ_meta(df, "layer", NA_character_))
      md <- as.character(occ_meta(df, "mode", NA_character_))
      bias <- report_occ_is_biased(id, lyr, df)

      data.frame(id = id,
                 set = k,
                 layer = lyr,
                 mode = md,
                 n_occ = as.numeric(occ_meta(df, "n_occ", NA)),
                 n = nrow(df),
                 seed = as.numeric(occ_meta(df, "seed", NA)),
                 sampling = as.character(occ_meta(df, "sampling", NA)),
                 strict = as.logical(occ_meta(df, "strict", NA)),
                 mask = as.character(occ_meta(df, "mask", "none")),
                 truncate = as.logical(occ_meta(df, "truncate", NA)),
                 effect = as.character(occ_meta(df, "effect", NA)),
                 biased = bias,
                 source_obj = report_occ_source_obj(id, bias, md),
                 stringsAsFactors = FALSE)
    })

    parts <- Filter(Negate(is.null), parts)
    if(length(parts) == 0) return(NULL)

    do.call(rbind, parts)
  })

  rows <- Filter(Negate(is.null), rows)
  if(length(rows) == 0) return(empty)

  do.call(rbind, rows)
}


#' Occurrence sets grouped by the call that produced them
#'
#' occ_set_signature() already keys a set on every parameter that defines it
#' and does not include the ellipsoid, so two sets sharing a key were made in
#' the same run with the same settings and differ only in which niche they came
#' from. That is the app's own identity rule rather than a similarity test, so
#' a loop over the group is honest. The source object is folded into the key,
#' since two niches with identical settings can still sit in different emitted
#' prediction lists.
#'
#' Sets whose source object could not be resolved are dropped: their chunk
#' would index a list that does not exist.
#'
#' @returns The table from report_occ_sets() with a group column added, ordered
#' by group.
#'
#' @noRd
report_occ_groups <- function(){

  sets <- report_occ_sets()
  if(nrow(sets) == 0) return(cbind(sets, data.frame(group = integer(0))))

  keep <- identical(session_data$input_mode, "virtual") | !is.na(sets$source_obj)
  sets <- sets[keep, , drop = FALSE]
  if(nrow(sets) == 0) return(cbind(sets, data.frame(group = integer(0))))

  full <- paste(sets$set, sets$source_obj, sep = "||")
  sets$group <- match(full, unique(full))

  sets[order(sets$group), , drop = FALSE]
}


# GENERATE, CODE ----------------------------------------------------------

#' Code that reloads a sampling mask
#'
#' The mask is an uploaded file and cannot be located from the app, so the
#' script names it and leaves the folder for the user to set. data_dir is
#' declared here only when no earlier chunk declared it.
#'
#' Reads session_data from the enclosing server environment.
#'
#' @param mask The recorded mask file name.
#'
#' @returns A character vector of code lines, empty when no mask was used.
#'
#' @noRd
report_occ_mask_code <- function(mask){

  if(is.na(mask) || identical(as.character(mask), "none")){
    return(character(0))
  }

  head <- if(session_data$input_mode %in% c("example", "virtual")){
    c("# Set this to the folder holding the mask you uploaded to the app",
      "data_dir <- \".\"")
  } else {
    character(0)
  }

  c(head,
    paste0("sampling_mask <- terra::rast(file.path(data_dir, \"", mask, "\"))"),
    "")
}


#' Code that samples one group of occurrence sets
#'
#' Three samplers, chosen the way the app chooses them. Biased surfaces go to
#' sample_biased_data(), which takes no strategy and no method because the
#' surface values are the weights. Virtual sessions have no surface at all and
#' draw from the ellipsoid itself.
#'
#' @param objs Object names of the ellipsoids in the group.
#' @param row One row of report_occ_groups(), the shared settings.
#' @param out Name to give the resulting list.
#'
#' @returns A character vector of code lines.
#'
#' @noRd
report_occ_code <- function(objs, row, out){

  if(identical(row$mode, "virtual")){

    return(c(paste0(out, "_ellipsoids <- list(",
                    paste0(objs, " = ", objs, collapse = ", "), ")"),
             "",
             paste0(out, " <- lapply(", out, "_ellipsoids, function(e){"),
             "",
             "  d <- virtual_data(object = e,",
             paste0("                    n = ", report_num(row$n_occ), ","),
             paste0("                    truncate = ",
                    if(isTRUE(row$truncate)) "TRUE" else "FALSE", ","),
             paste0("                    effect = \"", row$effect, "\","),
             paste0("                    seed = ", report_num(row$seed), ")"),
             "",
             "  as.data.frame(d)",
             "})"))
  }

  mask_arg <- if(!is.na(row$mask) && !identical(as.character(row$mask), "none")){
    "sampling_mask"
  } else {
    NULL
  }

  strict_arg <- if(length(row$strict) == 1 && !is.na(row$strict)){
    if(isTRUE(row$strict)) "TRUE" else "FALSE"
  } else {
    "NULL"
  }

  args <- if(isTRUE(row$biased)){
    c(n_occ = report_num(row$n_occ),
      prediction = paste0(row$source_obj, "[[nm]]"),
      prediction_layer = paste0("\"", row$layer, "\""),
      if(!is.null(mask_arg)) c(sampling_mask = mask_arg),
      seed = report_num(row$seed),
      strict = strict_arg,
      verbose = "FALSE")
  } else {
    c(n_occ = report_num(row$n_occ),
      prediction = paste0(row$source_obj, "[[nm]]"),
      prediction_layer = paste0("\"", row$layer, "\""),
      sampling = paste0("\"", row$sampling, "\""),
      method = paste0("\"", report_occ_method(row$layer), "\""),
      if(!is.null(mask_arg)) c(sampling_mask = mask_arg),
      seed = report_num(row$seed),
      strict = strict_arg,
      verbose = "FALSE")
  }

  fn <- if(isTRUE(row$biased)) "sample_biased_data" else "sample_data"
  pad <- strrep(" ", nchar(fn) + 7)

  call <- c(paste0("  d <- ", fn, "(", names(args)[1], " = ", args[1], ","),
            paste0(pad, names(args)[-1], " = ", args[-1],
                   c(rep(",", length(args) - 2), ")")))

  c(paste0(out, "_targets <- ", report_chr(objs)),
    "",
    paste0(out, " <- lapply(setNames(", out, "_targets, ", out,
           "_targets), function(nm){"),
    "",
    call,
    "",
    "})")
}


# GENERATE, PROSE ---------------------------------------------------------

#' Sentence describing how sampling weights were formed
#'
#' One owner for the mapping, so the sentence and the emitted method argument
#' cannot describe different weights.
#'
#' @param row One row of report_occ_groups().
#'
#' @returns A single string.
#'
#' @noRd
report_occ_weight_prose <- function(row){

  if(isTRUE(row$biased)){
    return(paste0("The values of `", row$layer, "` were used directly as ",
                  "sampling weights, with no strategy applied. That surface ",
                  "is already the product of suitability and sampling effort, ",
                  "so a cell is drawn in proportion to both at once, and the ",
                  "resulting records carry the same uneven coverage a real ",
                  "dataset would."))
  }

  if(identical(row$mode, "virtual")){
    return(switch(as.character(row$effect),
                  "direct" = paste0("Points were weighted by the multivariate ",
                                    "normal density, so they concentrate near ",
                                    "the centroid and thin out toward the ",
                                    "boundary."),
                  "inverse" = paste0("Points were weighted by the complement ",
                                     "of that density, so they concentrate ",
                                     "toward the boundary and thin out near ",
                                     "the centroid."),
                  "uniform" = paste0("Every point inside the ellipsoid was ",
                                     "given equal weight, so they are spread ",
                                     "evenly through its volume."),
                  paste0("Points were drawn from the ellipsoid.")))
  }

  method <- report_occ_method(row$layer)

  if(identical(row$sampling, "random")){
    return(paste0("Every eligible cell was given equal weight, so the records ",
                  "say where the layer has a value at all rather than where it ",
                  "is high."))
  }

  if(method == "suitability"){
    if(identical(row$sampling, "centroid")){
      paste0("Weights were proportional to `", row$layer, "`, so cells with ",
             "higher values were more likely to be drawn and the records ",
             "concentrate where conditions sit closest to the centre of the ",
             "niche.")
    } else {
      paste0("Weights were proportional to 1 minus `", row$layer, "`, so cells ",
             "with lower values were more likely to be drawn and the records ",
             "concentrate near the boundary of the niche rather than its ",
             "centre.")
    }
  } else {
    if(identical(row$sampling, "centroid")){
      paste0("Weights were inversely proportional to `", row$layer,
             "`, a distance, so the smallest distances were the most likely to ",
             "be drawn and the records concentrate near the centre of the ",
             "niche.")
    } else {
      paste0("Weights were proportional to `", row$layer, "`, a distance, so ",
             "the largest distances were the most likely to be drawn and the ",
             "records concentrate near the boundary of the niche.")
    }
  }
}


#' Paragraphs introducing the generate step
#'
#' @param sets The table from report_occ_groups().
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_generate_prose <- function(sets){

  n_set <- nrow(sets)
  n_ell <- length(unique(sets$id))
  n_group <- length(unique(sets$group))
  virtual <- identical(session_data$input_mode, "virtual")

  out <- paste0(
    "A niche says which conditions are tolerable and a prediction says where ",
    "they occur. Neither is a dataset. This step draws point records from ",
    "what came before, so the result is an occurrence table whose origin is ",
    "known: the niche it came from, the rule that selected the points, and ",
    "the seed that makes the selection repeatable. That is what makes these ",
    "records useful as a benchmark, since anything fitted to them can be ",
    "compared against the niche that produced them.")

  out <- c(out, "",
           paste0(n_set, " occurrence set", if(n_set == 1) "" else "s",
                  " were drawn across ", n_ell, " niche",
                  if(n_ell == 1) "" else "s",
                  if(n_group > 1){
                    paste0(", in ", n_group,
                           " groups, since not every set was drawn the same way.")
                  } else {
                    "."
                  },
                  if(virtual){
                    paste0(" This session has no environmental layers, so the ",
                           "points are environmental values drawn from the ",
                           "ellipsoid itself and carry no coordinates.")
                  } else {
                    paste0(" Each set holds the coordinates of the cells that ",
                           "were drawn.")
                  }))

  c(out, "")
}


#' Paragraphs describing one group of occurrence sets
#'
#' Members of a group differ only in which niche they came from, so they are
#' listed rather than given a heading each. What varies between them is the
#' count, and only when a sampler returned fewer points than were asked for.
#'
#' @param row One row of report_occ_groups(), the shared settings.
#' @param members The rows of the group.
#' @param objs Object names of the ellipsoids in the group.
#' @param out Name of the list the sets are stored in.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_occ_group_prose <- function(row, members, objs, out){

  virtual <- identical(row$mode, "virtual")

  open <- paste0(
    report_num(row$n_occ), " point", if(row$n_occ == 1) "" else "s",
    " were drawn per niche",
    if(virtual){
      paste0(" directly from the ellipsoid, using a random seed of ",
             report_num(row$seed), ".")
    } else {
      paste0(" from `", row$layer, "`, using a random seed of ",
             report_num(row$seed), ".")
    },
    " The seed is what makes this repeatable: running the code below returns ",
    "the same points rather than a fresh draw.")

  weights <- report_occ_weight_prose(row)

  trunc_txt <- if(virtual && length(row$truncate) == 1 && !is.na(row$truncate)){
    if(isTRUE(row$truncate)){
      paste0("Points were constrained to fall inside the confidence boundary, ",
             "so none sit outside the niche the ellipsoid describes.")
    } else {
      paste0("Points were not constrained to the confidence boundary, so some ",
             "sit outside the niche the ellipsoid describes. They are drawn ",
             "from the distribution the ellipsoid was built from rather than ",
             "from its interior.")
    }
  } else {
    ""
  }

  strict_txt <- if(!virtual && length(row$strict) == 1 && !is.na(row$strict)){
    if(isTRUE(row$strict)){
      paste0("Cells with no value and cells equal to zero were removed before ",
             "sampling, so every record falls somewhere the layer treats as ",
             "inside the niche.")
    } else {
      paste0("Cells equal to zero were kept as candidates, so records can fall ",
             "outside the niche boundary wherever the layer still carries a ",
             "value there.")
    }
  } else {
    ""
  }

  mask_txt <- if(!virtual && !is.na(row$mask) &&
                 !identical(as.character(row$mask), "none")){
    paste0("Sampling was restricted to the area covered by `", row$mask,
           "`, so the records describe that area rather than the full study ",
           "extent.")
  } else {
    ""
  }

  short <- members[!is.na(members$n_occ) & members$n < members$n_occ, ,
                   drop = FALSE]

  counts <- if(nrow(short) > 0){
    c("", paste0("- ", vapply(seq_len(nrow(short)), function(i){
      paste0(session_data$ellipsoid_list[[short$id[i]]]$ell_name, " returned ",
             report_num(short$n[i]), " of the ", report_num(short$n_occ[i]),
             " points requested")
    }, character(1))), "")
  } else {
    character(0)
  }

  tail_txt <- paste0("These records are in `", out, "`, named after the niche ",
                     "each set came from, so one table is pulled out with, for ",
                     "example, `", out, "[[\"", objs[1], "\"]]`.")

  body <- paste(c(open, weights, trunc_txt, strict_txt, mask_txt)[
    nzchar(c(open, weights, trunc_txt, strict_txt, mask_txt))], collapse = " ")

  c(body, "", counts, tail_txt, "")
}


# GENERATE, SECTION -------------------------------------------------------

#' Rmd source for the Generate section
#'
#' Reads session_data from the enclosing server environment.
#'
#' @returns A character vector of Rmd lines.
#'
#' @noRd
report_generate_section <- function(){

  sets <- report_occ_groups()
  if(nrow(sets) == 0) return(character(0))

  objs_all <- report_ell_objects()
  groups <- unique(sets$group)
  multi <- length(groups) > 1

  blocks <- unlist(lapply(groups, function(g){

    members <- sets[sets$group == g, , drop = FALSE]
    row <- members[1, , drop = FALSE]
    objs <- objs_all[members$id]
    out <- if(multi) paste0("occurrences_", g) else "occurrences"

    open <- if(multi){
      c("***", "",
        paste0("### Occurrence set ", g), "")
    } else {
      character(0)
    }

    c(open,
      report_chunk(c(report_occ_mask_code(row$mask),
                     report_occ_code(objs, row, out)),
                   label = paste0("generate-", g),
                   eval = report_occ_can_eval(row$biased, row$mask, row$mode)),
      report_occ_group_prose(row, members, objs, out))
  }))

  c("## Generating occurrence records", "",
    report_generate_prose(sets),
    blocks)
}

# CITATION ----------------------------------------------------------------

#' Citation block for the report
#'
#' Built from citation() rather than written out, so a change to inst/CITATION
#' or DESCRIPTION propagates without editing this file.
#'
#' @returns A character vector of markdown lines.
#'
#' @noRd
report_citation <- function(){

  cit <- tryCatch(format(utils::citation("nicheR"), style = "text"),
                  error = function(e) NULL)

  c("## Citation", "",
    "This analysis was produced with the nicheR R package. If you use these",
    "results in a publication, please cite:", "",
    if(!is.null(cit)){
      paste0("> ", unlist(strsplit(cit, "\n")))
    } else {
      "> Castaneda-Guzman M., Hughes C., Paansri P., Cobos M. E. (2026). nicheR."
    },
    "",
    "You can regenerate this citation at any time with `citation(\"nicheR\")`.",
    "")
}
