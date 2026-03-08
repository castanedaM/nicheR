#' Save a nicheR object to disk
#'
#' @description
#' A wrapper around \code{\link[base]{saveRDS}} to save \code{nicheR_ellipsoid}
#' or \code{nicheR_community} objects to a file. Includes a safety check
#' for overwriting existing files and ensures the file extension is .rds.
#'
#' @param object A \code{nicheR_ellipsoid} or \code{nicheR_community} object.
#' @param file Character. The connection or name of the file where the
#'   object will be saved (usually ending in ".rds").
#' @param overwrite Logical. If \code{TRUE}, an existing file at the
#'   specified path will be replaced. Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[base]{saveRDS}}.
#'
#' @return
#' No return value. Saves the object to the specified path.
#'
#' @export

save_nicheR <- function(object, file, overwrite = FALSE, ...) {

  if (missing(object)) {
    stop("Argument 'object' is required.")
  }
  if (!inherits(object, c("nicheR_ellipsoid", "nicheR_community"))) {
    warning("The object is not a 'nicheR' class.")
  }
  if (missing(file)) {
    stop("Argument 'file' is required.")
  }
  if (!grepl("\\.rds$", tolower(file))) {
    file <- paste0(file, ".rds")
  }
  if (file.exists(file) && !overwrite) {
    stop("File '", file, "' already exists. Set 'overwrite' = TRUE.")
  }

  saveRDS(object = object, file = file, ...)
}


#' Read a nicheR object from disk
#'
#' @description
#' A wrapper around \code{\link[base]{readRDS}} to load saved nicheR
#' objects back into the R environment.
#'
#' @param file Character. The path to the file to be read.
#'
#' @return
#' The saved \code{nicheR_ellipsoid} or \code{nicheR_community} object.
#'
#' @seealso
#' \code{\link{save_nicheR}}
#'
#' @export

read_nicheR <- function(file) {

  if (missing(file)) {
    stop("Argument 'file' is required.")
  }
  if (!file.exists(file)) {
    stop("The file '", file, "' does not exist.")
  }

  obj <- readRDS(file)

  if (!inherits(obj, c("nicheR_ellipsoid", "nicheR_community"))) {
    warning("The loaded object does not appear to be a 'nicheR' class.")
  }

  return(obj)
}