#' Create a data.frame of ranges from a data.frame of values
#'
#' @param data a data.frame of at least one column that contains variable values.
#' @param var_names a character vector of variable names to be used. If NULL (default),
#' all variables are used.
#'
#' @return a data.frame with the minimum and maximum values of each variable.
#'
#' @export
#'
#' @examples
#' data <- data.frame(var1 = runif(10), var2 = runif(10), var3 = runif(10))
#' ranges_from_data(data, var_names = c("var1", "var2"))
ranges_from_data <- function(data, var_names = NULL) {
  if (!is.data.frame(data) || ncol(data) < 1) {
    stop("data must be a data.frame with at least one column")
  }

  if (!is.null(var_names)) {
    if (!all(var_names %in% colnames(data))) {
      stop("all var_names must be present in the data column names")
    }
    data <- data[, var_names, drop = FALSE]
  }

  ranges <- sapply(data, range, na.rm = TRUE)

  ranges_df <- as.data.frame(ranges)
  rownames(ranges_df) <- c("min", "max")

  return(ranges_df)
}
