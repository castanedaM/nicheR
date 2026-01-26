#' Computes ranges from a data.frame of values
#'
#' @param data a data.frame of at least two columns. Each column should
#' contain numeric values.
#'
#' @return
#' A data.frame with the minimum and maximum values of each variable.
#'
#' @export
#'
#' @examples
#' data <- data.frame(var1 = runif(10), var2 = runif(10), var3 = runif(10))
#' ranges_from_data(data)

ranges_from_data <- function(data) {
  if (!is.data.frame(data) || ncol(data) < 1) {
    stop("'data' must be a 'data.frame' with at least one column.")
  }

  ranges <- sapply(data, range, na.rm = TRUE)

  ranges_df <- as.data.frame(ranges)
  rownames(ranges_df) <- c("min", "max")

  return(ranges_df)
}
