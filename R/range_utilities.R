#' Compute variable ranges from data or statistics with optional expansion
#'
#' These functions compute the minimum and maximum values for variables,
#' either directly from a dataset or based on normal distribution parameters,
#' with the ability to expand the resulting ranges by a percentage.
#'
#' @usage
#' ranges_from_data(data, expand_min = NULL, expand_max = NULL)
#' ranges_from_stats(mean, sd, cl = 95, expand_min = NULL, expand_max = NULL)
#'
#' @param data A data.frame of at least two columns. Each column should
#' contain numeric values.
#' @param mean A named numeric vector of mean values for each variable.
#' @param sd A named numeric vector of standard deviation values for each
#' variable. Names must match those in `mean`.
#' @param cl A numeric value indicating the confidence level (default 95).
#' @param expand_min A named vector or list of percentages (e.g., 10 for 10%)
#' defining how much to expand the minimum value of specific variables.
#' @param expand_max A named vector or list of percentages defining how
#' much to expand the maximum value of specific variables.
#'
#' @return
#' A data.frame with the (potentially expanded) minimum and maximum
#' values of each variable.
#'
#' @name range_utilities
NULL

#' @rdname range_utilities
#' @export
#' @examples
#' # From data
#' df <- data.frame(var1 = c(0, 10), var2 = c(50, 100))
#' ranges_from_data(df, expand_min = list(var1 = 10),
#'                  expand_max = list(var2 = 20))

ranges_from_data <- function(data, expand_min = NULL, expand_max = NULL) {
  # Initial structural check
  if (!is.data.frame(data) || ncol(data) < 2) {
    stop("'data' must be a 'data.frame' with at least two columns.")
  }

  # Error checking: Verify that expansion names match at least one data column
  col_names <- colnames(data)

  if (!is.null(expand_min)) {
    if (!any(names(expand_min) %in% col_names)) {
      stop("None of the names in 'expand_min' match the names in 'data'.")
    }
  }
  
  if (!is.null(expand_max)) {
    if (!any(names(expand_max) %in% col_names)) {
      stop("None of the names in 'expand_max' match the names in 'data'.")
    }
  }

  # Compute ranges
  ranges <- sapply(data, range, na.rm = TRUE)
  
  for (col_name in colnames(ranges)) {
    current_min <- ranges[1, col_name]
    current_max <- ranges[2, col_name]
    abs_range <- current_max - current_min
    
    if (!is.null(expand_min) && col_name %in% names(expand_min)) {
      ranges[1, col_name] <- current_min - 
        (abs_range * expand_min[[col_name]] / 100)
    }
    if (!is.null(expand_max) && col_name %in% names(expand_max)) {
      ranges[2, col_name] <- current_max + 
        (abs_range * expand_max[[col_name]] / 100)
    }
  }

  ranges_df <- as.data.frame(ranges)
  rownames(ranges_df) <- c("min", "max")
  return(ranges_df)
}

#' @rdname range_utilities
#' @export
#' @examples
#' # From statistics
#' m <- c(var1 = 10, var2 = 100)
#' s <- c(var1 = 2, var2 = 15)
#' ranges_from_stats(mean = m, sd = s, cl = 95,
#'                   expand_min = list(var1 = 10))

ranges_from_stats <- function(mean, sd, cl = 95, expand_min = NULL,
                              expand_max = NULL) {
 # Check that mean and sd are numeric
  if (!is.numeric(mean) || !is.numeric(sd)) {
    stop("'mean' and 'sd' must be numeric vectors.")
  }

  # Check that mean and sd share names and length
  if (!all(names(mean) %in% names(sd)) || length(mean) != length(sd)) {
    stop("'mean' and 'sd' must be named vectors sharing the same length and names.")
  }

  # Check that at least one name in expand arguments matches names in mean/sd
  var_names <- names(mean)
  
  if (!is.null(expand_min)) {
    if (!any(names(expand_min) %in% var_names)) {
      stop("None of the names in 'expand_min' match the names in 'mean'.")
    }
  }
  
  if (!is.null(expand_max)) {
    if (!any(names(expand_max) %in% var_names)) {
      stop("None of the names in 'expand_max' match the names in 'mean'.")
    }
  }

  # Compute ranges
  alpha <- 1 - (cl / 100)
  p_lower <- alpha / 2
  p_upper <- 1 - (alpha / 2)

  mins <- qnorm(p_lower, mean = mean, sd = sd)
  maxs <- qnorm(p_upper, mean = mean, sd = sd)
  
  var_names <- names(mean)
  ranges <- matrix(c(mins, maxs), nrow = 2, byrow = TRUE)
  colnames(ranges) <- var_names
  
  for (var in var_names) {
    current_min <- ranges[1, var]
    current_max <- ranges[2, var]
    abs_range <- current_max - current_min
    
    if (!is.null(expand_min) && var %in% names(expand_min)) {
      ranges[1, var] <- current_min - (abs_range * expand_min[[var]] / 100)
    }
    if (!is.null(expand_max) && var %in% names(expand_max)) {
      ranges[2, var] <- current_max + (abs_range * expand_max[[var]] / 100)
    }
  }

  ranges_df <- as.data.frame(ranges)
  rownames(ranges_df) <- c("min", "max")
  return(ranges_df)
}