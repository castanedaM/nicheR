#' Pipe operator
#'
#' This re-exports the pipe operator from \code{magrittr} to make it available
#' for use in your package's functions and for package users.
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @export
#' @importFrom dplyr filter select mutate arrange group_by summarise
#' @importFrom magrittr %>%
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL
