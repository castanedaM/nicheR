#' Launch the nicheR Shiny application
#'
#' Opens a graphical interface to the nicheR workflow: building ellipsoids
#' from environmental ranges, predicting suitability across a study area,
#' preparing and applying sampling bias, and generating virtual data.
#'
#' The app is an interface to the package functions, not a separate
#' implementation. Every control maps onto a call that can be made directly
#' from the console, and sessions can be saved and reloaded as `.rds` files.
#' Environmental data can be supplied as rasters or as a table, or defined
#' directly in environmental space with no geographic component.
#'
#' The app runs in the current R session and blocks the console until it is
#' closed.
#'
#' @usage run_gui( )
#'
#' @returns Invisibly returns `NULL`. Called for its side effect of running
#' the application.
#'
#' @details Uploaded rasters and tables are held in memory for the duration of
#' the session and are not written to disk unless a session or a figure is
#' explicitly saved. For a walkthrough of the interface, see
#' `vignette("shiny_app_vignette", package = "nicheR")`.
#'
#' @seealso [build_ellipsoid()], [predict()], [prepare_bias()],
#' [sample_data()]
#'
#' @examples
#' \dontrun{
#' run_gui()
#' }
#'
#' @export
run_gui <- function( ){

  requireNamespace("shiny")

  app_dir <- system.file("shiny", package = "nicheR")

  if(app_dir == ""){
    stop("Could not find the shiny directory. Try re-installing nicheR.",
         call. = FALSE)
  }

  shiny::runApp(app_dir)

  invisible(NULL)
}
