#' Size aware conversion of raster stacks to data.frame for NicheR
#'
#' Converts a raster stack to a \code{data.frame} with \code{x} and \code{y}
#' coordinates and one column per raster layer. For objects with an estimated
#' in memory size below a threshold, it wraps \code{terra::as.data.frame()}.
#' For larger objects, it uses a chunked SQLite workflow to avoid loading the
#' entire raster into memory. A single temporary cache file is used to reuse
#' the most recent large raster conversion within the current R session.
#'
#' @param raster_stack A multi layer raster object. Supports
#'   \code{terra::SpatRaster} or \code{raster::Raster*}. Must have valid
#'   georeferencing so that \code{xmin}, \code{xmax}, \code{ymin}, and
#'   row to y conversions are defined.
#' @param out_filename Character base name (without extension) for on disk
#'   outputs. Files \code{<out_filename>.sqlite} and \code{<out_filename>.RData}
#'   are written into the current working directory when the chunked workflow
#'   is used or when \code{write_rdata = TRUE}.
#' @param chunk_height Positive integer. Number of raster rows per chunk to
#'   process and append to the SQL table. Smaller values reduce memory usage
#'   but increase the number of database writes.
#' @param keep_sql Logical. If \code{TRUE}, the temporary SQLite file is kept
#'   on disk after processing. If \code{FALSE} it is deleted at the end.
#' @param write_rdata Logical. If \code{TRUE}, saves the resulting
#'   \code{data.frame} as \code{<out_filename>.RData} (object name
#'   \code{env_df}) in addition to returning it.
#' @param verbose Logical. If \code{TRUE}, prints progress messages and file
#'   paths.
#' @param size_threshold_mb Numeric. Estimated size in megabytes above which
#'   the function switches from in memory conversion to the chunked SQLite
#'   workflow.
#' @param use_cache Logical. If \code{TRUE} and a NicheR temporary cache file
#'   exists in \code{tempdir()}, the function loads \code{env_df} from that
#'   file instead of recomputing the large raster conversion. The cache is
#'   always overwritten by the most recent large raster call.
#'
#' @return A \code{data.frame} with columns \code{x}, \code{y}, and one column
#'   per raster layer (matching \code{names(raster_stack)}). If
#'   \code{write_rdata = TRUE}, the same object is also saved as
#'   \code{<out_filename>.RData} (object name \code{env_df}). When the chunked
#'   workflow is used and \code{keep_sql = TRUE}, \code{<out_filename>.sqlite}
#'   is left on disk with table \code{raster_table}. Independently, the most
#'   recent large raster conversion is stored in a temporary cache file inside
#'   \code{tempdir()}.
#'
#' @export
as.data.frame.nicheR <- function(raster_stack,
                                 out_filename      = "NicheR_db",
                                 chunk_height      = 500,
                                 keep_sql          = FALSE,
                                 write_rdata       = FALSE,
                                 verbose           = FALSE,
                                 size_threshold_mb = 5000,
                                 use_cache         = FALSE) {

  gc()

  if (missing(raster_stack)) {
    stop("raster_stack is required")
  }


  # Estimate size in MB (assuming 4-byte numerics)
  ncell_total <- ncell(raster_stack)
  nlyr_total  <- nlyr(raster_stack)
  est_size_mb <- ncell_total * nlyr_total * 4 / 1024^2

  if (verbose) {
    message("Estimated raster size ~ ", round(est_size_mb, 2), " MB (",
            ncell_total, " cells x ", nlyr_total, " layers)")
  }

  ## Small case in memory
  if (est_size_mb <= size_threshold_mb) {
    if (verbose) message("Using in-memory terra::as.data.frame (below threshold).")

    env_df <- terra::as.data.frame(raster_stack, xy = TRUE, na.rm = TRUE)

    if (write_rdata) {
      rdata_file <- paste0(out_filename, ".RData")
      if (verbose) message("Saving RData file: ", rdata_file)
      save(env_df, file = rdata_file)
    }

    return(env_df)
  }


  # Internal temp cache file for large-raster conversion
  temp_cache_file <- file.path(tempdir(), "NicheR_temp_env_df.RData")

  # Cache check: if temp RData exists, just load and return
  if (use_cache && file.exists(temp_cache_file)) {
    if (verbose) message("Found NicheR temp cache, loading: ", temp_cache_file)
    loaded <- load(temp_cache_file)  # should create env_df

    if (!"env_df" %in% loaded) {
      stop("Temp cache file '", temp_cache_file,
           "' did not contain an object named 'env_df'.")
    }

    return(env_df)
  }


  ## Large case: chunked SQLite workflow
  if (verbose) {
    message("Using chunked SQLite workflow (above threshold of ",
            size_threshold_mb, " MB).")
  }

  sqlite_file <- paste0(out_filename, ".sqlite")
  rdata_file  <- paste0(out_filename, ".RData")

  conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = sqlite_file)
  on.exit({
    try(DBI::dbDisconnect(conn), silent = TRUE)
  }, add = TRUE)

  # Create empty table
  empty_df <- as.data.frame(matrix(ncol = (nlyr_total + 2), nrow = 0))
  colnames(empty_df) <- c("x", "y", names(raster_stack))
  DBI::dbWriteTable(conn, "raster_table", empty_df, overwrite = TRUE)

  total_rows <- nrow(raster_stack)

  for (start_row in seq(1, total_rows, by = chunk_height)) {
    end_row <- min(start_row + chunk_height - 1, total_rows)

    if (verbose) {
      message("Processing rows ", start_row, " to ", end_row, " of ", total_rows)
    }

    # row -> y coordinates
    y_top <- terra::yFromRow(raster_stack, start_row)

    if (end_row == total_rows) {
      y_bottom <- terra::ymin(raster_stack)
    } else {
      y_bottom <- terra::yFromRow(raster_stack, end_row + 1)
    }

    e <- terra::ext(terra::xmin(raster_stack),
                    terra::xmax(raster_stack),
                    y_bottom, y_top)

    r_block  <- terra::crop(raster_stack, e)
    df_block <- terra::as.data.frame(r_block, xy = TRUE, na.rm = TRUE)

    if (nrow(df_block) > 0) {
      DBI::dbWriteTable(conn, "raster_table", df_block, append = TRUE)
    }
  }

  env_df <- DBI::dbReadTable(conn, "raster_table")

  if (write_rdata) {
    if (verbose) message("Saving RData file: ", rdata_file)
    save(env_df, file = rdata_file)
  }

  if (isFALSE(keep_sql)) {
    if (verbose) message("Deleting SQLite file: ", sqlite_file)
    try(file.remove(sqlite_file), silent = TRUE)
  } else if (verbose) {
    message("SQLite file kept: ", sqlite_file)
  }

  if (verbose) {
    message("Processing complete. Data frame has ",
            nrow(env_df), " rows and ", ncol(env_df), " columns")
  }

  # Always update internal temp cache for large rasters
  if (verbose) message("Updating NicheR temp cache: ", temp_cache_file)
  save(env_df, file = temp_cache_file)

  gc()

  return(env_df)
}
