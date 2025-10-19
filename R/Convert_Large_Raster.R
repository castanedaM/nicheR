library(RSQLite)
library(DBI)


convert_large_raster <- function(raster_stack,
                                 out_filename = "raster_db",
                                 chunk_height = 500,
                                 keep_sql = FALSE,
                                 verbose = TRUE) {

  # Create database connection
  sqlite_file <- paste0(out_filename, ".sqlite")
  rdata_file <- paste0(out_filename, ".RData")

  conn <- dbConnect(SQLite(), dbname = sqlite_file)

  # Create empty table
  empty_df <- as.data.frame(matrix(ncol = (nlyr(raster_stack) + 2), nrow = 0))
  colnames(empty_df) <- c("x", "y", names(raster_stack))
  dbWriteTable(conn, "raster_table", empty_df, overwrite = TRUE)

  # Get total rows for chunking
  total_rows <- nrow(raster_stack)

  # Loop over chunks
  for (start_row in seq(1, total_rows, by = chunk_height)) {
    end_row <- min(start_row + chunk_height - 1, total_rows)
    if(verbose) message("Processing rows ", start_row, " to ", end_row, " of ", total_rows)

    # Convert row numbers to y-coordinates
    y_top <- yFromRow(raster_stack, start_row)

    # Handle last chunk - use ymin instead of going beyond raster extent
    if (end_row == total_rows) {
      y_bottom <- ymin(raster_stack)
    } else {
      y_bottom <- yFromRow(raster_stack, end_row + 1)
    }

    # Build extent in spatial coordinates
    e <- ext(xmin(raster_stack), xmax(raster_stack), y_bottom, y_top)

    # Crop raster for this block
    r_block <- crop(raster_stack, e)

    # Convert to dataframe
    df_block <- as.data.frame(r_block, xy = TRUE, na.rm = TRUE)

    # Append to database
    if (nrow(df_block) > 0) {
      dbWriteTable(conn, "raster_table", df_block, append = TRUE)
    }
  }

  # Read final table from database
  env_df <- dbReadTable(conn, "raster_table")

  # Disconnect from database
  dbDisconnect(conn)

  # Save as RData
  if(verbose) message("Saving RData file: ", rdata_file)
  save(env_df, file = rdata_file)

  # Delete SQL file if requested
  if (!keep_sql) {
    if(verbose) message("Deleting SQLite file: ", sqlite_file)
    file.remove(sqlite_file)
  } else {
    if(verbose) message("SQLite file kept: ", sqlite_file)
  }

  if(verbose) message("Processing complete. Data frame has ", nrow(env_df), " rows and ", ncol(env_df), " columns")

  return(invisible(env_df))
}
