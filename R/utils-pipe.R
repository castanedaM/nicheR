#' Internal imports for NicheR
#'
#' Centralized roxygen imports to keep NAMESPACE tidy.
#' Nothing in this file is exported.
#'
#' @keywords internal
#' @name NicheR-internal-imports
NULL

# ---- Core plotting ----
#' @import ggplot2
NULL

# dplyr verbs used unqualified
#' @importFrom dplyr mutate filter row_number
NULL

# pipe and pronoun
#' @importFrom magrittr %>%
#' @importFrom rlang .data
NULL

# palettes, arranging, basemap, sf
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr ggarrange
#' @importFrom sf st_as_sf
NULL

# terra/raster utilities (keep comments on their own line!)
#' @importFrom terra ncell as.data.frame ext nlyr xmin xmax ymin ymax crop yFromRow
#' @importFrom raster ncell extent
NULL

# Database and SQL utilities used in convert_large_raster()
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable dbReadTable
#' @importFrom RSQLite SQLite
NULL
