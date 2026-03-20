#' Background environmental data for examples
#'
#' A dataset containing geographic coordinates and two bioclimatic variables
#' used as a background point cloud for generating and testing niche ellipsoids.
#'
#' @format A data frame with 12,396 rows and 4 variables:
#' \describe{
#'   \item{x}{Longitude in decimal degrees (WGS84).}
#'   \item{y}{Latitude in decimal degrees (WGS84).}
#'   \item{bio_1}{Annual Mean Temperature (°C).}
#'   \item{bio_5}{Max Temperature of Warmest Month (°C).}
#'   \item{bio_6}{Min Temperature of Coldest Month (°C).}
#'   \item{bio_7}{Temperature Annual Range (bio_5 - bio_6) (°C).}
#'   \item{bio_12}{Annual Precipitation (mm).}
#'   \item{bio_13}{Precipitation of Wettest Month (mm).}
#'   \item{bio_14}{Precipitation of Driest Month (mm).}
#'   \item{bio_15}{Precipitation Seasonality (Coefficient of Variation).}
#' }
#'
#' @details
#' This dataset represents an irregular point cloud typical of environmental
#' background data used in ecological niche modeling (ENM). It is primarily
#' used in the `nicheR` package to provide the environmental space for
#' functions like \code{\link{conserved_ellipses}}.
#'
#' @source \url{https://www.worldclim.org}
#' @examples
#' data(back_data)
#' head(back_data)
"back_data"



#' Reference ellipse for virtual community examples
#'
#' A pre-calculated \code{nicheR_ellipsoid} object representing a hypothetical
#' species niche based on Annual Mean Temperature (bio_1) and Annual
#' Precipitation (bio_12).
#'
#' @format An object of class \code{nicheR_ellipsoid} (which is a \code{list})
#' with 13 elements:
#' \describe{
#'   \item{dimensions}{Integer. Number of dimensions (2).}
#'   \item{var_names}{Character vector. Names of variables (\code{"bio_1"}, \code{"bio_12"}).}
#'   \item{centroid}{Named numeric vector. The center of the niche (\eqn{\mu}).}
#'   \item{cov_matrix}{Matrix. The \eqn{2 \times 2} covariance matrix (\eqn{\Sigma}).}
#'   \item{Sigma_inv}{Matrix. The precision matrix (inverse covariance).}
#'   \item{chol_Sigma}{Matrix. Cholesky decomposition of the covariance.}
#'   \item{eigen}{List. Eigenvectors and eigenvalues of the covariance.}
#'   \item{cl}{Numeric. Confidence level used (e.g., 0.99).}
#'   \item{chi2_cutoff}{Numeric. The chi-square quantile for the given \code{cl}.}
#'   \item{semi_axes_lengths}{Numeric vector. Radii of the ellipsoid axes.}
#'   \item{axes_coordinates}{List. Vertices (endpoints) for each ellipsoid axis.}
#'   \item{volume}{Numeric. The hyper-volume of the ellipsoid.}
#'   \item{cov_limits}{List. Axis-aligned minimum and maximum limits.}
#' }
#'
#' @details
#' This object serves as a template for testing community simulation functions
#' like \code{\link{conserved_ellipses}}. It was generated using
#' \code{\link{build_ellipsoid}} with a centroid at (23.75, 1750) and
#' specific covariance structures to reflect a typical temperature-precipitation
#' relationship.
#'
#' @examples
#' data(ref_ellipse)
#' print(ref_ellipse)
#'
#' # Access the volume
#' ref_ellipse$volume
"ref_ellipse"



#' Bioclimatic variables for part of the Americas
#'
#' A \code{SpatRaster} containing 8 bioclimatic variables representing
#' present-day climatic conditions for an area that covers parts of South
#' and North America. Variables were obtained at a 10 arc-minute resolution.
#' Sourced from WorldClim 2.1: \url{https://worldclim.org/data/worldclim21.html}
#'
#' @format A \code{SpatRaster} with 8 layers:
#' \describe{
#'   \item{bio_1}{Annual Mean Temperature}
#'   \item{bio_5}{Max Temperature of Warmest Month}
#'   \item{bio_6}{Min Temperature of Coldest Month}
#'   \item{bio_7}{Temperature Annual Range (bio_5 - bio_6)}
#'   \item{bio_12}{Annual Precipitation}
#'   \item{bio_13}{Precipitation of Wettest Month}
#'   \item{bio_14}{Precipitation of Driest Month}
#'   \item{bio_15}{Precipitation Seasonality (Coefficient of Variation)}
#' }
#'
#' @name ma_bios
#'
#' @return No return value. Used with function \code{\link[terra]{rast}} to
#' load the GeoTIFF file from the package's \code{inst/extdata} folder.
#'
#' @examples
#' ma_bios <- terra::rast(system.file("extdata", "ma_bios.tif",
#'                                    package = "nicheR"))
#'
#' terra::plot(ma_bios[[1]])
NULL

