#' Background environmental data for examples
#'
#' A dataset containing geographic coordinates and two bioclimatic variables
#' used as a background point cloud for generating and testing niche ellipsoids.
#'
#' @format A data frame with 12,396 rows and 4 variables:
#' \describe{
#'   \item{x}{Longitude in decimal degrees (WGS84).}
#'   \item{y}{Latitude in decimal degrees (WGS84).}
#'   \item{bio1}{Annual Mean Temperature (units: degrees Celsius).}
#'   \item{bio12}{Annual Precipitation (units: mm).}
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



#' Reference niche ellipsoid for examples
#'
#' A pre-calculated \code{nicheR_ellipsoid} object representing a hypothetical
#' species niche based on Annual Mean Temperature (bio1) and Annual
#' Precipitation (bio12).
#'
#' @format An object of class \code{nicheR_ellipsoid} (which is a \code{list})
#' with 13 elements:
#' \describe{
#'   \item{dimensions}{Integer. Number of dimensions (2).}
#'   \item{var_names}{Character vector. Names of variables (\code{"bio1"}, \code{"bio12"}).}
#'   \item{centroid}{Named numeric vector. The center of the niche (\eqn{\mu}).}
#'   \item{cov_matrix}{Matrix. The \eqn{2 \times 2} covariance matrix (\eqn{\Sigma}).}
#'   \item{Sigma_inv}{Matrix. The precision matrix (inverse covariance).}
#'   \item{chol_Sigma}{Matrix. Cholesky decomposition of the covariance.}
#'   \item{eigen}{List. Eigenvectors and eigenvalues of the covariance.}
#'   \item{cl}{Numeric. Confidence level used (e.g., 0.95).}
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
#' \code{\link{build_ellipsoid}} with a centroid at (17.5, 1000) and
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