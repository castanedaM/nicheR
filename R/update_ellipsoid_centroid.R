#' Move the centroid of an ellipsoidal niche
#'
#' Translates a `nicheR_ellipsoid` through environmental space, leaving its
#' shape, size, and orientation unchanged.
#'
#' The covariance matrix is not touched, so the semi-axis lengths, the niche
#' volume, and the chi-square cutoff are identical before and after. Only the
#' position changes. The ranges are shifted by the same amount as the centroid,
#' so they continue to describe the extent of the ellipsoid rather than the
#' extent it was originally built from.
#'
#' Use this to ask what a species with the same tolerance breadth would look
#' like centred on different conditions, for example a congener occupying a
#' cooler part of the same environmental space. To change the breadth or the
#' orientation instead, rebuild with [build_ellipsoid()] or rotate with
#' [update_ellipsoid_covariance()].
#'
#' Nothing prevents moving a centroid outside the range of the available
#' environmental data. The resulting ellipsoid is valid, but predicting it onto
#' that data may return little or no suitable area.
#'
#' @usage update_ellipsoid_centroid(ell, new_centroid, verbose = FALSE)
#'
#' @param ell A `nicheR_ellipsoid` object, as returned by [build_ellipsoid()].
#' @param new_centroid A named numeric vector giving the new centroid position.
#' Names must include every entry of `ell$var_names`; extra names are ignored
#' and the order does not matter.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default is
#' `FALSE`.
#'
#' @returns A `nicheR_ellipsoid` object with `centroid` and `ranges` shifted by
#' the difference between `new_centroid` and the original centroid. All other
#' components are unchanged.
#'
#' @seealso [build_ellipsoid()], [update_ellipsoid_covariance()]
#'
#' @examples
#' range_df <- data.frame(bio_1 = c(22, 28),
#'                        bio_12 = c(1000, 3500))
#' ell <- build_ellipsoid(range = range_df)
#' ell$centroid
#' ell$ranges
#' rownames(ell$ranges) <- c("min", "max")
#'
#' # Same niche breadth, shifted to cooler and wetter conditions
#' ell_shifted <- update_ellipsoid_centroid(
#'   ell,
#'   new_centroid = c(bio_1 = 19, bio_12 = 2800)
#' )
#'
#' ell_shifted$centroid
#' ell_shifted$ranges
#'
#' # The shape is untouched
#' identical(ell$covariance_matrix, ell_shifted$covariance_matrix)
#'
#' # Both niches in the same environmental space
#' ma_bios <- terra::rast(
#'   system.file("extdata/ma_bios.tif", package = "nicheR"))
#' back_df <- as.data.frame(ma_bios[[c("bio_1", "bio_12")]])
#'
#' plot_ellipsoid(object = ell, background = back_df)
#' add_ellipsoid(object = ell_shifted, col = "red")
#'
#' @export
update_ellipsoid_centroid <- function(ell, new_centroid, verbose = FALSE){

  if(!inherits(ell, "nicheR_ellipsoid")){
    stop("ell must be a nicheR_ellipsoid object.")
  }

  vars <- ell$var_names

  if(is.null(names(new_centroid)) || !all(vars %in% names(new_centroid))){
    stop("new_centroid must be a named numeric vector with names matching ",
         "ell$var_names: ", paste(vars, collapse = ", "))
  }

  new_centroid <- new_centroid[vars]

  if(any(!is.finite(new_centroid))){
    stop("new_centroid contains non-finite values.")
  }

  if(verbose) message("Starting: moving ellipsoid centroid...")

  old_centroid <- ell$centroid
  delta <- new_centroid - old_centroid

  # Shift ranges by the same delta so the ellipsoid translates rigidly. Indexed
  # by column, since ranges built outside the app may not carry row names, and
  # sorted so the row names describe the values rather than assert an order.
  new_ranges <- ell$ranges

  for(v in vars){
    new_ranges[, v] <- sort(new_ranges[, v] + delta[[v]])
  }

  rownames(new_ranges) <- c("min", "max")

  ell$axes_coordinates <- lapply(ell$axes_coordinates,
                                    function(a){
                                      sweep(a, 2L, delta, "+")
                                    })

  ell$centroid <- new_centroid
  ell$ranges <- new_ranges

  if(verbose) message("Done: centroid moved.")

  ell
}
