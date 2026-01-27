#' Sample occurrence points from a suitable environment pool
#'
#' Sample occurrence points from a precomputed pool of suitable environments.
#' Sampling can be uniform, biased toward the niche center, or biased toward
#' the niche edge. Weights can be computed from Mahalanobis distance
#' (\code{dist_sq}) or from a multivariate normal (MVN) density.
#'
#' @export
get_sample_occ <- function(n_occ,
                           suitable_env,
                           method = c("random", "center", "edge"),
                           sampling = c("mahalanobis", "mvn"),
                           niche = NULL,
                           bias_surface = NULL,
                           seed = NULL,
                           verbose = TRUE) {
}

