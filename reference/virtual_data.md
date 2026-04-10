# Generate data based on a ellipsoidal niche

Simulates `n` random points from a multivariate normal distribution
defined by the centroid and covariance matrix of a `nicheR_ellipsoid`
object.

## Usage

``` r
virtual_data(object, n = 100, truncate = FALSE, effect = "direct", seed = 1)
```

## Arguments

- object:

  A `nicheR_ellipsoid` object containing at least `centroid` and
  `cov_matrix`.

- n:

  Integer. The number of virtual points to generate. Default = 100.

- truncate:

  Logical. If `TRUE` (default), points are constrained within the
  confidence limit (`cl`) defined in the object.

- effect:

  Character. The distribution pattern of points. `"direct"` (default)
  creates a concentration near the centroid. `"inverse"` creates higher
  density towards the edges. `"uniform"` distributes points evenly
  throughout the ellipsoid volume. Note: `"inverse"` and `"uniform"`
  require `truncate = TRUE`.

- seed:

  Integer. Random seed for reproducibility. Default = 1. Set to `NULL`
  for no seeding.

## Value

A matrix with `n` rows and columns corresponding to the environmental
variables (dimensions) of the input `object`.

## Details

When `truncate = FALSE`, the function generates points from a standard
multivariate normal distribution defined by the ellipsoid's centroid and
covariance matrix, without any constraints on their location. The
function uses eigen-decomposition to transform standard normal variables
into the coordinate system defined by the ellipsoid's covariance
structure.

When `truncate = TRUE`, the function generates candidate points
uniformly distributed within a bounding box (hyper-cube) defined by the
ellipsoid's `axes_coordinates`. Points falling outside the ellipsoid
(where Mahalanobis distance \\Md \>\\ `chi2_cutoff`) are removed.

From this filtered pool, `n` points are selected using weighted random
sampling without replacement. The weights are determined by the `effect`
argument:

- `"direct"`: Weights are proportional to the multivariate normal
  density (\\\exp(-0.5 \times Md)\\), clustering points near the
  centroid.

- `"inverse"`: Weights are proportional to the complement of the normal
  density (\\1 - \exp(-0.5 \times Md)\\), pushing points toward the
  edges.

- `"uniform"`: All points within the ellipsoid have equal weight,
  resulting in a uniform spatial distribution.
