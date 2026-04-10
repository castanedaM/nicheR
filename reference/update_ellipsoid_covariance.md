# Update covariances in a nicheR ellipsoid and recompute metrics

Updates one or more off-diagonal covariance values in a
`nicheR_ellipsoid` object and recomputes all ellipsoid metrics
(centroid, semi-axes, volume, etc.) from the new covariance matrix. This
allows iterative niche shaping by adjusting the rotation and correlation
structure of the ellipsoid without rebuilding it from scratch.

## Usage

``` r
update_ellipsoid_covariance(object, covariance, tol = 1e-06, verbose = TRUE)
```

## Arguments

- object:

  A `nicheR_ellipsoid` object, typically created with
  [`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md).

- covariance:

  Either a single numeric value applied to all off-diagonal elements, or
  a named numeric vector where names identify the variable pair in the
  format `"var1-var2"` (e.g., `c("bio_1-bio_12" = 0.3)`).

- tol:

  Small positive number used as tolerance when computing safe covariance
  limits for positive definiteness. Default is `1e-6`.

- verbose:

  Logical. If `TRUE` (default), prints progress messages.

## Value

A `nicheR_ellipsoid` object with the updated covariance matrix and
recomputed ellipsoid metrics. An additional element
`cov_limits_remaining` is attached, giving the remaining safe covariance
limits for any variable pairs that still have zero covariance.

## Details

Covariance values control the orientation and correlation structure of
the ellipsoid in environmental space. Setting a positive covariance
between two variables tilts the ellipsoid so that high values of one
variable tend to co-occur with high values of the other. Negative
covariance tilts it in the opposite direction.

The updated covariance matrix must remain positive definite — if the
requested value would violate this, use
[`covariance_limits`](https://castanedam.github.io/nicheR/reference/covariance_limits.md)
to find the safe range before calling this function.

## See also

[`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md)
to create the initial ellipsoid,
[`covariance_limits`](https://castanedam.github.io/nicheR/reference/covariance_limits.md)
to find valid covariance ranges before updating.
