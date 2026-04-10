# Predict suitability and Mahalanobis distance from a nicheR ellipsoid

Computes Mahalanobis distance and optional suitability values from a
`nicheR_ellipsoid` for either (1) environmental samples provided as a
`data.frame` or `matrix`, (2) a `SpatRaster` stack of predictors, or (3)
virtual samples drawn in environmental space when `newdata = NULL`.

## Usage

``` r
# S3 method for class 'nicheR_ellipsoid'
predict(
  object,
  newdata,
  adjust_truncation_level = NULL,
  include_suitability = TRUE,
  suitability_truncated = FALSE,
  include_mahalanobis = TRUE,
  mahalanobis_truncated = FALSE,
  keep_data = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  A `nicheR_ellipsoid` object produced by
  [`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md).

- newdata:

  Environmental predictors. One of:

  - A `SpatRaster` (or legacy `raster` classes, coerced automatically).

  - A `data.frame`, `tibble`, or `matrix` with columns named to match
    `object$var_names`.

- adjust_truncation_level:

  Optional numeric confidence level in `(0, 1)` to override `object$cl`
  when computing truncated outputs. Default is `NULL` (uses the level
  stored in `object`).

- include_suitability:

  Logical. If `TRUE` (default), returns suitability values (\\\exp(-0.5
  D^2)\\).

- suitability_truncated:

  Logical. If `TRUE`, returns a truncated suitability layer where values
  outside the chi-square contour are set to `0`. Default is `FALSE`.

- include_mahalanobis:

  Logical. If `TRUE` (default), returns Mahalanobis distance (\\D^2\\).

- mahalanobis_truncated:

  Logical. If `TRUE`, returns a truncated Mahalanobis layer where values
  outside the chi-square contour are set to `NA`. Default is `FALSE`.

- keep_data:

  Logical or `NULL`. If `TRUE`, includes the original predictors in the
  output. Default is `NULL`: `FALSE` for `SpatRaster` input, `TRUE` for
  tabular input.

- verbose:

  Logical. If `TRUE` (default), prints progress messages.

## Value

If `newdata` is a `SpatRaster`, returns a `SpatRaster` with the
requested prediction layers (and optionally the original predictor
layers if `keep_data = TRUE`).

If `newdata` is tabular, returns a `data.frame` of class
`"nicheR_prediction"` with the requested prediction columns (and
optionally the original predictor columns if `keep_data = TRUE`).

## Details

Suitability is computed as \\\exp(-0.5 D^2)\\, where \\D^2\\ is the
squared Mahalanobis distance from the niche centroid. Truncated outputs
use a chi-square cutoff based on the ellipsoid confidence level (`cl`).

For tabular inputs, coordinate columns (e.g., `x`, `y`, `lon`, `lat`)
are detected and retained when `keep_data = TRUE`. Extra non-predictor
columns are ignored.
