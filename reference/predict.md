# Predict suitability and Mahalanobis distance from a nicheR ellipsoid

Computes Mahalanobis distance and suitability values deriving from a
`nicheR_ellipsoid` or `nicheR_community` object, for `newdata` provided
as a `data.frame`, `matrix`, or `SpatRaster`.

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
  verbose = TRUE,
  ...
)

# S3 method for class 'nicheR_community'
predict(object, newdata, prediction = "Mahalanobis", verbose = TRUE, ...)
```

## Arguments

- object:

  An object of the classes `"nicheR_ellipsoid"` or `"nicheR_community"`.

- newdata:

  Environmental predictors. One of:

  - A `SpatRaster` (or legacy `raster` classes, coerced automatically).

  - A `data.frame` or `matrix` with columns named to match
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

- ...:

  Additional arguments, not currently used.

- prediction:

  Character. The type of prediction to return. One of: `"Mahalanobis"`
  (default), `"suitability"`, `"Mahalanobis_trunc"`, or
  `"suitability_trunc"`.

## Value

For `nicheR_ellipsoid` objects, if `newdata` is a `SpatRaster`, returns
a `SpatRaster` with the requested prediction as layers (and optionally
the original predictor layers if `keep_data = TRUE`). If `newdata` is
tabular, returns a `data.frame` with the requested predictions as
columns (by default returns the original predictors as columns).

For `nicheR_community` objects, if `newdata` is a `SpatRaster`, returns
a `SpatRaster` where each layer represents predictions for each ellipse.
If `newdata` is a `data.frame`, returns a `data.frame` with the original
data plus one prediction column per ellipse.

## Details

Suitability is computed as \\\exp(-0.5 D^2)\\, where \\D^2\\ is the
squared Mahalanobis distance from the niche centroid. Truncated outputs
use a chi-square cutoff based on the ellipsoid confidence level (`cl`).

For tabular inputs, coordinate columns (e.g., `x`, `y`, `lon`, `lat`)
are detected and retained when `keep_data = TRUE`. Extra non-predictor
columns are ignored.
