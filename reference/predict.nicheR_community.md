# Predict method for a nicheR Community

Iterates through all ellipses in a `nicheR_community` and calculates
predictions (suitability or Mahalanobis distance) using
[`predict.nicheR_ellipsoid`](https://castanedam.github.io/nicheR/reference/predict.nicheR_ellipsoid.md).

## Usage

``` r
# S3 method for class 'nicheR_community'
predict(object, newdata, prediction = "Mahalanobis", verbose = TRUE)
```

## Arguments

- object:

  A `nicheR_community` object.

- newdata:

  A `SpatRaster`, `data.frame`, or `matrix` containing at least the
  environmental variables used to create the reference ellipsoid in the
  community.

- prediction:

  Character. The type of prediction to return. One of: `"Mahalanobis"`
  (default), `"suitability"`, `"Mahalanobis_trunc"`, or
  `"suitability_trunc"`.

- verbose:

  Logical. If `TRUE`, prints progress messages. Default is `TRUE`.

## Value

If `newdata` is a `SpatRaster`, returns a `SpatRaster` where each layer
represents one ellipse. If `newdata` is a `data.frame`, returns a
`data.frame` with the original data plus one column per ellipse.
