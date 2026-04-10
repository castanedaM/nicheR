# Sample occurrence data from a prediction surface

Samples `n_occ` virtual occurrence points from a suitability or
Mahalanobis distance prediction surface. Supports centroid, edge, and
random sampling strategies, and accepts both raster (`SpatRaster`) and
data frame inputs.

## Usage

``` r
sample_data(
  n_occ,
  prediction,
  prediction_layer = NULL,
  sampling = "centroid",
  method = "suitability",
  sampling_mask = NULL,
  seed = 1,
  strict = NULL,
  verbose = TRUE
)
```

## Arguments

- n_occ:

  Integer. Number of occurrence points to sample.

- prediction:

  A `SpatRaster` or data frame containing the prediction surface to
  sample from.

- prediction_layer:

  Character. Name of the layer or column to use as the prediction
  values. Required when `prediction` contains multiple layers or
  columns.

- sampling:

  Character. Sampling strategy. One of `"centroid"` (default), `"edge"`,
  or `"random"`. Controls where within the niche points are
  preferentially drawn from.

- method:

  Character. Weighting method. One of `"suitability"` (default) or
  `"mahalanobis"`. Must match the type of values in `prediction_layer`:
  suitability values must be in `[0, 1]`, Mahalanobis values must be
  non-negative.

- sampling_mask:

  A `SpatRaster` or `SpatVector` used to restrict sampling to a
  geographic area. Only supported when `prediction` is a `SpatRaster`.

- seed:

  Integer. Random seed for reproducibility. Default is `1`.

- strict:

  Logical or `NULL`. If `TRUE`, removes `NA` and zero-valued cells
  before sampling (recommended with truncated prediction layers). If
  `NULL` (default), auto-detected from the layer name and the proportion
  of zeros and `NA`s in the prediction values.

- verbose:

  Logical. If `TRUE` (default), prints progress messages.

## Value

A data frame of sampled occurrence points with the same columns as the
input `prediction` (minus the internal `pred` column). If `prediction`
is a `SpatRaster`, the output includes `x` and `y` coordinate columns.

## Details

The `sampling` and `method` arguments interact to define the probability
weights used when drawing points:

- `sampling = "centroid"`, `method = "suitability"`: weights
  proportional to suitability — higher near the niche center.

- `sampling = "edge"`, `method = "suitability"`: weights proportional to
  \\1 - \text{suitability}\\ — higher near the niche boundary.

- `sampling = "centroid"`, `method = "mahalanobis"`: weights inversely
  proportional to Mahalanobis distance — higher near the centroid.

- `sampling = "edge"`, `method = "mahalanobis"`: weights proportional to
  Mahalanobis distance — higher near the boundary.

- `sampling = "random"`: equal weights regardless of method.

When `strict = NULL`, the function auto-detects truncation by checking
whether the layer name contains `"trunc"` or whether the proportion of
zeros or `NA`s exceeds 25%.
