# Sample occurrence data from a bias-weighted prediction surface

Samples `n_occ` virtual occurrence points using the bias-weighted
prediction values directly as sampling probabilities. Unlike
[`sample_data()`](https://castanedam.github.io/nicheR/reference/sample_data.md),
there is no sampling strategy argument — the prediction layer values
themselves define where points are drawn from, making this function
suited for simulating realistically biased occurrence records.

## Usage

``` r
sample_biased_data(
  n_occ,
  prediction,
  prediction_layer = NULL,
  sampling_mask = NULL,
  seed = 1,
  verbose = TRUE,
  strict = NULL
)
```

## Arguments

- n_occ:

  Integer. Number of occurrence points to sample.

- prediction:

  A `SpatRaster` or data frame containing the bias-weighted prediction
  surface to sample from.

- prediction_layer:

  Character. Name of the layer or column to use as sampling weights.
  Required when `prediction` contains multiple layers or columns.

- sampling_mask:

  A `SpatRaster` or `SpatVector` used to restrict sampling to a
  geographic area. Only supported when `prediction` is a `SpatRaster`.

- seed:

  Integer. Random seed for reproducibility. Default is `1`.

- verbose:

  Logical. If `TRUE` (default), prints progress messages.

- strict:

  Logical or `NULL`. If `TRUE`, removes `NA` and zero-valued cells
  before sampling. If `NULL` (default), auto-detected from the layer
  name and the proportion of zeros and `NA`s in the prediction values.

## Value

A data frame of sampled occurrence points with the same columns as the
input `prediction` (minus the internal `pred` column). If `prediction`
is a `SpatRaster`, the output includes `x` and `y` coordinate columns.

## Details

Prediction values are used directly as sampling weights, so they must be
non-negative. Higher values correspond to higher sampling probability,
reflecting areas of greater bias (e.g., higher detectability or observer
effort). This is in contrast to
[`sample_data()`](https://castanedam.github.io/nicheR/reference/sample_data.md),
which transforms prediction values according to a `sampling` and
`method` argument.

Auto-detection of `strict` follows the same logic as
[`sample_data()`](https://castanedam.github.io/nicheR/reference/sample_data.md):
it is set to `TRUE` if the layer name contains `"trunc"` or if the
proportion of zeros or `NA`s exceeds 25%.

## See also

[`sample_data`](https://castanedam.github.io/nicheR/reference/sample_data.md)
for unbiased sampling with explicit strategy and method control,
[`apply_bias`](https://castanedam.github.io/nicheR/reference/apply_bias.md)
for generating the bias-weighted prediction surface used as input here.
