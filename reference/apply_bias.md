# Apply sampling bias to suitability surfaces

Applies a prepared composite sampling bias surface to a suitability
raster by multiplication. The bias surface is aligned to the suitability
grid when needed and the result is cropped and masked to the suitability
domain. The output is a product of suitability and bias and is therefore
no longer interpretable as a probability.

## Usage

``` r
apply_bias(
  prepared_bias,
  prediction,
  prediction_layer = NULL,
  effect_direction = "direct",
  verbose = TRUE
)
```

## Arguments

- prepared_bias:

  A single-layer `SpatRaster` composite bias surface, or the list output
  from
  [`prepare_bias`](https://castanedam.github.io/nicheR/reference/prepare_bias.md)
  containing a `composite_surface` element.

- prediction:

  A `SpatRaster` containing one or more suitability layers with values
  in `[0, 1]`.

- prediction_layer:

  Character. Name of the layer to extract from `prediction` when it
  contains multiple layers. If `NULL` (default) and `prediction` has a
  single layer, that layer is used.

- effect_direction:

  Character. How the bias surface is applied to the suitability layer.
  `"direct"` (default) multiplies suitability by the bias directly —
  higher bias increases sampling probability. `"inverse"` multiplies by
  \\1 - \text{bias}\\ — higher bias decreases sampling probability.

- verbose:

  Logical. If `TRUE` (default), prints progress messages.

## Value

A named list of class `"nicheR_biased_surface"` containing:

- One `SpatRaster` per input suitability layer, named
  `"<layer>_biased"`. The raster layer name includes the applied
  direction (e.g., `"suitability_biased_direct"`).

- `combination_formula`: a character string describing the operation
  applied (e.g., `"suitability * bias"` or `"suitability * (1-bias)"`).

## Details

The function performs the following steps:

1.  Extracts the composite bias surface from `prepared_bias`.

2.  Verifies both bias and suitability values are within `[0, 1]`.

3.  Aligns the bias surface to the suitability grid if geometries
    differ, using
    [`terra::resample()`](https://rspatial.github.io/terra/reference/resample.html)
    with nearest-neighbor interpolation.

4.  Multiplies suitability by the (possibly inverted) bias surface.

5.  Crops and masks the output to the suitability domain.

## See also

[`prepare_bias`](https://castanedam.github.io/nicheR/reference/prepare_bias.md)
to build the composite bias surface,
[`sample_biased_data`](https://castanedam.github.io/nicheR/reference/sample_biased_data.md)
to sample occurrences from the output.
