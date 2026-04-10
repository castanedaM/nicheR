# Prepare sampling bias surfaces

Standardizes and combines one or more bias layers into a composite bias
surface for use in biased occurrence sampling. Each layer is min-max
normalized to `[0, 1]` and assigned a directional effect (`"direct"` or
`"inverse"`) before being multiplied together into a single composite
surface.

## Usage

``` r
prepare_bias(
  bias_surface,
  effect_direction = c("direct", "inverse"),
  template_layer = NULL,
  include_composite = TRUE,
  include_processed_layers = FALSE,
  mask_na = FALSE,
  verbose = TRUE
)
```

## Arguments

- bias_surface:

  A `SpatRaster` (single or multi-layer) or a list of `SpatRaster`
  objects representing the raw bias layers.

- effect_direction:

  Character vector. Direction of effect for each bias layer. Each
  element must be `"direct"` (higher values increase sampling
  probability) or `"inverse"` (higher values decrease sampling
  probability). Length 1 recycles to all layers. Must otherwise match
  the number of bias layers.

- template_layer:

  Optional `SpatRaster`. Reference layer used to align all bias layers
  (resolution, extent, CRS). If `NULL` (default), the finest-resolution
  bias layer is used as the template.

- include_composite:

  Logical. If `TRUE` (default), includes the composite bias surface in
  the output.

- include_processed_layers:

  Logical. If `TRUE`, includes the standardized individual layers (after
  directional transformation) in the output. Default is `FALSE`.

- mask_na:

  Logical. Controls how `NA` values are handled when combining layers.
  If `TRUE`, uses the intersection of layer extents — any pixel with an
  `NA` in any layer becomes `NA` in the composite. If `FALSE` (default),
  uses the union of extents and ignores `NA`s where other layers have
  valid values.

- verbose:

  Logical. If `TRUE` (default), prints progress messages.

## Value

A named list with some or all of the following elements depending on the
`include_composite` and `include_processed_layers` arguments:

- `composite_surface`: A `SpatRaster` with the combined bias surface,
  named `"standarized_composite_bias_surface"`.

- `processed_layers`: A multi-layer `SpatRaster` with the standardized
  and direction-transformed individual layers.

- `combination_formula`: A character string showing the formula used to
  combine layers (e.g., `"sp_richness * (1-nighttime)"`).

## Details

Each bias layer is processed as follows:

1.  Resampled and cropped to the template grid if needed.

2.  Min-max standardized to `[0, 1]`.

3.  Transformed by direction: `"direct"` layers are kept as-is;
    `"inverse"` layers are replaced by \\1 - x\\.

The composite surface is the product of all transformed layers. If only
one layer is provided, it is returned as-is.

If both `include_composite` and `include_processed_layers` are `FALSE`,
the function defaults to `include_composite = TRUE` with a warning.

## See also

[`apply_bias`](https://castanedam.github.io/nicheR/reference/apply_bias.md)
to apply the prepared bias surface to a suitability prediction,
[`sample_biased_data`](https://castanedam.github.io/nicheR/reference/sample_biased_data.md)
to sample occurrences from the resulting bias-weighted surface.
