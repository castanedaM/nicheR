# nicheR 0.2.0

## Breaking changes

* `virtual_data()`: the `effect` argument is now `sampling`, and its values
  `"direct"` and `"inverse"` are now `"centroid"` and `"edge"`, the words
  `sample_data()` already uses. `"uniform"` is unchanged. Code that passes
  `effect` needs to be updated.

* `sample_data()`: the `sampling` value `"random"` is now `"uniform"`, so
  `sample_data()` and `virtual_data()` take the same three values. Code that
  passes `sampling = "random"` needs to be updated.

* `apply_bias()` now weights the prediction the same way `sample_data()`
  does. The `effect_direction` argument is removed and two arguments replace
  it: `sampling` (`"centroid"`, `"edge"`, or `"uniform"`) and `method`
  (`"suitability"` or `"mahalanobis"`). The weight multiplied by the bias
  surface is the suitability or 1 / Mahalanobis distance for `"centroid"`,
  1 - suitability or the Mahalanobis distance for `"edge"`, and 1 for
  `"uniform"`. The old `"inverse"` option, which reflected the layer as
  (max + min) - x, is removed. Code that passes `effect_direction` to
  `apply_bias()` needs to be updated. `prepare_bias()` keeps its own
  `effect_direction` argument.

    - A Mahalanobis distance layer now needs `method = "mahalanobis"`. With
      the old default it was multiplied by the bias as is, which favors the
      niche edge. The new default, `sampling = "centroid"`, favors the
      centroid.
    - Outputs are now named `<layer>_<sampling>_biased`, for example
      `suitability_centroid_biased`, for both the list element and the raster
      layer. They were `<layer>_biased` and `<layer>_biased_<direction>`.
    - Cells where the prediction is zero are set to `NA` in the output, the
      same as `strict = TRUE` in `sample_data()`, so data points are never
      drawn outside a truncated niche.

* `sample_data()` and `apply_bias()` now stop when `method` disagrees with
  the layer name, for example `method = "mahalanobis"` on a layer called
  `suitability_trunc`. They warn when `method = "mahalanobis"` is used on a
  layer whose name gives no hint and whose values all fall in [0, 1].

## New features

* `nicheR_gui()` launches a Shiny interface to the package workflow. The
  application covers building ellipsoids from environmental ranges, predicting
  suitability, preparing and applying sampling bias, and generating virtual
  data points. Sessions can be saved and reloaded, and several ellipsoids can
  be carried through the workflow at once for comparison.

* `update_ellipsoid_centroid()` translates an ellipsoid through environmental
  space, leaving its shape, size, and orientation unchanged. Use it to ask
  what a species with the same tolerance breadth would look like centred on
  different conditions.

## Documentation

* New vignette, `shiny_app`, documenting the interface.

# nicheR 0.1.0

* Initial CRAN submission.
