# Plot all pairwise 2D ellipsoid projections

Plots all pairwise two-dimensional slices of a `nicheR_ellipsoid` in a
multi-panel layout using
[`plot_ellipsoid()`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid.md).
When `background` or `prediction` is supplied, axis limits are computed
once from the global range of all variables and shared across every
panel, so projections are directly comparable without distortion from
per-panel rescaling.

## Usage

``` r
plot_ellipsoid_pairs(object, background = NULL, prediction = NULL, ...)
```

## Arguments

- object:

  A `nicheR_ellipsoid` object.

- background:

  Optional data frame or matrix of background points passed to each
  [`plot_ellipsoid()`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid.md)
  call. When provided, global axis limits are computed from the range of
  all variables in `background` combined with all pairwise ellipsoid
  boundaries.

- prediction:

  Optional data frame or matrix of prediction values passed to each
  [`plot_ellipsoid()`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid.md)
  call. Used when `background` is `NULL`. Global limits are computed
  from the range of all variables in `prediction`.

- ...:

  Additional graphical arguments passed to
  [`plot_ellipsoid()`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid.md).

## Value

Invisibly returns `NULL`.

## Details

Global limits are computed per variable across the full data and all
ellipsoid boundary projections, then passed to each panel via the
`fixed_lims` argument of
[`plot_ellipsoid()`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid.md).
This prevents individual panels from rescaling to their own data extent,
which would make niche widths appear identical across dimensions even
when they differ. If neither `background` nor `prediction` is provided,
each panel shows only the ellipsoid boundary and limits come from that
boundary alone, which is the intended behavior for a boundary-only view.
