# Add an ellipsoid boundary to an existing E-space plot

Draws the 2D boundary of a `nicheR_ellipsoid` object onto an existing
environmental space plot created with
[`plot_ellipsoid()`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid.md).
The boundary is computed as a cross-section of the ellipsoid at the
chosen pair of dimensions.

## Usage

``` r
add_ellipsoid(
  object,
  dim = c(1, 2),
  lty = 1,
  lwd = 1,
  col_ell = "#000000",
  alpha_ell = 1,
  cex_ell = 1,
  ...
)
```

## Arguments

- object:

  A `nicheR_ellipsoid` object.

- dim:

  Integer vector of length 2. Indices of the two dimensions to plot.
  Default is `c(1, 2)`.

- lty:

  Integer. Line type. Default is `1` (solid).

- lwd:

  Numeric. Line width. Default is `1`.

- col_ell:

  Character. Color of the ellipsoid boundary line. Default is
  `"#000000"` (black).

- alpha_ell:

  Numeric in `[0, 1]`. Transparency of the ellipsoid boundary line.
  Default is `1` (fully opaque).

- cex_ell:

  Numeric. Size scaling for the ellipsoid boundary. Default is `1`.

- ...:

  Additional arguments passed to
  [`lines`](https://rdrr.io/r/graphics/lines.html).

## Value

Called for its side effect of adding lines to the current plot. Returns
`NULL` invisibly.

## See also

[`plot_ellipsoid`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid.md),
[`add_data`](https://castanedam.github.io/nicheR/reference/add_data.md)
