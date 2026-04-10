# Compute safe axis limits covering data and ellipsoid boundary

Returns x and y ranges that span both a matrix of data points and a data
frame of ellipsoid boundary points, so the ellipsoid boundary is never
clipped by the data extent when passed to
[`plot()`](https://rspatial.github.io/terra/reference/plot.html).

## Usage

``` r
safe_lims(pts_xy, ell_xy)
```

## Arguments

- pts_xy:

  A data frame or matrix with at least two columns. The first column is
  used for x, the second for y.

- ell_xy:

  A data frame or matrix of ellipsoid boundary points with the same
  column structure as `pts_xy`.

## Value

A named list with elements `xlim` and `ylim`, each a numeric vector of
length 2.
