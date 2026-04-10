# Generate 2D ellipsoid boundary points for plotting

Computes ordered boundary points for a two-dimensional slice of a
`nicheR_ellipsoid`, suitable for plotting with
[`lines()`](https://rspatial.github.io/terra/reference/lines.html) or
`plot(type = "l")`. The boundary is derived from the selected covariance
submatrix and the stored chi-square cutoff.

## Usage

``` r
ellipsoid_boundary_2d(object, n_segments = 50, dim = c(1, 2))
```

## Arguments

- object:

  A `nicheR_ellipsoid` object.

- n_segments:

  Integer. Number of boundary points to generate (must be \>= 4).

- dim:

  Integer vector of length 2 indicating which dimensions (indices of the
  original variables) to use for the 2D slice.

## Value

A `data.frame` with `n_segments` ordered boundary points in the selected
dimensions.
