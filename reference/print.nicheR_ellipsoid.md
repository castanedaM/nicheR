# Print a nicheR Ellipsoid Object

Provides a concise summary of a `nicheR_ellipsoid` object created by
[`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md).
The printed output includes dimensionality, chi-square cutoff, centroid,
covariance matrix, principal semi-axis lengths, axis endpoints, and
ellipsoid volume.

## Usage

``` r
# S3 method for class 'nicheR_ellipsoid'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `nicheR_ellipsoid` object.

- digits:

  Integer. Number of decimal places used when printing numeric values.
  Default is 3.

- ...:

  Additional arguments.

## Value

The input object `x`, returned invisibly.

## Details

This is an S3 method for objects of class `"nicheR_ellipsoid"`. The
function formats and rounds key quantities for readability but does not
modify the underlying object.

## See also

[`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md)
