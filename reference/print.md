# Print method for nicheR objects

Provides a concise summary of `nicheR` objects.

## Usage

``` r
# S3 method for class 'nicheR_ellipsoid'
print(x, digits = 3, ...)

# S3 method for class 'nicheR_community'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of the classes `"nicheR_ellipsoid"` or `"nicheR_community"`.

- digits:

  Integer. Number of decimal places used when printing numeric values.
  Default is 3.

- ...:

  Additional arguments.

## Value

The input object `x`, returned invisibly.

## Details

The function formats and rounds key quantities for readability but does
not modify the underlying object.

## See also

[`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md)
and `generate_community`
