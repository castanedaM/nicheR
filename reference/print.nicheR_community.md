# Print a nicheR Community Object

Provides a structured summary of a `nicheR_community` object. Includes
the generation metadata, a summary of the reference ellipsoid, and
descriptive statistics (mean and standard deviation) for the centroids
and volumes of the generated community.

## Usage

``` r
# S3 method for class 'nicheR_community'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `nicheR_community` object.

- digits:

  Integer. Number of decimal places used when printing numeric values.
  Default is 3.

- ...:

  Additional arguments passed to the `nicheR_ellipsoid` print method.
