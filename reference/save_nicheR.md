# Save a nicheR object to disk

A wrapper around [`saveRDS`](https://rdrr.io/r/base/readRDS.html) to
save `nicheR_ellipsoid` or `nicheR_community` objects to a file.
Includes a safety check for overwriting existing files and ensures the
file extension is .rds.

## Usage

``` r
save_nicheR(object, file, overwrite = FALSE, ...)
```

## Arguments

- object:

  A `nicheR_ellipsoid` or `nicheR_community` object.

- file:

  Character. The connection or name of the file where the object will be
  saved (usually ending in ".rds").

- overwrite:

  Logical. If `TRUE`, an existing file at the specified path will be
  replaced. Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`saveRDS`](https://rdrr.io/r/base/readRDS.html).

## Value

No return value. Saves the object to the specified path.
