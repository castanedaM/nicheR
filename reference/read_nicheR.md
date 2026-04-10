# Read a nicheR object from disk

A wrapper around [`readRDS`](https://rdrr.io/r/base/readRDS.html) to
load saved nicheR objects back into the R environment.

## Usage

``` r
read_nicheR(file)
```

## Arguments

- file:

  Character. The path to the file to be read.

## Value

The saved `nicheR_ellipsoid` or `nicheR_community` object.

## See also

[`save_nicheR`](https://castanedam.github.io/nicheR/reference/save_nicheR.md)
