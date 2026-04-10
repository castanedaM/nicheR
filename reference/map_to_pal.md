# Map numeric values to palette indices

Rescales a numeric vector linearly from its observed range onto integer
indices in `[1, pal_len]` for use in palette lookups. When all values
are identical, returns the middle index for every element.

## Usage

``` r
map_to_pal(vals, pal_len)
```

## Arguments

- vals:

  Numeric vector of values to map.

- pal_len:

  Integer. Length of the target palette.

## Value

Integer vector of the same length as `vals`, with values in
`[1, pal_len]`.
