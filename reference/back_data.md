# Background environmental data for examples

A dataset containing geographic coordinates and two bioclimatic
variables used as a background point cloud for generating and testing
niche ellipsoids.

## Usage

``` r
back_data
```

## Format

A data frame with 12,396 rows and 4 variables:

- x:

  Longitude in decimal degrees (WGS84).

- y:

  Latitude in decimal degrees (WGS84).

- bio_1:

  Annual Mean Temperature (°C).

- bio_5:

  Max Temperature of Warmest Month (°C).

- bio_6:

  Min Temperature of Coldest Month (°C).

- bio_7:

  Temperature Annual Range (bio_5 - bio_6) (°C).

- bio_12:

  Annual Precipitation (mm).

- bio_13:

  Precipitation of Wettest Month (mm).

- bio_14:

  Precipitation of Driest Month (mm).

- bio_15:

  Precipitation Seasonality (Coefficient of Variation).

## Source

<https://www.worldclim.org>

## Details

This dataset represents an irregular point cloud typical of
environmental background data used in ecological niche modeling (ENM).
It is primarily used in the \`nicheR\` package to provide the
environmental space for functions like
[`conserved_ellipses`](https://castanedam.github.io/nicheR/reference/conserved_ellipses.md).

## Examples

``` r
data(back_data)
head(back_data)
#>           x        y    bio_1    bio_5   bio_6    bio_7 bio_12 bio_13 bio_14
#> 1 -99.91667 29.91667 18.16097 33.23550 0.86900 32.36650    680     84     26
#> 2 -99.75000 29.91667 18.06556 33.32575 0.65550 32.67025    703     87     28
#> 3 -99.58333 29.91667 17.95946 33.33925 0.44600 32.89325    725     92     31
#> 4 -99.41667 29.91667 18.01018 33.34200 0.62200 32.72000    734     95     33
#> 5 -99.25000 29.91667 18.14458 33.40400 0.99125 32.41275    748     97     34
#> 6 -99.08333 29.91667 18.36623 33.76550 1.02025 32.74525    771    101     36
#>     bio_15
#> 1 39.75968
#> 2 38.44158
#> 3 37.43598
#> 4 36.24147
#> 5 34.95365
#> 6 33.73626
```
