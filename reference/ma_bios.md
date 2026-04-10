# Bioclimatic variables for part of the Americas

A `SpatRaster` containing 8 bioclimatic variables representing
present-day climatic conditions for an area that covers parts of South
and North America. Variables were obtained at a 10 arc-minute
resolution. Sourced from WorldClim 2.1:
<https://worldclim.org/data/worldclim21.html>

## Format

A `SpatRaster` with 8 layers:

- bio_1:

  Annual Mean Temperature

- bio_5:

  Max Temperature of Warmest Month

- bio_6:

  Min Temperature of Coldest Month

- bio_7:

  Temperature Annual Range (bio_5 - bio_6)

- bio_12:

  Annual Precipitation

- bio_13:

  Precipitation of Wettest Month

- bio_14:

  Precipitation of Driest Month

- bio_15:

  Precipitation Seasonality (Coefficient of Variation)

## Value

No return value. Used with function
[`rast`](https://rspatial.github.io/terra/reference/rast.html) to load
the GeoTIFF file from the package's `inst/extdata` folder.

## Examples

``` r
ma_bios <- terra::rast(system.file("extdata", "ma_bios.tif",
                                   package = "nicheR"))

terra::plot(ma_bios[[1]])
```
