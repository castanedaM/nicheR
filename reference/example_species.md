# Example species occurrence data

Datasets containing geographic occurrence records for four example
species. These datasets are used to demonstrate and test functions
related to ecological niche modeling, such as building niche ellipsoids
and simulating virtual communities in the `nicheR` package.

## Usage

``` r
example_sp_1

example_sp_2

example_sp_3

example_sp_4
```

## Format

Data frames containing occurrence points. While the exact number of rows
varies per dataset, they typically share the following variables:

- longitude:

  Numeric. Longitude in decimal degrees (WGS84).

- latitude:

  Numeric. Latitude in decimal degrees (WGS84).

- species:

  Character (Optional). The identifier for the species.

An object of class `nicheR_ellipsoid` of length 13.

An object of class `nicheR_ellipsoid` of length 13.

An object of class `nicheR_ellipsoid` of length 14.

## Details

`example_sp_1`, `example_sp_2`, `example_sp_3`, and `example_sp_4`
represent typical presence-only occurrence records. They can be
extracted against background environmental data (like
[`back_data`](https://castanedam.github.io/nicheR/reference/back_data.md))
to fit `nicheR_ellipsoid` objects or test overlap functions.

## Examples

``` r
data(example_sp_1)
head(example_sp_1)
#> $dimensions
#> [1] 2
#> 
#> $var_names
#> [1] "bio_1"  "bio_12"
#> 
#> $centroid
#>  bio_1 bio_12 
#>     26   2375 
#> 
#> $cov_matrix
#>        bio_1   bio_12
#> bio_1      4    750.0
#> bio_12   750 293402.8
#> 
#> $Sigma_inv
#>              [,1]          [,2]
#> [1,]  0.480113636 -1.227273e-03
#> [2,] -0.001227273  6.545455e-06
#> 
#> $chol_Sigma
#>        bio_1  bio_12
#> bio_1      2 375.000
#> bio_12     0 390.868
#> 

data(example_sp_2)
head(example_sp_2)
#> $dimensions
#> [1] 2
#> 
#> $var_names
#> [1] "bio_1"  "bio_12"
#> 
#> $centroid
#>  bio_1 bio_12 
#>     18   2950 
#> 
#> $cov_matrix
#>        bio_1   bio_12
#> bio_1      4    500.0
#> bio_12   500 266944.4
#> 
#> $Sigma_inv
#>              [,1]          [,2]
#> [1,]  0.326426630 -6.114130e-04
#> [2,] -0.000611413  4.891304e-06
#> 
#> $chol_Sigma
#>        bio_1   bio_12
#> bio_1      2 250.0000
#> bio_12     0 452.1553
#> 
```
