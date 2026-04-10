# Update covariance values and calculate remaining safe limits

Updates a variance-covariance matrix with specific values and identifies
the safe limits for all remaining zero-covariance combinations.

## Usage

``` r
update_covariance(varcov_matrix, covariance = 0, tol = 1e-6)
```

## Arguments

- varcov_matrix:

  A square numerical matrix with variances on diagonal.

- covariance:

  A single value to apply to all off-diagonals, or a named vector where
  names follow the "dim1-dim2" format.

- tol:

  Small value to subtract from limits to ensure strict PD.

## Value

A list containing `up_mat` (updated matrix) and `limits` (data frame of
safe limits for remaining zero-cov pairs).

## Examples

``` r
# Setup
v_mat_3d <- matrix(c(10, 0, 0, 0, 5, 0, 0, 0, 7), 3, 3)
colnames(v_mat_3d) <- c("temp", "precip", "hum")

# 1. 2D Update: Change only the single available covariance
update_covariance(v_mat_3d[1:2, 1:2], covariance = c("temp-precip" = 3))
#> $updated_matrix
#>      temp precip
#> [1,]   10      3
#> [2,]    3      5
#> 
#> $remaining_limits
#> NULL
#> 

# 2. 3D Full Update: Set all combinations to 2.0
update_covariance(v_mat_3d, covariance = 2.0)
#> $updated_matrix
#>      temp precip hum
#> [1,]   10      2   2
#> [2,]    2      5   2
#> [3,]    2      2   7
#> 
#> $remaining_limits
#> NULL
#> 

# 3. 3D Partial: Define two, find limits for the remaining precip-hum
update_covariance(v_mat_3d, covariance = c("temp-precip" = -4, "temp-hum" = 5))
#> $updated_matrix
#>      temp precip hum
#> [1,]   10     -4   5
#> [2,]   -4      5   0
#> [3,]    5      0   7
#> 
#> $remaining_limits
#>                  min     max
#> precip-hum -5.911512 1.91152
#> 
```
