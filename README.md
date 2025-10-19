
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NicheR: An R package for ellipsoid-based virtual species and niche visualization in environmental (E-space) and geographic (G-space) space

<!-- badges: start -->

<!-- badges: end -->

The goal of **NicheR** is to provide a robust set of tools for
researchers and modelers to:

- Construct and define virtual niches using ellipsoid geometries in 2D
  or 3D environmental dimensions.  
- Identify and extract suitable environmental areas from various data
  sources (rasters, data frames).  
- Simulate species occurrence points within suitable habitats using
  flexible sampling strategies (random, center-biased, or
  edge-biased).  
- Visualize niche boundaries and simulated occurrences in both
  environmental space (E-space) and geographic space (G-space).

Inspired by the methodologies of **NicheA** and the **virtualspecies** R
package, **NicheR** streamlines niche conceptualization and data
generation for ecological studies.

## Authors

Mariana Castaneda-Guzman

Paanwaris Paansri

Connor Hughes

## Installation

You can install the development version of NicheR from GitHub as
follows:

``` r
# 1. Install devtools if you don't have it yet
if (!require("devtools")) {
 install.packages("devtools")
}

# 2. Install bean from GitHub
devtools::install_github("castanedaM/NicheR")
```

## Dependencies NicheR relies on the following R packages, which will be

installed automatically with the package: dplyr, ggplot2, ggpubr,
plotly, RColorBrewer, and terra.

### Example This is a basic example demonstrating the core functionality

of NicheR.

### 1. Load Data and Define Environmental Background We’ll load NicheR and

prepare some example environmental raster data. For this, we’ll use
example files that should be located in your package’s inst/extdata
directory.

``` r
library(NicheR)
library(terra)
#> terra 1.8.70
library(ggplot2)
library(ggpubr)
#> 
#> Attaching package: 'ggpubr'
#> The following object is masked from 'package:terra':
#> 
#>     rotate
library(rgl)
library(plotly)
#> 
#> Attaching package: 'plotly'
#> The following object is masked from 'package:ggplot2':
#> 
#>     last_plot
#> The following object is masked from 'package:stats':
#> 
#>     filter
#> The following object is masked from 'package:graphics':
#> 
#>     layout

# Access example environmental raster files from the package
# Make sure "Bio1.tif", "Bio4.tif", "Bio12.tif" exist in inst/extdata
bio_1 <- terra::rast(system.file("extdata", "Bio1.tif", package = "NicheR"))
bio_4 <- terra::rast(system.file("extdata", "Bio4.tif", package = "NicheR"))
bio_12 <- terra::rast(system.file("extdata", "Bio12.tif", package = "NicheR"))

# Stack and name the environmental layers
env_stack <- c(bio_1, bio_12, bio_4)
names(env_stack) <- c("mean_temp", "temp_seasonality", "annual_precip")

# Convert the raster stack to a data frame for E-space plotting
env_df <- as.data.frame(env_stack, xy = TRUE, na.rm = TRUE)
```

### 2. Build a 3D Ellipsoid Niche The build_ellipsoid() function creates

the core niche object by specifying its center, semi-axis lengths, and
rotation angles.

``` r
# Define niche parameters
niche_center <- c(mean_temp = 280, temp_seasonality = 1500, annual_precip = 2000)
niche_axes <- c(100, 500, 2500) # Semi-axis lengths (breadth)
niche_angles <- c(0, 0, pi/12) # Small rotation around Z-axis

# Create the ellipsoid niche object
my_niche <- build_ellps(
  center = niche_center,
  axes = niche_axes,
  # angles = niche_angles,
  n_points = 50 # Resolution for surface points
)
#> Warning in build_ellps(center = niche_center, axes = niche_axes, n_points =
#> 50): All axes will be applied the same angle specified: 0
```

### 3. Visualize the Niche in Environmental Space (E-space)

Function ‘plot_e_space()’ allows you to visualize e-sapce and the niche
ellipsoidal boundary.

It can generate both static 2D pairwise plots and interactive 3D plots.

#### 3.1. 2D Pairwise Plots

``` r
plot_e_space(
  env_bg = env_df,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  labels = c("Mean T (°C)", "Temp Seasonality", "Annual Prec (mm)"),
  niche = my_niche, 
  plot.3d = F
)
#> Sampling 10000 of 12412229 rows from 'env_bg' for plotting.
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

#### 3.2. 3D Interactive Plot

To render a 3D interactive plot within the same plotting function as 2D,
add the argument `plot.3d = TRUE`.

The code will then generate an interactive **plotly** plot.

``` r

plot_e_space(
  env_bg = env_df,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  labels = c("Mean T (°C)", "Temp Seasonality", "Annual Prec (mm)"),
  niche = my_niche,
  plot.3d = TRUE # Toggle 3D plotting by setting plot.3d = TRUE
)
```

<figure>
<img src="man/figures/3D.png" alt="3D plot Example" />
<figcaption aria-hidden="true">3D plot Example</figcaption>
</figure>

### 4. Extract Suitable Environments and Sample Occurrences

The `get_suitable_env()` function identifies all environmental points
that fall within the niche.  
It includes an argument `out` that can take the values `"spatial"`,
`"data.frame"`, or `"both"`.

- `"spatial"` returns a **SpatRaster**, allowing visualization of the
  suitable area in geographic space (G-space).  
- `"data.frame"` returns a table of points, better suited for
  visualization in environmental space (E-space).  
- `"both"` returns both formats simultaneously.

``` r
# Extract suitable environmental areas as a spatial raster
suitable_g_space <- get_suitable_env(
  niche = my_niche,
  env_bg = env_stack, # Use the raster stack for spatial output
  out = "spatial"
)

# Plot the suitable area on a map
plot_g_space(env_stack, show.suitable  = TRUE, niche = my_niche)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

**Note:** For 3D and 2D plotting in E-space, there is no need to run the
function beforehand.

Simply toggle the argument `show.pts.in` in `plot_e_space()`, and the
function will automatically calculate and plot the suitable
environments.

``` r
# Plot in 2D
plot_e_space(
  env_bg = env_df,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  niche = my_niche,
  show.pts.in = TRUE # Change this argument to TRUE to show suitable environments
)
#> Sampling 10000 of 12412229 rows from 'env_bg' for plotting.
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

Next, let’s sample some occurrences within the suitable area.  
The function `get_sample_occ()` draws a subset of points from within the
niche boundary and suitable environment.

``` r
# Sample 500 occurrences, biased towards the center of the niche
sampled_occ_center <- get_sample_occ(
  n_occ = 500,
  niche = my_niche,
  env_bg = env_stack, # Can use raster or data.frame here
  method = "center",
  seed = 42
)
#> Done sampling 500 occurrences

# Expand right margin so there's room for legend
plot_g_space(env_bg = env_stack, niche = my_niche,
             show.suitable = FALSE,
             occ_pts = sampled_occ_center)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Visualize Sampled Occurrences in Environmental Space You can integrate
these simulated occurrences back into your environmental space plots.

``` r
# Plot in 2D
plot_e_space(
  env_bg = env_df,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  niche = my_niche,
  occ_pts = sampled_occ_center,
  show.occ.density = TRUE # Change to remove density plots for occurrences points
)
#> Sampling 10000 of 12412229 rows from 'env_bg' for plotting.
#> Warning: Removed 1 row containing non-finite outside the scale range
#> (`stat_density()`).
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

Finally, we visualize all objects in one visualization.

``` r
# Plot in 2D
plot_e_space(
  env_bg = env_df,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  niche = my_niche,
  occ_pts = sampled_occ_center,
  show.pts.in = TRUE, # Highlight background points within niche
  show.occ.density = TRUE,
  palette = "default"# Show marginal density plots for occurrences
)
#> Sampling 10000 of 12412229 rows from 'env_bg' for plotting.
#> Warning: Removed 1 row containing non-finite outside the scale range
#> (`stat_density()`).
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

Try some of the other color palettes within the Package!!

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

Example: Palette 4

    #> Sampling 10000 of 12412229 rows from 'env_bg' for plotting.
    #> Warning: Removed 1 row containing non-finite outside the scale range
    #> (`stat_density()`).

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

Example: Palette 6

    #> Sampling 10000 of 12412229 rows from 'env_bg' for plotting.
    #> Warning: Removed 1 row containing non-finite outside the scale range
    #> (`stat_density()`).

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

Don’t like these palettes? Provide your own! Use the colors argument
with a list of hex codes. If colors are unnamed, they will be assigned
in order: bg (background), ellipsoid, centroid, tolerance (axis ranges),
suitable_env (environment inside ellipsoid), and occ (occurrence
records). Best practice is to use named values, for example:

``` r
plot_e_space(
  env_bg = env_df,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  niche = my_niche,
  show.pts.in = TRUE,
  occ_pts = sampled_occ_center,
  show.occ.density = TRUE,
  colors = list(bg = "lightgrey",
                ellipsoid = "lightblue", 
                centroid = "#FF6347", 
                tolerance = "orange",
                suitable_env = "#EEEE00",
                occ = "purple4")
  )
#> Sampling 10000 of 12412229 rows from 'env_bg' for plotting.
#> Warning: Removed 1 row containing non-finite outside the scale range
#> (`stat_density()`).
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

**Note**: for 3D plotting is better to provide the HEX color codes!

Want to save your plot in high resolution, publication-ready format?
Below is an example of how to save it as a TIFF with recommended sizing
for 2D plots.

``` r

tiff(filename = "NicheR_plot2D.tiff", width = 11, height = 8, 
     units = "in", res = 300, compression = "lzw")

plot_e_space(
  env_bg = env_df,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  labels = c("Mean T (°C)", "Temp Seasonality", "Annual Prec (mm)"),
  niche = my_niche,
  show.pts.in = TRUE, 
  occ_pts = sampled_occ_center,
  show.occ.density = TRUE,
  palette = "palette3"
)

dev.off()
```

### Contributing

We welcome contributions! If you have suggestions, bug reports, or would
like to contribute code, please feel free to open an issue or submit a
pull request on the NicheR GitHub repository.
