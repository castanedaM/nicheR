
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NicheR: Ellipsoid-Based Virtual Species and Niche Visualization in E-space and G-space

NicheR is an R package for building, visualizing, and sampling
**ellipsoid-based ecological niches** using environmental data.

It provides:

- Explicit ellipsoid construction in n-dimensional environmental
  space.  
- Efficient detection of suitable environments using Mahalanobis
  distance.  
- Sampling of virtual occurrences with flexible strategies and optional
  bias.  
- Visualization tools in environmental space (E-space) and geographic
  space (G-space).  
- A unified R-native workflow, with both high-level and low-level
  functions.

Inspired by the conceptual foundations of **NicheA** and the scripting
flexibility of the **virtualspecies** package, **nicheR** provides a
reproducible, scriptable ecosystem for virtual species simulation and
niche theory exploration.

------------------------------------------------------------------------

## Authors

- Mariana Castaneda-Guzman
- Connor Hughes  
- Paanwaris Paansri  
- Marlon Cobos

## Last Updated

01/21/2025

------------------------------------------------------------------------

## Installation

You can install the development version of NicheR from GitHub:

``` r
# 1. Install devtools if you do not have it yet
if (!require("devtools")) {
  install.packages("devtools")
}

# 2. Install NicheR from GitHub
devtools::install_github("castanedaM/NicheR")
```

**Dependencies**

NicheR relies on the following R packages (installed automatically when
you install NicheR):

- `terra` – raster / spatial data handling  
- `ggplot2`, `ggpubr` – plotting and layout  
- `plotly`, `rgl` – interactive and 3D visualization  
- `RColorBrewer` – color palettes  
- plus standard base/system packages

------------------------------------------------------------------------

# Quick Start: Full Workflow in One Step

The easiest way to use NicheR is through the high-level wrapper
`create_virtual_species()`, which:

1.  Builds an ellipsoid niche in E-space  
2.  Computes suitable environments  
3.  Samples occurrence records

------------------------------------------------------------------------

# Contributing

We welcome contributions!

If you have:

- suggestions,  
- bug reports, or  
- would like to contribute code,

please open an issue or submit a pull request on the **NicheR** GitHub
repository.

When contributing, please:

- keep functions documented with roxygen2,  
- include minimal reproducible examples, and  
- add tests whenever appropriate.

Thank you for helping improve NicheR!
