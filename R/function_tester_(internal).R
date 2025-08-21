# Title: Internal File to test functionalities of each of the functions
#
# Authors: Mariana Castaneda-Guzman,
#          Paanwaris Paansri, and
#          Connor Hughes
#
# Date Created: 07/17/2025
# Date Last Updated: 08/21/2025

# Description: This is the Demo use of the package 'virtualniche'. With
# adapatations of NicheA and R package 'virtualspecies'

# Packages ----------------------------------------------------------------

library(ggplot2)
library(terra)
library(dplyr)
library(ggpubr)
library(rnaturalearth)

# Environmental Predictors ------------------------------------------------

# Dummy data

# Based on WorldClim Data
bio_1 <- terra::rast("inst/extdata/Bio1.tif")
bio_4 <- terra::rast("inst/extdata/Bio4.tif")
bio_12 <- terra::rast("inst/extdata/Bio12.tif")

world_sf <- ne_countries(scale = "medium", returnclass = "sf")

# Sample this distributions assume normality
bio_stack <- c(bio_1, bio_12, bio_4)
names(bio_stack) <- c("env_x", "env_y", "env_z")

bio_stack <- mask(bio_stack, vect(world_sf))

bio_df <- as.data.frame(bio_stack, xy = TRUE)

df <- bio_df
names(df) <- c("x", "y", "env_x", "env_y", "env_z")

if(nrow(df) > 10000){
  set.seed(1234)
  df_plotting <- df[sample(1:nrow(df), size = 10000, replace = FALSE), ]
}

# Visualize Data ----------------------------------------------------------

ggplot(data = df_plotting) +
  geom_density(aes(x = env_x), fill = "lightblue", alpha = 0.5) +
  geom_density(aes(x = env_y), fill = "lightgreen", alpha = 0.5) +
  geom_density(aes(x = env_z), fill = "red", alpha = 0.5) +
  theme_bw()

# pairs(df_plotting[3:5])
# corrplot::corrplot(cor(df_plotting))


# Visualize E-space -------------------------------------------------------
source("R/plot_e_space.R")

# 2D plot
plot_e_space(df,
             x = "env_x", y = "env_y", z = "env_z",
             labels = c("BIO 1", "BIO 4", "BIO 12"))


plot_e_space(df,
             x = "env_x", y = "env_y", z = "env_z",
             labels = c("BIO 1", "BIO 4", "BIO 12"),
             plot.3d = TRUE)


# Create Ellipsoids -------------------------------------------------------

source("R/build_ellipsoid.R")

# 3D Ellipsoid #
center <- c(300, 1500, 2000) # center
axes <- c(100, 500, 2500) # offsets
angles <- c(0, 0, pi/12)
n_points <- 50 # more for plotting

# Build Ellipsoid
FN_1 <- build_ellipsoid(center = center,
                        axes = axes,
                        angles = angles,
                        n_points = n_points)

# Visual in 2D
plot_e_space(df,
             x = "env_x", y = "env_y", z = "env_z",
             niche = FN_1)

# Visual in 3D
plot_e_space(df,
             x = "env_x", y = "env_y", z = "env_z",
             niche = FN_1, plot.3d = TRUE)


# Extract Points from Ellipsoid -------------------------------------------

source("R/get_suitable_environment.R")

# Extract points
pts_in <- get_suitable_environment(niche = FN_1,
                                   env_bg = df[, c("env_x", "env_y", "env_z")],
                                   out = "data.frame")

ras_in <- get_suitable_environment(niche = FN_1,
                                   env_bg = bio_stack,
                                   out = "spatial")

plot(ras_in, col = rev(terrain.colors(5)))

both_in <- get_suitable_environment(niche = FN_1,
                                   env_bg = bio_stack,
                                   out = "both")


# Visual in 2D
# Function for plotting does it inside of it, just set show.pts.in = TRUE
plot_e_space(env_bg = df, x = "env_x", y = "env_y", z = "env_z",
             niche = FN_1, show.pts.in = TRUE)

# Visual in 3D
# Function for plotting does it inside of it, just set show.pts.in = TRUE
plot_e_space(env_bg = df, x = "env_x", y = "env_y", z = "env_z",
             niche = FN_1, show.pts.in = TRUE, plot.3d = TRUE)


# Sampling Occurrences ----------------------------------------------------

source("R/get_sample_occurrences.R")


sampled_pts <- get_sample_occurrences(n_occ = 100,
                                      niche = FN_1,
                                      env_bg = df,
                                      seed = 101)

sampled_pts_random <- get_sample_occurrences(n = 100,
                                          niche = FN_1,
                                          env_bg = bio_stack,
                                          seed = 101)

sampled_pts_center <- get_sample_occurrences(n = 100,
                                          niche = FN_1,
                                          env_bg = bio_stack,
                                          method = "center")

sampled_pts_edge <- get_sample_occurrences(n = 100,
                                        niche = FN_1,
                                        env_bg = bio_stack,
                                        method = "edge")

occ_x <- ggplot() +
  geom_density(data = sampled_pts_random,
               aes(x = env_x, fill = "Random"), alpha = 0.5) +
  geom_density(data = sampled_pts_center,
               aes(x = env_x, fill = "Center"), alpha = 0.5) +
  geom_density(data = sampled_pts_edge,
               aes(x = env_x, fill = "Edge"), alpha = 0.5) +
  geom_vline(aes(xintercept = FN_1$center[1], color = "Centroid"), linetype = 2) +
  scale_fill_manual(name = "Sampling strategy",
                    values = c("Random" = "lightblue",
                               "Center" = "lightgreen",
                               "Edge" = "red")) +
  scale_color_manual(name = "",
                     values = c("Centroid" = "black")) +
  theme_bw()

occ_y <- ggplot() +
  geom_density(data = sampled_pts_random,
               aes(x = env_y, fill = "Random"), alpha = 0.5) +
  geom_density(data = sampled_pts_center,
               aes(x = env_y, fill = "Center"), alpha = 0.5) +
  geom_density(data = sampled_pts_edge,
               aes(x = env_y, fill = "Edge"), alpha = 0.5) +
  geom_vline(aes(xintercept = FN_1$center[2], color = "Centroid"), linetype = 2) +
  scale_fill_manual(name = "Sampling strategy",
                    values = c("Random" = "lightblue",
                               "Center" = "lightgreen",
                               "Edge" = "red")) +
  scale_color_manual(name = "",
                     values = c("Centroid" = "black")) +
  theme_bw()

occ_z <- ggplot() +
  geom_density(data = sampled_pts_random,
               aes(x = env_z, fill = "Random"), alpha = 0.5) +
  geom_density(data = sampled_pts_center,
               aes(x = env_z, fill = "Center"), alpha = 0.5) +
  geom_density(data = sampled_pts_edge,
               aes(x = env_z, fill = "Edge"), alpha = 0.5) +
  geom_vline(aes(xintercept = FN_1$center[3], color = "Centroid"), linetype = 2) +
  scale_fill_manual(name = "Sampling strategy",
                    values = c("Random" = "lightblue",
                               "Center" = "lightgreen",
                               "Edge" = "red")) +
  scale_color_manual(name = "",
                     values = c("Centroid" = "black")) +
  theme_bw()

ggpubr::ggarrange(occ_x, occ_y, occ_z,
                  legend = "right", ncol = 1,
                  common.legend = TRUE)

ggplot(sampled_pts_center, aes(x=env_y, y=env_x) ) +
  geom_bin2d() +
  theme_bw()
ggplot(sampled_pts_center, aes(x=env_z, y=env_x) ) +
  geom_bin2d() +
  theme_bw()
ggplot(sampled_pts_center, aes(x=env_z, y=env_y) ) +
  geom_bin2d() +
  theme_bw()



# 2D plot
plot_e_space(df, x = "env_x", y = "env_y", z = "env_z",
             niche = FN_1,
             show.pts.in = FALSE,
             occ_pts = sampled_pts)

# Show density
plot_e_space(df, x = "env_x", y = "env_y", z = "env_z",
             niche = FN_1,
             show.pts.in = FALSE,
             occ_pts = sampled_pts_center,
             show.occ.density = TRUE)

# 3D Plot
plot_e_space(df, x = "env_x", y = "env_y", z = "env_z",
             niche = FN_1,
             show.pts.in = TRUE,
             occ_pts = sampled_pts, plot.3d = TRUE)


# Plot in geography -------------------------------------------------------


# 1. Plot the background layer
plot(ras_in,
     main = "Suitable Geographic Environemt",
     col  = rev(terrain.colors(50)))

# 2. Overlay the sampled points
points(sampled_pts$x, sampled_pts$y,
       pch   = 20,         # solid circle
       col   = "red",
       cex   = 0.4)        # point size multiplier










