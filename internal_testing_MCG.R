# This is use internally to test function much quicker and make sure we know all
# the parameters required

# Restart session and do not load any packaged besides NicheR. This ensures that
# we have loaded all package dependencies into the R package

library(NicheR)

# Data and Test Basic Niche -----------------------------------------------

# Load or read in the environmetal conditions
bio_1 <- terra::rast(system.file("extdata", "Bio1.tif", package = "NicheR"))
bio_4 <- terra::rast(system.file("extdata", "Bio4.tif", package = "NicheR"))
bio_12 <- terra::rast(system.file("extdata", "Bio12.tif", package = "NicheR"))

# Stack and name the environmental layers
env_stack_small <- c(bio_1, bio_12, bio_4)
names(env_stack_small) <- c("mean_temp", "temp_seasonality", "annual_precip")

# Define niche parameters
niche_center <- c(mean_temp = 200,
                  temp_seasonality = 1000,
                  annual_precip = 6500)
niche_axes <- c(75, 250, 2000) # Semi-axis lengths (breadth)

my_niche_small <- build_ellps(center = niche_center,
                              axes = niche_axes)



# Load or read in the environmental conditions
env_1 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_pr_sum_mm-mosaic.tif")
env_2 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_tmmn_min_C-mosaic.tif")
env_3 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_tmmx_max_C-mosaic.tif")

# Stack and name the environmental layers
env_stack_large <- c(env_1, env_2, env_3)
names(env_stack_large) <- c("pr_sum", "tmmn_min", "tmmx_max")

niche_center <- c(0.4, 25, 30)
niche_axes <- c(0.75, 0.75, 0.75) # Semi-axis lengths (breadth)

my_niche_large <- build_ellps(center = niche_center,
                              axes = niche_axes)

# TEST 1: Basic Plotting --------------------------------------------------

env_df_small <- as.data.frame.nicheR(raster_stack = env_stack_small)

# Plot E-space
plot_e_space(env_bg = env_stack_small) # should not fail

plot_e_space(env_bg = env_df_small)

plot_e_space(env_bg = env_df_small,
             plot.3d = TRUE)

plot_e_space(env_bg = env_df_small,
             x = "mean_temp",
             y = "temp_seasonality",
             z = "annual_precip",
             niche = my_niche_small)

plot_e_space(env_bg = env_df_small,
             niche = my_niche_small,
             plot.3d = TRUE)

# Plot G-space
plot_g_space(env_bg = env_df_small)

plot_g_space(env_bg = env_df_small, palette = "palette2")


# TEST 2: Suitability Function (Small Raster Stack) -----------------------

# IN: Stack (small)
# OUT: data frame
suitable_area_t1 <- get_suitable_env(niche = my_niche_small,
                                     env_bg = env_stack_small,
                                     out = "data.frame",
                                     distances = TRUE)

# IN: Stack (small)
# OUT: spatial
suitable_area_t2 <- get_suitable_env(niche = my_niche_small,
                                     env_bg = env_stack_small,
                                     out = "spatial",
                                     distances = TRUE)


# IN: Stack (small)
# OUT: both
suitable_area_t3 <- get_suitable_env(niche = my_niche_small,
                                     env_bg = env_stack_small,
                                     out = "both",
                                     distances = TRUE)

# IN: Data Frame (small)
# OUT: both
suitable_area_t4 <- get_suitable_env(niche = my_niche_small,
                                     env_bg = env_df_small,
                                     out = "both",
                                     distances = TRUE)


# TEST 3: Suitable Plotting -----------------------------------------------

plot_e_space(niche = my_niche_small,
             suitable_env = suitable_area_t4,
             plot.3d = FALSE)

plot_g_space(suitable_env = suitable_area_t1)

plot_g_space(env_bg = env_stack_small,
             suitable_env = suitable_area_t1)

plot_g_space(env_bg = env_df_small,
             suitable_env = suitable_area_t1,
             surface = "suit")

plot_g_space(env_bg = env_df_small,
             suitable_env = suitable_area_t2$dist_sq,
             surface = "dist")


# Test 4: Bias Surface Function -------------------------------------------
# set_bias_surface <- function(bias_surface,
#                              bias_dir = 1,
#                              suitable_env = NULL,
#                              out.bias = c("biased", "standardized", "both"),
#                              verbose = TRUE)

sp_rich <- terra::rast(system.file("extdata", "sr_sp_host_0.05.tif", package = "NicheR"))
nighttime <- terra::rast(system.file("extdata", "nighttime.tif", package = "NicheR"))


# IN: bias surface
# OUT: standardized
set_bias_surface_t1 <- set_bias_surface(bias_surface = c(sp_rich, nighttime),
                                        bias_dir = c(1, -1),
                                        out = "standardized")
# needs to fail
set_bias_surface_t2 <- set_bias_surface(bias_surface = c(sp_rich, nighttime),
                                        bias_dir = c(1, -1),
                                        out = "both")

# IN: bias surface, and suitable_area_t2
# OUT: standardized
set_bias_surface_t3 <- set_bias_surface(bias_surface = c(sp_rich, nighttime),
                                        bias_dir = c(1, -1),
                                        suitable_env = suitable_area_t2,
                                        out = "both")


# IN: bias surface, and suitable_area_t3
# OUT: standardized
set_bias_surface_t4 <- set_bias_surface(bias_surface = c(sp_rich, nighttime),
                                        bias_dir = c(1, -1),
                                        suitable_env = suitable_area_t3,
                                        out = "biased")


# Test 5: Sample Occurence Function ---------------------------------------

# get_sample_occ <- function(n_occ,
#                            suitable_env,
#                            method = c("random", "center", "edge"),
#                            bias_surface = NULL,
#                            seed = NULL,
#                            verbose = TRUE)


# IN: n_occ = 100, random
# OUT: DF

get_sample_occ_t1 <- get_sample_occ(n_occ = 100,
                                    suitable_env = suitable_area_t1)

plot_g_space(occ_pts = get_sample_occ_t1,
             suitable_env = suitable_area_t1, surface = "suit")

plot_e_space(env_bg = env_df_small,
             suitable_env = suitable_area_t1,
             occ_pts = get_sample_occ_t1,
             niche = my_niche_small)

# Bug in plotting: suitbale area does not plot if niche is not give, it should
# not matter, also it takes a long time so suitable are should be thinned too
# And points are SMALL

# IN: n_occ = 100, center
# OUT: DF
get_sample_occ_t2 <- get_sample_occ(n_occ = 100,
                                    suitable_env = suitable_area_t1,
                                    method = "center")

plot_g_space(occ_pts = get_sample_occ_t2)


# TEST if bias works
get_sample_occ_t3 <- get_sample_occ(n_occ = 100,
                                    suitable_env = suitable_area_t2$dist_sq,
                                    method = "edge", seed = 123)

get_sample_occ_t4 <- get_sample_occ(n_occ = 100,
                                    suitable_env = suitable_area_t3,
                                    method = "edge", seed = 123)

get_sample_occ_t5 <- get_sample_occ(n_occ = 100,
                                    suitable_env = suitable_area_t3,
                                    bias_surface = set_bias_surface_t3$pooled_bias_sp,
                                    method = "edge", seed = 123)
# Tiny bug to read the bias surface right

plot_g_space(occ_pts = get_sample_occ_t3)
plot_g_space(occ_pts = get_sample_occ_t4) #above should be the same
plot_g_space(occ_pts = get_sample_occ_t5)


plot_e_space(env_bg = env_df_small, niche = my_niche_small,
             x = "mean_temp",
             y = "temp_seasonality",
             z = "annual_precip",
             labels = c("Mean T (°C)", "Temp Seasonality", "Annual Prec (mm)"),
             palette = "default",
             occ_pts = get_sample_occ_t1)



# TEST #: Large raster ----------------------------------------------------

# Will take long, did all at once below
# ptm <- proc.time()
# env_df_db <- as.data.frame.nicheR(raster_stack = env_stack_large,
verbose = TRUE, use_cache = TRUE)
# proc.time() - ptm

gc()

plot_e_space(env_bg = env_df_db, niche = my_niche_large,
             x = "pr_sum", y = "tmmn_min", z = "tmmx_max")


# This should fail
suitbale_area <- get_suitable_env(niche = my_niche_large,
                                  env_bg = env_stack_large,
                                  out = "both")


ptm <- proc.time()
suitbale_area <- get_suitable_env(niche = my_niche_large,
                                  env_bg = env_df_db,
                                  out = "both")
proc.time() - ptm

gc()

plot_e_space(env_bg = env_df_db,
             niche = my_niche_large,
             suitable_env = suitbale_area,
             x = "pr_sum", y = "tmmn_min", z = "tmmx_max")




# TEST #: Create virtual species (Small) -----------------------------------------------

# Load or read in the environmetal conditions
bio_1 <- terra::rast(system.file("extdata", "Bio1.tif", package = "NicheR"))
bio_4 <- terra::rast(system.file("extdata", "Bio4.tif", package = "NicheR"))
bio_12 <- terra::rast(system.file("extdata", "Bio12.tif", package = "NicheR"))

# Stack and name the environmental layers
env_stack_small <- c(bio_1, bio_12, bio_4)
names(env_stack_small) <- c("mean_temp", "temp_seasonality", "annual_precip")

# Define niche parameters
niche_center <- c(mean_temp = 200,
                  temp_seasonality = 1000,
                  annual_precip = 6500)
niche_axes <- c(75, 250, 2000)

# Bias
sp_rich <- terra::rast(system.file("extdata", "sr_sp_host_0.05.tif", package = "NicheR"))
nighttime <- terra::rast(system.file("extdata", "nighttime.tif", package = "NicheR"))


vs_small <- create_virtual_species(env_bg = env_stack_small,
                                   center = niche_center,
                                   axes = niche_axes,
                                   n_occ = 75,
                                   # bias_surface = c(sp_rich, nighttime),
                                   # bias_dir = c(1, -1),
                                   out.suit = "both",
                                   out.bias = "both",
                                   distances = TRUE,
                                   out.file = TRUE)
plot_e_space(vs = vs_small)
plot_e_space(env_bg  = env_stack_small, vs = vs_small)

plot_g_space(vs = vs_small, surface = "both")

plot_g_space(env_bg = env_stack_small)

plot_g_space(env_bg = nr_get_env(vs),
             # occ_pts = nr_get_occ(vs),
             suitable_env = nr_get_dist_sq(vs),
             surface = "dist")


# Explore virtual species object
vs_small
vs_small$niche
head(vs_small$occurrences)
vs_small$suitability
vs_small$call_args
vs_small$save_path
vs_small$bias_surface

# reload later
# readRDS(vs_1$save_path)


# TEST #: Create Virtual Species (Large) ----------------------------------

# Load or read in the environmental conditions
env_1 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_pr_sum_mm-mosaic.tif")
env_2 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_tmmn_min_C-mosaic.tif")
env_3 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_tmmx_max_C-mosaic.tif")

# Stack and name the environmental layers
env_stack_large <- c(env_1, env_2, env_3)
names(env_stack_large) <- c("pr_sum", "tmmn_min", "tmmx_max")

niche_center <- c(0.4, 25, 30)
niche_axes <- c(0.75, 0.75, 0.75) # Semi-axis lengths (breadth)


# Bias
sp_rich <- terra::rast(system.file("extdata", "sr_sp_host_0.05.tif", package = "NicheR"))
nighttime <- terra::rast(system.file("extdata", "nighttime.tif", package = "NicheR"))

ptm <- proc.time()
env_df_db <- as.data.frame.nicheR(raster_stack = env_stack_large,
                                  verbose = TRUE, use_cache = TRUE)
proc.time() - ptm

gc()

ptm <- proc.time()
vs_large <- create_virtual_species(env_bg = env_df_db,
                                   center = niche_center,
                                   axes = niche_axes,
                                   out.suit = "both",
                                   n_occ = 75,
                                   bias_surface = c(sp_rich, nighttime),
                                   bias_dir = c(1, -1),
                                   out.bias = "both",
                                   distances = TRUE,
                                   out.file = FALSE)
proc.time() - ptm

plot_e_space(env_bg = env_df_db, vs = vs_large)

plot_g_space(env_bg = env_df_db, vs = vs_large)




