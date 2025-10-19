# This is use internally to test function much quicker and make sure we know all
# the parameters required

# Restart session and do not load any packaged besides NicheR. This ensures that
# we have loaded all package dependencies into the R package

library(NicheR)

# 1. Data Preparation -----------------------------------------------------

# Load or read in the environmetal conditions
bio_1 <- terra::rast(system.file("extdata", "Bio1.tif", package = "NicheR"))
bio_4 <- terra::rast(system.file("extdata", "Bio4.tif", package = "NicheR"))
bio_12 <- terra::rast(system.file("extdata", "Bio12.tif", package = "NicheR"))

# Stack and name the environmental layers
env_stack <- c(bio_1, bio_12, bio_4)
names(env_stack) <- c("mean_temp", "temp_seasonality", "annual_precip")

# Convert the raster stack to a data frame (not nesseary but can save you some
# time in the latters steps)
ptm <- proc.time()
env_df <- as.data.frame(env_stack, xy = TRUE, na.rm = TRUE)
proc.time() - ptm

# If file small this takes longer
ptm <- proc.time()
env_df_db <- convert_large_raster(raster_stack = env_stack,
                                  out_filename = "env_stack_bio_1_4_12",
                                  chunk_height = 1000,
                                  verbose = TRUE)
proc.time() - ptm


# 2. Specify the Fundamental Niche ----------------------------------------

# Define niche parameters
niche_center <- c(mean_temp = 280,
                  temp_seasonality = 1500,
                  annual_precip = 2000)
niche_axes <- c(100, 500, 2500) # Semi-axis lengths (breadth)
niche_angles <- c(0, 0, pi/12) # Small rotation around Z-axis


# 3. Create virtual species -----------------------------------------------


vs_1 <- create_virtual_species(env_bg = env_stack,
                               center = niche_center,
                               axes = niche_axes,
                               angles = niche_angles,
                               n_occ = 100,
                               out = "both",
                               distances = TRUE,
                               out.file = FALSE)

vs_1
vs_1$niche
head(vs_1$occurrences)
vs_1$suitability
vs_1$call_args
vs_1$save_path


# reload later
# readRDS(vs_1$save_path)



# 4. Visualize virtual species object -------------------------------------


# Visualize the niche
plot_e_space(env_bg = env_df,
             x = "mean_temp",
             y = "temp_seasonality",
             z = "annual_precip",
             labels = c("Mean T (°C)", "Temp Seasonality", "Annual Prec (mm)"),
             niche = vs_1$niche)

plot_g_space(env_bg = vs_1$call_args$env_bg,
             niche = vs_1$niche,
             show.suitable = FALSE,
             show.distance = TRUE,
             occ_pts = vs_1$occurrences)


# 1. Our data -------------------------------------------------------------


# Load or read in the environmetal conditions
env_1 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_pr_sum_mm-mosaic.tif")
env_2 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_tmmn_min_C-mosaic.tif")
env_3 <- terra::rast("../../30yrNormals/30yrNormals/ERA5_30yr_Normal_tmmx_max_C-mosaic.tif")

# Stack and name the environmental layers
env_stack <- c(env_1, env_2, env_3)
names(env_stack) <- c("pr_sum", "tmmn_min", "tmmx_max")

# If file small this takes longer
ptm <- proc.time()
env_df_db <- convert_large_raster(raster_stack = env_stack,
                                  out_filename = "env_stack_30yrNormals",
                                  chunk_height = 1000,
                                  verbose = TRUE)
proc.time() - ptm

head(env_df_db)


# 2. Specify the Fundamental Niche ----------------------------------------

# Define niche parameters
niche_center <- c(0.2, 25, 29)
niche_axes <- c(0.05, 0.05, 0.05) # Semi-axis lengths (breadth)


# 3. Create virtual species -----------------------------------------------

vs_1 <- create_virtual_species(env_bg = env_df_db,
                               center = niche_center,
                               axes = niche_axes,
                               n_occ = 100,
                               distances = TRUE,
                               out.file = FALSE)

vs_1

plot_e_space(env_bg = env_df_db,
             x = "pr_sum",
             y = "tmmn_min",
             z = "tmmx_max",
             niche = vs_1$niche)

