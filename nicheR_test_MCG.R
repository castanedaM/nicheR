
devtools::load_all()


# Data collection and preparation ----------------------------------------

bios <- terra::rast("inst/extdata/ma_bios.tif")
bios <- bios[[c("bio_1", "bio_12", "bio_15")]]

bios_df <- as.data.frame(bios, xy = TRUE)

# Bias layers
# For this example our bias_layers have the speciens richeness and nightitme
# species richeness as it increases is a good thing, but nighttime as it increrases is a bad thing, so the effct of nighttime has to be inverse
bias <- prepare_bias(bias_surface = terra::rast("inst/extdata/ma_biases.tif"),
                     effect_direction = c("direct", "inverse"),
                     include_processed_layers = TRUE)


# Here you can see how the standariszation happened and what inverse and not was applied
bias$combination_formula

# To see the inside of the bias object prepare bias created you can access it with opertion $ and the layer has a name composit bias. for more deatils go to the vigneette about bias (coming)
bias$composite_surface

# Visuale the resulting bias, this includes the standarization of the original bias layers and the coposite bias since the argument include_processed_layers was TRUE, by daulft composite is alwyas returned.
terra::plot(bias$processed_layers)
terra::plot(bias$composite_surface)


# Build Ellipsoid ---------------------------------------------------------
# Now we build our elliposid at it simplest we give it a range
range <- data.frame(bio_1 = c(20, 25),
                    bio_12 = c(1000, 1500),
                    bio_15 = c(60, 75))

ell <- build_ellipsoid(range = range)

# This is the information contained in that object
names(ell)

# It also has a pretty print
ell

# You can plot using function in nicheR designed in baser r to simplify the plotting
plot_ellipsoid(ell)
plot_ellipsoid_pairs(ell)

# To visualize with the bacgound you need to use the data frame of your bios, and it is alwyas recommended to sample the bacgorund for quicker plotting, this is set to NUll so if large it will take a long time. and change dim to see others.
plot_ellipsoid(ell,
               background = bios_df, dim = c(1,2),
               bg_sample = 5000)

plot_ellipsoid(ell,
               background = bios_df, dim = c(1,3),
               bg_sample = 5000)

# Or use the pairs
plot_ellipsoid_pairs(ell, background = bios_df, bg_sample = 5000, col_ell = "red")

# Predict Suitability -----------------------------------------------------

# Now lets predict suitbaility

# When RASTER data in
ell_predict_r <- predict(ell, #ellipsoid object
                         newdata = bios, #backgorund
                         # what layers to incude and if truncation aka only
                         # inside the elliposid is wanted as a return output,
                         # deafualt to FALSE
                         suitability_truncated   = TRUE,
                         mahalanobis_truncated   = TRUE,
                         keep_data = TRUE #if you want to keep your original data
)

# Plot all predictions, and sicne keep_data was true, this will aslo include the bios layers
terra::plot(ell_predict_r)


# DATA FRAME input
ell_predict_df <- predict(ell,
                          newdata = bios_df,
                          suitability_truncated   = TRUE)
head(ell_predict_df) # default suitability and mahalanobis, no truncation

ell_predict_df_r <- as.data.frame(ell_predict_r, xy = TRUE)

# For DF
set.seed(123)

# No TRUNCATED ZOOMED OUT
plot_ellipsoid(ell,
               background = ell_predict_df,
               pch = 20,
               bg_sample = 5000)
add_data(data = ell_predict_df, x = "bio_1", y = "bio_12",
         col_layer = "suitability",
         rev_pal = T,
         pal = terrain.colors(100),
         bg_sample = 10000,
         pch = 20)
add_ellipsoid(ell, col_ell = "red")

# TRUCATED ZOOMED IN
plot_ellipsoid(ell)
add_data(data = ell_predict_df, x = "bio_1", y = "bio_12",
         col_layer = "suitability_trunc", rev_pal = T, pch = 20)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


# Apply bias to prediction ------------------------------------------------

biased_predict_r <- apply_bias(prepared_bias = bias,
                               prediction = ell_predict_r,
                               prediction_layer = "suitability")

# Sample data -------------------------------------------------------------

# UNBIASED SAMPLING - Raster if we dont want to use the bias layers we use sample_data
sample_data_r <- sample_data(n_occ = 100,
                             prediction = ell_predict_r,
                             prediction_layer = "suitability",
                             sampling = "centroid",
                             method = "suitability",
                             strict = TRUE)
head(sample_data_r)
# To do: check for suitability prediction layer with mahalanobis method

# UNBIASED SAMPLING - Data frame
sample_data_df_cn <- sample_data(n_occ = 100,
                                 prediction = ell_predict_df,
                                 prediction_layer = "suitability",
                                 sampling = "centroid",
                                 method = "mahalanobis")
head(sample_data_df_cn)


plot_ellipsoid(ell, main = "Sampled Center Unbiased")
add_data(data = ell_predict_df,
         x = "bio_1", y = "bio_12",
         pts_col = "grey",
         pch = 20,
         bg_sample = 8000)
add_data(data = sample_data_df_cn,
         x = "bio_1", y = "bio_12",
         pts_col = "black",
         pch = 4,
         lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)

# And important difference is that in the privous example the sample values can still be outside the ellispois since the sampling was not restricted, to restrict set stric = true and use the trucated layer of suitbiality/ mahalanobis accordingly
sample_data_df_cn <- sample_data(n_occ = 100,
                                 prediction = ell_predict_df,
                                 prediction_layer = "suitability_trunc",
                                 sampling = "centroid",
                                 method = "mahalanobis",
                                 strict = TRUE)
head(sample_data_df_cn)


plot_ellipsoid(ell, main = "Sampled Center Unbiased")
add_data(data = ell_predict_df,
         x = "bio_1", y = "bio_12",
         pts_col = "grey",
         pch = 20,
         bg_sample = 8000)
add_data(data = sample_data_df_cn,
         x = "bio_1", y = "bio_12",
         pts_col = "black",
         pch = 4,
         lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)

# BIASED SAMPLING - Raster
sample_data_b <-
  sample_biased_data(n_occ = 100,
                     prediction = biased_predict_r,
                     prediction_layer = "suitability_biased_direct")

# Notice how even with raster input the ourput is a data frame, this are not spactial points.
head(sample_data_b)


# UNBIASED SAMPLING - Data frame - Edge, notice that if layer is truncated your message will also say strict beacome true based on detection
sample_data_df_ed <- sample_data(n_occ = 100,
                                 prediction = ell_predict_df,
                                 prediction_layer = "suitability_trunc",
                                 sampling = "edge",
                                 method = "suitability")

plot_ellipsoid(ell,
               main = "Sampled edge Unbiased",
               # background = bios_df,
               # bg_sample = 10000,
               pch = 20)
add_data(data = ell_predict_df,
         x = "bio_1", y = "bio_12",
         pts_col = "grey",
         pch = 20,
         bg_sample = 8000)
add_data(data = sample_data_df_ed,
         x = "bio_1", y = "bio_12",
         pts_col = "black",
         pch = 4,
         lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


# PLOT FOR CENTER AND EDGE - normal data
par(mfrow = c(2,1))
plot_ellipsoid(ell,
               main = "Sampled center Unbiased")
add_data(data = ell_predict_df,
         x = "bio_1", y = "bio_12",
         pts_col = "grey",
         pch = 20,
         bg_sample = 8000)
add_data(data = sample_data_df_ed,
         x = "bio_1", y = "bio_12",
         pts_col = "black",
         pch = 4,
         lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)



plot_ellipsoid(ell,
               main = "Sampled edge Unbiased",
               pch = 20)
add_data(data = ell_predict_df,
         x = "bio_1", y = "bio_12",
         pts_col = "grey",
         pch = 20,
         bg_sample = 8000)
add_data(data = sample_data_df_ed,
         x = "bio_1", y = "bio_12",
         pts_col = "black",
         pch = 4,
         lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


# All of this can also be done purely with virtual data.. coming soon example

