# nicheR_example_application.R ----
# Reproduces the Example Application section of the nicheR manuscript.
# Runs on data that ship with the package. Imports nicheR and terra only.
# Nothing here is fitted. The true niche is known throughout and is used only to
# describe where each sampling design put its records.

# library(nicheR)
library(terra)

out_dir <- "example_application"
dir.create(out_dir, showWarnings = FALSE)

# ENVIRONMENT ----
# Two bioclim layers from the package: annual mean temperature and annual precipitation.
bios_path <- system.file("extdata", "ma_bios.tif", package = "nicheR")
if(bios_path == ""){ bios_path <- system.file("ma_bios.tif", package = "nicheR") }
env <- rast(bios_path)
env <- env[[c("bio_1", "bio_12")]]

# NICHE ----
# Build one niche from ranges, tilt it with a covariance, then translate it.
# Everything after this point uses the translated niche, n2.
rng <- rbind(c(14, 375), c(26, 3625))
rownames(rng) <- c("min", "max")
colnames(rng) <- c("bio_1", "bio_12")

n1 <- build_ellipsoid(ranges = rng, cl = 0.99)
n11 <- update_ellipsoid_covariance(n1, covariance = c("bio_1-bio_12" = 750))
n2 <- update_ellipsoid_centroid(n11, new_centroid = c(bio_1 = 16, bio_12 = 1500))

cat("built centroid:", n11$centroid, "\n")
cat("moved centroid:", n2$centroid, "\n")
cat("true correlation:", round(cov2cor(n2$cov_matrix)[1, 2], 3), "\n")

# Niche SDs, for the table caption so a reader can scale the differences.
cat("niche SDs:", round(sqrt(diag(n2$cov_matrix)), 1), "\n")

# PREDICT ----
# pred is the working surface. pred_built exists only for the map panel, so the
# reader can see what update_ellipsoid_centroid() did. Delete both if not wanted.
env_df <- as.data.frame(env, xy = TRUE)

pred <- predict(n2, env,  keep_data = TRUE,
                mahalanobis_truncated  = TRUE,
                suitability_truncated = TRUE,
                include_suitability = FALSE, include_mahalanobis = FALSE)

cat("prediction layers:", names(pred), "\n")

suit <- pred[["suitability_trunc"]]
inside_sui <- suit > 0
terra::plot(inside_sui, main = "inside sui")

maha <- pred[["Mahalanobis_trunc"]]
inside_mah <- suit > 0
terra::plot(inside_mah, main = "inside mah")

cat("cells inside the niche sui:", sum(terra::values(inside_sui), na.rm = TRUE), "\n")
cat("cells inside the niche mah:", sum(terra::values(inside_sui), na.rm = TRUE), "\n")

# SAMPLING DESIGNS ----
# strict = TRUE is required: without it the edge design puts its heaviest
# weights on every zero cell outside the niche.
n_occ <- 100

# SUITABILITY
occ_centroid_sui <- sample_data(n_occ = n_occ, prediction = pred,
                            prediction_layer = "suitability_trunc", sampling = "centroid",
                            method = "suitability", seed = 1, strict = TRUE)

occ_random_sui <- sample_data(n_occ = n_occ, prediction = pred,
                          prediction_layer = "suitability_trunc", sampling = "random",
                          method = "suitability", seed = 2, strict = TRUE)

occ_edge_sui <- sample_data(n_occ = n_occ, prediction = pred,
                        prediction_layer = "suitability_trunc", sampling = "edge",
                        method = "suitability", seed = 3, strict = TRUE)

# MAHALANOBIS
occ_centroid_mah <- sample_data(n_occ = n_occ, prediction = pred,
                            prediction_layer = "Mahalanobis_trunc", sampling = "centroid",
                            method = "mahalanobis", seed = 1, strict = TRUE)

occ_random_mah <- sample_data(n_occ = n_occ, prediction = pred,
                          prediction_layer = "Mahalanobis_trunc", sampling = "random",
                          method = "mahalanobis", seed = 2, strict = TRUE)

occ_edge_mah <- sample_data(n_occ = n_occ, prediction = pred,
                        prediction_layer = "Mahalanobis_trunc", sampling = "edge",
                        method = "mahalanobis", seed = 3, strict = TRUE)

# Different Draws with same seed
head(occ_centroid_sui)
head(occ_centroid_mah)

# BIAS SURFACES ----
# One composite from the package survey proxies, one synthetic and extreme.
bias_path <- system.file("extdata", "ma_biases.tif", package = "nicheR")
if(bias_path == ""){ bias_path <- system.file("ma_biases.tif", package = "nicheR") }
biases <- rast(bias_path)

effort <- prepare_bias(biases, effect_direction = "direct", floor = 0.05)

w_effort_sui <- apply_bias(effort, pred, prediction_layer = "suitability_trunc")[[1]]

w_effort_mah <- apply_bias(effort, pred, effect_direction = "inverse",
                           prediction_layer = "Mahalanobis_trunc")[[1]]


# Temperature standardized across the cells the niche occupies. The 0.05 floor
# makes the gradient exactly twentyfold, coldest to warmest of those cells.
t_in <- terra::mask(env[["bio_1"]], inside_sui, maskvalues = c(0, NA))
t_mm <- terra::minmax(t_in)
steep_rast <- (t_in - t_mm[1]) / (t_mm[2] - t_mm[1])
steep <- prepare_bias(steep_rast, effect_direction = "direct", floor = 0.05)

w_steep_sui <- apply_bias(steep, pred, prediction_layer = "suitability_trunc")[[1]]

w_steep_mah <- apply_bias(steep, pred, prediction_layer = "Mahalanobis_trunc")[[1]]

# The weight surface is suitability times bias, so dividing suitability back out
# gives the bias composite on the niche's own cells. Numbers to quote in the text.
b <- terra::values(w_effort_sui) / terra::values(suit)
b <- b[is.finite(b) & b > 0]
cat("effort gradient over niche cells:", round(max(b) / min(b), 1), "\n")

b <- terra::values(w_steep_sui) / terra::values(suit)
b <- b[is.finite(b) & b > 0]
cat("steep gradient over niche cells:", round(max(b) / min(b), 1), "\n")


# BIASED DRAWS ----
# sample_biased_data() treats the supplied surface values as the weights directly.
occ_effort_sui <- sample_biased_data(n_occ = n_occ, prediction = w_effort_sui,
                                 seed = 5, strict = TRUE)
occ_steep_sui <- sample_biased_data(n_occ = n_occ, prediction = w_steep_sui,
                                seed = 6, strict = TRUE)

occ_effort_mah <- sample_biased_data(n_occ = n_occ, prediction = w_effort_mah,
                                 seed = 5, strict = TRUE)
occ_steep_mah <- sample_biased_data(n_occ = n_occ, prediction = w_steep_mah,
                                seed = 6, strict = TRUE)

# TABLE ----
# All suitable cells are the reference set: what the designs are drawing from.
pred_df_sui <- as.data.frame(pred[["suitability_trunc"]], xy = TRUE)
pred_df_sui <- pred_df_sui[pred_df_sui$suitability_trunc != 0, ]

pred_df_mah <- as.data.frame(pred[["Mahalanobis_trunc"]], xy = TRUE)
pred_df_mah <- pred_df_mah[pred_df_mah$Mahalanobis_trunc != 0, ]


draws <- list(
  "all suitable cells" = list(pred_df_sui, pred_df_mah),
  "centroid" = list(occ_centroid_sui,occ_centroid_mah),
  "random" = list(occ_random_sui, occ_random_mah),
  "edge" = list(occ_edge_sui, occ_edge_mah),
  "effort bias" = list(occ_effort_sui, occ_effort_mah),
  "steep bias" = list(occ_steep_sui, occ_steep_mah)
)

# One helper only: every draw needs the same columns off the prediction, and
# sample_data() returns coordinates rather than layer values.
get_cells <- function(occ){
  if(is.list(occ) && !is.data.frame(occ)){ occ <- occ[[1]] }
  occ <- as.data.frame(occ)
  j <- grep("^(x|y|lon|long|longitude|lat|latitude)$", names(occ), ignore.case = TRUE)[1:2]
  terra::extract(pred, as.matrix(occ[, j]))
}

# TABLE ----
methods <- c("suitability", "mahalanobis")

tab <- data.frame(
  design = "niche (true)",
  method = NA_character_,
  n = NA_integer_,
  mean_bio_1 = round(as.numeric(n2$centroid)[1], 2),
  mean_bio_12 = round(as.numeric(n2$centroid)[2], 0),
  mean_suit = 1,
  mean_mahal = 0
)

for(i in seq_along(draws)){
  for(j in seq_along(methods)){
    e <- get_cells(draws[[i]][[j]])
    tab <- rbind(tab, data.frame(
      design = names(draws)[i],
      method = methods[j],
      n = nrow(e),
      mean_bio_1 = round(mean(e$bio_1), 2),
      mean_bio_12 = round(mean(e$bio_12), 0),
      mean_suit = round(mean(e$suitability_trunc), 3),
      mean_mahal = round(mean(e$Mahalanobis_trunc), 3)
    ))
  }
}

# Shifts in niche SD units, so temperature and precipitation are comparable.
sds <- sqrt(diag(n2$cov_matrix))

ref <- c("suitability" = "centroid", "mahalanobis" = "edge")
tab$shift_bio_1 <- NA
tab$shift_bio_12 <- NA
for(i in which(grepl("bias", tab$design))){
  k <- which(tab$design == ref[[tab$method[i]]] & tab$method == tab$method[i])
  tab$shift_bio_1[i] <- round((tab$mean_bio_1[i] - tab$mean_bio_1[k]) / sds[1], 2)
  tab$shift_bio_12[i] <- round((tab$mean_bio_12[i] - tab$mean_bio_12[k]) / sds[2], 2)
}

tab

write.csv(tab, file.path(out_dir, "table_1.csv"), row.names = FALSE)


# FIGURE: BIAS ----
oldpar <- par(no.readonly = TRUE) # Save current settings
par(oldpar)                       # Restore original settings

xlab <- "Annual Mean Temperature"
ylab <- "Annual Precipitation"
bg <- as.data.frame(env, xy = FALSE)

n2_col <- "#1A1A1A"
n1_col <- "#0072B2"
n11_col <- "#E69F00"

par(mar = c(5.1, 4.1, 1, 2.1))
plot_ellipsoid(object = n2,
               background = bg,
               pch = 20, cex_bg = 0.6,
               col_bg = "lightgrey",
               xlab = xlab,
               ylab = ylab, lwd = 3, col_ell = "#1A1A1A")
add_ellipsoid(object = n1, lty = 2, col_ell = "#0072B2", lwd = 3)
add_ellipsoid(object = n11, lty = 3, col_ell = "#E69F00", lwd = 3)
add_data(data = as.data.frame(t(n2$centroid)),
         x = "bio_1", y = "bio_12", pts_col = "#1A1A1A", pch = 16, cex = 1.5)
add_data(data = as.data.frame(t(n1$centroid)),
         x = "bio_1", y = "bio_12", pts_col = "#0072B2", pch = 15, cex = 1.5)
add_data(data = as.data.frame(t(n11$centroid)),
         x = "bio_1", y = "bio_12", pts_col = "#E69F00", pch = 17, cex = 1.5)
text(x = 12, y = 2000, "C", col = n2_col, cex = 1.5)
text(x = 16, y = 3600, "A", col = n1_col, cex = 1.5)
text(x = 25, y = 3950, "B", col = n11_col,cex = 1.5)


par(oldpar)

pred2 <- predict(n2, env,  keep_data = FALSE,
                mahalanobis_truncated  = FALSE,
                suitability_truncated = TRUE,
                include_suitability = FALSE,
                include_mahalanobis = FALSE)

pred1 <- predict(n1, env,  keep_data = FALSE,
                mahalanobis_truncated  = FALSE,
                suitability_truncated = TRUE,
                include_suitability = FALSE,
                include_mahalanobis = FALSE)

pred11 <- predict(n11, env,  keep_data = FALSE,
                mahalanobis_truncated  = FALSE,
                suitability_truncated = TRUE,
                include_suitability = FALSE,
                include_mahalanobis = FALSE)

pred2_bin <- terra::subst(pred2 > 0, TRUE, 1)
pred1_bin <- terra::subst(pred1 > 0, TRUE, 1)
pred11_bin <- terra::subst(pred11 > 0, TRUE, 1)

terra::plot(pred1_bin,
            las = 1,
            mar = c(3.1, 2.1, 1.0, 1.0),
            col = c("lightgrey", n1_col),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE),
            plg = list(legend = c("Out", "In"), bty = "o",
                       x = -66.5, y = 29.4, y.intersp = 1.8,
                       text.width = strwidth("Out") * 1.6))
text(-96, 8, "A", cex = 2, col = n1_col)

terra::plot(pred11_bin,
            las = 1,
            mar = c(3.1, 2.1, 1.0, 1.0),
            col = c("lightgrey", n11_col),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE),
            plg = list(legend = c("Out", "In"), bty = "o",
                       x = -66.5, y = 29.4, y.intersp = 1.8,
                       text.width = strwidth("Out") * 1.8))
text(-96, 8, "B", cex = 2, col = n11_col)

terra::plot(pred2_bin,
            las = 1,
            mar = c(3.1, 2.1, 1.0, 1.0),
            col = c("lightgrey", n2_col),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE),
            plg = list(legend = c("Out", "In"), bty = "o",
                       x = -66.5, y = 29.4, y.intersp = 1.8,
                       text.width = strwidth("Out") * 1.8))
text(-96, 8, "C", cex = 2, col = n2_col)



terra::plot(env[[1]], main = xlab, font.main = 1,
            las = 1,
            mar = c(3.1, 2.1, 1.5, 1.0),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE))
text(-95, 8, "X-axis", cex = 1, col = n2_col)

terra::plot(env[[2]], main = ylab,
            font.main = 1,
            las = 1,
            mar = c(3.1, 2.1, 1.5, 1.0),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE))
text(-95, 8, "Y-axis", cex = 1, col = n2_col)

par(oldpar)


# BIAS

bg_mask <- terra::subst(env[[1]] > 0, TRUE, 1)

par(oldpar)

terra::plot(effort$composite_surface, main = "Samplig Effort Bias",
            font.main = 1,
            las = 1,
            mar = c(1.5, 1.0, 1.5, 1.0),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE))

terra::plot(bg_mask, col = "lightgrey", legend = FALSE,
            main = "Steep Temperature Bias",
            font.main = 1,
            las = 1,
            mar = c(1.5, 1.0, 1.5, 1.0),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE))
terra::plot(steep$composite_surface, add = TRUE)

par(oldpar)

bias_stack <- c(env, effort$composite_surface, steep$composite_surface)
names(bias_stack)[3:4] <- c("effort", "steep")
bias_stack_df <- as.data.frame(bias_stack)

par(mar = c(5.1, 4.1, 1, 1))
plot_ellipsoid(object = n2,
               prediction = bias_stack_df,
               col_layer = "effort",
               pch = 20, cex_bg = 0.7,
               col_bg = "lightgrey",
               xlab = xlab,
               ylab = ylab, lwd = 3, col_ell = "#1A1A1A")

plot_ellipsoid(object = n2,
               prediction = bias_stack_df,
               col_layer = "steep",
               pch = 20, cex_bg = 0.7,
               col_bg = "lightgrey",
               xlab = xlab,
               ylab = ylab, lwd = 3, col_ell = "#1A1A1A")
par(oldpar)




# Prediction and biased

pred <- predict(n2, env,  keep_data = TRUE,
                 mahalanobis_truncated  = TRUE,
                 suitability_truncated = TRUE,
                 include_suitability = FALSE,
                 include_mahalanobis = FALSE)
pred_df <- as.data.frame(pred)
head(pred_df)


pred_df <- as.data.frame(pred)
pred_effort_df <- as.data.frame(c(env, w_effort_sui))
pred_steep_df <- as.data.frame(c(env, w_steep_sui))



terra::plot(pred$suitability_trunc, main = "Suitability",
            font.main = 1,
            las = 1,
            mar = c(1.5, 1.0, 1.5, 1.0),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE))

terra::plot(w_effort_sui$suitability_trunc_biased_direct, main = "Effort",
            font.main = 1,
            las = 1,
            mar = c(1.5, 1.0, 1.5, 1.0),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE))

terra::plot(w_steep_sui$suitability_trunc_biased_direct, main = "Steep",
            font.main = 1,
            las = 1,
            mar = c(1.5, 1.0, 1.5, 1.0),
            pax = list(side = c(1,2), cex.axis = 1, retro = TRUE))

par(oldpar)
par(mar = c(5.1, 4.1, 1, 1))
plot_ellipsoid(object = n2,
               prediction = pred_df,
               col_layer = "suitability_trunc",
               pch = 20, cex_bg = 0.7,
               col_bg = "lightgrey",
               xlab = xlab,
               ylab = ylab, lwd = 3, col_ell = "#1A1A1A")

plot_ellipsoid(object = n2,
               prediction = pred_effort_df,
               col_layer = "suitability_trunc_biased_direct",
               pch = 20, cex_bg = 0.7,
               col_bg = "lightgrey",
               xlab = xlab,
               ylab = ylab, lwd = 3, col_ell = "#1A1A1A")

plot_ellipsoid(object = n2,
               prediction = pred_steep_df,
               col_layer = "suitability_trunc_biased_direct",
               pch = 20, cex_bg = 0.7,
               col_bg = "lightgrey",
               xlab = xlab,
               ylab = ylab, lwd = 3, col_ell = "#1A1A1A")
par(oldpar)






par(oldpar)
par(mar = c(5.1, 4.1, 1, 1))
plot_ellipsoid(object = n2,
               prediction = pred_df,
               col_layer = "Mahalanobis_trunc",
               pch = 20,
               col_bg = "lightgrey",
               xlab = xlab,
               ylab = ylab, lwd = 3, col_ell = "#1A1A1A")
par(oldpar)





# Each bias surface against the unbiased design it modifies, with the shift drawn
bg <- as.data.frame(env, xy = FALSE)
sds <- sqrt(diag(n2$cov_matrix))

ref <- c(suitability = "centroid", mahalanobis = "edge")
bias_rows <- c("effort bias", "steep bias")
bias_col <- c("effort bias" = "#1B9E77", "steep bias" = "#D95F02")

op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(4, 0, 2, 0))
for(i in seq_along(bias_rows)){
  for(j in seq_along(methods)){
    k <- bias_rows[i]
    r <- get_cells(draws[[ref[[methods[j]]]]][[j]])
    e <- get_cells(draws[[k]][[j]])

    plot_ellipsoid(n2, background = bg, col_ell = "black",
                   col_bg = "grey", pch = ".", main = "", xlab = "BIO1", ylab = "BIO12")
    add_data(data = r, x = "bio_1", y = "bio_12",
             pts_col = "grey45", pch = 1, cex = 0.7)
    add_data(data = e, x = "bio_1", y = "bio_12",
             pts_col = bias_col[[k]], pch = 16, cex = 0.7)

    # Arrow from the reference mean to the biased mean, labeled in niche SDs.
    m0 <- c(mean(r$bio_1), mean(r$bio_12))
    m1 <- c(mean(e$bio_1), mean(e$bio_12))
    points(m0[1], m0[2], pch = 20, col = "grey20", cex = 2)
    arrows(m0[1], m0[2], m1[1], m1[2], length = 0.08, lwd = 2, col = "black")
    shift <- round((m1 - m0) / sds, 2)
    mtext(paste0("shift in BIO1: ", shift[1], ", BIO12: ", shift[2]),
          side = 3, line = -1.2, adj = 0.97, cex = 0.7)

    if(i == 1){ mtext(methods[j], side = 3, line = 0.5, font = 2) }
  }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", horiz = TRUE, bty = "n", cex = 0.9,
       legend = c("unbiased reference", names(bias_col), "mean shift"),
       col = c("grey45", bias_col, "black"), pch = c(1, 16, 16, NA),
       lty = c(NA, NA, NA, 1), lwd = c(NA, NA, NA, 2))
par(op)



# SAVE ----
# The full specification of the simulated truth, so the section is reproducible.
save_nicheR(n1, file.path(out_dir, "niche_built.rds"))  # CHECK: argument names
save_nicheR(n2, file.path(out_dir, "niche_moved.rds"))
saveRDS(draws[-1], file.path(out_dir, "draws.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
