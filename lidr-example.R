library(lidR)
library(proceduralnames) # used in the source functions.
library(terra)
library(sf)
source("lidar-metrics.R")


# load an example dataset from the lidR package
LASfile <- system.file("extdata", "MixedConifer.laz", package = "lidR")
trees <- readLAS(LASfile, filter = "-drop_z_below 0")

plot(trees)

# This uses the `compute_pixel_metrics` which is a wrapper around
# `lidR::pixel_metrics`. It generates 40 rastersied metrics from the
# point cloud data.
gridded_metrics <- compute_pixel_metrics(trees, res = 2)
gm <- terra::rast(gridded_metrics)
print(names(gm))
plot(gm[[1:10]], col = hcl.colors(100, "viridis"))


# make some fake plots:
tree_crs <- st_crs(trees)
trees_bbox <- lidR::st_bbox(trees) |>
  sf::st_as_sfc()
trees_cent <- st_centroid(trees_bbox)
plot_max_ext <- (trees_bbox - trees_cent) * 0.75 + trees_cent
fake_plots <- sf::st_sample(x = plot_max_ext, size = 9, type = "regular") |>
  sf::st_buffer(5) |>
  sf::st_set_crs(tree_crs)

plot(gm[[2]], col = hcl.colors(100, "viridis"))
plot(fake_plots, add = TRUE, col = "#cc4040c9")

# now to extract metrics for these plots we use `compute`
plot_metrics_sf <- compute_plot_metrics(trees, geometry = fake_plots)

# or if you just want relative height metrics at 1 percentile intervals you
# can use:

rh_mets_sf <- lidR::plot_metrics(
  trees,
  func = rh_preds(Z),
  geometry = fake_plots
)
print(rh_mets_sf)
# the rh_preds function is just a reduced/modified version of the lidar_preds
# these examples should demonstrate how to build any custom functions that
# can be used with `lidR::plot_metrics` or `lidR::pixel_metrics`
