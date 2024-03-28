density_preds <- function(Z, percentiles = seq(0.1, 0.9, 0.1)) {
  z_range <- range(Z)
  z_range_diff <- z_range[[2]] - z_range[[1]]

  lapply(
    z_range_diff * percentiles,
    \(h) mean(Z > (h + z_range[[1]]))
  ) |>
    stats::setNames(paste0("d", (percentiles * 100))) |>
    list2DF()
}

height_preds <- function(Z, percentiles = c(seq(0.1, 0.9, 0.1), 0.95, 0.99)) {
  stats::quantile(Z, percentiles) |>
    t() |>
    stats::setNames(paste0("h", (percentiles * 100))) |>
    as.data.frame()
}

l_moment_preds <- function(Z) {
  if (length(unique(Z)) == 1) {
    out <- data.frame(
      L2 = 0,
      L3 = 0,
      L4 = 0,
      L_cv = 0,
      L_skew = 0,
      L_kurt = 0
    )
  } else {
    l_moments <- lmomco::lmom.ub(Z)
    out <- data.frame(
      L2 = l_moments$L2,
      L3 = l_moments$L3,
      L4 = l_moments$L4,
      L_cv = l_moments$LCV,
      L_skew = l_moments$TAU3,
      L_kurt = l_moments$TAU4
    )
  }
  out
}

#' Standard LiDAR metrics
#'
#' This function computes a set of metrics from LiDAR point clouds. For metric
#' definitions and references, see Table 1 of Mahoney et al (2022).
#'
#' @param Z Return height for each point
#' @param ReturnNumber Return number of each point
#' @param min,max The minimum and maximum valid Z value. Returns outside of
#' this range will be discarded before metric computation.
#'
#' @references
#' Michael J Mahoney, Lucas K Johnson, Eddie Bevilacqua & Colin M Beier
#' (2022) Filtering ground noise from LiDAR returns produces inferior models of
#' forest aboveground biomass in heterogenous landscapes, GIScience & Remote
#' Sensing, 59:1, 1266-1280, DOI: 10.1080/15481603.2022.2103069
#'
#' @return A `data.frame` with 40 columns, containing various LiDAR metrics.
#'
#' @examples
#' lidar_preds(rnorm(1000, 3), sample(1:3, 1000, TRUE, c(0.7, 0.2, 0.1)))
#'
#' @export
lidar_preds <- function(Z, ReturnNumber, min = 0, max = Inf) {
  density_pred_template <- lapply(
    paste0("d", seq(0.1, 0.9, 0.1) * 100),
    \(x) stats::setNames(data.frame(NA_real_), x)
  ) |>
    do.call(what = cbind)

  height_pred_template <- lapply(
    paste0("h", c(seq(0.1, 0.9, 0.1), 0.95, 0.99) * 100),
    \(x) stats::setNames(data.frame(NA_real_), x)
  ) |>
    do.call(what = cbind)

  out <- cbind(
    n = NA_integer_,
    zmean = NA_real_,
    max = NA_real_,
    min = NA_real_,
    quad_mean = NA_real_,
    cv = NA_real_,
    z_kurt = NA_real_,
    z_skew = NA_real_,
    L2 = NA_real_,
    L3 = NA_real_,
    L4 = NA_real_,
    L_cv = NA_real_,
    L_skew = NA_real_,
    L_kurt = NA_real_,
    height_pred_template,
    density_pred_template,
    cancov = NA_real_,
    quad_mean_c = NA_real_,
    zmean_c = NA_real_,
    cv_c = NA_real_,
    hvol = NA_real_,
    rpc1 = NA_real_
  )

  include <- Z >= min & Z <= max

  if (sum(include) < 10) {
    return(out)
  }

  Z <- Z[include]
  ReturnNumber <- ReturnNumber[include]

  out[["n"]] <- length(Z)
  out[["zmean"]] <- mean(Z)
  out[["max"]] <- max(Z)
  out[["min"]] <- min(Z)
  out[["quad_mean"]] <- sqrt(mean(Z^2))
  out[["cv"]] <- stats::sd(Z) / mean(Z)
  out[["z_kurt"]] <- e1071::kurtosis(Z, type = 2)
  out[["z_skew"]] <- e1071::skewness(Z, type = 2)
  l_moments <- l_moment_preds(Z)
  out[names(l_moments)] <- l_moments
  out[names(height_pred_template)] <- height_preds(Z)
  out[names(density_pred_template)] <- density_preds(Z)
  out[["cancov"]] <- mean(Z > 2)

  if (sum(Z > 2.5) > 1) {
    out[["quad_mean_c"]] <- sqrt(mean((Z[Z > 2.5])^2))
    out[["zmean_c"]] <- mean(Z[Z > 2.5])
    out[["cv_c"]] <- stats::sd(Z[Z > 2.5]) / out$zmean_c
  } else {
    out[["quad_mean_c"]] <- 0
    out[["zmean_c"]] <- 0
    out[["cv_c"]] <- 0
  }

  out$hvol <- out$cancov * out$zmean
  out$rpc1 <- mean(ReturnNumber == 1)
  out
}

#' Relative height metrics function
#' For use with lidR::plot_metrics or lidR::pixel_metrics functions. Assumes
#' the point cloud has been normalised to ground level.
#' @param Z Return height for each point
#' @param min,max The minimum and maximum valid Z value. Returns outside of
#' this range will be discarded before metric computation.
#' @return A `data.frame` with 100 columns, containing relative height metrics
#' at 1 percentile intervals.
#' @export
rh_preds <- function(Z, min = 0, max = Inf) {
  out <- lapply(
    paste0("rh", seq(0.01, 1, 0.01) * 100),
    \(x) stats::setNames(data.frame(NA_real_), x)
  ) |>
    do.call(what = cbind)
  include <- Z >= min & Z <= max

  if (sum(include) < 10) {
    return(out)
  }

  Z <- Z[include]
  Z <- Z[include]
  out[names(out)] <- height_preds(Z, percentiles = seq(0.01, 1, 0.01))
  return(out)
}

#' Compute pixel-level metrics from LiDAR data
#'
#' These are incredibly thin wrappers around [lidR::pixel_metrics()]
#' and [lidR::plot_metrics()] with different default argument values.
#'
#' @inheritParams lidR::pixel_metrics
#' @inheritParams lidR::plot_metrics
#' @param output_filename randomly generated filename for the output raster.
#' Where las is a LAScatalog, the processing is carried out in a temporary
#' directory and the output vrt is resaved to the output_filename in order to
#' preserve file naming and reduce clutter.
#'
#' @returns The name of the output file.
#'
#' @examplesIf rlang::is_installed("lidR")
#' system.file("extdata/MixedConifer.laz", package = "lidR") |>
#'   lidR::readLAS() |>
#'   compute_pixel_metrics(output_filename = tempfile(fileext = ".tif"))
#'
#' @name compute_pixel_metrics
#' @export
compute_pixel_metrics <- function(las,
                                  func = lidar_preds(Z, ReturnNumber),
                                  res = 30,
                                  start = c(0, 0),
                                  ...,
                                  output_filename = random_tif(),
                                  overwrite = TRUE) {
  UseMethod("compute_pixel_metrics")
}

#' @exportS3Method
compute_pixel_metrics.LAScatalog <- function(
    las,
    func = ~ lidar_preds(Z, ReturnNumber),
    res = 30,
    start = c(0, 0),
    ...,
    output_filename = random_tif(),
    overwrite = TRUE) {
  temp_dir <- file.path(
    tempdir(),
    proceduralnames::make_english_names(1),
    "pixel_metrics_{ID}"
  )

  lidR::opt_output_files(las) <- temp_dir
  lidR::opt_merge(las) <- TRUE
  pm <- lidR::pixel_metrics(las, func,
    res = res, start = start, ...
  )

  terra::writeRaster(pm, output_filename, overwrite = overwrite)
  return(output_filename)
}

#' @exportS3Method
compute_pixel_metrics.LAS <- function(las,
                                      func = lidar_preds(Z, ReturnNumber),
                                      res = 30,
                                      start = c(0, 0),
                                      ...,
                                      output_filename = random_tif(),
                                      overwrite = TRUE) {
  if (output_filename == tools::file_path_sans_ext(output_filename)) {
    output_filename <- paste0(output_filename, ".tif")
  }
  out <- lidR::pixel_metrics(las, func, res = res, start = start, ...)
  terra::writeRaster(out, output_filename, overwrite = overwrite)
  return(output_filename)
}


#' @rdname compute_pixel_metrics
#' @export
compute_plot_metrics <- function(
    las,
    func = lidar_preds(Z, ReturnNumber),
    geometry,
    ...,
    radius) {
  lidR::plot_metrics(las, func, geometry = geometry, ..., radius = radius)
}


random_tif <- function() {
  paste0(proceduralnames::make_english_names(1), ".tif")
}
