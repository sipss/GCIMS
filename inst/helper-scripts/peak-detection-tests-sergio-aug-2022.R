library(BiocParallel)
library(ggplot2)
library(GCIMS)


show_progress_bar <- interactive() && is.null(getOption("knitr.in.progress"))
# disable parallellization: (Useful for better error reporting)
register(SerialParam(progressbar = show_progress_bar), default = TRUE)

saveRDS(ket1, "Ketones1_for_peak_detection.rds")
ket1 <- readRDS("Ketones1_for_peak_detection.rds")

dt <- dtime(ket1)
rt <- rtime(ket1)
int_mat <- intensity(ket1)

dt_length_pts <- GCIMS:::units_to_points(length_phys = 0.14, step_phys = dt[2] - dt[1], must_odd = TRUE)
rt_length_pts <- GCIMS:::units_to_points(length_phys = 3, step_phys = rt[2] - rt[1], must_odd = TRUE)


peak_list_ket1 <- peak_detection(
    drift_time = dt,
    retention_time = rt,
    int_mat = int_mat,
    noise_level = 2,
    verbose = TRUE,
    dt_length_pts = dt_length_pts,
    rt_length_pts = rt_length_pts,
    iou_overlap_threshold = 0.2
)

plt <- plotRaw(ket1)
plt <- add_peaklist_rect(plt, peak_list_ket1, color_by = "PeakID")
plt


## Step by step

the_rip <- find_rip(int_mat, verbose = TRUE, retention_time = rt, drift_time = dt)
dt_minpeakdistance_pts <- estimate_minpeakdistance(the_rip$rip, verbose = TRUE, drift_time = dt)


# Compute the 2nd derivative for both axes
deriv2 <- compute_second_deriv(
  int_mat,
  dt_length_pts = dt_length_pts,
  rt_length_pts = rt_length_pts,
  dt_order = 2L,
  rt_order = 2L
)
drt <- -deriv2$drt
ddt <- -deriv2$ddt
daux <- drt + ddt


minmax <- range(intmat)

cubic_root_trans <- scales::trans_new(
  name = "cubic_root",
  transform = function(x) sign(x)*abs(x)^(1/3),
  inverse = function(x) sign(x)*abs(x)^3,
  breaks = function(x, n = 5) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(numeric())
    }
    rng <- range(x)
    rng <- sign(rng)*abs(rng)^(1/3)
    out <- labeling::extended(rng[1], rng[2], n)
    out <- sign(out)*abs(out)^3
    out
  }
)

mat_to_plot <- int_mat
mat_to_plot_trans <- cubic_root_trans$transform(mat_to_plot)
nr <- mat_to_nativeRaster(mat_to_plot_trans, COLORMAP_VIRIDIS_256_A_m1)
idx <- dt_rt_range_normalization(dt, rt, dt_range = NULL, rt_range = NULL)
minmax <- range(mat_to_plot)

# The geom_rect is fake and it is only used to force the fill legend to appear
gplt <- ggplot() +
  geom_rect(
    xmin = idx$dt_ms_min, xmax = idx$dt_ms_min,
    ymin = idx$rt_s_min, ymax = idx$rt_s_min,
    ggplot2::aes(fill = .data$x),
    data = data.frame(x = NA_real_)
  ) +
  ggplot2::annotation_raster(
    nr,
    xmin = idx$dt_ms_min, xmax = idx$dt_ms_max,
    ymin = idx$rt_s_min, ymax = idx$rt_s_max
  ) +
  ggplot2::scale_fill_viridis_c( # This has to match with the COLORMAP above
    direction = -1,
    option = "A",
    limits = minmax,
    na.value = "#00000000",
    trans = cubic_root_trans
  ) +
  ggplot2::lims(
    x = c(idx$dt_ms_min, idx$dt_ms_max),
    y = c(idx$rt_s_min, idx$rt_s_max)
  ) +
  ggplot2::labs(
    x = "Drift time (ms)",
    y = "Retention time (s)",
    fill = "Intensity (a.u.)"
  ) +
  ggplot2::theme_minimal()
gplt

region <- find_region_without_peaks(int_mat, half_min_size = c(10, 10), noise_quantile = 0.25)
sigmaNoise_drt <- stats::sd(drt[region$row_min:region$row_max, region$col_min:region$col_max])
sigmaNoise_ddt <- stats::sd(ddt[region$row_min:region$row_max, region$col_min:region$col_max])
sigmaNoise_daux <- stats::sd(daux[region$row_min:region$row_max, region$col_min:region$col_max])
noise_level <- 4

peaks_zeros <- find_all_peaks_and_zero_crossings(
  ddt, drt,
  drift_time = dt,
  retention_time = rt,
  min_peak_height_ddt = noise_level*sigmaNoise_ddt,
  min_peak_height_drt = noise_level*sigmaNoise_drt,
  dt_minpeakdistance_pts = dt_minpeakdistance_pts,
  rt_minpeakdistance_pts = 1,
  verbose = FALSE
)

if (0) {
  pts <- purrr::imap_dfr(
    peaks_zeros$peaksrt,
    function(x, rt_apex_idx) data.frame(rt_apex_idx = rt_apex_idx, dt_apex_idx = x[,1])
  ) %>%
    dplyr::mutate(
      rt_apex_s = rt[rt_apex_idx],
      dt_apex_ms = dt[dt_apex_idx]
    )

  pts2 <- purrr::imap_dfr(
    peaks_zeros$peaksdt,
    function(x, dt_apex_idx) {
      x <- x[,1]
      data.frame(dt_apex_idx = rep(dt_apex_idx, length(x)), rt_apex_idx = x)
    }
  ) %>%
    dplyr::mutate(
      rt_apex_s = rt[rt_apex_idx],
      dt_apex_ms = dt[dt_apex_idx]
    )



  gplt +
    geom_point(data = pts, mapping = aes(x = dt_apex_ms, y = rt_apex_s), color = "green", shape = "|") +
    geom_point(data = pts2, mapping = aes(x = dt_apex_ms, y = rt_apex_s), color = "green", shape = "-") +
    lims(x = c(13, 14.5), y = c( 200, 400))
}


rois <- intersect_peaks_in_both_directions(peaks_zeros$peaksrt, peaks_zeros$zeros_rt, peaks_zeros$peaksdt, peaks_zeros$zeros_dt, the_rip, spec_length = nrow(int_mat), num_spec = ncol(int_mat))

# Merging algorithm:
ROIs_overlap <- merge_overlapping_rois(rois, int_mat, iou_overlap_threshold = 0.9)

roi_center_of_mass <- compute_center_of_mass(ROIs = ROIs_overlap, int_mat = int_mat)
# Convert into a data frame / peak list:
peak_list <- rois_to_peaklist(
  ROIs = ROIs_overlap,
  roi_center_of_mass = roi_center_of_mass,
  drift_time = dt,
  retention_time = rt
)

GCIMS:::add_peaklist_rect(gplt, peak_list)
