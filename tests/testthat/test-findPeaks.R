gaussian_2d_signal <- function(dt, rt, centers, heights, dt_sd, rt_sd) {
  int_mat <- matrix(0, nrow = length(dt), ncol = length(rt))
  for (i in seq_along(centers$dt)) {
    int_mat <- int_mat + outer(
      exp(-(dt - centers$dt[i])^2 / (2 * dt_sd^2)),
      exp(-(rt - centers$rt[i])^2 / (2 * rt_sd^2))
    ) * heights[i]
  }
  int_mat
}

# MassSpecWavelet::peakDetectionCWT() defaults to nearbyPeak = TRUE, which sets
# excludeBoundariesSize = 75 regardless of signal length/units. Both axes need
# to comfortably exceed 150 points, or no peak (wherever placed) can ever pass
# that filter. dt/rt peak-width ranges must also produce wavelet scales whose
# max exceeds peakDetectionCWT's default ridgeLength = 24, or real peaks get
# silently dropped.
dt_axis <- seq(0, 3, by = 0.01) # 301 points
rt_axis <- seq(0, 150, by = 0.5) # 301 points
common_args <- list(
  dt_length_ms = 0.07, rt_length_s = 3,
  dt_peakwidth_range_ms = c(0.05, 0.5),
  rt_peakwidth_range_s = c(2, 20)
)

test_that("findPeaksImpl detects a single well-separated 2D peak", {
  int_mat <- gaussian_2d_signal(
    dt_axis, rt_axis,
    centers = list(dt = 1, rt = 50), heights = 100, dt_sd = 0.05, rt_sd = 1.5
  )

  out <- rlang::exec(
    findPeaksImpl,
    drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
    !!!common_args
  )
  peaks <- out$peak_list

  expect_equal(nrow(peaks), 1L)
  expect_equal(peaks$dt_apex_ms, 1, tolerance = 0.01)
  expect_equal(peaks$rt_apex_s, 50, tolerance = 0.5)
  expect_equal(peaks$int_apex_au, 100, tolerance = 1e-6)
})

test_that("findPeaksImpl detects two well-separated 2D peaks at the right positions and heights", {
  int_mat <- gaussian_2d_signal(
    dt_axis, rt_axis,
    centers = list(dt = c(1, 2), rt = c(50, 100)),
    heights = c(100, 60), dt_sd = 0.05, rt_sd = 1.5
  )

  out <- rlang::exec(
    findPeaksImpl,
    drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
    !!!common_args
  )
  peaks <- out$peak_list

  expect_equal(nrow(peaks), 2L)
  ord <- order(peaks$rt_apex_s)
  expect_equal(peaks$dt_apex_ms[ord], c(1, 2), tolerance = 0.01)
  expect_equal(peaks$rt_apex_s[ord], c(50, 100), tolerance = 0.5)
  expect_equal(peaks$int_apex_au[ord], c(100, 60), tolerance = 1e-6)
})

test_that("findPeaksImpl merges two 2D peaks closer than the overlap threshold into one", {
  int_mat <- gaussian_2d_signal(
    dt_axis, rt_axis,
    centers = list(dt = c(1.5, 1.5), rt = c(70, 73)),
    heights = c(100, 80), dt_sd = 0.05, rt_sd = 1.5
  )

  out <- rlang::exec(
    findPeaksImpl,
    drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
    !!!common_args, iou_overlap_threshold = 0.2
  )

  expect_equal(nrow(out$peak_list), 1L)
})

test_that("findPeaksImpl returns an empty peak list for pure noise", {
  set.seed(42)
  int_mat <- matrix(
    rnorm(length(dt_axis) * length(rt_axis), mean = 100, sd = 0.5),
    nrow = length(dt_axis), ncol = length(rt_axis)
  )

  out <- rlang::exec(
    findPeaksImpl,
    drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
    !!!common_args
  )

  expect_equal(nrow(out$peak_list), 0L)
})

test_that("findPeaksImpl returns an empty peak list for a flat/constant matrix", {
  int_mat <- matrix(5, nrow = length(dt_axis), ncol = length(rt_axis))

  out <- rlang::exec(
    findPeaksImpl,
    drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
    !!!common_args
  )

  expect_equal(nrow(out$peak_list), 0L)
})

test_that("findPeaksImpl's peak_list has the documented schema", {
  int_mat <- gaussian_2d_signal(
    dt_axis, rt_axis,
    centers = list(dt = 1, rt = 50), heights = 100, dt_sd = 0.05, rt_sd = 1.5
  )

  out <- rlang::exec(
    findPeaksImpl,
    drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
    !!!common_args
  )
  peaks <- out$peak_list

  expect_named(
    peaks,
    c(
      "PeakID", "dt_apex_ms", "rt_apex_s", "int_apex_au", "ddt_apex_au", "drt_apex_au",
      "dt_min_ms", "dt_max_ms", "rt_min_s", "rt_max_s", "dt_cm_ms", "rt_cm_s",
      "dt_apex_idx", "rt_apex_idx", "dt_min_idx", "dt_max_idx", "rt_min_idx", "rt_max_idx",
      "dt_cm_idx", "rt_cm_idx"
    )
  )
  expect_true(all(peaks$dt_min_idx <= peaks$dt_apex_idx & peaks$dt_apex_idx <= peaks$dt_max_idx))
  expect_true(all(peaks$rt_min_idx <= peaks$rt_apex_idx & peaks$rt_apex_idx <= peaks$rt_max_idx))
})

test_that("findPeaksImpl with debug_idx also returns debug_info", {
  int_mat <- gaussian_2d_signal(
    dt_axis, rt_axis,
    centers = list(dt = 1, rt = 50), heights = 100, dt_sd = 0.05, rt_sd = 1.5
  )

  out_nodebug <- rlang::exec(
    findPeaksImpl,
    drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
    !!!common_args
  )
  out_debug <- rlang::exec(
    findPeaksImpl,
    drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
    !!!common_args, debug_idx = list(dt = 1, rt = 1)
  )

  expect_null(out_nodebug$debug_info)
  expect_named(out_debug$debug_info, c("rt", "dt"))
})

test_that("findPeaksImpl with verbose = TRUE reports the scales used", {
  int_mat <- gaussian_2d_signal(
    dt_axis, rt_axis,
    centers = list(dt = 1, rt = 50), heights = 100, dt_sd = 0.05, rt_sd = 1.5
  )

  expect_message(
    rlang::exec(
      findPeaksImpl,
      drift_time = dt_axis, retention_time = rt_axis, int_mat = int_mat,
      !!!common_args, verbose = TRUE
    ),
    "Using the following scales"
  )
})

test_that("intersectionOverUnion (2D) computes the expected ratios", {
  expect_equal(
    intersectionOverUnion(
      list(dt_idx_min = 0, dt_idx_max = 10, rt_idx_min = 0, rt_idx_max = 10),
      list(dt_idx_min = 20, dt_idx_max = 30, rt_idx_min = 0, rt_idx_max = 10)
    ),
    0
  )
  # Partial overlap: intersection 5x5=25, union 100+100-25=175
  expect_equal(
    intersectionOverUnion(
      list(dt_idx_min = 0, dt_idx_max = 10, rt_idx_min = 0, rt_idx_max = 10),
      list(dt_idx_min = 5, dt_idx_max = 15, rt_idx_min = 5, rt_idx_max = 15)
    ),
    25 / 175
  )
  # Full containment: intersection 10x10=100, union 400+100-100=400
  expect_equal(
    intersectionOverUnion(
      list(dt_idx_min = 0, dt_idx_max = 20, rt_idx_min = 0, rt_idx_max = 20),
      list(dt_idx_min = 5, dt_idx_max = 15, rt_idx_min = 5, rt_idx_max = 15)
    ),
    0.25
  )
  # Zero-area ROI gives an undetermined (zero) iou:
  expect_equal(
    intersectionOverUnion(
      list(dt_idx_min = 5, dt_idx_max = 5, rt_idx_min = 0, rt_idx_max = 10),
      list(dt_idx_min = 0, dt_idx_max = 10, rt_idx_min = 0, rt_idx_max = 10)
    ),
    0
  )
})

test_that("compute_center_of_mass (2D) weights the center of mass by intensity", {
  patch_mat <- matrix(0, nrow = 5, ncol = 5)
  patch_mat[3, 3] <- 10 # single, symmetric spike at the center
  rois <- tibble::tibble(dt_idx_min = 1, dt_idx_max = 5, rt_idx_min = 1, rt_idx_max = 5)

  out <- compute_center_of_mass(rois, patch_mat)

  expect_equal(out$dt_cm_idx, 3L)
  expect_equal(out$rt_cm_idx, 3L)
})

test_that("merge_overlapping_rois (2D) merges ROIs above the overlap threshold, keeping the taller apex", {
  rois <- tibble::tibble(
    dt_idx_apex = c(5, 6), rt_idx_apex = c(5, 6),
    dt_idx_min = c(0, 3), dt_idx_max = c(10, 12),
    rt_idx_min = c(0, 3), rt_idx_max = c(10, 12),
    int_apex_au = c(50, 80)
  )

  merged <- merge_overlapping_rois(rois, int_mat = matrix(0, 1, 1), iou_overlap_threshold = 0.2)

  expect_equal(nrow(merged), 1L)
  expect_equal(merged$dt_idx_apex, 6)
  expect_equal(merged$rt_idx_apex, 6)
  expect_equal(merged$dt_idx_min, 0)
  expect_equal(merged$dt_idx_max, 12)
  expect_equal(merged$int_apex_au, 80)
})

test_that("merge_overlapping_rois (2D) leaves non-overlapping ROIs untouched", {
  rois <- tibble::tibble(
    dt_idx_apex = c(5, 50), rt_idx_apex = c(5, 50),
    dt_idx_min = c(0, 45), dt_idx_max = c(10, 55),
    rt_idx_min = c(0, 45), rt_idx_max = c(10, 55),
    int_apex_au = c(50, 80)
  )

  merged <- merge_overlapping_rois(rois, int_mat = matrix(0, 1, 1), iou_overlap_threshold = 0.2)

  expect_equal(nrow(merged), 2L)
})

test_that("merge_overlapping_rois (2D) skips merging when iou_overlap_threshold > 1", {
  rois <- tibble::tibble(
    dt_idx_apex = c(5, 6), rt_idx_apex = c(5, 6),
    dt_idx_min = c(0, 3), dt_idx_max = c(10, 12),
    rt_idx_min = c(0, 3), rt_idx_max = c(10, 12),
    int_apex_au = c(50, 80)
  )

  merged <- merge_overlapping_rois(rois, int_mat = matrix(0, 1, 1), iou_overlap_threshold = 1.1)

  expect_equal(nrow(merged), 2L)
})

test_that("findZeroCrossings finds sign changes in the requested direction", {
  x <- c(-2, -1, 1, 2, 1, -1, -2)

  expect_equal(findZeroCrossings(x, "both"), c(3, 5))
  expect_equal(findZeroCrossings(x, "up"), 3)
  expect_equal(findZeroCrossings(x, "down"), 5)
})

test_that("interleave_peaks_with_zeros pairs each peak with its surrounding zero crossings", {
  peak_idx <- c(5, 15)
  zero_idx <- c(2, 8, 12, 18)

  out <- interleave_peaks_with_zeros(peak_idx, zero_idx, signal_idx = 3L)

  expect_equal(out$signal_idx, c(3L, 3L))
  expect_equal(out$idx_apex, c(5, 15))
  expect_equal(out$idx_min, c(2, 12))
  expect_equal(out$idx_max, c(8, 18))
})

test_that("interleave_peaks_with_zeros returns an empty tibble when there are no peaks or no zeros", {
  out <- interleave_peaks_with_zeros(integer(0), c(1, 2, 3), signal_idx = 1L)

  expect_equal(nrow(out), 0L)
  expect_named(out, c("signal_idx", "idx_apex", "idx_min", "idx_max"))
})
