gaussian_signal <- function(x, centers, heights, sd) {
  y <- rep(0, length(x))
  for (i in seq_along(centers)) {
    y <- y + heights[i] * exp(-(x - centers[i])^2 / (2 * sd^2))
  }
  y
}

test_that("findPeaksImpl1D detects a single well-separated peak", {
  x <- seq(0, 20, by = 0.01)
  y <- gaussian_signal(x, centers = 10, heights = 100, sd = 0.3)

  out <- findPeaksImpl1D(x, y, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))
  peaks <- out$peak_list

  expect_equal(nrow(peaks), 1L)
  expect_equal(peaks$apex, 10, tolerance = 0.02)
  expect_equal(peaks$int_apex_au, 100, tolerance = 1e-6)
  expect_true(peaks$min < peaks$apex && peaks$apex < peaks$max)
})

test_that("findPeaksImpl1D detects two well-separated peaks at the right positions and heights", {
  x <- seq(0, 20, by = 0.01)
  y <- gaussian_signal(x, centers = c(5, 15), heights = c(100, 60), sd = 0.3)

  out <- findPeaksImpl1D(x, y, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))
  peaks <- out$peak_list

  expect_equal(nrow(peaks), 2L)
  expect_equal(sort(peaks$apex), c(5, 15), tolerance = 0.02)
  expect_equal(peaks$int_apex_au[order(peaks$apex)], c(100, 60), tolerance = 1e-6)
  # Regions of interest should not overlap:
  expect_true(peaks$max[peaks$apex == min(peaks$apex)] <= peaks$min[peaks$apex == max(peaks$apex)])
})

test_that("findPeaksImpl1D merges two peaks closer than the overlap threshold into one", {
  x <- seq(0, 20, by = 0.01)
  y <- gaussian_signal(x, centers = c(10, 10.4), heights = c(100, 80), sd = 0.3)

  out <- findPeaksImpl1D(
    x, y,
    length_in_xunits = 0.07,
    peakwidth_range_xunits = c(0.3, 3),
    iou_overlap_threshold = 0.2
  )

  expect_equal(nrow(out$peak_list), 1L)
})

test_that("findPeaksImpl1D returns an empty peak list for pure noise", {
  set.seed(42)
  x <- seq(0, 20, by = 0.01)
  y <- rnorm(length(x), mean = 0, sd = 0.01)

  out <- findPeaksImpl1D(x, y, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))

  expect_equal(nrow(out$peak_list), 0L)
})

test_that("findPeaksImpl1D returns an empty peak list for a flat/constant signal", {
  x <- seq(0, 20, by = 0.01)
  y <- rep(5, length(x))

  out <- findPeaksImpl1D(x, y, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))

  expect_equal(nrow(out$peak_list), 0L)
})

test_that("findPeaksImpl1D's peak_list has the documented schema", {
  x <- seq(0, 20, by = 0.01)
  y <- gaussian_signal(x, centers = 10, heights = 100, sd = 0.3)

  out <- findPeaksImpl1D(x, y, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))

  expect_named(
    out$peak_list,
    c(
      "PeakID", "apex", "int_apex_au", "deriv_apex_au", "min", "max", "cm",
      "apex_idx", "min_idx", "max_idx", "cm_idx"
    )
  )
  expect_true(is.character(out$peak_list$PeakID))
  expect_true(all(out$peak_list$min_idx <= out$peak_list$apex_idx))
  expect_true(all(out$peak_list$apex_idx <= out$peak_list$max_idx))
})

test_that("findPeaksImpl1D with debug = TRUE also returns debug_info", {
  x <- seq(0, 20, by = 0.01)
  y <- gaussian_signal(x, centers = 10, heights = 100, sd = 0.3)

  out_nodebug <- findPeaksImpl1D(x, y, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))
  out_debug <- findPeaksImpl1D(
    x, y,
    length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3), debug = TRUE
  )

  expect_null(out_nodebug$debug_info)
  expect_false(is.null(out_debug$debug_info))
})

test_that("findPeaksImpl1D with verbose = TRUE reports the scales used", {
  x <- seq(0, 20, by = 0.01)
  y <- gaussian_signal(x, centers = 10, heights = 100, sd = 0.3)

  expect_message(
    findPeaksImpl1D(
      x, y,
      length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3), verbose = TRUE
    ),
    "Using the following scales"
  )
})

test_that("intersectionOverUnion1D computes the expected ratios", {
  # No overlap:
  expect_equal(
    intersectionOverUnion1D(list(idx_min = 0, idx_max = 10), list(idx_min = 20, idx_max = 30)),
    0
  )
  # Touching boundaries count as no overlap:
  expect_equal(
    intersectionOverUnion1D(list(idx_min = 0, idx_max = 10), list(idx_min = 10, idx_max = 20)),
    0
  )
  # Partial overlap: intersection 5, union 15
  expect_equal(
    intersectionOverUnion1D(list(idx_min = 0, idx_max = 10), list(idx_min = 5, idx_max = 15)),
    1 / 3
  )
  # Full containment: intersection 10, union 20
  expect_equal(
    intersectionOverUnion1D(list(idx_min = 0, idx_max = 20), list(idx_min = 5, idx_max = 15)),
    0.5
  )
  # Zero-area ROI gives an undetermined (zero) iou:
  expect_equal(
    intersectionOverUnion1D(list(idx_min = 5, idx_max = 5), list(idx_min = 0, idx_max = 10)),
    0
  )
})

test_that("compute_center_of_mass_1D weights the center of mass by intensity", {
  # A symmetric patch centers exactly on its middle point:
  y_sym <- c(1, 2, 3, 2, 1)
  rois_sym <- tibble::tibble(idx_min = 1, idx_max = 5)
  expect_equal(compute_center_of_mass_1D(rois_sym, y_sym)$cm_idx, 3L)

  # An asymmetric patch shifts the center of mass towards the heavier side:
  y_asym <- c(1, 1, 5, 1, 1, 1, 1)
  rois_asym <- tibble::tibble(idx_min = 1, idx_max = 7)
  expect_equal(compute_center_of_mass_1D(rois_asym, y_asym)$cm_idx, 4L)
})

test_that("merge_overlapping_rois_1D merges ROIs above the overlap threshold, keeping the taller apex", {
  rois <- tibble::tibble(
    idx_apex = c(5, 12),
    idx_min = c(0, 5),
    idx_max = c(10, 15),
    int_apex_au = c(50, 80)
  )

  merged <- merge_overlapping_rois_1D(rois, y = numeric(0), iou_overlap_threshold = 0.2)

  expect_equal(nrow(merged), 1L)
  expect_equal(merged$idx_apex, 12)
  expect_equal(merged$idx_min, 0)
  expect_equal(merged$idx_max, 15)
  expect_equal(merged$int_apex_au, 80)
})

test_that("merge_overlapping_rois_1D leaves non-overlapping ROIs untouched", {
  rois <- tibble::tibble(
    idx_apex = c(5, 25),
    idx_min = c(0, 20),
    idx_max = c(10, 30),
    int_apex_au = c(50, 80)
  )

  merged <- merge_overlapping_rois_1D(rois, y = numeric(0), iou_overlap_threshold = 0.2)

  expect_equal(nrow(merged), 2L)
  expect_equal(merged$idx_apex, c(5, 25))
})

test_that("merge_overlapping_rois_1D skips merging when iou_overlap_threshold > 1", {
  rois <- tibble::tibble(
    idx_apex = c(5, 12),
    idx_min = c(0, 5),
    idx_max = c(10, 15),
    int_apex_au = c(50, 80)
  )

  merged <- merge_overlapping_rois_1D(rois, y = numeric(0), iou_overlap_threshold = 1.1)

  expect_equal(nrow(merged), 2L)
})
