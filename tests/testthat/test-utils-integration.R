gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

test_that("find_half_max_boundaries locates the half-max crossings of a single peak", {
  x <- c(0, 1, 2, 4, 8, 10, 8, 4, 2, 1, 0)

  out <- find_half_max_boundaries(x)

  expect_equal(out$apex_idx, 6L)
  expect_equal(out$left_bound_idx, 5L)
  expect_equal(out$right_bound_idx, 7L)
  expect_equal(out$half_width_left, 1L)
  expect_equal(out$half_width_right, 1L)
  expect_equal(out$fwhm, 2L)
  expect_equal(out$fwhm_sym, 2)
})

test_that("find_half_max_boundaries returns NaN fwhm_sym for a flat vector with no real peak", {
  out <- find_half_max_boundaries(rep(5, 50))

  expect_true(is.na(out$fwhm_sym))
})

test_that("find_rip locates an isolated RIP peak and clips its boundaries to the FWHM rule", {
  dt_idx <- 1:200
  tis_shape <- gauss(dt_idx, center = 100, height = 1000, sd = 10)
  int_mat <- matrix(0, nrow = 200, ncol = 20)
  rt_apex_col <- 12
  for (j in 1:20) {
    int_mat[, j] <- tis_shape * if (j == rt_apex_col) 1 else 0.6
  }

  out <- find_rip(int_mat)

  # No internal valleys exist in an isolated Gaussian, so the boundaries come
  # from clipping at the apex +/- 3*fwhm_sym of the TIS profile itself:
  fwhm_sym <- find_half_max_boundaries(rowSums(int_mat))$fwhm_sym
  expect_equal(out$dt_idx_apex, 100L)
  expect_equal(out$rt_idx_apex, rt_apex_col)
  expect_equal(out$dt_idx_start, 100 - 3 * fwhm_sym)
  expect_equal(out$dt_idx_end, 100 + 3 * fwhm_sym)
  expect_equal(length(out$rip), out$dt_idx_end - out$dt_idx_start + 1)
})

test_that("find_rip prefers a nearby valley over the FWHM clip when the valley is tighter", {
  dt_idx <- 1:200
  # Two small satellite peaks close enough to the main RIP to create real
  # valleys well inside the +/- 3*fwhm_sym clip radius:
  tis_shape <- gauss(dt_idx, center = 100, height = 1000, sd = 10) +
    gauss(dt_idx, center = 130, height = 200, sd = 5) +
    gauss(dt_idx, center = 65, height = 150, sd = 5)
  int_mat <- matrix(tis_shape, nrow = 200, ncol = 5)

  out <- find_rip(int_mat)

  expect_equal(out$dt_idx_apex, 100L)
  expect_equal(out$dt_idx_start, 74)
  expect_equal(out$dt_idx_end, 122)
})

test_that("find_rip aborts with a clear message when the RIP width can't be located", {
  int_mat_flat <- matrix(5, nrow = 50, ncol = 10)

  expect_error(
    suppressWarnings(find_rip(int_mat_flat)),
    "could not locate the RIP width"
  )
})

test_that("find_rip with verbose = TRUE reports the RIP position in physical units", {
  dt_idx <- 1:200
  tis_shape <- gauss(dt_idx, center = 100, height = 1000, sd = 10)
  int_mat <- matrix(tis_shape, nrow = 200, ncol = 5)
  dt <- dt_idx * 0.1
  rt <- seq_len(5)

  expect_message(
    find_rip(int_mat, verbose = TRUE, retention_time = rt, drift_time = dt),
    "RIP was detected"
  )
})

test_that("find_regions_rip_saturated flags retention-time regions where the RIP is depleted", {
  dt_idx <- 1:200
  tis_shape <- gauss(dt_idx, center = 100, height = 1000, sd = 10)
  int_mat <- matrix(0, nrow = 200, ncol = 20)
  for (j in 1:20) int_mat[, j] <- tis_shape
  int_mat[, 8:10] <- int_mat[, 8:10] * 0.05 # a saturated/depleted RIP window

  rt <- seq(0, by = 1, length.out = 20)
  out <- find_regions_rip_saturated(int_mat, rip_saturation_threshold = 0.1, retention_time = rt)

  expect_equal(dim(out), c(1L, 2L))
  expect_equal(unname(out[1, ]), c(rt[8], rt[10]))
})

test_that("find_regions_rip_saturated returns an empty matrix when nothing is saturated", {
  dt_idx <- 1:200
  tis_shape <- gauss(dt_idx, center = 100, height = 1000, sd = 10)
  int_mat <- matrix(0, nrow = 200, ncol = 20)
  for (j in 1:20) int_mat[, j] <- tis_shape

  rt <- seq(0, by = 1, length.out = 20)
  out <- find_regions_rip_saturated(int_mat, rip_saturation_threshold = 0.1, retention_time = rt)

  expect_equal(nrow(out), 0L)
})
