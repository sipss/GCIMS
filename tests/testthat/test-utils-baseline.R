gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

test_that("compute_baseline connects per-region local minima with linear interpolation", {
  aux <- matrix(c(5, 1, 8, 6, 2, 9, 7, 0, 3), ncol = 1)
  x <- seq_len(nrow(aux))

  bl <- compute_baseline(aux, x, number_of_regions = 3, region_size = 3)

  # Region minima are at (x=2, y=1), (x=5, y=2), (x=8, y=0); connected by
  # linear interpolation (and extrapolated at the edges):
  expected <- signal::interp1(x = c(2, 5, 8), y = c(1, 2, 0), xi = 1:9, method = "linear", extrap = TRUE)
  expect_equal(as.numeric(bl), expected)
})

test_that("compute_baseline handles a shorter last region", {
  aux <- matrix(c(5, 1, 8, 6, 2), ncol = 1) # 5 points, region_size 3 -> regions [1-3], [4-5]

  bl <- compute_baseline(aux, x = 1:5, number_of_regions = 2, region_size = 3)

  # Region minima at (x=2, y=1) and (x=5, y=2):
  expected <- signal::interp1(x = c(2, 5), y = c(1, 2), xi = 1:5, method = "linear", extrap = TRUE)
  expect_equal(as.numeric(bl), expected)
})

test_that("compute_baseline processes each column independently", {
  aux <- matrix(c(5, 1, 8, 10, 2, 16, 6, 9, 3, 12, 18, 6), nrow = 6, ncol = 2)

  bl <- compute_baseline(aux, x = 1:6, number_of_regions = 2, region_size = 3)

  # Column 1: region minima at (x=2, y=1) and (x=5, y=2)
  expect_equal(bl[, 1], signal::interp1(x = c(2, 5), y = c(1, 2), xi = 1:6, method = "linear", extrap = TRUE))
  # Column 2: region minima at (x=3, y=3) and (x=6, y=6)
  expect_equal(bl[, 2], signal::interp1(x = c(3, 6), y = c(3, 6), xi = 1:6, method = "linear", extrap = TRUE))
})

test_that("estimate_baseline_td estimates a baseline close to a known flat offset", {
  y <- gauss(1:200, 50, 100, 5)
  aux <- matrix(y + 10, ncol = 1) # a peak riding on a flat baseline of 10

  bl <- estimate_baseline_td(aux, sig_mult = 12, peak_fwhm_pts = 12)

  expect_equal(dim(bl), dim(aux))
  # Far from the peak, the baseline should track the flat offset:
  expect_equal(bl[1, 1], 10, tolerance = 1)
})

test_that("estimate_baseline_td auto-detects peak_fwhm_pts from the RIP position when not given", {
  y <- gauss(1:200, 50, 100, 5)
  aux <- matrix(y, ncol = 1)

  bl <- estimate_baseline_td(aux, sig_mult = 12)

  expect_equal(dim(bl), dim(aux))
  expect_lt(bl[50, 1], 10) # near the peak apex, baseline should stay low
})

test_that("estimate_baseline_tr applies compute_baseline along the retention-time axis", {
  aux <- matrix(1:12, nrow = 3, ncol = 4)

  bl <- estimate_baseline_tr(aux, region_size = 2)

  expected <- t(compute_baseline(t(aux), x = seq_len(ncol(aux)), number_of_regions = ceiling(ncol(aux) / 2), region_size = 2))
  expect_identical(bl, expected)
  expect_equal(dim(bl), dim(aux))
})
