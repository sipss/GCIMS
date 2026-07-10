gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

make_sample_with_baseline <- function() {
  dt <- seq(0, 4, by = 0.02) # 201 pts
  rt <- seq(0, 50, by = 0.5) # 101 pts
  int_mat <- matrix(50, nrow = length(dt), ncol = length(rt)) + # flat baseline offset
    outer(gauss(dt, 2, 500, 0.1), gauss(rt, 25, 1, 3)) # a peak
  GCIMSSample(drift_time = dt, retention_time = rt, data = int_mat)
}

test_that("baseline() errors before estimateBaseline() has been run", {
  s <- make_sample_with_baseline()

  expect_error(baseline(s), "estimateBaseline")
  expect_null(baseline(s, .error_if_missing = FALSE))
})

test_that("estimateBaseline() with remove = TRUE estimates and subtracts a flat baseline", {
  s <- make_sample_with_baseline()

  out <- estimateBaseline(s, dt_peak_fwhm_ms = 0.3, dt_region_multiplier = 6, rt_length_s = 10, remove = TRUE)
  bl <- baseline(out)

  expect_equal(dim(bl), dim(intensity(s)))
  # Far from the peak, the estimated baseline tracks the known flat offset (50)
  # and the residual intensity is close to zero:
  expect_equal(bl[1, 1], 50, tolerance = 5)
  expect_equal(intensity(out)[1, 1], 0, tolerance = 5)
  # The peak itself survives baseline removal:
  dt <- dtime(s)
  rt <- rtime(s)
  expect_gt(intensity(out)[which(dt == 2), which(rt == 25)], 400)
})

test_that("estimateBaseline() with remove = FALSE sets the baseline but leaves the data untouched", {
  s <- make_sample_with_baseline()

  out <- estimateBaseline(s, dt_peak_fwhm_ms = 0.3, dt_region_multiplier = 6, rt_length_s = 10, remove = FALSE)

  expect_identical(intensity(out), intensity(s))
  expect_false(is.null(baseline(out, .error_if_missing = FALSE)))
})

test_that("baseline() subsets by dt_range/rt_range and sets dimnames", {
  s <- make_sample_with_baseline()
  out <- estimateBaseline(s, dt_peak_fwhm_ms = 0.3, dt_region_multiplier = 6, rt_length_s = 10)

  bl_sub <- baseline(out, dt_range = c(0, 1), rt_range = c(0, 10))

  expect_equal(dim(bl_sub), c(51L, 21L))
  expect_named(dimnames(bl_sub), c("dt_ms", "rt_s"))
})

test_that("baseline<- rejects a value with the wrong dimensions", {
  s <- make_sample_with_baseline()

  expect_error(
    {
      baseline(s) <- matrix(1, 2, 2)
    },
    "same length as the intensity"
  )
})

test_that("baseline<- accepts a matrix matching the sample's dimensions", {
  s <- make_sample_with_baseline()
  bl <- matrix(0, nrow = nrow(intensity(s)), ncol = ncol(intensity(s)))

  baseline(s) <- bl

  expect_equal(baseline(s), bl, ignore_attr = TRUE)
})
