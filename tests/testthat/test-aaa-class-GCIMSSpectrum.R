test_that("GCIMSSpectrum() constructs a valid object", {
  sp <- GCIMSSpectrum(drift_time = 1:5, intensity = c(1, 2, 5, 2, 1))

  expect_s4_class(sp, "GCIMSSpectrum")
  expect_equal(dtime(sp), 1:5)
})

test_that("GCIMSSpectrum() rejects unknown arguments, listing the valid ones", {
  expect_error(
    GCIMSSpectrum(drift_time = 1:5, intensity = 1:5, bogus_arg = "x"),
    "bogus_arg"
  )
})

test_that("GCIMSSpectrum() requires drift_time and intensity to have the same length", {
  expect_error(GCIMSSpectrum(drift_time = 1:5, intensity = 1:3))
})

test_that("GCIMSSpectrum() requires baseline to be NULL or the same length as drift_time", {
  expect_error(GCIMSSpectrum(drift_time = 1:5, intensity = 1:5, baseline = 1:3))

  sp <- GCIMSSpectrum(drift_time = 1:5, intensity = c(1, 2, 5, 2, 1), baseline = c(0, 0, 0, 0, 0))
  expect_equal(baseline(sp, .error_if_missing = FALSE), c(0, 0, 0, 0, 0), ignore_attr = TRUE)
})
