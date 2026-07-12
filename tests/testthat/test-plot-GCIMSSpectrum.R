test_that("plot() on a GCIMSSpectrum returns a ggplot with data, title and labels from the object", {
  sp <- GCIMSSpectrum(drift_time = 1:5, intensity = c(1, 2, 5, 2, 1), description = "Spec1")

  p <- plot(sp)

  expect_s3_class(p, "ggplot")
  expect_equal(p$data$x, dtime(sp))
  expect_equal(p$data$y, intensity(sp), ignore_attr = TRUE)
  expect_equal(p$labels$title, "Spec1")
  expect_equal(p$labels$x, "Drift time (ms)")
  expect_equal(p$labels$y, "Intensity (a.u.)")
})

test_that("plot() on a GCIMSSpectrum builds a single-value retention time subtitle", {
  sp <- GCIMSSpectrum(
    drift_time = 1:5, intensity = c(1, 2, 5, 2, 1),
    retention_time_idx = 3L, retention_time_s = 3.5
  )

  p <- plot(sp)

  expect_equal(as.character(p$labels$subtitle), "Retention time 3.5 s")
})

test_that("plot() on a GCIMSSpectrum builds a range retention time subtitle for multiple values", {
  sp <- GCIMSSpectrum(
    drift_time = 1:5, intensity = c(1, 2, 5, 2, 1),
    retention_time_idx = c(3L, 4L), retention_time_s = c(3, 4)
  )

  p <- plot(sp)

  expect_equal(as.character(p$labels$subtitle), "Retention time 3 - 4 s")
})

test_that("plot() on a GCIMSSpectrum has no subtitle when retention_time_s is empty", {
  sp <- GCIMSSpectrum(
    drift_time = 1:5, intensity = c(1, 2, 5, 2, 1),
    retention_time_idx = integer(0), retention_time_s = numeric(0)
  )

  p <- plot(sp)

  expect_null(p$labels$subtitle)
})
