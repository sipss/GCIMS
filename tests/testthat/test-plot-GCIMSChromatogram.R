test_that("plot() on a GCIMSChromatogram returns a ggplot with data, title and labels from the object", {
  ch <- GCIMSChromatogram(retention_time = 1:5, intensity = c(1, 2, 5, 2, 1), description = "Chrom1")

  p <- plot(ch)

  expect_s3_class(p, "ggplot")
  expect_equal(p$data$x, rtime(ch))
  expect_equal(p$data$y, intensity(ch), ignore_attr = TRUE)
  expect_equal(p$labels$title, "Chrom1")
  expect_equal(p$labels$x, "Retention time (s)")
  expect_equal(p$labels$y, "Intensity (a.u.)")
})

test_that("plot() on a GCIMSChromatogram builds a single-value drift time subtitle", {
  ch <- GCIMSChromatogram(
    retention_time = 1:5, intensity = c(1, 2, 5, 2, 1),
    drift_time_idx = 3L, drift_time_ms = 3.5
  )

  p <- plot(ch)

  expect_equal(as.character(p$labels$subtitle), "Drift time 3.5 ms")
})

test_that("plot() on a GCIMSChromatogram builds a range drift time subtitle for multiple values", {
  ch <- GCIMSChromatogram(
    retention_time = 1:5, intensity = c(1, 2, 5, 2, 1),
    drift_time_idx = c(3L, 4L), drift_time_ms = c(3, 4)
  )

  p <- plot(ch)

  expect_equal(as.character(p$labels$subtitle), "Drift time 3 - 4 ms")
})

test_that("plot() on a GCIMSChromatogram has no subtitle when drift_time_ms is empty", {
  ch <- GCIMSChromatogram(
    retention_time = 1:5, intensity = c(1, 2, 5, 2, 1),
    drift_time_idx = integer(0), drift_time_ms = numeric(0)
  )

  p <- plot(ch)

  expect_null(p$labels$subtitle)
})
