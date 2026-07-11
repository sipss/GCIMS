test_that("GCIMSChromatogram() sets the mandatory slots", {
  ch <- GCIMSChromatogram(retention_time = 1:3, intensity = c(10, 20, 30))

  expect_s4_class(ch, "GCIMSChromatogram")
  expect_equal(rtime(ch), 1:3)
  expect_equal(intensity(ch), c(10, 20, 30), ignore_attr = TRUE)
  expect_null(baseline(ch, .error_if_missing = FALSE))
})

test_that("GCIMSChromatogram() forwards baseline, peaks and peaks_debug_info to the object", {
  # Regression test: the constructor used to hardcode baseline/peaks/peaks_debug_info
  # to NULL in the methods::new() call, silently discarding these arguments.
  peaks_df <- S4Vectors::DataFrame(dt_apex_ms = 5, rt_apex_s = 10)

  ch <- GCIMSChromatogram(
    retention_time = 1:3,
    intensity = c(10, 20, 30),
    baseline = c(1, 2, 3),
    peaks = peaks_df,
    peaks_debug_info = list(foo = "bar")
  )

  expect_equal(baseline(ch), c(`1` = 1, `2` = 2, `3` = 3))
  expect_equal(nrow(peaks(ch)), 1L)
  expect_equal(peaks(ch)$dt_apex_ms, 5)
  expect_equal(ch@peaks_debug_info, list(foo = "bar"))
})

test_that("GCIMSChromatogram() errors when retention_time and intensity lengths differ", {
  expect_error(
    GCIMSChromatogram(retention_time = 1:3, intensity = c(1, 2)),
    "length"
  )
})

test_that("GCIMSChromatogram() errors when baseline length does not match", {
  expect_error(
    GCIMSChromatogram(retention_time = 1:3, intensity = c(1, 2, 3), baseline = c(1, 2)),
    "length"
  )
})
