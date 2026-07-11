gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

make_sample <- function() {
  dt <- seq(0, 4, by = 0.02)
  rt <- seq(0, 50, by = 0.5)
  int_mat <- matrix(50, nrow = length(dt), ncol = length(rt)) +
    outer(gauss(dt, 2, 500, 0.1), gauss(rt, 25, 1, 3))
  GCIMSSample(drift_time = dt, retention_time = rt, data = int_mat)
}

test_that("as.data.frame() on a GCIMSSpectrum includes drift_time_ms and intensity, and baseline if set", {
  s <- make_sample()
  sp <- getSpectrum(s)

  out_no_basel <- as.data.frame(sp)
  expect_named(out_no_basel, c("drift_time_ms", "intensity"))
  expect_equal(out_no_basel$drift_time_ms, dtime(sp))
  expect_equal(out_no_basel$intensity, intensity(sp), ignore_attr = TRUE)

  sp_b <- estimateBaseline(sp, dt_peak_fwhm_ms = 0.3, dt_region_multiplier = 6)
  out_basel <- as.data.frame(sp_b)
  expect_named(out_basel, c("drift_time_ms", "intensity", "baseline"))
})

test_that("as.data.frame() on a GCIMSSpectrum uses given row.names", {
  s <- make_sample()
  sp <- getSpectrum(s)

  out <- as.data.frame(sp, row.names = paste0("r", seq_along(dtime(sp))))

  expect_equal(rownames(out)[1:3], c("r1", "r2", "r3"))
})

test_that("as.data.frame() on a GCIMSChromatogram includes retention_time_s and intensity, and baseline if set", {
  s <- make_sample()
  ch <- getChromatogram(s)

  out_no_basel <- as.data.frame(ch)
  expect_named(out_no_basel, c("retention_time_s", "intensity"))
  expect_equal(out_no_basel$retention_time_s, rtime(ch))

  ch_b <- estimateBaseline(ch, rt_length_s = 10)
  out_basel <- as.data.frame(ch_b)
  expect_named(out_basel, c("retention_time_s", "intensity", "baseline"))
})

test_that("as.data.frame() on a GCIMSChromatogram uses given row.names", {
  s <- make_sample()
  ch <- getChromatogram(s)

  out <- as.data.frame(ch, row.names = paste0("r", seq_along(rtime(ch))))

  expect_equal(rownames(out)[1:3], c("r1", "r2", "r3"))
})
