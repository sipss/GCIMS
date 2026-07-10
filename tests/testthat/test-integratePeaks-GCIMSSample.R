gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

dt_axis <- seq(0, 4, by = 0.02) # 201 points
rt_axis <- seq(0, 100, by = 0.5) # 201 points

# A sample with a RIP row (dominant TIS row, near the low end of drift time)
# that's depleted/"saturated" over retention-time indices 90-100, plus two
# constant-intensity peak blocks elsewhere in drift time: one inside a
# non-saturated retention-time window, one inside the saturated window.
make_integration_sample <- function() {
  int_mat <- matrix(0.1, nrow = length(dt_axis), ncol = length(rt_axis))
  rip_profile <- gauss(seq_along(dt_axis), center = 20, height = 1000, sd = 5)
  for (j in seq_along(rt_axis)) {
    scale_j <- if (j >= 90 && j <= 100) 0.02 else 1
    int_mat[, j] <- int_mat[, j] + rip_profile * scale_j
  }
  int_mat[150:160, 40:50] <- int_mat[150:160, 40:50] + 50 # peak A: non-saturated window
  int_mat[150:160, 92:98] <- int_mat[150:160, 92:98] + 50 # peak B: saturated window
  s <- GCIMSSample(drift_time = dt_axis, retention_time = rt_axis, data = int_mat)
  description(s) <- "S1"
  s
}

make_peak_list <- function() {
  data.frame(
    SampleID = c("S1", "S1"),
    dt_min_ms = dt_axis[150], dt_max_ms = dt_axis[160],
    rt_min_s = c(rt_axis[40], rt_axis[92]),
    rt_max_s = c(rt_axis[50], rt_axis[98]),
    rt_apex_s = c(rt_axis[45], rt_axis[95]),
    rt_cm_s = c(rt_axis[45], rt_axis[95]),
    fixedsize_dt_min_ms = dt_axis[148], fixedsize_dt_max_ms = dt_axis[162],
    fixedsize_rt_min_s = c(rt_axis[38], rt_axis[90]),
    fixedsize_rt_max_s = c(rt_axis[52], rt_axis[100])
  )
}

test_that("integratePeaks computes Area from the ROI's physical dimensions", {
  s <- make_integration_sample()
  peak_list <- make_peak_list()

  out <- integratePeaks(s, peak_list, integration_size_method = "free_size")
  pl <- peaks(out)

  expected_area <- (dt_axis[160] - dt_axis[150]) * (rt_axis[50] - rt_axis[40])
  expect_equal(pl$Area[1], expected_area)
})

test_that("integratePeaks computes Volume as the sum of intensities in the ROI", {
  s <- make_integration_sample()
  peak_list <- make_peak_list()

  out <- integratePeaks(s, peak_list, integration_size_method = "free_size")
  pl <- peaks(out)

  patch <- intensity(s, dt_range = c(dt_axis[150], dt_axis[160]), rt_range = c(rt_axis[40], rt_axis[50]))
  expected_volume <- sum(patch) * (rt_axis[2] - rt_axis[1]) * (dt_axis[2] - dt_axis[1])
  expect_equal(pl$Volume[1], expected_volume)
})

test_that("integratePeaks with fixed_size uses the wider fixedsize_* bounds instead of the ROI's own bounds", {
  s <- make_integration_sample()
  peak_list <- make_peak_list()

  out_free <- integratePeaks(s, peak_list, integration_size_method = "free_size")
  out_fixed <- integratePeaks(s, peak_list, integration_size_method = "fixed_size")

  expect_gt(peaks(out_fixed)$Volume[1], peaks(out_free)$Volume[1])
})

test_that("integratePeaks computes a symmetric Asymmetry of 0 for a centered apex", {
  s <- make_integration_sample()
  peak_list <- make_peak_list()
  peak_list$rt_apex_s[1] <- (rt_axis[40] + rt_axis[50]) / 2 # dead center

  out <- integratePeaks(s, peak_list, integration_size_method = "free_size")

  expect_equal(peaks(out)$Asymmetry[1], 0)
})

test_that("integratePeaks computes a nonzero Asymmetry for an off-center apex", {
  s <- make_integration_sample()
  peak_list <- make_peak_list()
  peak_list$rt_apex_s[1] <- rt_axis[43]

  out <- integratePeaks(s, peak_list, integration_size_method = "free_size")

  rising <- rt_axis[43] - rt_axis[40]
  falling <- rt_axis[50] - rt_axis[43]
  expect_equal(peaks(out)$Asymmetry[1], round(falling / rising - 1, 2))
})

test_that("integratePeaks flags Saturation only for peaks inside a RIP-depleted retention-time window", {
  s <- make_integration_sample()
  peak_list <- make_peak_list()

  out <- integratePeaks(s, peak_list, integration_size_method = "free_size", rip_saturation_threshold = 0.1)
  pl <- peaks(out)

  expect_false(pl$Saturation[1]) # non-saturated window
  expect_true(pl$Saturation[2]) # saturated window
})

test_that("integratePeaks only integrates peaks belonging to the given sample", {
  s <- make_integration_sample()
  peak_list <- make_peak_list()
  peak_list$SampleID[2] <- "SomeOtherSample"

  out <- integratePeaks(s, peak_list, integration_size_method = "free_size")

  expect_equal(nrow(peaks(out)), 1L)
  expect_equal(peaks(out)$SampleID, "S1")
})
