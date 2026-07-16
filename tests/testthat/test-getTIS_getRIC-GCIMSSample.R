make_rip_sample <- function() {
  dt <- seq(1, 10, by = 0.5) # 19 pts, dt[10] == 5.5
  rt <- seq(1, 60, by = 1)
  data <- matrix(1, nrow = length(dt), ncol = length(rt))
  data[10, ] <- 1000 # unambiguous RIP row at dt = 5.5
  GCIMSSample(drift_time = dt, retention_time = rt, data = data, description = "s1")
}

test_that("plotTIS(GCIMSSample) plots the full spectrum by default, with an accurate subtitle", {
  s <- make_rip_sample()

  p <- plotTIS(s)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Drift time (ms)")
  expect_equal(p$labels$y, "TIS Intensity (a.u.)")
  expect_equal(p$labels$title, "s1")
  expect_match(p$labels$subtitle, "Retention time 1 - 60 s")
  expect_equal(unname(p$data$y), unname(getTIS(s)))
})

test_that("plotTIS(GCIMSSample) restricts the displayed drift time range without changing values", {
  s <- make_rip_sample()

  p_full <- plotTIS(s)
  p_sub <- plotTIS(s, dt_range = c(2, 8))

  expect_equal(range(p_sub$data$x), c(2, 8))
  # Values within the window are identical to the unrestricted plot's values:
  common <- merge(p_full$data, p_sub$data, by = "x")
  expect_equal(common$y.x, common$y.y)
})

test_that("plotRIC(GCIMSSample) has an accurate single-value drift time subtitle (the RIP position)", {
  s <- make_rip_sample()

  p <- plotRIC(s)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Retention time (s)")
  expect_equal(p$labels$y, "RIC Intensity (a.u.)")
  expect_equal(p$labels$title, "s1")
  expect_equal(p$labels$subtitle, "Drift time 5.5 ms")
})

test_that("plotRIC(GCIMSSample) restricts the displayed retention time range without changing values", {
  s <- make_rip_sample()

  p_full <- plotRIC(s)
  p_sub <- plotRIC(s, rt_range = c(10, 30))

  expect_equal(range(p_sub$data$x), c(10, 30))
  common <- merge(p_full$data, p_sub$data, by = "x")
  expect_equal(common$y.x, common$y.y)
})
