make_sample <- function() {
  GCIMSSample(drift_time = 1:8, retention_time = 1:7, data = matrix(1:56, nrow = 8))
}

test_that("filterRt() on a GCIMSSample subsets by retention time", {
  s <- make_sample()

  out <- filterRt(s, rt_range = c(2, 5))

  expect_equal(rtime(out), 2:5)
})

test_that("filterDt() on a GCIMSSample subsets by drift time", {
  s <- make_sample()

  out <- filterDt(s, dt_range = c(3, 6))

  expect_equal(dtime(out), 3:6)
})
