test_that("show() prints axis ranges, step and point counts for a multi-point sample", {
  s <- GCIMSSample(
    drift_time = seq(0, 10, by = 0.5),
    retention_time = seq(0, 20, by = 1),
    data = matrix(1, nrow = 21, ncol = 21)
  )

  out <- capture.output(print(s))

  expect_equal(
    out,
    c(
      "A GCIMS Sample",
      " with drift time from 0 to 10 ms (step: 0.5 ms, points: 21)",
      " with retention time from 0 to 20 s (step: 1 s, points: 21)"
    )
  )
})

test_that("show() handles a single-point axis without a step", {
  s <- GCIMSSample(drift_time = 5, retention_time = 10, data = matrix(1, 1, 1))

  out <- capture.output(print(s))

  expect_equal(
    out,
    c(
      "A GCIMS Sample",
      " with drift time from 5 to 5 ms (step: NaN ms, points: 1)",
      " with retention time from 10 to 10 s (step: NaN s, points: 1)"
    )
  )
})
