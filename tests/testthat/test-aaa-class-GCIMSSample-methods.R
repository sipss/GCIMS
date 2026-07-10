make_sample <- function() {
  GCIMSSample(drift_time = 1:2, retention_time = 1:3, data = matrix(1:6, nrow = 2, ncol = 3))
}

test_that("initialize() errors when extra arguments are unnamed", {
  expect_error(
    methods::new(
      "GCIMSSample",
      drift_time = 1:2, retention_time = 1:3, data = matrix(1:6, 2, 3),
      "nitrogen"
    ),
    "All extra arguments should be named"
  )
})

test_that("[.GCIMSSample with both indices missing returns the object unchanged", {
  x <- make_sample()
  expect_identical(x[], x)
})

test_that("[.GCIMSSample subsets by drift time (i) and retention time (j)", {
  x <- make_sample()

  by_dt <- x[1, ]
  expect_equal(dim(by_dt), c(1L, 3L))
  expect_equal(dtime(by_dt), 1)

  by_rt <- x[, 1:2]
  expect_equal(dim(by_rt), c(2L, 2L))
  expect_equal(rtime(by_rt), 1:2)
})

test_that("[.GCIMSSample validates logical index lengths", {
  x <- make_sample()

  expect_error(
    x[c(TRUE, FALSE, TRUE), ],
    "Incorrect number of logical values provided to subset drift time"
  )
  expect_error(
    x[, c(TRUE, FALSE)],
    "Incorrect number of logical values provided to subset retention time"
  )
})

test_that("dim.GCIMSSample() returns the dimensions of the data matrix", {
  x <- make_sample()
  expect_equal(dim(x), c(2L, 3L))
})

test_that("subset.GCIMSSample() also subsets the baseline when it is set", {
  x <- make_sample()
  baseline(x) <- matrix(0, nrow = 2, ncol = 3)

  xs <- subset.GCIMSSample(x, dt_idx = 1)

  expect_equal(dim(baseline(xs, .error_if_missing = FALSE)), c(1L, 3L))
})

test_that("getChromatogram() errors when the aggregate function returns the wrong length", {
  x <- make_sample()
  expect_error(
    getChromatogram(x, aggregate = function(m) 1),
    "aggregation outcome should be a vector of length 3"
  )
})

test_that("getChromatogram() propagates the sample's baseline, aggregated the same way as intensity", {
  x <- make_sample()
  baseline(x) <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)

  ch <- getChromatogram(x)

  expect_equal(baseline(ch), colSums(baseline(x)), ignore_attr = TRUE)
})

test_that("getSpectrum() errors when the aggregate function returns the wrong length", {
  x <- make_sample()
  expect_error(
    getSpectrum(x, aggregate = function(m) 1),
    "aggregation outcome should be a vector of length 2"
  )
})

test_that("getSpectrum() propagates the sample's baseline, aggregated the same way as intensity", {
  x <- make_sample()
  baseline(x) <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)

  sp <- getSpectrum(x)

  expect_equal(baseline(sp), rowSums(baseline(x)), ignore_attr = TRUE)
})
