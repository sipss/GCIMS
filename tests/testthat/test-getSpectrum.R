test_that("getSpectrum() works", {
  s1 <- GCIMSSample(
    drift_time = 1:8,
    retention_time = 1:7,
    data = matrix(1:56, nrow = 8)
  )
  expected <- rowSums(intensity(s1)[2:4, 2:5])
  sp1 <- getSpectrum(s1, dt_range = c(1.5, 4.5), rt_range = c(1.5, 5.5))
  implemented <- intensity(sp1)
  expect_equal(implemented, expected)

  expected <- apply(intensity(s1)[2:4, 2:5], 1, max)
  sp1 <- getSpectrum(s1, dt_range = c(1.5, 4.5), rt_range = c(1.5, 5.5), aggregate = function(x) apply(x, 1, max))
  implemented <- intensity(sp1)
  expect_equal(implemented, expected)
})
