test_that("dtime works on GCIMSSampleS3", {
  obj <- new_GCIMSSampleS3(
    drift_time = 1:3,
    retention_time = 1:5,
    intensity = matrix(1:15, nrow = 3, ncol = 5),
    metadata = list(Name = "Sample1")
  )
  expect_true(all(dtime(obj) == 1:3))
})
