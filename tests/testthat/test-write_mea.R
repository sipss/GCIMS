test_that("Write and read a mea file to get the same object", {
  obj <- GCIMSSample(
   drift_time = 0:1,
   retention_time = 0:2,
   data = matrix(1:6, nrow = 2, ncol = 3)
  )
  mea_file <- tempfile(fileext = ".mea")
  write_mea(obj, mea_file)
  obj2 <- read_mea(mea_file)
  expect_equal(dtime(obj2) , dtime(obj))
  expect_equal(rtime(obj2) , rtime(obj))
  expect_equal(intensity(obj2) , intensity(obj))
  file.remove(mea_file)
})


test_that("read, write, read mea file to get the same object", {
  obj <- read_mea(system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS"))
  mea_file <- tempfile(fileext = ".mea")
  write_mea(obj, mea_file)
  obj2 <- read_mea(mea_file)
  expect_equal(dtime(obj2) , dtime(obj))
  expect_equal(rtime(obj2) , rtime(obj))
  expect_equal(intensity(obj2) , intensity(obj))
  file.remove(mea_file)
})

