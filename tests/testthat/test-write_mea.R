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

test_that("write_mea validates the filename argument", {
  obj <- GCIMSSample(drift_time = 0:1, retention_time = 0:2, data = matrix(1:6, nrow = 2, ncol = 3))

  expect_error(write_mea(obj, 123), "filename should be a string")
  expect_error(write_mea(obj, c("a", "b")), "filename should be of length one")
})

test_that("write_mea gzip-compresses when given a .mea.gz filename, and appends .mea otherwise", {
  obj <- GCIMSSample(drift_time = 0:1, retention_time = 0:2, data = matrix(1:6, nrow = 2, ncol = 3))

  gz_file <- tempfile(fileext = ".mea.gz")
  write_mea(obj, gz_file)
  expect_true(file.exists(gz_file))
  obj_gz <- read_mea(gz_file)
  expect_equal(intensity(obj_gz), intensity(obj))
  file.remove(gz_file)

  no_ext_file <- tempfile()
  write_mea(obj, no_ext_file)
  expect_true(file.exists(paste0(no_ext_file, ".mea")))
  file.remove(paste0(no_ext_file, ".mea"))
})

