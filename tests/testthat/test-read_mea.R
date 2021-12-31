test_that("read_mea does not crash", {
  mea_file <- system.file("extdata/sample_formats/small.mea.gz", package = "GCIMS")
  sample <- read_mea(mea_file)
  expect_s4_class(sample, "GCIMSSample")
})
