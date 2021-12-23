test_that("read_mea does not crash", {
  mea_file <- system.file("extdata/sample_formats/211108_153700.mea.gz", package = "GCIMS")
  sample <- read_mea(mea_file)
})
