make_dataset <- function() {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  annot <- data.frame(SampleID = "Sample1", FileName = basename(sample_file))
  GCIMSDataset$new(annot, base_dir = dirname(sample_file), on_ram = TRUE)
}

test_that("dtime()/rtime() realize a never-realized dataset to compute the reference axes", {
  ds <- make_dataset()

  dt <- dtime(ds)
  rt <- rtime(ds)

  expect_true(length(dt) > 0)
  expect_true(length(rt) > 0)
  expect_false(ds$hasDelayedOps())
})

test_that("dtime()/rtime() with sample= return that sample's own axis", {
  ds <- make_dataset()

  dt_sample <- dtime(ds, sample = 1)
  rt_sample <- rtime(ds, sample = "Sample1")

  s <- ds$getSample(1)
  expect_equal(dt_sample, dtime(s))
  expect_equal(rt_sample, rtime(s))
})
