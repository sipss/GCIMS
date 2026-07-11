make_dataset <- function() {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  annot <- data.frame(SampleID = "Sample1", FileName = basename(sample_file))
  GCIMSDataset$new(annot, base_dir = dirname(sample_file), on_ram = TRUE)
}

test_that("filterRt() queues a retention-time filter applied on realize", {
  ds <- make_dataset()

  filterRt(ds, rt_range = c(5, 50))
  s <- ds$getSample(1)

  expect_true(min(rtime(s)) >= 5)
  expect_true(max(rtime(s)) <= 50)
})

test_that("filterDt() queues a drift-time filter applied on realize", {
  ds <- make_dataset()

  filterDt(ds, dt_range = c(5, 10))
  s <- ds$getSample(1)

  expect_true(min(dtime(s)) >= 5)
  expect_true(max(dtime(s)) <= 10)
})
