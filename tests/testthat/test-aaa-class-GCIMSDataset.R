make_ram_dataset <- function() {
  s1 <- GCIMSSample(1:2, 1:2, matrix(1, nrow = 2, ncol = 2))
  s2 <- GCIMSSample(1:2, 1:2, matrix(2, nrow = 2, ncol = 2))
  GCIMSDataset$new_from_list(samples = list(a = s1, b = s2), on_ram = TRUE, scratch_dir = NULL)
}

test_that("GCIMSDataset can be constructed from files on disk", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  pData <- data.frame(SampleID = "s1", FileName = basename(sample_file))

  ds <- GCIMSDataset$new(pData = pData, base_dir = dirname(sample_file), on_ram = TRUE)

  expect_equal(ds$sampleNames, "s1")
  expect_s4_class(ds$getSample(1), "GCIMSSample")
})

test_that("GCIMSDataset construction from files errors clearly when a file is missing", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  pData <- data.frame(SampleID = "s1", FileName = "does_not_exist.mea")

  expect_error(
    GCIMSDataset$new(pData = pData, base_dir = dirname(sample_file), on_ram = TRUE),
    "does_not_exist.mea"
  )
})

test_that("dataset$sampleNames<- renames samples, updates pData, and validates its input", {
  ds <- make_ram_dataset()
  ds$realize() # sampleNames<- requires no pending operations

  ds$sampleNames <- c("x", "y")

  expect_equal(ds$sampleNames, c("x", "y"))
  expect_equal(as.character(ds$pData$SampleID), c("x", "y"))

  expect_error(
    {
      ds$sampleNames <- "only_one"
    },
    "number of sample names"
  )
  expect_error(
    {
      ds$sampleNames <- c("z", "z")
    },
    "unique and not missing"
  )
})

test_that("is_on_disk distinguishes RAM datasets from disk datasets", {
  ds <- make_ram_dataset()
  expect_false(ds$is_on_disk())
})

test_that("copy() on a RAM dataset returns an independent GCIMSDataset with the same content", {
  ds <- make_ram_dataset()

  ds2 <- ds$copy()

  expect_false(identical(ds, ds2))
  expect_s3_class(ds2, "GCIMSDataset")
  expect_equal(ds2$sampleNames, ds$sampleNames)
})

test_that("updateScratchDir() and getCurrentDir() warn and no-op on a RAM dataset", {
  ds <- make_ram_dataset()

  expect_warning(ds$updateScratchDir(tempdir()), "has no effect")
  expect_warning(result <- ds$getCurrentDir(), "only makes sense for on-disk datasets")
  expect_null(result)
})

test_that("print() describes the dataset without erroring", {
  ds <- make_ram_dataset()

  out <- capture.output(print(ds))

  expect_true(any(grepl("GCIMSDataset", out)))
  expect_true(any(grepl("Stored on RAM", out)))
})

test_that("subset(inplace = FALSE) returns a copy, leaving the original dataset untouched", {
  ds <- make_ram_dataset()

  subset_copy <- ds$subset("a", inplace = FALSE)

  expect_equal(ds$sampleNames, c("a", "b"))
  expect_equal(subset_copy$sampleNames, "a")
})

test_that("subset(inplace = TRUE) modifies the dataset itself", {
  ds <- make_ram_dataset()

  ds$subset("a", inplace = TRUE)

  expect_equal(ds$sampleNames, "a")
})

test_that("subset() accepts a numeric vector of indices with more than one element", {
  # Regression test: sample_name_or_number_to_both() used to hard-error on
  # any numeric index vector of length > 1 (an `||` vs `|` bug), so this is
  # a common, previously-broken usage pattern.
  s1 <- GCIMSSample(1:2, 1:2, matrix(1, nrow = 2, ncol = 2))
  s2 <- GCIMSSample(1:2, 1:2, matrix(2, nrow = 2, ncol = 2))
  s3 <- GCIMSSample(1:2, 1:2, matrix(3, nrow = 2, ncol = 2))
  ds <- GCIMSDataset$new_from_list(samples = list(a = s1, b = s2, c = s3), on_ram = TRUE, scratch_dir = NULL)

  ds$subset(c(1, 3), inplace = TRUE)

  expect_equal(ds$sampleNames, c("a", "c"))
})

test_that("subset() also filters the peaks slot down to the remaining samples", {
  ds <- make_ram_dataset()
  ds$peaks <- data.frame(SampleID = c("a", "a", "b"), PeakID = c("1", "2", "1"))

  ds$subset("a", inplace = TRUE)

  expect_equal(ds$peaks$SampleID, c("a", "a"))
})
