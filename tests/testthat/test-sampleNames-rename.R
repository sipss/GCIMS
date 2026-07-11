# Regression tests for a bug where dataset$sampleNames <- newNames could
# crash the next getSample()/realize() with "description should be a
# string". Root cause: the delayed "setSampleNamesAsDescription" operation
# that keeps each sample's internal description in sync with its
# dataset-level name is matched against the sample's identity at the moment
# it actually runs during realize() -- but if other operations were still
# pending (including the one queued automatically at construction time),
# that identity could already have moved past what an earlier-queued copy
# of this operation expected.
#
# Fixed by requiring the dataset to have no pending operations before a
# rename is allowed, so there is always at most one such operation queued,
# unambiguously matching the sample's current identity. See
# aaa-class-GCIMSDataset.R.
#
# One narrow exception: if the *only* pending operation is itself an
# earlier, not-yet-realized rename, it is safe to just replace it, since
# nothing else was queued in between that could have observed it. See
# DelayedDatasetBase$dropSolePendingOp().

make_ram_pair <- function() {
  s1 <- GCIMSSample(drift_time = 1:2, retention_time = 1:2, data = matrix(1, nrow = 2, ncol = 2))
  s2 <- GCIMSSample(drift_time = 1:2, retention_time = 1:2, data = matrix(2, nrow = 2, ncol = 2))
  GCIMSDataset_fromList(list(a = s1, b = s2))
}

test_that("renaming a dataset with pending operations aborts with a clear message", {
  ds <- make_ram_pair()

  expect_true(ds$hasDelayedOps()) # construction always queues some ops
  expect_error(
    {
      ds$sampleNames <- c("x", "y")
    },
    "pending operations"
  )
})

test_that("renaming a RAM dataset works correctly right after an explicit realize()", {
  ds <- make_ram_pair()
  ds$realize()

  ds$sampleNames <- c("x", "y")

  expect_equal(ds$sampleNames, c("x", "y"))
  s <- ds$getSample("x")
  expect_equal(intensity(s)[1, 1], 1)
  expect_equal(description(s), "x")
})

test_that("renaming a RAM dataset works correctly after getSample() has implicitly realized it", {
  ds <- make_ram_pair()
  ds$getSample(1)

  ds$sampleNames <- c("x", "y")

  s <- ds$getSample("x")
  expect_equal(description(s), "x")
})

test_that("renaming an on-disk dataset with pending operations also aborts", {
  scratch <- tempfile("gcims_scratch_")
  dir.create(scratch)
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  pData <- data.frame(SampleID = "s1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(pData = pData, base_dir = dirname(sample_file), on_ram = FALSE, scratch_dir = scratch)

  expect_error(
    {
      ds$sampleNames <- "renamed"
    },
    "pending operations"
  )
})

test_that("renaming an on-disk dataset works correctly after realize()", {
  scratch <- tempfile("gcims_scratch_")
  dir.create(scratch)
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  pData <- data.frame(SampleID = "s1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(pData = pData, base_dir = dirname(sample_file), on_ram = FALSE, scratch_dir = scratch)
  ds$getSample(1)

  ds$sampleNames <- "renamed"

  s <- ds$getSample("renamed")
  expect_equal(description(s), "renamed")
})

test_that("swapping two sample names resolves each sample to the correct data and description", {
  ds <- make_ram_pair()
  ds$realize()

  ds$sampleNames <- c("b", "a")

  sb <- ds$getSample("b")
  sa <- ds$getSample("a")
  expect_equal(intensity(sb)[1, 1], 1)
  expect_equal(description(sb), "b")
  expect_equal(intensity(sa)[1, 1], 2)
  expect_equal(description(sa), "a")
})

test_that("a second rename replaces the first when nothing else is pending in between (RAM)", {
  ds <- make_ram_pair()
  ds$realize()

  ds$sampleNames <- c("x", "y")
  ds$sampleNames <- c("p", "q") # no realize() needed in between

  expect_equal(ds$sampleNames, c("p", "q"))
  sp <- ds$getSample("p")
  sq <- ds$getSample("q")
  expect_equal(intensity(sp)[1, 1], 1)
  expect_equal(description(sp), "p")
  expect_equal(intensity(sq)[1, 1], 2)
  expect_equal(description(sq), "q")
})

test_that("a second rename replaces the first when nothing else is pending in between (disk)", {
  scratch <- tempfile("gcims_scratch_")
  dir.create(scratch)
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  pData <- data.frame(SampleID = "s1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(pData = pData, base_dir = dirname(sample_file), on_ram = FALSE, scratch_dir = scratch)
  ds$getSample(1)

  ds$sampleNames <- "x"
  ds$sampleNames <- "p"

  s <- ds$getSample("p")
  expect_equal(description(s), "p")
})

test_that("renaming still aborts if another operation was queued between two renames", {
  ds <- make_ram_pair()
  ds$realize()
  ds$sampleNames <- c("x", "y")
  ds$appendDelayedOp(DelayedOperation(name = "noop", fun = function(sample) sample))

  expect_error(
    {
      ds$sampleNames <- c("p", "q")
    },
    "pending operations"
  )
})

test_that("dropSolePendingOp() only drops when it is the sole pending operation with that name", {
  dd <- DelayedDatasetRAM$new(samples = list(a = 1))

  # nothing pending: nothing to drop
  expect_false(dd$dropSolePendingOp("foo"))

  dd$appendDelayedOp(DelayedOperation(name = "foo", fun = function(x) x))
  expect_true(dd$dropSolePendingOp("foo"))
  expect_false(dd$hasDelayedOps())

  # two operations pending: refuses, even if one of them matches
  dd$appendDelayedOp(DelayedOperation(name = "foo", fun = function(x) x))
  dd$appendDelayedOp(DelayedOperation(name = "bar", fun = function(x) x))
  expect_false(dd$dropSolePendingOp("foo"))
  expect_true(dd$hasDelayedOps())
})
