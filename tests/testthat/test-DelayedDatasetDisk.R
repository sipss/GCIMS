new_scratch_dir <- function() {
  dir <- tempfile("gcims_disk_test_")
  dir.create(dir)
  dir
}

make_disk_dataset <- function(n = 2, scratch_dir = new_scratch_dir(), keep_intermediate = FALSE) {
  samples <- lapply(seq_len(n), function(i) GCIMSSample(1:2, 1:2, matrix(i, 2, 2)))
  names(samples) <- letters[seq_len(n)]
  DelayedDatasetDisk$new(
    samples = samples,
    scratch_dir = scratch_dir,
    sample_class = "GCIMSSample",
    keep_intermediate = keep_intermediate
  )
}

test_that("DelayedDatasetDisk dumps a list of samples to disk on construction", {
  ds <- make_disk_dataset()

  expect_equal(ds$sampleNames, c("a", "b"))
  files <- list.files(ds$getCurrentDir())
  expect_setequal(files, c("sample_a.rds", "sample_b.rds"))

  s <- ds$getSample("a", dataset = list())
  expect_s4_class(s, "GCIMSSample")
})

test_that("GCIMSDataset(on_ram = FALSE) constructs a disk-backed dataset from files", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  pData <- data.frame(SampleID = "s1", FileName = basename(sample_file))
  scratch <- new_scratch_dir()

  ds <- GCIMSDataset$new(pData = pData, base_dir = dirname(sample_file), on_ram = FALSE, scratch_dir = scratch)

  expect_true(ds$is_on_disk())
  expect_null(ds$getCurrentDir()) # nothing realized yet

  s <- ds$getSample(1)
  expect_s4_class(s, "GCIMSSample")
  expect_false(is.null(ds$getCurrentDir()))
})

test_that("sample_rds_basenames builds sample_<name>.rds paths, with or without a prefix", {
  expect_equal(sample_rds_basenames(c("a", "b")), c("sample_a.rds", "sample_b.rds"))
  expect_equal(sample_rds_basenames("a", prefix_dir = "/some/dir"), "/some/dir/sample_a.rds")
})

test_that("checkSampleFiles reports missing sample files without aborting by default", {
  ds <- make_disk_dataset()
  unlink(file.path(ds$getCurrentDir(), "sample_a.rds"))

  missing <- ds$checkSampleFiles()

  expect_named(missing, "a")
})

test_that("checkSampleFiles(on_error = 'abort') aborts with a clear message when files are missing", {
  ds <- make_disk_dataset()
  unlink(file.path(ds$getCurrentDir(), "sample_a.rds"))

  expect_error(ds$checkSampleFiles(on_error = "abort"), "missing")
})

test_that("checkSampleFiles passes silently when all files exist", {
  ds <- make_disk_dataset()

  expect_equal(ds$checkSampleFiles(), character(0))
})

test_that("scratchDir is read-only and directs the user to updateScratchDir()", {
  ds <- make_disk_dataset()

  expect_true(dir.exists(ds$scratchDir))
  expect_error(
    {
      ds$scratchDir <- "/tmp/foo"
    },
    "updateScratchDir"
  )
})

test_that("updateScratchDir copies sample files to the new scratch directory", {
  ds <- make_disk_dataset()
  new_scratch <- new_scratch_dir()

  ds$updateScratchDir(new_scratch, dataset = list())

  expect_equal(ds$scratchDir, new_scratch)
  expect_true(file.exists(file.path(ds$getCurrentDir(), "sample_a.rds")))
})

test_that("the sampleNames setter renames sample files on disk", {
  ds <- make_disk_dataset()
  old_dir <- ds$getCurrentDir()

  ds$sampleNames <- c("x", "y")

  expect_equal(ds$sampleNames, c("x", "y"))
  expect_setequal(list.files(old_dir), c("sample_x.rds", "sample_y.rds"))
})

test_that("the sampleNames setter validates length and uniqueness", {
  ds <- make_disk_dataset()

  expect_error(
    {
      ds$sampleNames <- "only_one"
    },
    "length"
  )
  expect_error(
    {
      ds$sampleNames <- c("z", "z")
    },
    "unique and not missing"
  )
})

test_that("subset() keeps only the requested samples and deletes the others' files", {
  ds <- make_disk_dataset()
  dir <- ds$getCurrentDir()

  ds$subset("a")

  expect_equal(ds$sampleNames, "a")
  expect_true(file.exists(file.path(dir, "sample_a.rds")))
  expect_false(file.exists(file.path(dir, "sample_b.rds")))
})

test_that("subset() refuses to run with pending operations", {
  ds <- make_disk_dataset(n = 1)
  ds$appendDelayedOp(DelayedOperation(name = "op", fun = function(x) x))

  expect_error(ds$subset("a"), "pending operations")
})

test_that("realize() runs the queued operation, moves samples to a new directory, and cleans up the old one when keep_intermediate = FALSE", {
  ds <- make_disk_dataset(n = 1, keep_intermediate = FALSE)
  dir_before <- ds$getCurrentDir()
  ds$appendDelayedOp(DelayedOperation(
    name = "double-intensity",
    fun = function(sample) {
      intensity(sample) <- intensity(sample) * 2
      sample
    }
  ))

  ds$realize(dataset = list())

  dir_after <- ds$getCurrentDir()
  expect_false(identical(dir_before, dir_after))
  expect_false(dir.exists(dir_before))
  s <- ds$getSample("a", dataset = list())
  expect_true(all(intensity(s) == 2))
})

test_that("realize() keeps the previous directory when keep_intermediate = TRUE", {
  ds <- make_disk_dataset(n = 1, keep_intermediate = TRUE)
  dir_before <- ds$getCurrentDir()
  ds$appendDelayedOp(DelayedOperation(name = "noop", fun = function(x) x))

  ds$realize(dataset = list())

  expect_true(dir.exists(dir_before))
})
