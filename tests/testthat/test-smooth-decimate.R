gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

make_sample <- function() {
  dt <- seq(0, 4, by = 0.02) # 201 pts
  rt <- seq(0, 50, by = 0.5) # 101 pts
  int_mat <- matrix(50, nrow = length(dt), ncol = length(rt)) +
    outer(gauss(dt, 2, 500, 0.1), gauss(rt, 25, 1, 3))
  GCIMSSample(drift_time = dt, retention_time = rt, data = int_mat)
}

test_that("smooth() on a GCIMSSample applies a Savitzky-Golay filter without changing dimensions", {
  s <- make_sample()

  out <- smooth(s, dt_length_ms = 0.1, rt_length_s = 2)

  expect_equal(dim(out), dim(s))
  expect_false(identical(intensity(out), intensity(s)))
})

test_that("smooth() on a GCIMSChromatogram filters the intensity vector in place", {
  s <- make_sample()
  ch <- getChromatogram(s)

  out <- smooth(ch, rt_length_s = 2)

  expect_equal(length(intensity(out)), length(intensity(ch)))
  expect_false(identical(intensity(out), intensity(ch)))
})

test_that("smooth() on a GCIMSSpectrum filters the intensity vector in place", {
  s <- make_sample()
  sp <- getSpectrum(s)

  out <- smooth(sp, dt_length_ms = 0.1)

  expect_equal(length(intensity(out)), length(intensity(sp)))
  expect_false(identical(intensity(out), intensity(sp)))
})

test_that("smooth() on a GCIMSDataset queues a delayed operation applied on realize", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  annot <- data.frame(SampleID = "Sample1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(annot, base_dir = dirname(sample_file), on_ram = TRUE)
  before <- intensity(ds$getSample(1))

  smooth(ds, dt_length_ms = 0.1, rt_length_s = 2)
  after <- intensity(ds$getSample(1))

  expect_equal(dim(after), dim(before))
  expect_false(identical(after, before))
})

test_that("decimate() on a GCIMSSample keeps one every n points", {
  s <- make_sample()

  out <- decimate(s, dt_factor = 2, rt_factor = 2)

  expect_equal(dim(out), c(length(seq(1, 201, by = 2)), length(seq(1, 101, by = 2))))
})

test_that("decimate() with factors of 1 is a no-op", {
  s <- make_sample()

  out <- decimate(s, dt_factor = 1L, rt_factor = 1L)

  expect_identical(out, s)
})

test_that("decimate() rejects non-positive factors", {
  s <- make_sample()

  expect_error(decimate(s, dt_factor = 0), "positive integers")
  expect_error(decimate(s, rt_factor = -1), "positive integers")
})

test_that("decimate() on a GCIMSDataset queues a delayed operation applied on realize", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  annot <- data.frame(SampleID = "Sample1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(annot, base_dir = dirname(sample_file), on_ram = TRUE)
  dims_before <- dim(ds$getSample(1))

  decimate(ds, dt_factor = 2, rt_factor = 2)
  dims_after <- dim(ds$getSample(1))

  expect_true(dims_after[1] < dims_before[1])
  expect_true(dims_after[2] < dims_before[2])
})

test_that("decimate() on a GCIMSDataset with factors of 1 is a no-op that returns the dataset unchanged", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  annot <- data.frame(SampleID = "Sample1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(annot, base_dir = dirname(sample_file), on_ram = TRUE)
  ds$realize()
  had_pending <- ds$hasDelayedOps()

  out <- decimate(ds, dt_factor = 1L, rt_factor = 1L)

  expect_identical(out, ds)
  expect_equal(ds$hasDelayedOps(), had_pending)
})

test_that("decimate() on a GCIMSDataset rejects non-positive factors", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  annot <- data.frame(SampleID = "Sample1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(annot, base_dir = dirname(sample_file), on_ram = TRUE)

  expect_error(decimate(ds, dt_factor = 0), "positive integers")
})
