gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

make_sample <- function() {
  dt <- seq(0, 4, by = 0.02) # 201 pts
  rt <- seq(0, 50, by = 0.5) # 101 pts
  int_mat <- matrix(50, nrow = length(dt), ncol = length(rt)) +
    outer(gauss(dt, 2, 500, 0.1), gauss(rt, 25, 1, 3))
  GCIMSSample(drift_time = dt, retention_time = rt, data = int_mat)
}

test_that("estimateBaseline() on a GCIMSChromatogram estimates a flat baseline near the known offset", {
  s <- make_sample()
  ch <- getChromatogram(s)

  out <- estimateBaseline(ch, rt_length_s = 10)

  expect_equal(length(baseline(out)), length(intensity(out)))
  expect_equal(baseline(out)[1], intensity(out)[1], tolerance = 5)
})

test_that("baseline() on a GCIMSChromatogram errors before estimateBaseline() has been run", {
  s <- make_sample()
  ch <- getChromatogram(s)

  expect_error(baseline(ch), "estimateBaseline")
  expect_null(baseline(ch, .error_if_missing = FALSE))
})

test_that("baseline<- on a GCIMSChromatogram rejects a value with the wrong length", {
  s <- make_sample()
  ch <- getChromatogram(s)

  expect_error(
    {
      baseline(ch) <- 1:3
    },
    "same length as the intensity"
  )
})

test_that("estimateBaseline() on a GCIMSSpectrum estimates a flat baseline near the known offset", {
  s <- make_sample()
  sp <- getSpectrum(s)

  out <- estimateBaseline(sp, dt_peak_fwhm_ms = 0.3, dt_region_multiplier = 6)

  expect_equal(length(baseline(out)), length(intensity(out)))
  expect_equal(baseline(out)[1], intensity(out)[1], tolerance = 5)
})

test_that("baseline() on a GCIMSSpectrum errors before estimateBaseline() has been run", {
  s <- make_sample()
  sp <- getSpectrum(s)

  expect_error(baseline(sp), "estimateBaseline")
  expect_null(baseline(sp, .error_if_missing = FALSE))
})

test_that("baseline<- on a GCIMSSpectrum rejects a value with the wrong length", {
  s <- make_sample()
  sp <- getSpectrum(s)

  expect_error(
    {
      baseline(sp) <- 1:3
    },
    "same length as the intensity"
  )
})

test_that("estimateBaseline() on a GCIMSDataset queues a delayed operation applied on realize", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  annot <- data.frame(SampleID = "Sample1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(annot, base_dir = dirname(sample_file), on_ram = TRUE)

  estimateBaseline(ds, dt_peak_fwhm_ms = 0.3, dt_region_multiplier = 6, rt_length_s = 10)
  s <- ds$getSample(1)

  expect_false(is.null(baseline(s, .error_if_missing = FALSE)))
})
