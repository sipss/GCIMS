gaussian_signal <- function(x, centers, heights, sd) {
  y <- rep(0, length(x))
  for (i in seq_along(centers)) y <- y + heights[i] * exp(-(x - centers[i])^2 / (2 * sd^2))
  y
}

make_spectrum <- function() {
  dt <- seq(0, 20, by = 0.01)
  y <- gaussian_signal(dt, centers = 10, heights = 100, sd = 0.3)
  GCIMSSpectrum(drift_time = dt, intensity = y, description = "S1")
}

make_chromatogram <- function() {
  rt <- seq(0, 20, by = 0.01)
  y <- gaussian_signal(rt, centers = 10, heights = 100, sd = 0.3)
  GCIMSChromatogram(retention_time = rt, intensity = y, description = "S1")
}

test_that("description()/description<- work on a GCIMSChromatogram and validate their input", {
  ch <- make_chromatogram()

  description(ch) <- "myChrom"

  expect_equal(description(ch), "myChrom")
  expect_error(
    {
      description(ch) <- 5
    },
    "description should be a string"
  )
})

test_that("description()/description<- work on a GCIMSSpectrum and validate their input", {
  sp <- make_spectrum()

  description(sp) <- "mySpec"

  expect_equal(description(sp), "mySpec")
  expect_error(
    {
      description(sp) <- 5
    },
    "description should be a string"
  )
})

test_that("peaks() on a GCIMSSpectrum errors before findPeaks() has been run", {
  sp <- make_spectrum()

  expect_error(peaks(sp), "findPeaks")
})

test_that("findPeaks() on a GCIMSSpectrum detects the peak and stamps it with the sample's description", {
  sp <- make_spectrum()

  out <- findPeaks(sp, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))
  p <- peaks(out)

  expect_equal(nrow(p), 1L)
  expect_equal(p$SampleID, "S1")
  expect_equal(p$apex, 10, tolerance = 0.02)
})

test_that("peaks<- on a GCIMSSpectrum with an empty data frame gives an empty peak list", {
  sp <- make_spectrum()
  out <- findPeaks(sp, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))
  p <- peaks(out)

  peaks(out) <- p[0, ]

  expect_null(peaks(out))
})

test_that("findPeaks() on a GCIMSChromatogram detects the peak and stamps it with the sample's description", {
  ch <- make_chromatogram()

  out <- findPeaks(ch, length_in_xunits = 0.07, peakwidth_range_xunits = c(0.3, 3))
  p <- peaks(out)

  expect_equal(nrow(p), 1L)
  expect_equal(p$SampleID, "S1")
  expect_equal(p$apex, 10, tolerance = 0.02)
})

test_that("findPeaks() on a GCIMSDataset queues peak detection applied on realize, aggregating into the peak table", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  annot <- data.frame(SampleID = "Sample1", FileName = basename(sample_file))
  ds <- GCIMSDataset$new(annot, base_dir = dirname(sample_file), on_ram = TRUE)

  findPeaks(ds)
  p <- peaks(ds)

  expect_true(nrow(p) > 0)
  expect_true(all(p$SampleID == "Sample1"))
})
