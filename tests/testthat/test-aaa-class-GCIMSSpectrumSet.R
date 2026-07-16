make_spec <- function(sample_id) {
  GCIMSSpectrum(
    drift_time = 1:5,
    intensity = (1:5) * 10,
    description = sample_id
  )
}

test_that("GCIMSSpectrumSet stores spectra and pData", {
  spectra <- list(s1 = make_spec("s1"), s2 = make_spec("s2"))
  pd <- data.frame(SampleID = c("s1", "s2"), Group = c("A", "B"))

  ss <- GCIMSSpectrumSet(spectra = spectra, pData = pd)

  expect_s4_class(ss, "GCIMSSpectrumSet")
  expect_equal(sampleNames(ss), c("s1", "s2"))
  expect_equal(length(ss), 2L)
  expect_equal(pData(ss)$Group, c("A", "B"))
  expect_identical(ss[["s1"]], spectra$s1)
})

test_that("GCIMSSpectrumSet can be built with no pData", {
  spectra <- list(s1 = make_spec("s1"))

  ss <- GCIMSSpectrumSet(spectra = spectra)

  expect_null(pData(ss))
  expect_equal(sampleNames(ss), "s1")
})

test_that("GCIMSSpectrumSet can be empty", {
  ss <- GCIMSSpectrumSet()

  expect_equal(length(ss), 0L)
  expect_equal(sampleNames(ss), character(0))
})

test_that("GCIMSSpectrumSet rejects an unnamed list of spectra", {
  spectra <- list(make_spec("s1"), make_spec("s2"))

  expect_error(GCIMSSpectrumSet(spectra = spectra), "named list")
})

test_that("GCIMSSpectrumSet rejects elements that are not GCIMSSpectrum objects", {
  spectra <- list(s1 = make_spec("s1"), s2 = "not a spectrum")

  expect_error(GCIMSSpectrumSet(spectra = spectra), "GCIMSSpectrum objects")
})

test_that("GCIMSSpectrumSet rejects pData whose SampleID doesn't match the spectra names", {
  spectra <- list(s1 = make_spec("s1"), s2 = make_spec("s2"))

  # Different sample entirely:
  pd_mismatched <- data.frame(SampleID = c("s1", "other"), Group = c("A", "B"))
  expect_error(GCIMSSpectrumSet(spectra = spectra, pData = pd_mismatched), "does not match")

  # Missing a sample:
  pd_missing <- data.frame(SampleID = "s1", Group = "A")
  expect_error(GCIMSSpectrumSet(spectra = spectra, pData = pd_missing), "does not match")

  # No SampleID column at all:
  pd_no_sampleid <- data.frame(Group = c("A", "B"))
  expect_error(GCIMSSpectrumSet(spectra = spectra, pData = pd_no_sampleid), "SampleID column")
})

test_that("sampleNames() follows pData's row order, even when it differs from the spectra list order", {
  spectra <- list(s1 = make_spec("s1"), s2 = make_spec("s2"))
  pd_reversed <- data.frame(SampleID = c("s2", "s1"), Group = c("B", "A"))

  ss <- GCIMSSpectrumSet(spectra = spectra, pData = pd_reversed)

  expect_equal(sampleNames(ss), c("s2", "s1"))
  # [[ is unaffected by pData order, it always looks up by name:
  expect_identical(ss[["s1"]], spectra$s1)
  expect_identical(ss[["s2"]], spectra$s2)
})
