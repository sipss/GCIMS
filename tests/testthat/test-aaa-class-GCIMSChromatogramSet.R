make_chrom <- function(sample_id) {
  GCIMSChromatogram(
    retention_time = 1:5,
    intensity = (1:5) * 10,
    description = sample_id
  )
}

test_that("GCIMSChromatogramSet stores chromatograms and pData", {
  chroms <- list(s1 = make_chrom("s1"), s2 = make_chrom("s2"))
  pd <- data.frame(SampleID = c("s1", "s2"), Group = c("A", "B"))

  cs <- GCIMSChromatogramSet(chromatograms = chroms, pData = pd)

  expect_s4_class(cs, "GCIMSChromatogramSet")
  expect_equal(sampleNames(cs), c("s1", "s2"))
  expect_equal(length(cs), 2L)
  expect_equal(pData(cs)$Group, c("A", "B"))
  expect_identical(cs[["s1"]], chroms$s1)
})

test_that("GCIMSChromatogramSet can be built with no pData", {
  chroms <- list(s1 = make_chrom("s1"))

  cs <- GCIMSChromatogramSet(chromatograms = chroms)

  expect_null(pData(cs))
  expect_equal(sampleNames(cs), "s1")
})

test_that("GCIMSChromatogramSet can be empty", {
  cs <- GCIMSChromatogramSet()

  expect_equal(length(cs), 0L)
  expect_equal(sampleNames(cs), character(0))
})

test_that("GCIMSChromatogramSet rejects an unnamed list of chromatograms", {
  chroms <- list(make_chrom("s1"), make_chrom("s2"))

  expect_error(GCIMSChromatogramSet(chromatograms = chroms), "named list")
})

test_that("GCIMSChromatogramSet rejects elements that are not GCIMSChromatogram objects", {
  chroms <- list(s1 = make_chrom("s1"), s2 = "not a chromatogram")

  expect_error(GCIMSChromatogramSet(chromatograms = chroms), "GCIMSChromatogram objects")
})
