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

test_that("GCIMSChromatogramSet rejects pData whose SampleID doesn't match the chromatogram names", {
  chroms <- list(s1 = make_chrom("s1"), s2 = make_chrom("s2"))

  # Different sample entirely:
  pd_mismatched <- data.frame(SampleID = c("s1", "other"), Group = c("A", "B"))
  expect_error(GCIMSChromatogramSet(chromatograms = chroms, pData = pd_mismatched), "does not match")

  # Missing a sample:
  pd_missing <- data.frame(SampleID = "s1", Group = "A")
  expect_error(GCIMSChromatogramSet(chromatograms = chroms, pData = pd_missing), "does not match")

  # No SampleID column at all:
  pd_no_sampleid <- data.frame(Group = c("A", "B"))
  expect_error(GCIMSChromatogramSet(chromatograms = chroms, pData = pd_no_sampleid), "SampleID column")
})

test_that("sampleNames() follows pData's row order, even when it differs from the chromatograms list order", {
  chroms <- list(s1 = make_chrom("s1"), s2 = make_chrom("s2"))
  pd_reversed <- data.frame(SampleID = c("s2", "s1"), Group = c("B", "A"))

  cs <- GCIMSChromatogramSet(chromatograms = chroms, pData = pd_reversed)

  expect_equal(sampleNames(cs), c("s2", "s1"))
  # [[ is unaffected by pData order, it always looks up by name:
  expect_identical(cs[["s1"]], chroms$s1)
  expect_identical(cs[["s2"]], chroms$s2)
})
