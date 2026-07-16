make_dataset_mismatched_axes <- function() {
  # s1 and s2 deliberately have retention/drift time axes of different
  # length and range, to prove that getChromatogram()/getSpectrum() on a
  # GCIMSDataset never interpolates samples onto a common grid: each
  # extracted chromatogram/spectrum keeps its own sample's native axis.
  s1 <- GCIMSSample(
    drift_time = seq(1, 5, by = 1),
    retention_time = seq(1, 10, by = 1),
    data = matrix(seq_len(5 * 10), nrow = 5, ncol = 10)
  )
  s2 <- GCIMSSample(
    drift_time = seq(1, 5, by = 1),
    retention_time = seq(1, 20, by = 2),
    data = matrix(seq_len(5 * 10), nrow = 5, ncol = 10)
  )
  pd <- data.frame(SampleID = c("s1", "s2"), Group = c("A", "B"))
  GCIMSDataset$new_from_list(samples = list(s1 = s1, s2 = s2), pData = pd, on_ram = TRUE, scratch_dir = NULL)
}

test_that("getChromatogram(dataset) returns one GCIMSChromatogram per sample, each on its own retention time axis", {
  ds <- make_dataset_mismatched_axes()

  cs <- getChromatogram(ds, dt_range = c(2, 4))

  expect_s4_class(cs, "GCIMSChromatogramSet")
  expect_equal(sampleNames(cs), c("s1", "s2"))
  expect_equal(rtime(cs[["s1"]]), 1:10)
  expect_equal(rtime(cs[["s2"]]), seq(1, 20, by = 2))
  expect_equal(unname(intensity(cs[["s1"]])), unname(getChromatogram(ds$getSample("s1"), dt_range = c(2, 4))@intensity))
})

test_that("getChromatogram(dataset) attaches pData(dataset)", {
  ds <- make_dataset_mismatched_axes()

  cs <- getChromatogram(ds)

  expect_equal(pData(cs)$SampleID, c("s1", "s2"))
  expect_equal(pData(cs)$Group, c("A", "B"))
})

test_that("getSpectrum(dataset) returns one GCIMSSpectrum per sample, each on its own drift time axis", {
  ds <- make_dataset_mismatched_axes()

  ss <- getSpectrum(ds, rt_range = c(2, 5))

  expect_s4_class(ss, "GCIMSSpectrumSet")
  expect_equal(sampleNames(ss), c("s1", "s2"))
  expect_equal(dtime(ss[["s1"]]), 1:5)
  expect_equal(dtime(ss[["s2"]]), 1:5)
})

test_that("plot(GCIMSChromatogramSet) colors by SampleID by default and combines mismatched axes without error", {
  ds <- make_dataset_mismatched_axes()
  cs <- getChromatogram(ds)

  p <- plot(cs)

  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(unique(p$data$SampleID)), c("s1", "s2"))
  expect_equal(range(p$data$retention_time_s[p$data$SampleID == "s2"]), c(1, 19))
})

test_that("plot(GCIMSChromatogramSet) attributes each line to the right sample even when pData's row order differs from the chromatograms list order", {
  chroms <- list(
    s1 = GCIMSChromatogram(retention_time = 1:3, intensity = c(100, 100, 100)),
    s2 = GCIMSChromatogram(retention_time = 1:3, intensity = c(1, 1, 1))
  )
  pd_reversed <- data.frame(SampleID = c("s2", "s1"), Group = c("B", "A"))
  cs <- GCIMSChromatogramSet(chromatograms = chroms, pData = pd_reversed)

  p <- plot(cs)

  expect_equal(unique(p$data$intensity[p$data$SampleID == "s1"]), 100)
  expect_equal(unique(p$data$intensity[p$data$SampleID == "s2"]), 1)
})

test_that("plot(GCIMSChromatogramSet, color_by=) colors by an arbitrary pData column", {
  ds <- make_dataset_mismatched_axes()
  cs <- getChromatogram(ds)

  p <- plot(cs, color_by = "Group")

  expect_equal(sort(unique(p$data$Group)), c("A", "B"))
})

test_that("plot(GCIMSChromatogramSet, color_by=) errors clearly for an unknown column", {
  ds <- make_dataset_mismatched_axes()
  cs <- getChromatogram(ds)

  expect_error(plot(cs, color_by = "NotAColumn"), "not a column")
})

test_that("plot(GCIMSSpectrumSet) colors by SampleID by default", {
  ds <- make_dataset_mismatched_axes()
  ss <- getSpectrum(ds)

  p <- plot(ss)

  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(unique(p$data$SampleID)), c("s1", "s2"))
})

test_that("plot() errors clearly on an empty set", {
  cs <- GCIMSChromatogramSet()
  expect_error(plot(cs), "empty")

  ss <- GCIMSSpectrumSet()
  expect_error(plot(ss), "empty")
})
