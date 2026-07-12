gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

make_dataset <- function() {
  dt <- seq(0, 4, by = 0.02) # 201 pts
  rt <- seq(0, 50, by = 0.5) # 101 pts
  make_s <- function(peak_dt, peak_rt) {
    int_mat <- matrix(50, nrow = length(dt), ncol = length(rt)) +
      outer(gauss(dt, peak_dt, 500, 0.1), gauss(rt, peak_rt, 1, 3))
    GCIMSSample(drift_time = dt, retention_time = rt, data = int_mat)
  }
  s1 <- make_s(2, 25)
  s2 <- make_s(2.2, 25)
  GCIMSDataset$new_from_list(samples = list(s1 = s1, s2 = s2), on_ram = TRUE, scratch_dir = NULL)
}

test_that("getTIS()/getRIC() extract per-sample TIS and RIC matrices, keyed by SampleID", {
  # Registering SerialParam runs the per-sample delayed operation
  # (.extract_RIC_and_TIS_fun_extract) in this same process, so covr can
  # instrument it -- it would otherwise run invisibly inside a forked worker.
  old_bpparam <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(old_bpparam))

  ds <- make_dataset()

  tis <- getTIS(ds)
  ric <- getRIC(ds)

  expect_equal(dim(tis), c(2L, 201L))
  expect_equal(dim(ric), c(2L, 101L))
  expect_equal(dimnames(tis)[[1]], c("s1", "s2"))
})

test_that("plotTIS() returns a ggplot with one line per sample over the full drift time range by default", {
  ds <- make_dataset()

  p <- plotTIS(ds)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Drift time (ms)")
  expect_equal(p$labels$y, "TIS Intensity (a.u.)")
  expect_setequal(as.character(unique(p$layers[[1]]$data$SampleID)), c("s1", "s2"))
  expect_equal(range(p$layers[[1]]$data$drift_time_ms), c(0, 4))
})

test_that("plotTIS() restricts to dt_range and a single sample when given", {
  ds <- make_dataset()

  p <- plotTIS(ds, dt_range = c(1, 3), sample = "s1")

  expect_equal(as.character(unique(p$layers[[1]]$data$SampleID)), "s1")
  expect_equal(range(p$layers[[1]]$data$drift_time_ms), c(1, 3))
})

test_that("plotRIC() returns a ggplot with one line per sample over the full retention time range by default", {
  ds <- make_dataset()

  p <- plotRIC(ds)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Retention time (s)")
  expect_equal(p$labels$y, "RIC Intensity (a.u.)")
  expect_setequal(as.character(unique(p$layers[[1]]$data$SampleID)), c("s1", "s2"))
})

test_that("plotRIC() restricts to rt_range and a single sample when given", {
  ds <- make_dataset()

  p <- plotRIC(ds, rt_range = c(10, 30), sample = 1)

  expect_equal(as.character(unique(p$layers[[1]]$data$SampleID)), "s1")
  expect_equal(range(p$layers[[1]]$data$retention_time_s), c(10, 30))
})
