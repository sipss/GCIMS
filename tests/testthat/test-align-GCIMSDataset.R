gauss <- function(rt, center, height, sd) height * exp(-(rt - center)^2 / (2 * sd^2))

# A GCIMSSample whose single informative drift-time row (the RIP row) has
# both an injection-point spike (early in rt, for align_ip) and an
# analyte-elution dip (later in rt, for the ptw retention-time warp),
# mirroring how getRIC()'s max(row)-row inversion turns a raw dip into a
# RIC peak and a raw spike into a RIC minimum.
make_rip_sample <- function(dt, rt, rip_idx, injection_rt, analyte_rt) {
  data <- matrix(1, nrow = length(dt), ncol = length(rt))
  data[rip_idx, ] <- 1000 + gauss(rt, injection_rt, 2000, 1.5) - gauss(rt, analyte_rt, 900, 3)
  GCIMSSample(drift_time = dt, retention_time = rt, data = data)
}

make_dataset <- function() {
  dt <- seq(1, 10, by = 0.5)
  rt <- 1:80
  s1 <- make_rip_sample(dt, rt, rip_idx = 10, injection_rt = 5, analyte_rt = 40)
  s2 <- make_rip_sample(dt, rt, rip_idx = 8, injection_rt = 5, analyte_rt = 34)
  s3 <- make_rip_sample(dt, rt, rip_idx = 12, injection_rt = 5, analyte_rt = 46)
  GCIMSDataset$new_from_list(
    samples = list(s1 = s1, s2 = s2, s3 = s3),
    on_ram = TRUE,
    scratch_dir = NULL
  )
}

test_that("align() on a GCIMSDataset aligns retention time and drift time across samples", {
  ds <- make_dataset()

  align(ds, method_rt = "ptw", align_dt = TRUE, align_ip = TRUE)
  ds$realize()

  samples <- lapply(ds$sampleNames, ds$getSample)
  for (s in samples) {
    expect_false(anyNA(intensity(s)))
  }

  # align_ip assigns the very same rt_ref vector to every sample:
  rt_lists <- lapply(samples, rtime)
  expect_true(all(vapply(rt_lists[-1], identical, logical(1), rt_lists[[1]])))

  # Drift-time correction was actually computed (samples had different RIP
  # positions), not left at the trivial default:
  expect_false(all(ds$align$dt_kcorr == 1))
  expect_named(ds$align$dt_kcorr, ds$sampleNames)
})

test_that("align() honors an explicit reference_sample_idx", {
  ds <- make_dataset()

  align(ds, method_rt = "ptw", align_dt = TRUE, align_ip = TRUE, reference_sample_idx = 2)
  ds$realize()

  samples <- lapply(ds$sampleNames, ds$getSample)
  for (s in samples) {
    expect_false(anyNA(intensity(s)))
  }
  rt_lists <- lapply(samples, rtime)
  expect_true(all(vapply(rt_lists[-1], identical, logical(1), rt_lists[[1]])))
})

test_that("align() with align_dt = FALSE, align_ip = FALSE skips prealign without crashing", {
  # Registering SerialParam (instead of the default forking MulticoreParam)
  # runs the per-sample delayed operation in this same process, so covr can
  # instrument .align_fun_extract()'s dt_kcorr fallback below -- it would
  # otherwise run invisibly inside a forked worker.
  old_bpparam <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(old_bpparam))

  ds <- make_dataset()
  dt_before <- dtime(ds$getSample(1))

  align(ds, method_rt = "ptw", align_dt = FALSE, align_ip = FALSE)
  ds$realize()

  expect_identical(dtime(ds$getSample(1)), dt_before)
  # .align_fun_extract() falls back to dt_kcorr = 1 when prealign() (and
  # therefore alignDt()) never ran, since no drift-time correction was
  # ever recorded on the sample:
  expect_true(all(ds$align$dt_kcorr == 1))
})

test_that("align(method_rt = 'pow') on a GCIMSDataset fails fast with a clear message", {
  skip_if(requireNamespace("pow", quietly = TRUE), "pow is installed; not testing the missing-dependency guard")

  ds <- make_dataset()

  expect_error(
    align(ds, method_rt = "pow"),
    "pow.*is not on CRAN/Bioconductor"
  )
})

test_that("alignPlots() returns rt/dt correction plots and a dt_kcorr plot with the expected data", {
  ds <- make_dataset()
  align(ds, method_rt = "ptw", align_dt = TRUE, align_ip = TRUE)
  ds$realize()

  plots <- alignPlots(ds)

  expect_named(plots, c("rt_diff_plot", "dt_diff_plot", "dt_kcorr_plot"))
  for (p in plots) expect_s3_class(p, "ggplot")

  expect_equal(plots$rt_diff_plot$labels$x, "Retention time (s)")
  expect_equal(plots$rt_diff_plot$labels$y, "Ret. time correction (s)")
  expect_setequal(unique(as.character(plots$rt_diff_plot$data$SampleID)), ds$sampleNames)

  expect_equal(plots$dt_diff_plot$labels$x, "Drift time (ms)")
  expect_equal(plots$dt_diff_plot$labels$y, "Drift time correction (ms)")
  expect_setequal(unique(as.character(plots$dt_diff_plot$data$SampleID)), ds$sampleNames)

  expect_equal(plots$dt_kcorr_plot$labels$x, "SampleID")
  expect_equal(as.character(plots$dt_kcorr_plot$data$x), ds$sampleNames)
  expect_equal(plots$dt_kcorr_plot$data$y, unname(ds$align$dt_kcorr))
})
