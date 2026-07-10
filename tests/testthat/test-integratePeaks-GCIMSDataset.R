test_that("integratePeaks on a GCIMSDataset integrates each sample and combines them with a UniqueID", {
  dt <- seq(0, 4, by = 0.02)
  rt <- seq(0, 100, by = 0.5)
  make_sample <- function() {
    int_mat <- matrix(0.1, nrow = length(dt), ncol = length(rt))
    int_mat[150:160, 40:50] <- int_mat[150:160, 40:50] + 50
    GCIMSSample(drift_time = dt, retention_time = rt, data = int_mat)
  }
  ds <- GCIMSDataset$new_from_list(
    samples = list(s1 = make_sample(), s2 = make_sample()),
    on_ram = TRUE, scratch_dir = NULL
  )

  peak_list <- data.frame(
    SampleID = c("s1", "s2"),
    PeakID = c("1", "1"),
    dt_min_ms = dt[150], dt_max_ms = dt[160],
    rt_min_s = rt[40], rt_max_s = rt[50],
    rt_apex_s = rt[45], rt_cm_s = rt[45],
    fixedsize_dt_min_ms = dt[148], fixedsize_dt_max_ms = dt[162],
    fixedsize_rt_min_s = rt[38], fixedsize_rt_max_s = rt[52]
  )

  integratePeaks(ds, peak_list, integration_size_method = "free_size")
  ds$realize()
  pl <- peaks(ds)

  expect_equal(nrow(pl), 2L)
  expect_equal(sort(pl$UniqueID), c("s1/1", "s2/1"))
  expect_equal(pl$Volume[pl$SampleID == "s1"], pl$Volume[pl$SampleID == "s2"])
})
