make_impute_fixtures <- function() {
  cluster_stats <- data.frame(
    cluster = c("Cluster1", "Cluster2"),
    dt_min_ms = c(8, 10),
    dt_max_ms = c(9, 12),
    rt_min_s = c(120, 300),
    rt_max_s = c(128, 320)
  )

  peak_table <- matrix(NA_real_, nrow = 2, ncol = 2)
  rownames(peak_table) <- c("Sample1", "Sample2")
  colnames(peak_table) <- c("Cluster1", "Cluster2")
  peak_table["Sample1", "Cluster2"] <- 9.5
  peak_table["Sample2", "Cluster1"] <- 3.6

  dt <- seq(from = 0, to = 13, by = 0.1)
  rt <- seq(from = 0, to = 350, by = 1)
  set.seed(42)
  s1_intensity <- matrix(
    rnorm(length(dt) * length(rt), sd = 0.1),
    nrow = length(dt),
    ncol = length(rt)
  )
  s2_intensity <- matrix(
    rnorm(length(dt) * length(rt), sd = 0.1),
    nrow = length(dt),
    ncol = length(rt)
  )
  # Plateau in the Sample1/Cluster1 region so the imputed value stands out
  s1_intensity[dt > 8.25 & dt < 8.75, rt > 122 & rt < 126] <- 1

  s1 <- GCIMSSample(drift_time = dt, retention_time = rt, data = s1_intensity)
  s2 <- GCIMSSample(drift_time = dt, retention_time = rt, data = s2_intensity)
  dataset <- GCIMSDataset_fromList(list(Sample1 = s1, Sample2 = s2))

  list(peak_table = peak_table, cluster_stats = cluster_stats, dataset = dataset)
}

test_that("imputePeakTable() fills NA cells by integrating the cluster region and leaves existing values untouched", {
  fx <- make_impute_fixtures()

  out <- imputePeakTable(fx$peak_table, fx$dataset, fx$cluster_stats)

  expect_equal(out["Sample1", "Cluster2"], 9.5)
  expect_equal(out["Sample2", "Cluster1"], 3.6)
  expect_equal(out["Sample1", "Cluster1"], 1.476875, tolerance = 1e-6)
  expect_equal(out["Sample2", "Cluster2"], -0.2335684, tolerance = 1e-6)
})

test_that("imputePeakTable() preserves dimensions and dimnames", {
  fx <- make_impute_fixtures()

  out <- imputePeakTable(fx$peak_table, fx$dataset, fx$cluster_stats)

  expect_equal(dim(out), dim(fx$peak_table))
  expect_equal(dimnames(out), dimnames(fx$peak_table))
})

test_that("imputePeakTable() is a no-op when the peak table has no missing values", {
  fx <- make_impute_fixtures()
  pt_full <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("Sample1", "Sample2"), c("Cluster1", "Cluster2"))
  )

  out <- imputePeakTable(pt_full, fx$dataset, fx$cluster_stats)

  expect_identical(out, pt_full)
})

test_that("imputePeakTable() imputes multiple missing values in the same sample row", {
  fx <- make_impute_fixtures()
  pt_all_na <- fx$peak_table
  pt_all_na[, ] <- NA_real_

  out <- imputePeakTable(pt_all_na, fx$dataset, fx$cluster_stats)

  expect_false(anyNA(out))
  expect_equal(out["Sample1", "Cluster1"], 1.476875, tolerance = 1e-6)
  expect_equal(out["Sample1", "Cluster2"], 0.1385071728, tolerance = 1e-6)
  expect_equal(out["Sample2", "Cluster1"], 0.1079987948, tolerance = 1e-6)
  expect_equal(out["Sample2", "Cluster2"], -0.2335684497, tolerance = 1e-6)
})
