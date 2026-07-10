test_that("peakTable pivots a clustered peak list into a sample x cluster matrix", {
  pl <- data.frame(
    SampleID = c("S1", "S1", "S2", "S2"),
    cluster = c("Cluster1", "Cluster2", "Cluster1", "Cluster2"),
    Volume = c(10, 20, 8, 18)
  )

  out <- peakTable(pl)

  expect_equal(
    out$peak_table_matrix,
    matrix(
      c(10, 8, 20, 18),
      nrow = 2, ncol = 2,
      dimnames = list(c("S1", "S2"), c("Cluster1", "Cluster2"))
    )
  )
})

test_that("peakTable excludes peaks that were never assigned to a cluster", {
  pl <- data.frame(
    SampleID = c("S1", "S1", "S2"),
    cluster = c("Cluster1", NA, "Cluster1"),
    Volume = c(10, 99, 8)
  )

  out <- peakTable(pl)

  expect_equal(colnames(out$peak_table_matrix), "Cluster1")
  expect_equal(unname(out$peak_table_matrix[, "Cluster1"]), c(10, 8))
})

test_that("peakTable requires a Volume column", {
  expect_error(
    peakTable(data.frame(SampleID = "S1", cluster = "C1")),
    "Volume"
  )
})

test_that("peakTable errors when two peaks from the same sample land in the same cluster and no aggregate function is given", {
  pl <- data.frame(SampleID = c("S1", "S1"), cluster = c("C1", "C1"), Volume = c(5, 7))

  expect_error(
    peakTable(pl),
    "more than one peak assigned to the same cluster"
  )
})

test_that("peakTable resolves conflicting peaks with aggregate_conflicting_peaks", {
  pl <- data.frame(SampleID = c("S1", "S1"), cluster = c("C1", "C1"), Volume = c(5, 7))

  out <- peakTable(pl, aggregate_conflicting_peaks = max)

  expect_equal(unname(out$peak_table_matrix[, "C1"]), 7)
  expect_equal(out$peak_table_duplicity$S1, 2L)
})

test_that("omit_times drops ROIs at the given drift times", {
  peak_list <- data.frame(
    rt_apex_s = c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6),
    dt_apex_ms = c(2, 4, 6, 4, 8, 4, 10, 4, 4, 12)
  )

  out <- omit_times(peak_list, dt_time_2_omit = 4)

  expect_equal(nrow(out), 5L)
  expect_false(any(out$dt_apex_ms == 4))
})

test_that("omit_times drops ROIs at the given retention times", {
  peak_list <- data.frame(
    rt_apex_s = c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6),
    dt_apex_ms = c(2, 4, 6, 4, 8, 4, 10, 4, 4, 12)
  )

  out <- omit_times(peak_list, rt_time_2_omit = c(3, 5))

  expect_equal(nrow(out), 6L)
  expect_false(any(out$rt_apex_s %in% c(3, 5)))
})

test_that("omit_times with no times to omit returns the peak list unchanged", {
  peak_list <- data.frame(
    rt_apex_s = c(1, 2, 3),
    dt_apex_ms = c(2, 4, 6)
  )

  out <- omit_times(peak_list)

  expect_equal(out, peak_list)
})
