make_peaks_df <- function(unique_id, sample_id, dt_apex, rt_apex, dt_half = 0.05, rt_half = 2) {
  data.frame(
    UniqueID = unique_id,
    SampleID = sample_id,
    dt_apex_ms = dt_apex,
    rt_apex_s = rt_apex,
    dt_min_ms = dt_apex - dt_half,
    dt_max_ms = dt_apex + dt_half,
    rt_min_s = rt_apex - rt_half,
    rt_max_s = rt_apex + rt_half
  )
}

test_that("clusterPeaks groups the same analyte across samples into one cluster", {
  peaks <- make_peaks_df(
    unique_id = c("S1/A", "S1/B", "S1/C", "S2/A", "S2/B", "S2/C"),
    sample_id = c("S1", "S1", "S1", "S2", "S2", "S2"),
    dt_apex = c(1, 1.5, 2, 1.02, 1.48, 2.03),
    rt_apex = c(30, 60, 90, 31, 59, 91)
  )

  out <- clusterPeaks(peaks)
  pl <- out$peak_list_clustered

  # Each analyte's two peaks (one per sample) share a cluster, and different
  # analytes get different clusters:
  cluster_of <- function(uid) pl$cluster[pl$UniqueID == uid]
  expect_equal(cluster_of("S1/A"), cluster_of("S2/A"))
  expect_equal(cluster_of("S1/B"), cluster_of("S2/B"))
  expect_equal(cluster_of("S1/C"), cluster_of("S2/C"))
  expect_equal(length(unique(pl$cluster)), 3L)
})

test_that("clusterPeaks names clusters Cluster1, Cluster2, ... ordered by dt then rt apex", {
  peaks <- make_peaks_df(
    unique_id = c("S1/A", "S2/A", "S1/B", "S2/B"),
    sample_id = c("S1", "S2", "S1", "S2"),
    dt_apex = c(2, 2.01, 1, 1.01), # cluster around dt=2 defined first in the data, but should be named 2nd
    rt_apex = c(90, 91, 30, 31)
  )

  out <- clusterPeaks(peaks)
  cluster_stats <- out$cluster_stats

  expect_equal(cluster_stats$cluster[order(cluster_stats$dt_apex_ms)], c("Cluster1", "Cluster2"))
})

test_that("clusterPeaks computes fixedsize_* bounds shared by all peaks in a cluster", {
  peaks <- make_peaks_df(
    unique_id = c("S1/A", "S2/A"),
    sample_id = c("S1", "S2"),
    dt_apex = c(1, 1.02),
    rt_apex = c(30, 31),
    dt_half = c(0.05, 0.07),
    rt_half = c(2, 3)
  )

  out <- clusterPeaks(peaks)
  pl <- out$peak_list_clustered

  # The fixedsize window is a constant WIDTH (derived from the cluster's
  # median apex-to-boundary distances) centered on each peak's own apex, so
  # the width is shared across the cluster even though the absolute bounds
  # differ (each peak has a slightly different apex position):
  widths_dt <- pl$fixedsize_dt_max_ms - pl$fixedsize_dt_min_ms
  widths_rt <- pl$fixedsize_rt_max_s - pl$fixedsize_rt_min_s
  expect_equal(widths_dt[pl$UniqueID == "S1/A"], widths_dt[pl$UniqueID == "S2/A"])
  expect_equal(widths_rt[pl$UniqueID == "S1/A"], widths_rt[pl$UniqueID == "S2/A"])
})

test_that("clusterPeaks excludes peaks whose only merge exceeds the cluster spread ('breaks')", {
  peaks <- make_peaks_df(
    unique_id = c("S1/A", "S2/A", "S1/X", "S2/X"),
    sample_id = c("S1", "S2", "S1", "S2"),
    dt_apex = c(1, 1.01, 5, 5.01),
    rt_apex = c(30, 31, 30, 30 + 25) # 25s apart > default rt_cluster_spread_s = 20
  )

  out <- clusterPeaks(peaks)
  pl <- out$peak_list_clustered

  expect_setequal(pl$UniqueID, c("S1/A", "S2/A"))
})

test_that("clusterPeaks returns gracefully (not crash) when no peaks form a valid cluster", {
  peaks <- make_peaks_df(
    unique_id = c("S1/X", "S2/X"),
    sample_id = c("S1", "S2"),
    dt_apex = c(1, 1.01),
    rt_apex = c(30, 30 + 25) # exceeds rt_cluster_spread_s, and it's the only pair
  )

  out <- suppressWarnings(clusterPeaks(peaks))

  expect_equal(nrow(out$peak_list_clustered), 0L)
  expect_equal(nrow(out$cluster_stats), 0L)
})

test_that("clusterPeaks: same-sample peaks CAN still end up in the same cluster if physically close (known limitation)", {
  # The same-sample distance inflation only delays when same-sample peaks
  # merge in the dendrogram; the final cluster assignment is decided purely
  # by the physical dt/rt spread of the largest still-valid merge. If two
  # peaks from the same sample are close enough that even their (forced-last)
  # merge stays under the spread threshold, they get assigned the same
  # cluster anyway, despite the code comment's intent that this "never"
  # happens. This test documents today's actual behavior.
  peaks <- make_peaks_df(
    unique_id = c("S1/A", "S1/B", "S2/A"),
    sample_id = c("S1", "S1", "S2"),
    dt_apex = c(1, 1.001, 1.0005),
    rt_apex = c(30, 30.1, 30.05)
  )

  out <- clusterPeaks(peaks)
  pl <- out$peak_list_clustered

  cluster_of <- function(uid) pl$cluster[pl$UniqueID == uid]
  expect_equal(cluster_of("S1/A"), cluster_of("S1/B"))
})

test_that("clusterPeaks errors on an unsupported clustering method", {
  peaks <- make_peaks_df(c("S1/A", "S2/A"), c("S1", "S2"), c(1, 1.02), c(30, 31))

  expect_error(
    clusterPeaks(peaks, clustering = list(method = "bogus")),
    "Unsupported clustering method"
  )
})

test_that("clusterPeaks supports kmedoids clustering with an explicit number of clusters", {
  skip_if_not_installed("cluster")
  peaks <- make_peaks_df(
    unique_id = c("S1/A", "S1/B", "S2/A", "S2/B"),
    sample_id = c("S1", "S1", "S2", "S2"),
    dt_apex = c(1, 1.5, 1.02, 1.48),
    rt_apex = c(30, 60, 31, 59)
  )

  out <- clusterPeaks(peaks, clustering = list(method = "kmedoids", Nclusters = 2))

  expect_equal(length(unique(out$peak_list_clustered$cluster)), 2L)
})

test_that("clusterPeaks supports kmedoids with Nclusters = 'max_peaks_sample'", {
  skip_if_not_installed("cluster")
  peaks <- make_peaks_df(
    unique_id = c("S1/A", "S1/B", "S2/A", "S2/B"),
    sample_id = c("S1", "S1", "S2", "S2"),
    dt_apex = c(1, 1.5, 1.02, 1.48),
    rt_apex = c(30, 60, 31, 59)
  )

  out <- clusterPeaks(peaks, clustering = list(method = "kmedoids", Nclusters = "max_peaks_sample"))

  expect_equal(length(unique(out$peak_list_clustered$cluster)), 2L)
})

test_that("clusterPeaks errors when kmedoids Nclusters is neither numeric nor 'max_peaks_sample'", {
  skip_if_not_installed("cluster")
  peaks <- make_peaks_df(c("S1/A", "S2/A"), c("S1", "S2"), c(1, 1.02), c(30, 31))

  expect_error(
    clusterPeaks(peaks, clustering = list(method = "kmedoids", Nclusters = "bogus")),
    "Nclusters"
  )
})

test_that("set_peak_distances_within_groups sets within-group distances and leaves the rest untouched", {
  m <- matrix(
    c(0, 1, 2, 3, 1, 0, 4, 5, 2, 4, 0, 6, 3, 5, 6, 0),
    nrow = 4,
    dimnames = list(c("a", "b", "c", "d"), c("a", "b", "c", "d"))
  )
  d <- stats::as.dist(m)

  out <- as.matrix(set_peak_distances_within_groups(d, peak_groups = list(c("a", "b")), value = 99))

  expect_equal(out["a", "b"], 99)
  expect_equal(out["b", "a"], 99)
  expect_equal(out["a", "a"], 0)
  expect_equal(out["c", "d"], 6) # untouched
})

test_that("peak2peak_distance computes euclidean distances", {
  pm <- matrix(c(0, 0, 3, 4), nrow = 2, byrow = TRUE, dimnames = list(c("p1", "p2"), c("x", "y")))

  d <- peak2peak_distance(pm, distance_method = "euclidean")

  expect_equal(as.numeric(d), 5)
})

test_that("peak2peak_distance errors on an unsupported distance method", {
  pm <- matrix(c(0, 0, 3, 4), nrow = 2, byrow = TRUE, dimnames = list(c("p1", "p2"), c("x", "y")))

  expect_error(
    peak2peak_distance(pm, distance_method = "bogus"),
    "Unsupported distance"
  )
})
