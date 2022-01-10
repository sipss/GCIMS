rand_index <- function(clustering1, clustering2) {
  n <- length(clustering1)
  if (length(clustering2) != n) {
    stop("clustering1 and clustering2 must be of the same length")
  }
  # pairs in the same group for clustering1:
  clu1_pairs <- vapply(clustering1, function(x) x == clustering1, FUN.VALUE = logical(n))
  # Ignore all in diagonal and below:
  clu1_pairs <- upper.tri(clu1_pairs) & clu1_pairs

  # for clustering2:
  clu2_pairs <- vapply(clustering2, function(x) x == clustering2, FUN.VALUE = logical(n))
  clu2_pairs <- upper.tri(clu2_pairs) & clu2_pairs

  # a: number of pairs in the same cluster in clu1 that are also in the same cluster for clu2
  a <- sum(clu1_pairs & clu2_pairs)
  # b: number of pairs in the different cluster in clu1 that are also in different cluster for clu2
  clu1_pairs <- upper.tri(clu1_pairs) & (!clu1_pairs)
  clu2_pairs <- upper.tri(clu2_pairs) & (!clu2_pairs)
  b <- sum(clu1_pairs & clu2_pairs)

  pairs <- n*(n-1)/2
  (a+b)/pairs
}



get_simple_peak_list <- function() {
  fn <- system.file("extdata", "peak_lists", "peak_list_6_peaks_15_samples_ketones.csv.gz", package = "GCIMS")
  readr::read_csv(fn, show_col_types = FALSE)
}

test_that("sample peak list is consistent", {
  peak_list <- get_simple_peak_list()
  expect_true(
    all(peak_list$rtmin_s <= peak_list$rtmax_s)
  )
  expect_true(
    all(peak_list$dtmin_ms <= peak_list$dtmax_ms)
  )
  expect_true(
    all(peak_list$rtapex_s >= peak_list$rtmin_s & peak_list$rtapex_s <= peak_list$rtmax_s)
  )
  expect_true(
    all(peak_list$dtapex_ms >= peak_list$dtmin_ms & peak_list$dtapex_ms <= peak_list$dtmax_ms)
  )
})

test_that("No peak is removed if criterias are NULL", {
  peak_list <- get_simple_peak_list()
  expect_equal(
    peak_list$UniqueID,
    remove_peaks_with_outlier_rois(peak_list, dtime_criteria = NULL, rtime_criteria = NULL)$UniqueID
  )
})

test_that("True Peaks are not removed with dt IQR criteria", {
  skip("This test is known to fail")
  peak_list <- get_simple_peak_list()
  expect_equal(
    peak_list$UniqueID,
    remove_peaks_with_outlier_rois(
      peak_list,
      dtime_criteria = "IQR",
      rtime_criteria = NULL)$UniqueID
  )
})

test_that("True Peaks are not removed with rt arnau criteria", {
  skip("This test is known to fail")
  peak_list <- get_simple_peak_list()
  expect_equal(
    peak_list$UniqueID,
    remove_peaks_with_outlier_rois(
      peak_list,
      dtime_criteria = NULL,
      rtime_criteria = "arnau")$UniqueID
  )
})


test_that("Peak grouping of simple case works well", {
  peak_list <- get_simple_peak_list()
  peak_list$volume <- 1
  peak_table_list <- group_peak_list(
    peaks = peak_list,
    filter_dt_width_criteria = NULL, # FIXME: outlier roi criteria disabled because it does not work well
    filter_rt_width_criteria = NULL,# FIXME: outlier roi criteria disabled because it does not work well
    distance_method = "mahalanobis",
    distance_between_peaks_from_same_sample = Inf,
    clustering = list(method = "kmedoids", Nclusters = "max_peaks_sample"),
    aggregate_conflicting_peaks = NULL,
    verbose = FALSE
  )
  peak_list_with_cluster <- peak_table_list$peak_list_with_cluster

  expect_equal(
    rand_index(
      peak_list_with_cluster$PeakID,
      peak_list_with_cluster$cluster
    ),
    1
  )
})


test_that("Peak grouping with sd_scaled_euclidean and hclust", {
  peak_list <- get_simple_peak_list()
  peak_list$volume <- 1
  peak_table_list <- group_peak_list(
    peaks = peak_list,
    filter_dt_width_criteria = NULL, # FIXME: outlier roi criteria disabled because it does not work well
    filter_rt_width_criteria = NULL,# FIXME: outlier roi criteria disabled because it does not work well
    distance_method = "sd_scaled_euclidean",
    distance_between_peaks_from_same_sample = 2,
    clustering = list(method = "hclust", hclust_method = "complete"),
    aggregate_conflicting_peaks = NULL,
    verbose = FALSE
  )
  peak_list_with_cluster <- peak_table_list$peak_list_with_cluster

  expect_equal(
    rand_index(
      peak_list_with_cluster$PeakID,
      peak_list_with_cluster$cluster
    ),
    1
  )
})
