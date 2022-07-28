#' Peak grouping function, exposing a lot of options useful for benchmarking
#'
#' @param peaks A data frame with at least the following columns:
#'  - "UniqueID" A unique ID for each peak
#'  - "SampleID" The sample ID the peak belongs to
#'  - "dt_apex_ms", "rt_apex_s" The peak positions
#'  - "dt_max_ms", "dt_min_ms", "rt_max_s", "rt_min_s" (for filtering outlier peaks based on their size)
#' @param filter_dt_width_criteria,filter_rt_width_criteria A character with the method for outlier detection.
#'   - "IQR": Remove peaks with widths more than 1.5 interquartile ranges above upper quartile or
#'     below the lower quartile.
#'   - "arnau": FIXME Adhoc method from Arnau, where he removes peaks with widths above mean+4iqr or below median-0.75iqr
#'   - "none": Do not remove peaks based on their drift time width or retention time height
#' @param distance_method A string. One of the distance methods from [stats::dist], "sd_scaled_euclidean" or "mahalanobis"
#' @param distance_between_peaks_from_same_sample The distance between two peaks from the same sample will be set to `distance_between_peaks_from_same_sample*max(distance_matrix)`
#' @param clustering A named list with "method" and the supported method, as well as further options.
#'   For `method = "kmedoids"`, you must provide `Nclusters`, with either the number of clusters
#'   to use in the kmedoids algorithm ([cluster::pam]) or the string `"max_peaks_sample"` to use the maximum number of
#'   detected peaks per sample.
#'
#'   For `method = "hclust"`, you can provide `hclust_method`, with the `method` passed to [stats::hclust].
#' @param verbose logical, to control printing in the function
#'
#' @return A list with :
#' - peak_table: A peak table that includes peak position, median peak minimum/maximum retention and drift times and the peak Volume for each sample
#' - peak_table_duplicity: How many Volume values have been aggregated. Should be 1 for each sample/peak
#' - extra_clustering_info: Arbitrary clustering extra information, that depends on the clustering method
#' @export
#'
#' @examples
#' peak_list <- data.frame(
#'   UniqueID = c("P1", "P2", "P3", "P4"),
#'   SampleID = c("S1", "S1", "S2", "S2"),
#'   dt_apex_ms = c(7, 10, 7.1, 10.2),
#'   rt_apex_s = c(30, 250, 33, 247),
#'   dt_min_ms = c(6.5, 9.4, 6.6, 9.7),
#'   dt_max_ms = c(7.7, 10.8, 7.6, 11.1),
#'   rt_min_s = c(27, 246, 30, 245),
#'   rt_max_s = c(36, 260, 37, 255)
#' )
#' peak_table_list <- group_peak_list(
#'   peaks = peak_list,
#'   filter_dt_width_criteria = NULL,
#'   filter_rt_width_criteria = NULL,
#'   distance_method = "mahalanobis",
#'   distance_between_peaks_from_same_sample = Inf,
#'   clustering = list(method = "kmedoids", Nclusters = "max_peaks_sample"),
#'   verbose = FALSE
#' )
group_peak_list <- function(
  peaks,
  filter_dt_width_criteria = "IQR",
  filter_rt_width_criteria = "arnau",
  distance_method = "mahalanobis",
  distance_between_peaks_from_same_sample = Inf,
  clustering = list(method = "kmedoids", Nclusters = "max_peaks_sample"),
  verbose = FALSE
) {
  # 0. Warn if peaks with NA positions, and remove them
  peaks_with_na <- stats::complete.cases(peaks[,c("UniqueID", "SampleID", "dt_apex_ms", "rt_apex_s")])
  if (!all(peaks_with_na)) {
    rlang::warn("Some peaks in samples have wrong indexes leading to NA positions")
    print(peaks[!peaks_with_na,])
    peaks <- peaks[peaks_with_na,]
  }

  # 1. Filter peaks with weird width or height
  peaks <- remove_peaks_with_outlier_rois(
    peaks,
    dtime_criteria = filter_dt_width_criteria,
    rtime_criteria = filter_rt_width_criteria,
    verbose = verbose
  )

  # Compute the peak to peak distance:
  peak_matrix <- as.matrix(peaks[,c("dt_apex_ms", "rt_apex_s")])
  rownames(peak_matrix) <- peaks$UniqueID

  peak2peak_dist <- peak2peak_distance(
    peak_matrix = peak_matrix,
    distance_method = distance_method
  )

  # Set distances from pairs of peaks belonging to the same sample to Inf,
  # so they are never in the same cluster
  peakuids_by_sample <- peaks %>%
    dplyr::select(dplyr::all_of(c("SampleID", "UniqueID"))) %>%
    dplyr::group_by(.data$SampleID) %>%
    dplyr::summarize(UniqueIDs = list(.data$UniqueID)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.data$UniqueIDs)

  peak2peak_dist <- set_peak_distances_within_groups(
    dist_matrix = peak2peak_dist,
    peak_groups = peakuids_by_sample,
    value = distance_between_peaks_from_same_sample*max(peak2peak_dist)
  )

  extra_clustering_info <- list()
  if (clustering$method == "kmedoids") {
    if (clustering$Nclusters == "max_peaks_sample") {
      N_clusters <- max(purrr::map_int(peakuids_by_sample, length))
    } else if (is.numeric(clustering$Nclusters)) {
      N_clusters <- clustering$Nclusters
    } else {
      stop("When clustering$method is kmedoids, clustering$Nclusters must be an integer or the string 'max_peaks_sample'")
    }
    cluster <- cluster::pam(x = peak2peak_dist, k = N_clusters)
    peaks$cluster <- cluster$clustering
    extra_clustering_info$cluster_result <- cluster
    extra_clustering_info$silhouette <- cluster::silhouette(cluster, dist = peak2peak_dist)
  } else if (clustering$method == "hclust") {
    hclust_method <- ifelse(is.null(clustering$hclust_method), "complete", clustering$hclust_method)
    cluster <- stats::hclust(d = peak2peak_dist, method = hclust_method)
    # FIXME: Implement something more robust to estimate num_clusters or height to cut:
    num_clusters <- which.max(-diff(sort(cluster$height, decreasing = TRUE)))+1
    peaks$cluster <- stats::cutree(cluster, k = num_clusters)
    extra_clustering_info$cluster <- cluster
    extra_clustering_info$num_clusters <- num_clusters
  } else {
    stop(sprintf("Unsupported clustering method %s", clustering$method))
  }

  median_roi_per_cluster <- peaks %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarise(
      dplyr::across(
        c(dplyr::starts_with("dt"), dplyr::starts_with("rt")),
        stats::median
      )
    ) %>%
    dplyr::ungroup()


  list(
    peak_list_clustered = peaks,
    clusters_median_roi = median_roi_per_cluster,
    dist = peak2peak_dist,
    extra_clustering_info = extra_clustering_info
  )
}

#' Build a peak table
#'
#' @param peak_list_clustered The output of [group_peak_list]
#' @param aggregate_conflicting_peaks `NULL` or a function. When we build the peak table, with peaks in rows, samples in
#'  columns, `peak_table[i,j]` is the volume of the peak from sample `j` in cluster `i`. If the clustering process
#'  clusters together two peaks form the same sample, those peaks will conflict in the peak table. `NULL` will error
#'  in that case, another function will be applied on the conflicting volumes (e.g `mean` or `max` would be reasonable options)
#'
#'
#' @return A list with the peak table and the ROI duplicity information
#' @export
#' @examples
#' pl <- data.frame(
#'   SampleID = c("S1", "S1", "S2", "S2"),
#'   cluster = c(1, 2, 1, 2),
#'   Volume = c(10, 20, 8, 18)
#' )
#' build_peak_table(pl)
build_peak_table <- function(peak_list_clustered, aggregate_conflicting_peaks = NULL) {
  if (!"Volume" %in% colnames(peak_list_clustered)) {
    rlang::abort("Please compute a 'Volume' column in peak_list_clustered")
  }

  peak_table_duplicity <- peak_list_clustered %>%
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "Volume"))) %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of("SampleID"),
      values_from = dplyr::all_of("Volume"),
      values_fn = length
    )

  peak_table <- peak_list_clustered %>%
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "Volume"))) %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of("SampleID"),
      values_from = dplyr::all_of("Volume"),
      values_fn = aggregate_conflicting_peaks
    )
  # Missing values still need to be filled
  list(
    peak_table = peak_table,
    peak_table_duplicity = peak_table_duplicity
  )
}


#' @noRd
#' @param peaks A data frame with one peak per row and at least the following columns:
#'   - UniqueID A unique peak name
#'   - dt_max_ms, dt_min_ms
#'   - rt_max_s, rt_min_s
#' @param dtime_criteria,rtime_criteria A character with the method for outlier detection.
#'   - "IQR": Remove peaks more than 1.5 interquartile ranges above upper quartile or
#'     below the lower quartile.
#'   - "arnau": FIXME Adhoc method from Arnau, where he removes peaks above mean+4iqr or below median-0.75iqr
#'   - "none": Do not remove peaks
#' @param verbose logical. if `TRUE`, a message is printed with the peaks excluded
#' FIXME: These methods remove peaks in a simple well behaved dataset. Check the tests.
remove_peaks_with_outlier_rois <- function(
  peaks,
  dtime_criteria = "IQR",
  rtime_criteria = "arnau",
  verbose = FALSE
  ) {

  # 1. Filter peaks with weird width or height
  roi_sizes <- tibble::tibble(
    UniqueID = peaks$UniqueID,
    dt_width_ms = peaks$dt_max_ms - peaks$dt_min_ms,
    rt_width_s = peaks$rt_max_s - peaks$rt_min_s
  )

  peaks_to_exclude <- c()

  if (is.null(dtime_criteria)) {
    # do nothing
  } else if (dtime_criteria == "IQR") {
    quartiles_dt <- stats::quantile(roi_sizes$dt_width_ms)
    iqr_dt <- quartiles_dt["75%"] - quartiles_dt["25%"]
    lower_bound_iqr_dt <- quartiles_dt["25%"] - 1.5*iqr_dt
    higher_bound_iqr_dt <- quartiles_dt["75%"] + 1.5*iqr_dt
    peaks_to_exclude <- c(
      peaks_to_exclude,
      roi_sizes$UniqueID[
        roi_sizes$dt_width_ms < lower_bound_iqr_dt | roi_sizes$dt_width_ms > higher_bound_iqr_dt
      ]
    )
  } else {
    rlang::abort(sprintf("Unknown dtime_criteria: %s", dtime_criteria))
  }


  if (is.null(rtime_criteria)) {
    # do nothing
  } else if (rtime_criteria == "IQR") {
    quartiles_rt <- stats::quantile(roi_sizes$rt_width_s)
    iqr_rt <- quartiles_rt["75%"] - quartiles_rt["25%"]
    lower_bound_iqr_rt <- quartiles_rt["25%"] - 1.5*iqr_rt
    higher_bound_iqr_rt <- quartiles_rt["75%"] + 1.5*iqr_rt
    peaks_to_exclude <- c(
      peaks_to_exclude,
      roi_sizes$UniqueID[
        roi_sizes$rt_width_s < lower_bound_iqr_rt | roi_sizes$rt_width_s > higher_bound_iqr_rt
      ]
    )
  } else if (rtime_criteria == "arnau") { # FIXME: This method is very adhoc and hard to justify
    quartiles_rt <- stats::quantile(roi_sizes$rt_width_s)
    median_rt <- quartiles_rt["50%"]
    iqr_rt <- quartiles_rt["75%"] - quartiles_rt["25%"]
    lower_bound_rt <- mean(roi_sizes$rt_width_s) - 0.75*iqr_rt
    higher_bound_rt <- median_rt + 4*iqr_rt
    peaks_to_exclude <- c(
      peaks_to_exclude,
      roi_sizes$UniqueID[
        roi_sizes$rt_width_s < lower_bound_rt | roi_sizes$rt_width_s > higher_bound_rt
      ]
    )
  }else {
    rlang::abort(sprintf("Unknown rtime_criteria: %s", rtime_criteria))
  }
  if (verbose) {
    message(sprintf("Excluding %d/%d peaks", length(peaks_to_exclude), nrow(peaks)))
  }
  dplyr::filter(peaks, ! .data$UniqueID %in% peaks_to_exclude)
}


#' Override peak distances to infinity
#'
#' This function receives a distance matrix and a list of peak groups. Each group
#' consists of peaks that should not be grouped as the same peak (for instance because
#' they belong to the same sample). For each group, we set the distance between
#' all its peaks to infinity.
#'
#' @noRd
#'
#' @param dist_matrix A square matrix, where `dist_matrix[i,j]` is the distance
#'  from peak `i` to peak `j`. The matrix must have as row names and column names
#'  unique peak names.
#' @param peak_groups A list, where each element is a character vector with peak names
#' @param value `Inf` by default, but you could set to any other value
#'
#' @return An object of class "dist". See [stats::dist].
#'
set_peak_distances_within_groups <- function(dist_matrix, peak_groups, value = Inf) {
  # Set distances from pairs of peaks belonging to the same sample to Inf,
  # so they are never in the same cluster
  dist_matrix <- as.matrix(dist_matrix)
  for (peak_ids in peak_groups) {
    for (peak_i in peak_ids) {
      dist_matrix[peak_i, peak_ids] <- value
      dist_matrix[peak_ids, peak_i] <- value
      dist_matrix[peak_i, peak_i] <- 0
    }
  }
  stats::as.dist(dist_matrix)
}


# Mahalanobis distance:
# https://stats.stackexchange.com/a/81710/62083
mahalanobis_distance <- function(x) {
  covmat <- stats::cov(x)
  dec <- chol(covmat)
  tmp <- forwardsolve(t(dec), t(x))
  colnames(tmp) <- rownames(x)
  stats::dist(t(tmp))
}


peak2peak_distance <- function(peak_matrix, distance_method = "mahalanobis") {
  STATS_METHODS <- c("euclidean", "maximum", "manhattan", "canberra",
                     "binary", "minkowski")
  if (distance_method %in% STATS_METHODS) {
    peak2peak_dist <- stats::dist(peak_matrix, method = distance_method)
  } else if (distance_method == "sd_scaled_euclidean") {
    peak_matrix_scaled <- scale(peak_matrix, center = FALSE, scale = TRUE)
    peak2peak_dist <- stats::dist(peak_matrix_scaled, method = "euclidean")
  } else if (distance_method == "mahalanobis") {
    peak2peak_dist <- mahalanobis_distance(peak_matrix)
  } else {
    stop(sprintf("Unsupported distance %s", distance_method))
  }
  peak2peak_dist
}
