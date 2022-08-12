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
#' @examples
#' \donttest{
#' dir_in <- system.file("extdata", package = "GCIMS")
#' peak_list <- readRDS(file.path(dir_in, "peak_list.rds"))
#'
#' peak_clustering  <- group_peak_list(
#'   peaks = peak_list,
#'   filter_dt_width_criteria = NULL,
#'   filter_rt_width_criteria = NULL,
#'   distance_method = "mahalanobis",
#'   distance_between_peaks_from_same_sample = Inf,
#'   clustering = list(method = "kmedoids", Nclusters = "max_peaks_sample"),
#'   verbose = FALSE
#' )
#' }
group_peak_list <- function(
  peaks,
  filter_dt_width_criteria = "IQR",
  filter_rt_width_criteria = "arnau",
  distance_method = "mahalanobis",
  distance_between_peaks_from_same_sample = 100,
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
    cluster <- cluster::pam(x = peak2peak_dist, k = N_clusters, pamonce = 3)
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

  # Turn numeric peak clusters into IDs
  if (is.numeric(peaks$cluster)) {
    ndigits_print <- paste0("Cluster%0", nchar(max(peaks$cluster)), "d")
    peaks$cluster <- sprintf(ndigits_print, peaks$cluster)
  }

  cluster_stats <- peaks %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(
          c(
            "dt_apex_ms", "dt_min_ms", "dt_max_ms", "dt_cm_ms",
            "rt_apex_s", "rt_min_s", "rt_max_s", "rt_cm_s"
          )
        ),
        stats::median
      ),
      dplyr::across(
        dplyr::all_of(c("dt_min_idx", "rt_min_idx")),
        ~ floor(stats::median(.))
      ),
      dplyr::across(
        dplyr::all_of(c("dt_max_idx", "rt_max_idx")),
        ~ ceiling(stats::median(.))
      ),
      dplyr::across(
        dplyr::all_of(c("dt_apex_idx", "dt_cm_idx", "rt_apex_idx", "rt_cm_idx")),
        ~ round(stats::median(.))
      ),
    ) %>%
    dplyr::ungroup()

  list(
    peak_list_clustered = peaks,
    cluster_stats = cluster_stats,
    dist = peak2peak_dist,
    extra_clustering_info = extra_clustering_info
  )
}

#' Build a peak table
#'
#' @param peak_list_clustered The output of [gcims_figures_of_merit].
#' @param aggregate_conflicting_peaks `NULL` or a function. What to do, in case two peaks from the same sample
#' have been assigned to the same cluster. If `NULL`, throw an error. If `mean`, `max` or any other function,
#' we will summarize all the conflicting volumes into that number (e.g. "take the maximum of the peaks")
#'
#' @return A list with the peak table and the ROI duplicity information. The peak table
#' is a data frame with
#' @export
#' @examples
#' \donttest{
#'
#' # Create your peak table from scratch:
#' pl <- data.frame(
#'   SampleID = c("S1", "S1", "S2", "S2"),
#'   cluster = c(1, 2, 1, 2),
#'   Volume = c(10, 20, 8, 18)
#' )
#' build_peak_table(pl)
#'
#' # Create a peak table from the output of the function gcims_figure_of_merit()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' peak_list_fom <- readRDS(file.path(dir_in, "peak_list_fom.rds"))
#' peak_table <- build_peak_table(peak_list_fom, aggregate_conflicting_peaks = max)
#'
#' peak_table$peak_table_matrix
#' }
#'
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

  peak_table_mat <- peak_table %>%
    tidyr::pivot_longer(cols = -1, names_to = "SampleID", values_to = "Volume") %>%
    tidyr::pivot_wider(names_from = "cluster", values_from = "Volume") %>%
    tibble::column_to_rownames("SampleID") %>%
    as.matrix()


  # Missing values still need to be filled
  list(
    peak_table = peak_table,
    peak_table_matrix = peak_table_mat,
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



#' Identify ROIs our of the Cluster
#'
#' This function receives a clustering object and identifies, for each cluster,
#' if there is one ROI that do not belong to the cluster.
#' @noRd
#'
#' @param clustering_object The result of the clustering function.
#' @return A clustering object where, in the peak_list, some samples are identified
#' to check on them the clustering process.
#' @examples
#' peak_clustering$peak_list_clustered <- data.frame(
#' SampleID = c("S1", "S2", "S3", "S4"),
#' PeakID = c(1, 2, 1, 3),
#' dt_apex_ms = c(7, 7.2, 7.1, 10.2),
#' rt_apex_s = c(30, 33, 35, 255),
#' dt_min_ms = c(6.5, 6.9, 6.8, 9.9),
#' dt_max_ms = c(7.7, 8.3, 8.2, 12.1),
#' rt_min_s = c(27, 28, 30, 250),
#' rt_max_s = c(36, 36, 37, 265),
#' dt_cm_ms = c(7, 7.2, 7.2, 10.2),
#' rt_cm_s = c(30, 33, 35, 256),
#' dt_apex_idx = c(140, 160, 180, 395),
#' rt_apex_idx = c(30, 33, 33, 247),
#' dt_min_idx = c(100, 110, 110, 350),
#' dt_max_idx = c(200, 200, 205, 450),
#' rt_min_idx = c(27, 28, 27, 255),
#' rt_max_idx = c(36, 36, 37, 267),
#' dt_cm_idx = c(140, 150, 160, 395),
#' rt_cm_idx = c(30, 33, 33, 247),
#' uniqueID = c("S1", "S2", "S3", "S4"),
#' cluster = c(1, 1, 1, 1))
#'
removing_outsiders <- function(clustering_object){
  clusters_infor <- clustering_object$peak_list_clustered
  num_clusts <- unique(clusters_infor$cluster)
  threshold <- 0.6
  for(i in num_clusts){
    # 4 corrdinates
    ROIs2Fusion <- which(clusters_infor$cluster == i)
    ioverumatriz <- matrix(0, nrow = length(ROIs2Fusion), ncol = length(ROIs2Fusion))
    r <- 0
    for (j in ROIs2Fusion){
      r <- r + 1
      R1 <- clusters_infor[j, c(13:16)]
      c <- 0
      for (k in ROIs2Fusion){
        c <- c + 1
        R2 <- clusters_infor[k, c(13:16)]
        thrOverlap <- round(as.numeric(intersection_over_union(R1, R2)), digits = 2)
        ioverumatriz[r,c] <- thrOverlap
      }
    }
    median_overlaps <- apply(ioverumatriz, 1, medianvalues)
    samples_out_cluster <- which(median_overlaps <= threshold)
    samples_to_recluster <- ROIs2Fusion[samples_out_cluster]
    clusters_infor$cluster[samples_to_recluster] <- "Re-Cluster"
  }
  clustering_object$peak_list_clustered <- clusters_infor
  return(clustering_object)
}


medianvalues <- function(x) {
  if(length(x) > 1){
    sort(x)[round(length(x)/2)]
  } else {
    x
  }
}

intersection_over_union <- function(ROI1, ROI2){
  if (!(ROI1[1] >= ROI2[2] | ROI1[2] <= ROI2[1] | ROI1[3] >= ROI2[4] | ROI1[4] <= ROI2[3])){
    area1 <- (ROI1[2] - ROI1[1])*(ROI1[4] - ROI1[3])
    area2 <- (ROI2[2] - ROI2[1])*(ROI2[4] - ROI2[3])
    x_left <- max(ROI1[1], ROI2[1])
    y_top <- min(ROI1[4], ROI2[4])
    x_right <- min(ROI1[2], ROI2[2])
    y_bottom <- max(ROI1[3], ROI2[3])
    overlapping_area <- (x_right - x_left)*(y_top - y_bottom)
    if (area1 == overlapping_area | area2 == overlapping_area){
      p <- 1
    } else {
      p <- overlapping_area / (area1 + area2 - overlapping_area)
    }
  } else {
    p <- 0
  }
  return(p)
}



