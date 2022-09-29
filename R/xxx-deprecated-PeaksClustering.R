
peak_and_cluster_metrics_old <- function(peaks) {
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
  list(peaks = peaks, cluster_stats = cluster_stats)
}


#' Group peaks in clusters
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
#' @description Peak grouping function, exposing several options useful for benchmarking.
#'
#' @return A list with :
#' - peak_table: A peak table that includes peak position, median peak minimum/maximum retention and drift times and the peak Volume for each sample
#' - peak_table_duplicity: How many Volume values have been aggregated. Should be 1 for each sample/peak
#' - extra_clustering_info: Arbitrary clustering extra information, that depends on the clustering method
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
#'}
#' @export
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
    warn("Some peaks in samples have wrong indexes leading to NA positions")
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
    num_clusters <- clustering$num_clusters
    num_cluster_estimation <- NULL
    if (is.null(num_clusters)) {
      dt_ms_max_dist_thres <- clustering$dt_ms_max_dist_thres
      if (is.null(dt_ms_max_dist_thres)) {
        dt_ms_max_dist_thres <- signif(3*stats::median(peaks$dt_max_ms - peaks$dt_min_ms), digits = 2)
        if (verbose) {
          rlang::inform(c("i" = glue("The maximum distance between two peaks in the same cluster is of {dt_ms_max_dist_thres} ms")))
        }
      }
      rt_s_max_dist_thres <- clustering$rt_s_max_dist_thres
      if (is.null(rt_s_max_dist_thres)) {
        rt_s_max_dist_thres <- signif(3*stats::median(peaks$rt_max_s - peaks$rt_min_s), digits = 2)
        if (verbose) {
          rlang::inform(c("i" = glue("The maximum distance between two peaks in the same cluster is of {rt_s_max_dist_thres} s")))
        }
      }

      num_cluster_estimation <- estimate_num_clusters(
        peak_list = peaks,
        cluster = cluster,
        dt_ms_max_dist_thres = dt_ms_max_dist_thres,
        rt_s_max_dist_thres = rt_s_max_dist_thres
      )
      num_clusters <- num_cluster_estimation$num_clusters
    }
    peaks$cluster <- stats::cutree(cluster, k = num_clusters)
    extra_clustering_info$cluster <- cluster
    extra_clustering_info$num_clusters <- num_clusters
    extra_clustering_info$num_cluster_estimation <- num_cluster_estimation
  } else {
    stop(sprintf("Unsupported clustering method %s", clustering$method))
  }

  # Turn numeric peak clusters into IDs:
  if (is.numeric(peaks$cluster)) {
    ndigits_print <- paste0("Cluster%0", nchar(max(peaks$cluster)), "d")
    peaks$cluster <- sprintf(ndigits_print, peaks$cluster)
  }

  peak_cluster_stats <- peak_and_cluster_metrics_old(peaks)

  list(
    peak_list_clustered = peak_cluster_stats$peaks,
    cluster_stats = peak_cluster_stats$cluster_stats,
    dist = peak2peak_dist,
    extra_clustering_info = extra_clustering_info
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
    abort(sprintf("Unknown dtime_criteria: %s", dtime_criteria))
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
    abort(sprintf("Unknown rtime_criteria: %s", rtime_criteria))
  }
  if (verbose) {
    message(sprintf("Excluding %d/%d peaks", length(peaks_to_exclude), nrow(peaks)))
  }
  dplyr::filter(peaks, !.data$UniqueID %in% peaks_to_exclude)
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



