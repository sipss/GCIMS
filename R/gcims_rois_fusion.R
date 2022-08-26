#' ROIs Fusion
#'
#' For each ROI in `peak_list_clustered`, we define additional columns with ROI boundaries.
#' The center of the ROI is not changed, but the limits are taken so the size of the ROI is
#' the same for all ROIs in the same cluster, using the median size of the ROIs in the cluster.
#'
#'
#' @param peak_list_clustered The peak list clustered data frame from [group_peak_list()]
#' @param cluster_stats The cluster statistics data frame from [group_peak_list()]
#' @return The peak list clustered table with added `ref_roi_[dt|rt]_[min|max]_idx` columns and the cluster statistics
#' @export
#' @examples
#' \donttest{
#' dir_in <- system.file("extdata", package = "GCIMS")
#' clustering <- readRDS(file.path(dir_in, "peak_clustering.rds"))
#'
#' roi_fusion_out <- gcims_rois_fusion(clustering$peak_list_clustered, clustering$cluster_stats)
#' samplestocheck <- roi_fusion_out$peak_list_clustered[
#'   which(roi_fusion_out$peak_list_clustered$cluster == "Cluster01"),
#' ]
#' library(ggplot2)
#' ggplot() +
#'  geom_rect(
#'    data = samplestocheck,
#'    mapping = aes(
#'      xmin = dt_min_idx,
#'      xmax = dt_max_idx,
#'      ymin = rt_min_idx,
#'      ymax = rt_max_idx,
#'      color = SampleID,
#'      fill = SampleID
#'    ),
#'    alpha=0.5
#' ) +
#' labs(x = "Drift time (index)", "Retention time (index)", title = "Cluster01")
#' }
#'
gcims_rois_fusion <- function(peak_list_clustered, cluster_stats, drift_time = NULL, retention_time = NULL) {
  peak_list_clustered$ref_roi_dt_min_idx <- 0
  peak_list_clustered$ref_roi_dt_max_idx <- 0
  peak_list_clustered$ref_roi_rt_min_idx <- 0
  peak_list_clustered$ref_roi_rt_max_idx <- 0

  if (!is.null(drift_time)) {
    dt_allmax_idx <- length(drift_time)
  } else {
    dt_allmax_idx <- max(c(peak_list_clustered$dt_max_idx, cluster_stats$dt_max_idx))
  }

  if (!is.null(retention_time)) {
    rt_allmax_idx <- length(retention_time)
  } else {
    rt_allmax_idx <- max(c(peak_list_clustered$rt_max_idx, cluster_stats$rt_max_idx))
  }

  for (i in seq_len(nrow(peak_list_clustered))) {
    roi_prop <- as.list(peak_list_clustered[i,])
    cluster_id <- roi_prop$cluster
    cluster_prop <- as.list(cluster_stats[cluster_stats$cluster == cluster_id,])

    dt_cluster_half_length <- ceiling((cluster_prop$dt_max_idx - cluster_prop$dt_min_idx)/2)
    rt_cluster_half_length <- ceiling((cluster_prop$rt_max_idx - cluster_prop$rt_min_idx)/2)

    dt_roi_center_idx <- floor((roi_prop$dt_max_idx + roi_prop$dt_min_idx)/2)
    rt_roi_center_idx <- floor((roi_prop$rt_max_idx + roi_prop$rt_min_idx)/2)

    peak_list_clustered$ref_roi_dt_min_idx[i] <- max(1L,            dt_roi_center_idx - dt_cluster_half_length)
    peak_list_clustered$ref_roi_dt_max_idx[i] <- min(dt_allmax_idx, dt_roi_center_idx + dt_cluster_half_length)
    peak_list_clustered$ref_roi_rt_min_idx[i] <- max(1L,            rt_roi_center_idx - rt_cluster_half_length)
    peak_list_clustered$ref_roi_rt_max_idx[i] <- min(rt_allmax_idx, rt_roi_center_idx + rt_cluster_half_length)
  }

  if (!is.null(drift_time)) {
    peak_list_clustered$ref_roi_dt_min_ms <- drift_time[peak_list_clustered$ref_roi_dt_min_idx]
    peak_list_clustered$ref_roi_dt_max_ms <- drift_time[peak_list_clustered$ref_roi_dt_max_idx]
  }
  if (!is.null(retention_time)) {
    peak_list_clustered$ref_roi_rt_min_s <- retention_time[peak_list_clustered$ref_roi_rt_min_idx]
    peak_list_clustered$ref_roi_rt_max_s <- retention_time[peak_list_clustered$ref_roi_rt_max_idx]
  }



    # thrOverlap <- 0.8
    # for (j in ROIs2Fusion){
    #   R1 <- clusters_infor[j, c(13:18)]
    #   if (abs(overlapPercentage(R1, medianROI)) >= thrOverlap){
    #     clusters_infor$dt_min_idx[j] <- median_dtmin_clust
    #     clusters_infor$dt_max_idx[j] <- median_dtmax_clust
    #     clusters_infor$rt_min_idx[j] <- median_rtmin_clust
    #     clusters_infor$rt_max_idx[j] <- median_rtmax_clust
    #   } else {
    #      ## QUE HAGO SI NO PASA EL PORCENTAJE
    #     }
    # }

  list(
    peak_list_clustered = peak_list_clustered,
    cluster_stats = cluster_stats
  )
}

