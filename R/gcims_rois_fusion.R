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
#' peak_list <- data.frame(
#'   UniqueID = c("P1", "P2", "P3", "P4"),
#'   SampleID = c("S1", "S1", "S2", "S2"),
#'   dt_apex_idx = c(140, 340, 180, 380),
#'   dt_apex_ms = c(7, 10, 7.1, 10.2),
#'   rt_apex_idx = c(30, 250, 33, 247),
#'   rt_apex_s = c(30, 250, 33, 247),
#'   dt_min_idx = c(100, 300, 120, 320),
#'   dt_min_ms = c(6.5, 9.4, 6.6, 9.7),
#'   dt_max_idx = c(200, 400, 250, 440),
#'   dt_max_ms = c(7.7, 10.8, 7.6, 11.1),
#'   rt_min_idx = c(27, 246, 30, 245),
#'   rt_min_s = c(27, 246, 30, 245),
#'   rt_max_idx = c(36, 260, 37, 255),
#'   rt_max_s = c(36, 260, 37, 255),
#'   dt_cm_idx = c(140, 340, 180, 380),
#'   dt_cm_ms = c(7, 10, 7.1, 10.2),
#'   rt_cm_idx = c(30, 250, 33, 247),
#'   rt_cm_s = c(30, 250, 33, 247)
#' )
#' clustering <- group_peak_list(
#'   peaks = peak_list,
#'   filter_dt_width_criteria = NULL,
#'   filter_rt_width_criteria = NULL,
#'   distance_method = "mahalanobis",
#'   distance_between_peaks_from_same_sample = Inf,
#'   clustering = list(method = "kmedoids", Nclusters = "max_peaks_sample"),
#'   verbose = FALSE
#' )
#' roi_fusion_out <- gcims_rois_fusion(clustering$peak_list_clustered, clustering$cluster_stats)
#' samplestocheck <- roi_fusion_out$peak_list_clustered[which(roi_fusion_out$peak_list_clustered$cluster == "Cluster02"), ]
#' library(ggplot2)
#' ggplot() +
#'  geom_rect(
#'    data = samplestocheck,
#'    mapping = aes(
#'      xmin=dt_min_idx,
#'      xmax=dt_max_idx,
#'      ymin=rt_min_idx,
#'      ymax=rt_max_idx,
#'      color= SampleID
#'    ),
#'    alpha=0.5
#' )
#' invisible(file.remove(files))
gcims_rois_fusion <- function(peak_list_clustered, cluster_stats) {
  # FIXME: Provide these columns in physical units as well.
  # We will do this afterwards, when we have a dataset object that has
  # this information (to avoid reading the samples just to get their dtime and rtime)
  peak_list_clustered$ref_roi_dt_min_idx <- 0
  peak_list_clustered$ref_roi_dt_max_idx <- 0
  peak_list_clustered$ref_roi_rt_min_idx <- 0
  peak_list_clustered$ref_roi_rt_max_idx <- 0

  # Largest valid index. It should come from the sample/dataset, but I'm not going
  # to load all the samples just for that. We'll fix this later FIXME
  dt_allmax_idx <- max(c(peak_list_clustered$dt_max_idx, cluster_stats$dt_max_idx))
  rt_allmax_idx <- max(c(peak_list_clustered$rt_max_idx, cluster_stats$rt_max_idx))

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

