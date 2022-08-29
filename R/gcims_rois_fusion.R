#' Merge regions of interest
#'
#' @description For each ROI in `peak_list_clustered`, we define additional columns with ROI boundaries.
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
#' merge_rois_out <- gcims_merge_rois(clustering$peak_list_clustered, clustering$cluster_stats)
#' samplestocheck <- merge_rois_out$peak_list_clustered[
#'   which(merge_rois_out$peak_list_clustered$cluster == "Cluster01"),
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
gcims_merge_rois <- function(peak_list_clustered, cluster_stats) {
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

