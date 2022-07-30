#' ROIs Fusion
#' @param clustering_object   Objects obtained from the clustering function
#' containing a list of all the ROIs with tehir corresponding cluster.
#' @details `clustering_object` calculates de median coordinates of all the ROIs
#' of each cluster. All the ROis inside a cluster are resized according to their
#' geometrical center and the size of the median ROI.
#' @return A Set of S3 objects.
#' @family Utility functions
#' @export
#' @examples
#' # Example of ROIs FUsion
#' tt <- gcims_rois_fusion(peak_clustering)
#' samplestocheck <- tt$peak_list_clustered[which(tt$peak_list_clustered$cluster == 35), ]
#' ggplot() +
#'  geom_rect(data = samplestocheck, mapping = aes(xmin=dt_min_idx, xmax=dt_max_idx,
#'  ymin=rt_min_idx, ymax=rt_max_idx), color= as.factor(samplestocheck$SampleID), alpha=0.5)
#' invisible(file.remove(files))
gcims_rois_fusion <- function(clustering_object){
  clusters_infor <- clustering_object$peak_list_clustered
  clusters_infor$ref_roi_dt_min_idx <- rep(0, dim(clusters_infor)[1])
  clusters_infor$ref_roi_dt_max_idx <- rep(0, dim(clusters_infor)[1])
  clusters_infor$ref_roi_rt_min_idx <- rep(0, dim(clusters_infor)[1])
  clusters_infor$ref_roi_rt_max_idx <- rep(0, dim(clusters_infor)[1])
  clusters_infor$ref_roi_dt_apex <- rep(0, dim(clusters_infor)[1])
  clusters_infor$ref_roi_rt_apex <- rep(0, dim(clusters_infor)[1])
  num_clusts <- unique(clusters_infor$cluster)
  num_samps <- unique(clusters_infor$SampleID)
  for(i in num_clusts){
    # 4 corrdinates
    ROIs2Fusion <- which(clusters_infor$cluster == i)
    median_dtmin_clust <- floor(stats::median(clusters_infor$dt_min_idx[ROIs2Fusion]))
    median_dtmax_clust <- ceiling(stats::median(clusters_infor$dt_max_idx[ROIs2Fusion]))
    median_rtmin_clust <- floor(stats::median(clusters_infor$rt_min_idx[ROIs2Fusion]))
    median_rtmax_clust <- ceiling(stats::median(clusters_infor$rt_max_idx[ROIs2Fusion]))

    medianROI <- c(median_dtmin_clust, median_dtmax_clust, median_rtmin_clust, median_rtmax_clust)

    # Geometric Center
    median_dtapex_clust <- floor(medianROI[1] + (medianROI[2] - medianROI[1])/2)
    median_rtapex_clust <- floor(medianROI[3] + (medianROI[4] - medianROI[3])/2)

    clusters_infor$ref_roi_dt_min_idx[ROIs2Fusion] <- median_dtmin_clust
    clusters_infor$ref_roi_dt_max_idx[ROIs2Fusion] <- median_dtmax_clust
    clusters_infor$ref_roi_rt_min_idx[ROIs2Fusion] <- median_rtmin_clust
    clusters_infor$ref_roi_rt_max_idx[ROIs2Fusion] <- median_rtmax_clust
    clusters_infor$ref_roi_dt_apex[ROIs2Fusion] <- median_dtapex_clust
    clusters_infor$ref_roi_rt_apex[ROIs2Fusion] <- median_rtapex_clust

    for (j in ROIs2Fusion){
      ROI <- c(clusters_infor$dt_min_idx[j], clusters_infor$dt_max_idx[j],
               clusters_infor$rt_min_idx[j], clusters_infor$rt_max_idx[j])
      ROI_dtapex_clust <- floor(ROI[1] + (ROI[2] - ROI[1])/2)
      ROI_rtapex_clust <- floor(ROI[3] + (ROI[4] - ROI[3])/2)
      clusters_infor$dt_min_idx[j] <- ROI_dtapex_clust - (median_dtapex_clust - median_dtmin_clust)
      clusters_infor$dt_max_idx[j] <- ROI_dtapex_clust + (median_dtmax_clust - median_dtapex_clust)
      clusters_infor$rt_min_idx[j] <- ROI_rtapex_clust - (median_rtapex_clust - median_rtmin_clust)
      clusters_infor$rt_max_idx[j] <- ROI_rtapex_clust + (median_rtmax_clust - median_rtapex_clust)
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

  }
  clustering_object$peak_list_clustered <- clusters_infor
  return(clustering_object)
}

