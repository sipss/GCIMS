peaks <- peak_list
distance_method <- "sd_scaled_euclidean"
distance_between_peaks_from_same_sample <- 3

peak_matrix <- as.matrix(peaks[,c("dt_apex_ms", "rt_apex_s")])
rownames(peak_matrix) <- peaks$UniqueID

peak2peak_dist <- GCIMS:::peak2peak_distance(
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

peak2peak_dist <- GCIMS:::set_peak_distances_within_groups(
  dist_matrix = peak2peak_dist,
  peak_groups = peakuids_by_sample,
  value = distance_between_peaks_from_same_sample*max(peak2peak_dist)
)

N_clusters <- 450

cluster <- cluster::pam(x = peak2peak_dist, k = N_clusters, pamonce = 3)
peaks_out <- peaks
peaks_out$cluster <- cluster$clustering
cluster
silhouette <- cluster::silhouette(cluster, dist = peak2peak_dist)
cosa2 <- tibble::rownames_to_column(as.data.frame(as.matrix(silhouette)), "SampleID")

ggplot(cosa2) +
  geom_boxplot(aes(x = gsub("/.*", "", SampleID), y = sil_width)) +
  labs(x = "SampleID", y = "Silouhette width", title = sprintf("N clusters: %d", N_clusters)) +
  coord_flip()


peak_clustering_hclust2 <- clusterPeaks(
  peak_list,
  distance_method = "sd_scaled_euclidean",
  clustering = list(
    method = "hclust2",
    method = "complete",
    dt_ms_max_dist_thres = 0.58*3,
    rt_s_max_dist_thres = 28*3
  ),
  verbose = TRUE
)
