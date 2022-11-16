
peaks <- readRDS("peak_list_ketones.rds")
peaks <- readRDS("peak_list_urine.rds")
distance_method <- "euclidean"
dt_cluster_spread_ms <- 0.7
rt_cluster_spread_s <- 7
dt_cluster_spread_ms <- 0.1
rt_cluster_spread_s <- 20
distance_between_peaks_from_same_sample <- 100
verbose <- FALSE



# Compute the peak to peak distance:
peak_matrix <- as.matrix(peaks[,c("dt_apex_ms", "rt_apex_s")])
rownames(peak_matrix) <- peaks$UniqueID
peak_matrix_scaled <- peak_matrix
peak_matrix_scaled[,"dt_apex_ms"] <- peak_matrix_scaled[,"dt_apex_ms"]/dt_cluster_spread_ms
peak_matrix_scaled[,"rt_apex_s"] <- peak_matrix_scaled[,"rt_apex_s"]/rt_cluster_spread_s
colnames(peak_matrix_scaled) <- c("dt_apex_scaled", "rt_apex_scaled")

peak2peak_dist <- GCIMS:::peak2peak_distance(
  peak_matrix = peak_matrix_scaled,
  distance_method = distance_method
)

# Set distances from pairs of peaks belonging to the same sample to Inf,
# so they are never in the same cluster
peakuids_by_sample <- peaks |>
  dplyr::select(dplyr::all_of(c("SampleID", "UniqueID"))) |>
  dplyr::group_by(.data$SampleID) |>
  dplyr::summarize(UniqueIDs = list(.data$UniqueID)) |>
  dplyr::ungroup() |>
  dplyr::pull(.data$UniqueIDs)

peak2peak_dist <- GCIMS:::set_peak_distances_within_groups(
  dist_matrix = peak2peak_dist,
  peak_groups = peakuids_by_sample,
  value = distance_between_peaks_from_same_sample*max(peak2peak_dist)
)

extra_clustering_info <- list()

# Clustering:
cluster <- mdendro::linkage(peak2peak_dist, method = "complete")
cluster_hclust <- stats::hclust(d = peak2peak_dist, method = "complete")

cluster$merger[[1]]

# ind2pair <- function(ind, n) {
#   all_i <- integer(length(ind))
#   all_j <- integer(length(ind))
#   k <- seq_len(n - 1)
#   first_j <- 1 + n*(k - 1) - (k*(k - 1))/2
#   for (xi in seq_along(ind)) {
#     x <- ind[xi]
#     j <- which(x >= first_j)
#     j <- j[length(j)]
#     i <- (j + 1) + (x - first_j[j])
#     all_i[xi] <- i
#     all_j[xi] <- j
#   }
#   data.frame(ind = ind, i = all_i, j = all_j)
# }

peak2peak_mat <- as.matrix(peak2peak_dist)

merger_dist <- tibble::tibble(
  merger_id = seq_along(cluster$merger),
  dt_length_ms = Inf,
  rt_length_s = Inf,
  singletons = vector("list", length = length(cluster$merger))
)

get_merged_idx <- function(idx, merger) {
  singletons <- idx[idx < 0]
  clusters <- idx[idx > 0]
  singletons2 <- c()
  for (cl in clusters) {
    thiscl <- get_merged_idx(merger[[cl]], merger)
    singletons2 <- c(singletons2, thiscl)
  }
  abs(c(singletons, singletons2))
}

for (i in seq_len(nrow(merger_dist))) {
  peaks_merged <- get_merged_idx(cluster$merger[[i]], merger = cluster$merger)
  merger_dist$dt_length_ms[i] <- diff(range(peaks$dt_apex_ms[peaks_merged]))
  merger_dist$rt_length_s[i] <- diff(range(peaks$rt_apex_s[peaks_merged]))
  merger_dist$singletons[[i]] <- peaks_merged
}


merger_dist <- dplyr::mutate(
    merger_dist,
    dt_breaks = rt_length_s > rt_cluster_spread_s,
    rt_breaks = dt_length_ms > dt_cluster_spread_ms,
    breaks = dt_breaks | rt_breaks
  )

peaks_clustered <- data.frame(
  UniqueID = peaks$UniqueID,
  cluster = NA_integer_
)
for (i in seq_len(nrow(merger_dist))) {
  if (!merger_dist$breaks[i]) {
    peaks_clustered$cluster[merger_dist$singletons[[i]]] <- i
  }
}

peaks_clustered

peaks
saveRDS(object = peaks_clustered, "peaks_clustered_urine.rds")






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
# Clustering ends




# Turn numeric peak clusters into IDs:
if (is.numeric(peaks$cluster)) {
  ndigits_print <- paste0("Cluster%0", nchar(max(peaks$cluster)), "d")
  peaks$cluster <- sprintf(ndigits_print, peaks$cluster)
}

peak_cluster_stats <- peak_and_cluster_metrics(peaks)

list(
  peak_list_clustered = peak_cluster_stats$peaks,
  cluster_stats = peak_cluster_stats$cluster_stats,
  dist = peak2peak_dist,
  extra_clustering_info = extra_clustering_info
)
