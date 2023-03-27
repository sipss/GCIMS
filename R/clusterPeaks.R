#' Group peaks in clusters
#'
#' @param peaks A data frame with at least the following columns:
#'  - "UniqueID" A unique ID for each peak
#'  - "SampleID" The sample ID the peak belongs to
#'  - "dt_apex_ms", "rt_apex_s" The peak positions
#'  - "dt_max_ms", "dt_min_ms", "rt_max_s", "rt_min_s" (for filtering outlier peaks based on their size)
#' @param distance_method A string. One of the distance methods from [stats::dist], "sd_scaled_euclidean" or "mahalanobis"
#' @param distance_between_peaks_from_same_sample The distance between two peaks from the same sample will be set to `distance_between_peaks_from_same_sample*max(distance_matrix)`
#' @param clustering A named list with "method" and the supported method, as well as further options.
#'   For `method = "kmedoids"`, you must provide `Nclusters`, with either the number of clusters
#'   to use in the kmedoids algorithm ([cluster::pam]) or the string `"max_peaks_sample"` to use the maximum number of
#'   detected peaks per sample.
#'
#'   For `method = "hclust"`, you can provide `hclust_method`, with the `method` passed to [mdendro::linkage()].
#' @param verbose logical, to control printing in the function
#' @param ... Ignored. All other parameters beyond `peaks` should be named
#' @param dt_cluster_spread_ms,rt_cluster_spread_s The typical spread of the clusters. Used for scaling.
#' dimensions when computing distances. When `clustering$method` is `"hclust"`, these spreads are used to cut cluster sizes.
#' @description Peak grouping function, exposing several options useful for benchmarking.
#'
#' @return A list with :
#'  - `peak_list_clustered`: The peak list with a "cluster" column
#'  - `cluster_stats`: Cluster statistics (cluster size...)
#'  - `dist`: peak to peak distance object
#'  - `extra_clustering_info`: Arbitrary clustering extra information, that depends on the clustering method
#' @examples
#' peak_list_fn <- system.file("extdata", "peak_list.rds", package = "GCIMS")
#' peak_list <- readRDS(peak_list_fn)
#'
#' peak_clustering  <- clusterPeaks(peak_list)
#' @export
clusterPeaks <- function(
    peaks,
    ...,
    distance_method = "euclidean",
    dt_cluster_spread_ms = 0.1,
    rt_cluster_spread_s = 20,
    distance_between_peaks_from_same_sample = 100,
    clustering = list(method = "hclust"),
    verbose = FALSE
) {
  # Compute the peak to peak distance:
  peak_matrix <- as.matrix(peaks[,c("dt_apex_ms", "rt_apex_s")])
  rownames(peak_matrix) <- peaks$UniqueID
  peak_matrix_scaled <- peak_matrix
  peak_matrix_scaled[,"dt_apex_ms"] <- peak_matrix_scaled[,"dt_apex_ms"]/dt_cluster_spread_ms
  peak_matrix_scaled[,"rt_apex_s"] <- peak_matrix_scaled[,"rt_apex_s"]/rt_cluster_spread_s
  colnames(peak_matrix_scaled) <- c("dt_apex_scaled", "rt_apex_scaled")

  peak2peak_dist <- peak2peak_distance(
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

  peak2peak_dist <- set_peak_distances_within_groups(
    dist_matrix = peak2peak_dist,
    peak_groups = peakuids_by_sample,
    value = distance_between_peaks_from_same_sample*max(peak2peak_dist)
  )

  extra_clustering_info <- list()
  if (clustering$method == "kmedoids") {
    require_pkgs("cluster")
    if (clustering$Nclusters == "max_peaks_sample") {
      N_clusters <- max(purrr::map_int(peakuids_by_sample, length))
    } else if (is.numeric(clustering$Nclusters)) {
      N_clusters <- clustering$Nclusters
    } else {
      cli_abort("When clustering$method is kmedoids, clustering$Nclusters must be an integer or the string 'max_peaks_sample'")
    }
    cluster <- cluster::pam(x = peak2peak_dist, k = N_clusters, pamonce = 3)
    peaks$cluster <- cluster$clustering
    extra_clustering_info$cluster_result <- cluster
    extra_clustering_info$silhouette <- cluster::silhouette(cluster, dist = peak2peak_dist)
  } else if (clustering$method == "hclust") {
    hclust_method <- ifelse(is.null(clustering$hclust_method), "complete", clustering$hclust_method)
    cluster <- mdendro::linkage(peak2peak_dist, method = hclust_method)

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
      dt_breaks = .data$rt_length_s > rt_cluster_spread_s,
      rt_breaks = .data$dt_length_ms > dt_cluster_spread_ms,
      breaks = .data$dt_breaks | .data$rt_breaks
    )

    peaks$cluster <- NA_integer_
    for (i in seq_len(nrow(merger_dist))) {
      if (!merger_dist$breaks[i]) {
        peaks$cluster[merger_dist$singletons[[i]]] <- i
      }
    }
    # renumber the values to the 1:N
    peaks$cluster <- as.integer(as.factor(peaks$cluster))
    extra_clustering_info$cluster <- cluster
    extra_clustering_info$num_clusters <- length(unique(peaks$cluster))
    extra_clustering_info$merger_distances <- merger_dist
  } else {
    cli::cli_abort("Unsupported clustering method {clustering$method}")
  }

  # Turn numeric peak clusters into IDs:
  if (is.numeric(peaks$cluster)) {
    ndigits_print <- paste0("Cluster%0", nchar(max(peaks$cluster, na.rm = TRUE)), "d")
    peaks$cluster <- ifelse(
      is.na(peaks$cluster),
      NA_character_,
      sprintf(ndigits_print, peaks$cluster)
    )
  }

  peak_cluster_stats <- peak_and_cluster_metrics(peaks)

  list(
    peak_list_clustered = peak_cluster_stats$peaks,
    cluster_stats = peak_cluster_stats$cluster_stats,
    dist = peak2peak_dist,
    extra_clustering_info = extra_clustering_info
  )
}



peak_and_cluster_metrics <- function(peaks) {
  peaks <- peaks |>
    dplyr::mutate(
      dt_length_ms = .data$dt_max_ms - .data$dt_min_ms,
      rt_length_s = .data$rt_max_s - .data$rt_min_s,
      dt_center_ms = (.data$dt_max_ms + .data$dt_min_ms)/2,
      rt_center_s = (.data$rt_max_s + .data$rt_min_s)/2,
    )

  cluster_stats <- peaks |>
    dplyr::filter(!is.na(.data$cluster)) |>
    dplyr::mutate(
      dt_apex_to_min_ms = .data$dt_apex_ms - .data$dt_min_ms,
      dt_apex_to_max_ms = .data$dt_max_ms - .data$dt_apex_ms,
      rt_apex_to_min_s = .data$rt_apex_s - .data$rt_min_s,
      rt_apex_to_max_s = .data$rt_max_s - .data$rt_apex_s
    ) |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(
          c(
            "dt_apex_ms", "rt_apex_s",
            "dt_apex_to_min_ms", "dt_apex_to_max_ms",
            "rt_apex_to_min_s", "rt_apex_to_max_s"
          )
        ),
        stats::median
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      dt_min_ms = .data$dt_apex_ms - .data$dt_apex_to_min_ms,
      dt_max_ms = .data$dt_apex_ms + .data$dt_apex_to_max_ms,
      dt_length_ms = .data$dt_max_ms - .data$dt_min_ms,
      rt_min_s = .data$rt_apex_s - .data$rt_apex_to_min_s,
      rt_max_s = .data$rt_apex_s + .data$rt_apex_to_max_s,
      rt_length_s = .data$rt_max_s - .data$rt_min_s
    )

  peaks_fixed_size <- dplyr::left_join(
      dplyr::select(
        peaks,
        dplyr::all_of(c("UniqueID", "cluster", "dt_apex_ms", "rt_apex_s"))
      ),
      dplyr::select(
        cluster_stats,
        dplyr::all_of(c("cluster", "dt_apex_to_min_ms", "dt_apex_to_max_ms", "rt_apex_to_min_s", "rt_apex_to_max_s"))
      ),
      by = "cluster"
    ) |>
    dplyr::mutate(
      fixedsize_dt_min_ms = .data$dt_apex_ms - .data$dt_apex_to_min_ms,
      fixedsize_dt_max_ms = .data$dt_apex_ms + .data$dt_apex_to_max_ms,
      fixedsize_rt_min_s = .data$rt_apex_s - .data$rt_apex_to_min_s,
      fixedsize_rt_max_s = .data$rt_apex_s + .data$rt_apex_to_max_s,
    ) |>
    dplyr::select(dplyr::all_of(c("UniqueID", "fixedsize_dt_min_ms", "fixedsize_dt_max_ms", "fixedsize_rt_min_s", "fixedsize_rt_max_s")))
  peaks <- dplyr::left_join(
    peaks,
    peaks_fixed_size,
    by = "UniqueID"
  )
  list(peaks = peaks, cluster_stats = cluster_stats)
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
  } else if (distance_method == "mahalanobis") {
    peak2peak_dist <- mahalanobis_distance(peak_matrix)
  } else {
    cli_abort("Unsupported distance {distance_method}")
  }
  peak2peak_dist
}
