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

#'   For `method = "hclust_old"`, you can provide `hclust_method`, with the `method` passed to [stats::hclust()].
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
#' \donttest{
#' dir_in <- system.file("extdata", package = "GCIMS")
#' peak_list <- readRDS(file.path(dir_in, "peak_list.rds"))
#'
#' peak_clustering  <- clusterPeaks(
#'   peaks = peak_list,
#'   distance_method = "mahalanobis",
#'   distance_between_peaks_from_same_sample = Inf,
#'   clustering = list(method = "hclust"),
#'   verbose = FALSE
#' )
#'}
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
  } else if (clustering$method == "hclust_old") {
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
    stop(sprintf("Unsupported clustering method %s", clustering$method))
  }

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
}



peak_and_cluster_metrics <- function(peaks) {
  peaks <- peaks |>
    dplyr::mutate(
      dt_length_ms = .data$dt_max_ms - .data$rt_min_s,
      rt_length_s = .data$rt_max_s - .data$rt_min_s,
      dt_center_ms = (.data$dt_max_ms + .data$rt_min_s)/2,
      rt_center_s = (.data$rt_max_s + .data$rt_min_s)/2,
    )

  cluster_stats <- peaks |>
    dplyr::filter(!is.na(cluster)) |>
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
      fixedsize_dt_min_ms = .data$dt_apex_ms + .data$dt_apex_to_min_ms,
      fixedsize_dt_max_ms = .data$dt_apex_ms + .data$dt_apex_to_max_ms,
      fixedsize_rt_min_s = .data$rt_apex_s + .data$rt_apex_to_min_s,
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

#' @param num_clusters: A numeric vector with candidates for the number of clusters to choose
#' @param peak_list: A data frame with peaks, including "UniqueID", "dt_apex_ms", "rt_apex_s"
#' @param cluster: The outcome of the hierarchical clustering
#' @param max_dist_thresh_ppb, the maximum distance allowed within a cluster
#' @return A data frame with two columns: The given num_clusters and the maximum measured cluster size within
#'         clusters
#' @noRd
get_max_dist_ppb_for_num_clusters <- function(num_clusters, peak_list, cluster, dt_ms_max_dist_thres, rt_s_max_dist_thres) {
  peak_assignments <- stats::cutree(cluster, k = num_clusters)
  peak_assignments <- peak_assignments[peak_list$UniqueID, ]
  peak_list$cluster <- NULL
  dt_ms_max_dist <- numeric(length(num_clusters))
  rt_s_max_dist <- numeric(length(num_clusters))
  dt_break_in <- NULL
  rt_break_in <- NULL
  limiting_threshold <- "none"
  for (i in seq_len(ncol(peak_assignments))) {
    peak_list$cluster <- peak_assignments[,i]
    max_distances <- peak_list |>
      dplyr::group_by(.data$cluster)  |>
      dplyr::summarize(
        dt_max_dist_ms = max(.data$dt_apex_ms) - min(.data$dt_apex_ms),
        rt_max_dist_s = max(.data$rt_apex_s) - min(.data$rt_apex_s),
        .groups = "drop")

    dt_ms_max_dist[i] <- max(max_distances$dt_max_dist_ms)
    rt_s_max_dist[i] <- max(max_distances$rt_max_dist_s)

    if (!is.null(dt_ms_max_dist_thres) && is.null(dt_break_in) && dt_ms_max_dist[i] < dt_ms_max_dist_thres) {
      dt_break_in <- i + 10
    }
    if (!is.null(rt_s_max_dist_thres) && is.null(rt_break_in) && rt_s_max_dist[i] < rt_s_max_dist_thres) {
      rt_break_in <- i + 10
    }
    # Do least stringent to have more exploration
    if (is.null(dt_break_in) || is.null(rt_break_in)) {
      next
    }
    if (dt_break_in > rt_break_in) {
      limiting_threshold <- "rt"
    } else if (dt_break_in < rt_break_in) {
      limiting_threshold <- "dt"
    } else {
      limiting_threshold <- "both"
    }
    break
  }
  list(
    clust_dist = data.frame(
      num_clusters = num_clusters[seq_len(i)],
      dt_ms_max_dist = dt_ms_max_dist[seq_len(i)],
      rt_s_max_dist = rt_s_max_dist[seq_len(i)]
    ),
    limiting_threshold = limiting_threshold
  )
}

estimate_num_clusters <- function(peak_list, cluster, dt_ms_max_dist_thres, rt_s_max_dist_thres) {
  peaks_per_sample <- peak_list |>
    dplyr::group_by(.data$SampleID) |>
    dplyr::summarize(n = dplyr::n()) |>
    dplyr::pull("n")
  min_clusters_to_test <- max(peaks_per_sample)
  max_clusters_to_test <- sum(peaks_per_sample)
  if ((max_clusters_to_test - min_clusters_to_test) > 20) {
    num_clusters_coarse <- seq.int(from = max(peaks_per_sample), to = sum(peaks_per_sample), by = 10)
  } else {
    num_clusters_coarse <- seq.int(from = max(peaks_per_sample), to = sum(peaks_per_sample), by = 1)
  }
  clust_dist_and_limiting <- get_max_dist_ppb_for_num_clusters(
    num_clusters_coarse,
    peak_list,
    cluster,
    dt_ms_max_dist_thres,
    rt_s_max_dist_thres
  )
  clust_dist <- clust_dist_and_limiting$clust_dist
  limiting_threshold <- clust_dist_and_limiting$limiting_threshold
  # FIXME: FIX NEXT LINE
  num_clusters <- clust_dist$num_clusters[
    clust_dist$dt_ms_max_dist < dt_ms_max_dist_thres & clust_dist$rt_s_max_dist < rt_s_max_dist_thres
  ][1]
  # Refine:
  num_clusters_fine <- seq.int(
    from = max(min_clusters_to_test, num_clusters - 19),
    to = min(max_clusters_to_test, num_clusters + 11)
  )
  clust_dist2 <- get_max_dist_ppb_for_num_clusters(num_clusters_fine, peak_list, cluster, dt_ms_max_dist_thres = NULL, rt_s_max_dist_thres = NULL)
  # Combine:
  num_clusters_vs_max_distance <- dplyr::bind_rows(clust_dist, clust_dist2$clust_dist) |>
    dplyr::arrange(num_clusters) |>
    dplyr::distinct()
  dt_num_clusters <- num_clusters_vs_max_distance |>
    dplyr::filter(.data$dt_ms_max_dist < !!dt_ms_max_dist_thres) |>
    dplyr::pull("num_clusters")
  rt_num_clusters <- num_clusters_vs_max_distance |>
    dplyr::filter(.data$rt_s_max_dist < !!rt_s_max_dist_thres) |>
    dplyr::pull("num_clusters")
  num_clusters <- num_clusters_vs_max_distance |>
    dplyr::filter(.data$dt_ms_max_dist < !!dt_ms_max_dist_thres, .data$rt_s_max_dist < !!rt_s_max_dist_thres) |>
    dplyr::pull("num_clusters")
  dt_gplt <- ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$num_clusters,
        y = .data$dt_ms_max_dist
      ),
      data = num_clusters_vs_max_distance,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = dt_ms_max_dist_thres, color = "gray") +
    ggplot2::labs(x = "Number of clusters", y = "Max drift time distance within cluster (ms)")
  rt_gplt <- ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$num_clusters,
        y = .data$rt_s_max_dist
      ),
      data = num_clusters_vs_max_distance,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = rt_s_max_dist_thres, color = "gray") +
    ggplot2::labs(x = "Number of clusters", y = "Max retention time distance within cluster (s)")
  if (length(num_clusters) == 0) {
    rlang::abort(
      c(
        "Can't find a suitable number of clusters",
        "Probably the distance threshold is too small",
        "i" = "Please consider increasing the threshold of the maximum distance",
        "i" = glue("Current thresholds are dt_ms_max_dist_thres={dt_ms_max_dist_thres} ms and rt_s_max_dist_thres={rt_s_max_dist_thres} ss.)"),
        "i" = "Use `rlang::last_error()$dt_plot` and `rlang::last_error()$dt_plot` to see plots showing the maximum distance vs the number of clusters explored and guide you"
      ),
      dt_plot = dt_gplt,
      rt_plot = rt_gplt,
    )
  }
  num_clusters <- num_clusters[1]
  dt_gplt <- dt_gplt +
    ggplot2::geom_vline(xintercept = num_clusters, color = "red")
  rt_gplt <- rt_gplt +
    ggplot2::geom_vline(xintercept = num_clusters, color = "red")
  if (length(dt_num_clusters) > 0) {
    dt_gplt <- dt_gplt +
      ggplot2::geom_vline(xintercept = !!dt_num_clusters[1], color = "red", linetype = "dashed")
  }
  if (length(rt_num_clusters) > 0) {
    rt_gplt <- rt_gplt +
      ggplot2::geom_vline(xintercept = !!rt_num_clusters[1], color = "red", linetype = "dashed")
  }
  list(
    num_clusters = num_clusters,
    dt_num_clusters = utils::head(dt_num_clusters, 1),
    rt_num_clusters = utils::head(rt_num_clusters, 1),
    limiting_threshold = limiting_threshold,
    table = num_clusters_vs_max_distance,
    dt_ms_max_dist_thres = dt_ms_max_dist_thres,
    rt_s_max_dist_thres = rt_s_max_dist_thres,
    dt_plot = dt_gplt,
    rt_plot = rt_gplt
  )
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
    rlang::warn(
      message = c(
        "Deprecated distance metric",
        "i" = "The 'sd_scaled_euclidean' metric is now deprecated",
        "i" = "Please use the following arguments instead",
        "i" = '  distance_method="euclidean"',
        "i" = '  dt_cluster_spread_ms=<typical cluster spread in drift time in ms>',
        "i" = '  rt_cluster_spread_s=<typical cluster spread in retention time in s>',
        "i" = "Do not use the standard deviation of all peaks as the cluster spread, as it doesn't make sense"
      )
    )
    peak_matrix_scaled <- scale(peak_matrix, center = FALSE, scale = TRUE)
    peak2peak_dist <- stats::dist(peak_matrix_scaled, method = "euclidean")
  } else if (distance_method == "mahalanobis") {
    peak2peak_dist <- mahalanobis_distance(peak_matrix)
  } else {
    stop(sprintf("Unsupported distance %s", distance_method))
  }
  peak2peak_dist
}
