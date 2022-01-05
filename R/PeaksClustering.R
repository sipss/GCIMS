#' ROIs Clustering

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where peak table data file is
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples that are going to be included in the peak table.
#' @details `gcims_peaks_clustering`  does a clustering along the ROIs
#'   for a peak table creation. The Figures of merits of each ROI are also
#'   reported. In this table are included all samples in `samples`. Use this
#'   function if you are interested in obtaining a final peak table for future
#'   classification techniques.
#' @return A Set of S3 objects.
#' @family Peaks CLustering function
#' @export
#' @importFrom cluster pam
#' @importFrom plyr count
#' @importFrom Hotelling hotelling.test
#' @importFrom rlang .data
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 1:3
#'
#' # Example of ROIs Clustering for Peak Table Creation
#' # Need a proper dataset with peaks detected in the pkg
#' # (or maybe better a refactor the function arguments)
#' peak_table <- gcims_peaks_clustering(dir_in, dir_out, samples)
gcims_peaks_clustering <- function(dir_in, dir_out, samples){
  print(" ")
  print("  ///////////////////////////")
  print(" /    Clustering the ROIs  /")
  print("///////////////////////////")
  print(" ")


  samples_fn <- utils::setNames(
    object = file.path(dir_in, paste0("M", samples, ".rds")),
    nm = paste0("M", samples)
  )
  peaks <- rds_samples_to_peak_list(samples_fn)
  group_peak_list(peaks)
  group_peak_list(
    peaks = peaks,
    filter_dt_width_criteria = "IQR",
    filter_rt_width_criteria = "arnau",
    distance_method = "mahalanobis",
    distance_between_peaks_from_same_sample = Inf,
    clustering = list(method = "kmedoids", Nclusters = "max_peaks_sample"),
    aggregate_conflicting_peaks = mean,
    verbose = FALSE
  )
}

#' Peak grouping function, exposing a lot of options useful for benchmarking
#'
#' @param peaks A data frame with at least the following columns:
#'  - "UniqueID" A unique ID for each peak
#'  - "SampleID" The sample ID the peak belongs to
#'  - "dtapex_ms", "rtapex_s" The peak positions
#'  - "volume" The size of the peak
#'  - "dt_max_ms", "dt_min_ms", "rt_max_s", "rt_min_s" (for filtering outlier peaks based on their size)
#' @param filter_dt_width_criteria,filter_rt_width_criteria A character with the method for outlier detection.
#'   - "IQR": Remove peaks more than 1.5 interquartile ranges above upper quartile or
#'     below the lower quartile.
#'   - "arnau": FIXME Adhoc method from Arnau, where he removes peaks above mean+4iqr or below median-0.75iqr
#'   - "none": Do not remove peaks based on their drift time width or retention time height
#' @param distance_method A string. One of the distance methods from [stats::dist] or "mahalanobis"
#' @param distance_between_peaks_from_same_sample The distance between two peaks from the same sample will be set to `distance_between_peaks_from_same_sample*max(distance_matrix)`
#' @param clustering A named list with "method" and the supported method, as well as further options.
#'   For `method = "kmedoids"`, you must provide `Nclusters`, with either the number of clusters
#'   to use in the kmedoids algorithm or the string "max_peaks_sample" to use the maximum number of
#'   detected peaks per sample.
#'
#'   No other clustering method is currently supported
#' @param aggregate_conflicting_peaks `NULL` or a function. When we build the peak table, with peaks in rows, samples in
#'  columns, peak_table[i,j] is the volume of the peak from sample j in cluster i. If the clustering process
#'  clusters together two peaks form the same sample, those peaks will conflict in the peak table. `NULL` will error
#'  in that case, another function will be applied on the conflicting volumes (e.g `mean` or `max` would be reasonable options)
#'
#' @param verbose logical, to control printing in the function
#'
#' @return A list with :
#' - peak_table: A peak table that includes peak position, median peak minimum/maximum retention and drift times and the peak volume for each sample
#' - peak_table_duplicity: How many volume values have been aggregated. Should be 1 for each sample/peak
#' - extra_clustering_info: Arbitrary clustering extra information, that depends on the clustering method
#' @export
#'
#' @examples
#' # FIXME
group_peak_list <- function(
  peaks,
  filter_dt_width_criteria = "IQR",
  filter_rt_width_criteria = "arnau",
  distance_method = "mahalanobis",
  distance_between_peaks_from_same_sample = Inf,
  clustering = list(method = "kmedoids", Nclusters = "max_peaks_sample"),
  aggregate_conflicting_peaks = NULL,
  verbose = FALSE
) {
  # 0. Warn if peaks with NA positions, and remove them
  peaks_with_na <- complete.cases(peaks)
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

  peak_matrix <- as.matrix(peaks[,c("dtapex_ms", "rtapex_s")])
  rownames(peak_matrix) <- peaks$UniqueID

  STATS_METHODS <- c("euclidean", "maximum", "manhattan", "canberra",
                     "binary", "minkowski")
  if (distance_method %in% STATS_METHODS) {
    peak2peak_dist <- stats::dist(peak_matrix, method = distance_method)
  } else if (distance_method == "mahalanobis") {
    peak2peak_dist <- mahalanobis_distance(peak_matrix)
  } else {
    stop(sprintf("Unsupported distance %s", distance_method))
  }

  # Set distances from pairs of peaks belonging to the same sample to Inf,
  # so they are never in the same cluster
  peakuids_by_sample <- peaks %>%
    dplyr::select(dplyr::all_of(c("SampleID", "UniqueID"))) %>%
    dplyr::group_by(.data$SampleID) %>%
    dplyr::summarize(UniqueIDs = list(.data$UniqueID)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.data$UniqueIDs)

  peak2peak_distance <- set_peak_distances_within_groups(
    dist_matrix = peak2peak_dist,
    peak_groups = peakuids_by_sample,
    value = distance_between_peaks_from_same_sample*max(peak2peak_dist)
  )

  if (clustering$method == "kmedoids") {
    if (clustering$Nclusters == "max_peaks_sample") {
      N_clusters <- max(purrr::map_int(peakuids_by_sample, length))
    } else if (is.numeric(clustering$Nclusters)) {
      N_clusters <- clustering$Nclusters
    } else {
      stop("When clustering$method is kmedoids, clustering$Nclusters must be an integer or the string 'max_peaks_sample'")
    }
    cluster <- cluster::pam(x = peak2peak_distance, k = N_clusters)
    peaks$cluster <- cluster$clustering
    extra_clustering_info <- cluster
  } else {
    stop(sprintf("Unsupported clustering method %s", clustering$method))
  }

  median_roi_per_cluster <- peaks %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarise(
      dplyr::across(
        c(dplyr::starts_with("dt"), dplyr::starts_with("rt")),
        stats::median
      )
    ) %>%
    dplyr::ungroup()

  peak_table_duplicity <- peaks %>%
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "volume"))) %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of("SampleID"),
      values_from = dplyr::all_of("volume"),
      values_fn = length
    )

  peak_table <- peaks %>%
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "volume"))) %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of("SampleID"),
      values_from = dplyr::all_of("volume"),
      values_fn = aggregate_conflicting_peaks
    )

  peak_table_annotated <- dplyr::left_join(
    median_roi_per_cluster,
    peak_table,
    by = "cluster"
  )
  # Missing values still need to be filled
  list(
    peak_table = peak_table_annotated,
    peak_table_duplicity = peak_table_duplicity,
    peak_list_with_cluster = peaks,
    extra_clustering_info = extra_clustering_info
  )
}

rds_samples_to_peak_list <- function(samples) {
  pb <- progress::progress_bar$new(total = length(samples))

  peaks_df <- purrr::map_dfr(samples, function(sample_fn) {
    pb$tick()
    aux_list <- readRDS(sample_fn)

    parameters_df_i <- as.data.frame(t(aux_list$data$Parameters))
    rownames(parameters_df_i) <- sprintf("Peak%d", seq_len(nrow(parameters_df_i)))

    rois_i <- aux_list$data$ROIs
    rownames(rois_i) <- sprintf("Peak%d", seq_len(nrow(rois_i)))
    colnames(rois_i) <- c("dtmin_idx", "dtmax_idx", "rtmin_idx", "rtmax_idx")

    apex_i <- aux_list$data$Peaks
    colnames(apex_i) <- c("dtapex_idx", "rtapex_idx")
    rownames(rois_i) <- sprintf("Peak%d", seq_len(nrow(apex_i)))

    # Maybe instead of the Peaks, ROIs and Parameters, we could just
    # have one dataframe in the RDS:
    peaks_i <- tibble::rownames_to_column(
      cbind(apex_i, rois_i, parameters_df_i),
      var = "PeakID"
    )
    dt <- aux_list$data$drift_time
    rt <- aux_list$data$retention_time
    peaks_i <- peaks_i %>%
      dplyr::mutate(
        dplyr::across(dplyr::starts_with("dt"), ~ dt[.], .names = "{.col}DTDELETEME"),
        dplyr::across(dplyr::starts_with("rt"), ~ rt[.], .names = "{.col}RTDELETEME")
      ) %>%
      dplyr::rename_with(
        .fn = ~gsub("_idxDTDELETEME", "_ms", .),
        .cols = dplyr::ends_with("_idxDTDELETEME")
      ) %>%
      dplyr::rename_with(
        .fn = ~gsub("_idxRTDELETEME", "_s", .),
        .cols = dplyr::ends_with("_idxRTDELETEME")
      ) %>%
      dplyr::select(
        dplyr::all_of(
          c("PeakID", "dtapex_ms", "rtapex_s",
            "rtmin_s", "rtmax_s", "dtmin_ms", "dtmax_ms")
        ),
        dplyr::ends_with("_idx"),
        dplyr::everything()
      )
    peaks_i
  }, .id = "SampleID")

  # Add a unique peak ID, that combines Sample+Peak ids:
  peaks_df <- tibble::add_column(
    peaks_df,
    UniqueID = paste0(peaks_df$SampleID, peaks_df$PeakID),
    .before = 1
  )
  tibble::as_tibble(peaks_df)
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
    dt_width_ms = peaks$dtmax_ms - peaks$dtmin_ms,
    rt_width_s = peaks$rtmax_s - peaks$rtmin_s
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
