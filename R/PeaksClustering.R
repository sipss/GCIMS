#' ROIs Clustering

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where peak table data file is
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples that are going to be included in the peak table.
#' @details \code{gcims_peaks_clustering}  does a clustering along the ROIs
#'   for a peak table creation. The Figures of merits of each ROI are also
#'   reported. In this table are included all samples in \code{samples}. Use this
#'   function if you are interested in obtaining a final peak table for future
#'   classification techniques.
#' @return A Set of S3 objects.
#' @family Peaks CLustering function
#' @export
#' @importFrom cluster pam
#' @importFrom plyr count
#' @importFrom Hotelling hotelling.test
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 1:3
#'
#' # Example of ROIs Clustering for Peak Table Creation
#' gcims_peaks_clustering <- function(dir_in, dir_out, samples)
gcims_peaks_clustering <- function(dir_in, dir_out, samples){
  print(" ")
  print("  ///////////////////////////")
  print(" /    Clustering the ROIs  /")
  print("///////////////////////////")
  print(" ")

  load_samples <- function(dir_in, samples) {
    setwd(dir_in)
    num_samples <- length(samples)

    imgs <- vector("list", length = num_samples)
    peaks <- vector("list", length = num_samples)

    drift_times <- vector("list", length = num_samples)
    ret_times <-  vector("list", length = num_samples)

    sample_names <- paste0("M", samples)

    for (m in seq_along(samples)) {
      print(paste0("Sample ", m, " of ", num_samples))
      aux_list <- readRDS(paste0(sample_names[m], ".rds"))

      parameters_df_i <- as.data.frame(t(aux_list$data$Parameters))
      rownames(parameters_df_i) <- sprintf("Peak%d", seq_len(nrow(parameters_df_i)))

      rois_i <- aux_list$data$ROIs
      rownames(rois_i) <- sprintf("Peak%d", seq_len(nrow(rois_i)))
      colnames(rois_i) <- c("xmin", "xmax", "ymin", "ymax") # FIXME: Please check the order of the columns.

      apex_i <- aux_list$data$Peaks
      colnames(apex_i) <- c("xapex", "yapex")
      rownames(rois_i) <- sprintf("Peak%d", seq_len(nrow(apex_i)))

      # Maybe instead of the Peaks, ROIs and Parameters, we could just have one dataframe in the RDS:
      peaks_i <- cbind(apex_i, rois_i, parameters_df_i)
      imgs[[m]] <- as.matrix(aux_list$data$data_df)
      ret_times[[m]] <- aux_list$data$retention_time
      drift_times[[m]] <- aux_list$data$drift_time
      peaks[[m]] <- tibble::rownames_to_column(peaks_i, var = "PeakID")
    }
    names(imgs) <- sample_names
    names(peaks) <- sample_names
    peaks_df <- dplyr::bind_rows(peaks, .id = "SampleID")
    # Add a unique peak ID, that combines Sample+Peak ids:
    peaks_df <- tibble::add_column(
      peaks_df,
      UniqueID = paste0(peaks_df$SampleID, peaks_df$PeakID),
      .before = 1
    )
    list(
      imgs = imgs,
      ret_times = ret_times,
      drift_times = drift_times,
      peaks = tibble::as_tibble(peaks_df)
    )
  }
  loaded_samples <- load_samples(dir_in, samples)

  imgs <- loaded_samples$imgs
  peaks <- loaded_samples$peaks
  drift_times <- loaded_samples$drift_times
  ret_times <- loaded_samples$ret_times
  rm(loaded_samples, load_samples)

  # 1. Filter peaks with weird width or height
  roi_sizes <- tibble::tibble(
    UniqueID = peaks$UniqueID,
    width = peaks$xmax - peaks$xmin,
    height = peaks$ymax - peaks$ymin
  )
  quartiles_w <- stats::quantile(roi_sizes$width)
  lower_bound_iqr_w <- quartiles_w["75%"] - 1.5*(quartiles_w["75%"] - quartiles_w["25%"])
  higher_bound_iqr_w <- quartiles_w["75%"] + 1.5*(quartiles_w["75%"] - quartiles_w["25%"])
  quartiles_h <- stats::quantile(roi_sizes$height)
  lower_bound_iqr_h <- quartiles_h["75%"] - 1.5*(quartiles_h["75%"] - quartiles_h["25%"])
  higher_bound_iqr_h <- quartiles_h["75%"] + 1.5*(quartiles_h["75%"] - quartiles_h["25%"])

  peaks_to_exclude <- roi_sizes$UniqueID[
    roi_sizes$width < lower_bound_iqr_w | roi_sizes$width > higher_bound_iqr_w |
      roi_sizes$height < lower_bound_iqr_h | roi_sizes$height > higher_bound_iqr_h
  ]
  message(sprintf("Excluding %d/%d peaks", length(peaks_to_exclude), nrow(peaks)))
  peaks <- dplyr::filter(peaks, ! UniqueID %in% peaks_to_exclude)


  N_clusters <- max(by(peaks, peaks$SampleID, nrow))

  # Mahalanobis distance:
  # https://stats.stackexchange.com/a/81710/62083
  cholMaha <- function(X) {
    dec <- chol( stats::cov(X) )
    tmp <- forwardsolve(t(dec), t(X) )
    stats::dist(t(tmp))
  }
  peak2peak_D_mahal <- cholMaha(as.matrix(peaks[,c("xapex", "yapex")]))
  peak2peak_D_mahal_mat <- as.matrix(peak2peak_D_mahal)
  rownames(peak2peak_D_mahal_mat) <- peaks$UniqueID
  colnames(peak2peak_D_mahal_mat) <- peaks$UniqueID
  # Set distances from pairs of peaks belonging to the same sample to Inf,
  # so they are never in the same cluster
  max_dist <- max(peak2peak_D_mahal_mat)
  for (sampleid in unique(peaks$SampleID)) {
    peaks_i <- peaks$UniqueID[peaks$SampleID == sampleid]
    for (peak_i in peaks_i) {
      peak2peak_D_mahal_mat[peak_i, peaks_i] <- 10*nrow(peak2peak_D_mahal_mat)*max_dist
      peak2peak_D_mahal_mat[peaks_i, peak_i] <- 10*nrow(peak2peak_D_mahal_mat)*max_dist
      peak2peak_D_mahal_mat[peak_i, peak_i] <- 0
    }
  }
  peak2peak_D_mahal_inf_self_sample <- stats::as.dist(peak2peak_D_mahal_mat)
  rm(peak2peak_D_mahal_mat)

  cluster <- cluster::pam(x = peak2peak_D_mahal_inf_self_sample, k = N_clusters)

  peaks$cluster <- cluster$clustering

  ## Imputation and statistics

  indices_clusters <- unique(peaks$cluster)
  num_clusters <- length(indices_clusters)
  K <- length(indices_clusters)

  stats <- vector(mode = "list", length = num_clusters)

  median_roi_per_cluster <- peaks %>%
    dplyr::group_by(!!rlang::sym("clusters")) %>%
    dplyr::summarise(
      xmin = stats::median(!!rlang::sym("xmin")),
      xmax = stats::median(!!rlang::sym("xmax")),
      ymin = stats::median(!!rlang::sym("ymin")),
      ymax = stats::median(!!rlang::sym("ymax")),
    ) %>%
    dplyr::ungroup()

  peaks %>%
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "volume"))) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym("SampleID"), values_from = !!rlang::sym("volume"), values_fn = length)

  for (k in seq_len(num_clusters)) {
    n <- indices_clusters[k]
    pos_n <- idx_post[, 1] == n & idx_post[, 3] != 1
    samples_n <- idx_post[pos_n, 2]
    rest_samples_n <- setdiff((1:length(samples)), samples_n)
    ROI_median <- floor(apply(ROIs[pos_n, ], 2, stats::median))

    AsF_tmp <- rep(0, length(samples))
    AsF_tmp[samples_n] <- AsF[pos_n]
    AsF_tmp[rest_samples_n] <- stats::median(AsF[pos_n])

    volume_tmp <- rep(0, length(samples))
    volume_tmp[samples_n] <- volume[pos_n]

    I <- NULL
    if (length(rest_samples_n) >= 1){
      for (l in rest_samples_n){
        I[[l]] <- imgs[[l]][ROI_median[1]:ROI_median[2], ROI_median[3]:ROI_median[4]]
        volume_tmp[rest_samples_n] <- sum(apply(I[[l]], 1, sum)) #[1, dim(I[[l]], 3)]) # size que?
      }
    }

    ROI_tmp <- matrix(0, length(samples), 4)
    ROI_tmp[samples_n, ] <- as.matrix(ROIs[pos_n, ])
    ROI_tmp[rest_samples_n, ] <- matrix(rep(ROI_median, each = length(rest_samples_n)), nrow = length(rest_samples_n))

    peaks_tmp <- matrix(0, length(samples), 2)
    peaks_tmp[samples_n, ] <- as.matrix(datapoints[pos_n, ])
    peaks_tmp[rest_samples_n, ]  <- matrix(rep(floor(apply(datapoints[pos_n, ], 2, stats::median)), each = length(rest_samples_n)), nrow = length(rest_samples_n))

    stats[[k]]$Name <- paste("Cluster", k)
    stats[[k]]$AsF <- AsF_tmp
    stats[[k]]$volume <- volume_tmp

    stats[[k]]$AsF.mean  <- mean(AsF_tmp)
    stats[[k]]$volume_mean <- mean(volume_tmp)

    stats[[k]]$AsF_std <- sd(AsF_tmp)
    stats[[k]]$volume_std <- sd(volume_tmp)

    stats[[k]]$ROIs <- ROI_tmp
    stats[[k]]$peaks <- peaks_tmp

  }


  ## Make table

  peak_table <- matrix(0, (length(samples) * K), 7)

  acc <- 1
  for (i in (1:K)){
    c <- 1
    for (j in (acc : (acc + length(samples) - 1))){
      drift_time_idx <- stats[[i]]$peaks[c, 1]
      ret_time_idx <- stats[[i]]$peaks[c, 2]
      peak_table[j, 1] <- c
      peak_table[j, 2] <- i
      peak_table[j, 3] <- drift_time_idx
      peak_table[j, 4] <- ret_time_idx
      peak_table[j, 5] <- drift_times[[1]][drift_time_idx] # FIXME: ? Assuming homogeneus drift times across samples
      peak_table[j, 6] <- ret_times[[1]][ret_time_idx] # FIXME: ? Assuming homogeneous ret times across samples
      peak_table[j, 7] <- stats[[i]]$volume[c]

      c <- c + 1
    }
    acc <- acc + length(samples)
  }

  colnames(peak_table) <- c("Sample", "Cluster", "PicoX (indice)", "PicoY (indice)", "PicoX (tiempo)", "PicoY (tiempo)", "Intensity")
  setwd(dir_out)
  utils::write.csv(peak_table, "Peaktable.csv")
  saveRDS(stats, "ROIsParameters.rds")
}


