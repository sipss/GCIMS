#' Figures Of Merit Calculation
#'
#' Calculates the area, volume asymmetry and saturation of each peak roi in `peak_list`
#'
#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where files containing the Figures of Meir of each ROI are
#'   stored.
#' @param peak_list A data frame. The peak list where we will add the figures of merit.
#' @param cluster_stats A data frame with cluster statistics (ROI limits...)
#' @param integration_size Either "cluster_roi" or "individual_roi". When computing the volume,
#' the integration size can be set individually for each sample or use the reference cluster size given by `cluster_stats`.
#' @param rip_saturation_threshold A number. The fraction of the maximum RIP. If the RIP at the ROI is below that
#' fraction, we will consider the peak to be saturated.
#' @return The given peak list, with the added columns
#' @family Utility functions
#' @export
#' @examples
#' \donttest{
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' roi_fusion_out <- readRDS(file.path(dir_in, "peak_clustering.rds"))
#'
#' peak_list_fom <- gcims_figures_of_merit(
#'   dir_in = dir_in,
#'   dir_out = dir_out,
#'   peak_list = roi_fusion_out$peak_list_clustered,
#'   cluster_stats = roi_fusion_out$cluster_stats
#' )
#' head(peak_list_fom)
#'
#' files_fom <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files_fom))
#'}
gcims_figures_of_merit <- function(
    dir_in,
    dir_out,
    peak_list,
    cluster_stats,
    integration_size = c("cluster_roi", "individual_roi"),
    rip_saturation_threshold = 0.1
) {

  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  integration_size <- match.arg(integration_size)

  peak_list$Area <- 0
  peak_list$Volume <- 0
  peak_list$Asymmetry <- 0
  peak_list$Saturation <- 0


 # s = 0

  sample_names <- unique(peak_list$SampleID)

  for (sample_name in sample_names) {
    #s = s + 1

    # 1. Data load

    aux_string <- sample_name
    aux_list <- readRDS(file.path(dir_in, aux_string)) # Load RDS file
    aux <- as.matrix(aux_list$data$data_df)
    peak_list_rows <- which(peak_list$SampleID == sample_name)

    # 2. Find the retention time regions where the RIP is saturated:
    find_regions_rip_saturated <- function(aux) {
      the_rip <- find_rip(aux)
      # Search saturation regions
      rip_region <- aux[the_rip$dt_idx_start:the_rip$dt_idx_end, , drop = FALSE]
      rip_chrom <- rowSums(rip_region) / nrow(rip_region)
      rt_rip_saturated_indices <- which(rip_chrom <= rip_saturation_threshold * max(rip_chrom))
      if (length(rt_rip_saturated_indices) == 0L) {
        rt_saturated_regions <- matrix(nrow = 0L, ncol = 2L)
      } else {
        saturation_list <- split(rt_rip_saturated_indices, cumsum(c(1, diff(rt_rip_saturated_indices)) != 1))
        rt_saturated_regions <- matrix(0, nrow = length(saturation_list), ncol = 2)
        for (k in seq_along(saturation_list)) {
          rt_saturated_regions[k, 1L] <- min(saturation_list[[k]])
          rt_saturated_regions[k, 2L] <- max(saturation_list[[k]])
        }
      }
      colnames(rt_saturated_regions) <- c("begin", "end")
      rt_saturated_regions
    }
    rt_saturated_regions <- find_regions_rip_saturated(aux)

    for (peak_list_row in peak_list_rows) {
      roi_prop <- as.list(peak_list[peak_list_row,])
      cluster_id <- roi_prop$cluster
      cluster_prop <- as.list(cluster_stats[cluster_stats$cluster == cluster_id,])

      len_dt <- roi_prop$dt_max_idx - roi_prop$dt_min_idx
      len_rt <- roi_prop$rt_max_idx - roi_prop$rt_min_idx
      peak_list$Area[peak_list_row] <- len_dt * len_rt

      # ROI volume
      if (integration_size == "cluster_roi") {
        # Coordinates to integrate:
        # FIXME: The minmax stuff here should not be needed when fusion_rois
        # is able to verify the original index limits.
        dt_min_idx <- max(1L,        roi_prop$ref_roi_dt_min_idx)
        dt_max_idx <- min(nrow(aux), roi_prop$ref_roi_dt_max_idx)
        rt_min_idx <- max(1L,        roi_prop$ref_roi_rt_min_idx)
        rt_max_idx <- min(ncol(aux), roi_prop$ref_roi_rt_max_idx)
        patch <- aux[dt_min_idx:dt_max_idx,
                     rt_min_idx:rt_max_idx]
      } else if (integration_size == "individual_roi") {
        patch <- aux[roi_prop$dt_min_idx:roi_prop$dt_max_idx,
                     roi_prop$rt_min_idx:roi_prop$rt_max_idx]
      } else {
        stop("Invalid integration_limits parameter")
      }
      peak_list$Volume[peak_list_row] <- compute_integral2(patch)

      # ROI asymmetries
      rt_rising_length <- roi_prop$rt_apex_s - roi_prop$rt_min_s
      rt_falling_length <- roi_prop$rt_max_s - roi_prop$rt_apex_s
      peak_list$Asymmetry[peak_list_row] <- round(rt_falling_length/rt_rising_length - 1, digits = 2)


      # ROI saturation
      for (l in seq_len(nrow(rt_saturated_regions))) {
        if (rt_saturated_regions[l, 1L] < roi_prop$rt_cm_idx && rt_saturated_regions[l, 2L] > roi_prop$rt_cm_idx) {
          peak_list$Saturation[peak_list_row] <- TRUE
          break
        }
      }
    }
    peak_list_this_sample <- peak_list[peak_list_rows, , drop = FALSE]
    aux_list$data$Peaktable <- peak_list_this_sample
    utils::write.csv(
      peak_list_this_sample,
      file = file.path(dir_out, paste0("PeakTable", tools::file_path_sans_ext(sample_name), ".csv"))
    )
    saveRDS(aux_list, file = file.path(dir_out, sample_name))
  }
  return(peak_list)
}


#----------------------#
#   compute_integral2  #
#----------------------#

compute_integral2 <- function(data){

  # Set up dimensions and integral limits
  n <- dim(data)[1]
  m <- dim(data)[2]
  xa <- 1
  xb <- n
  ya <- 1
  yb <- m

  # Set up Gauss-Legendre Method
  cx <- pracma::gaussLegendre(n, xa, xb)
  x <- cx$x
  wx <- cx$w
  cy <- pracma::gaussLegendre(m, ya, yb)
  y <- cy$x
  wy <- cy$w

  # Compute the integral
  I <- 0
  for (i in 1:n) {
    for (j in 1:m) {
      I <- I + wx[i] * wy[j] * data[x[i], y[j]]
    }
  }
  return(I)
}





