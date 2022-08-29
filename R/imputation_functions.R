
#' Impute missing values
#'
#' Compute the volume of ROIS with missing values.
#'
#' @param peak_table The matrix to be imputed, with samples in rows (see example)
#' @param dir_in Directory with baseline removed samples
#' @param cluster_stats Cluster statistics (reference ROI limits)
#' @details `gcims_impute_missing()` uses the coordinates of the reference ROI for each
#' cluster and substitute the missing value by the volume in this region without
#' baseline.
#' @return A matrix with samples in rows, clusters in columns and volumes
#'  as values, without missing values
#' @family Imputation functions
#' @export
#' @examples
#' \donttest{
#' dir_in <- system.file("extdata", package = "GCIMS")
#' peak_table_obj <- readRDS(file.path(dir_in, "peak_table.rds"))
#' roi_merging_out <- readRDS(file.path(dir_in, "roi_fusion_out.rds"))
#' peak_table_to_impute <- peak_table_obj$peak_table_mat
#'
#' # Peak table before imputation
#' peak_table_to_impute
#'
#'
#' peak_table_imputed <- gcims_impute_missing(peak_table = peak_table_to_impute,
#'                                                dir_in = dir_in, # use bsln
#'                                                cluster_stats = roi_merging_out$cluster_stats
#'                                                )
#' # Peak table after imputation
#' peak_table_imputed
#' }
#'
gcims_impute_missing <- function(peak_table, dir_in, cluster_stats){

  #-------------#
  #     MAIN    #
  #-------------#



  num_samps <- nrow(peak_table)

  for(i in seq_len(num_samps)){
    sample_name <- rownames(peak_table)[i]
    missing_values <- which(is.na(peak_table[i,]))
    if (length(missing_values) >= 1) {
      aux_list <- readRDS(file.path(dir_in, sample_name)) # Load RDS file
      aux <- as.matrix(aux_list$data$data_df) # The data is in data_df
      for (j in missing_values) {
        #4 coordinates
        cluster_name <- colnames(peak_table)[j]
        cluster_info_row <- which(cluster_stats$cluster == cluster_name)
        dtmin_clust <- cluster_stats$dt_min_idx[cluster_info_row]
        dtmax_clust <- cluster_stats$dt_max_idx[cluster_info_row]
        rtmin_clust <- cluster_stats$rt_min_idx[cluster_info_row]
        rtmax_clust <- cluster_stats$rt_max_idx[cluster_info_row]
        patch <- aux[dtmin_clust:dtmax_clust, rtmin_clust:rtmax_clust]
        peak_table[i, j] <- round(compute_integral2(patch), digits = 0)
      }
    }
  }
  peak_table
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

