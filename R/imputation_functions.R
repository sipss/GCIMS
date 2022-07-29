#' Peak Table imputation using RF

#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param prop_samples    Mimimum proportion of hits per feature in the dataset to
#'                        consider a feature as valid.
#' @return A roi table for the whole dataset.
#' @family Imputation functions
#' @export
#' @importFrom utils capture.output
#' @importFrom tidyr spread
#' @importFrom missForest missForest
gcims_peak_imputation <- function(dir_in, dir_out, prop_samples){
  roi_cluster <- volume <- NULL

  print("Performing Peak Imputation")

  # 1) Read the roi table in long format
  roi_table_long <- readRDS(file.path(dir_in, "all_roi_df.rds"))

  # 2) Select only the interesting roi variables:
  roi_table_long <- roi_table_long[c("roi_cluster", "sample_id", "volume")]

  # 3) Don't include data from the reference sample
  roi_table_long <- roi_table_long[roi_table_long$sample_id != 0, ]

  # 4) Convert data fron long to wide format (this action generates missing values)
  roi_table_wide_mis  <- spread(roi_table_long , key = roi_cluster, value = volume)

  # 5) Remove features that appear less than minumum allowed proportion in the dataset
  cond_prop_samples <- round(nrow(roi_table_wide_mis) * prop_samples)
  reliable_columns <- colSums(!is.na(roi_table_wide_mis)) > cond_prop_samples

  # 6) Perform roi imputation using Random Forest
  roi_table_wide_mis <- roi_table_wide_mis[, reliable_columns]
  invisible(capture.output({roi_table_wide_imp <- missForest(roi_table_wide_mis)$ximp}))

  # 7) Save results
  saveRDS(roi_table_wide_imp, file = file.path(dir_out, "roi_table.rds"))
}




#' Missing Values Imputation

#' @param dir_in              Input directory. Where input data files are loaded
#'   from. It should be the one with the baseline removed.
#' @param peak_table_with_na  Peak table that contains missing values
#' @param peak_list_foms      Peak list containing all the figures of merit.
#' @details `gcims_missing_imputation` calculates the volume of the ROIs that
#' have missing values. It uses the coordinates of the reference ROI for each
#' cluster and subistitue the missing value by the volume in this region without
#' baseline.
#' @return A Set of S3 objects.
#' @family Imputation functions
#' @export
#' @examples
#' peak_table_imputed <- gcims_missing_imputation(dir_in = bslnr,
#' peak_table_with_na, peak_list_foms)
#' head(peak_table_imputed)
gcims_missing_imputation <- function(dir_in, peak_table_with_na, peak_list_foms){

  #-------------#
  #     MAIN    #
  #-------------#


  print(" ")
  print("  /////////////////////////////////")
  print(" /    Missing Values Imputation   /")
  print("///////////////////////////////////")
  print(" ")

  peaktable <- peak_table_with_na
  clusters_infor <- peak_list_foms
  num_clusts <- dim(peaktable)[1]
  num_samps <- dim(peaktable)[2] - 1
  setwd(dir_in)

  for(i in seq_len(num_samps)){
    missing_value <- which(is.na(peaktable[, i+1]) == TRUE)
    if (length(missing_value) >= 1){
      aux_list <- readRDS(paste0("M", i, ".rds")) # Load RDS file
      aux <- (as.matrix(aux_list$data$data_df)) # The data is in data_df
      for(j in missing_value) {
        #4 corrdinates
        dtmin_clust <- clusters_infor$ref_roi_dt_min_idx[j]
        dtmax_clust <- clusters_infor$ref_roi_dt_max_idx[j]
        rtmin_clust <- clusters_infor$ref_roi_rt_min_idx[j]
        rtmax_clust <- clusters_infor$ref_roi_rt_max_idx[j]
        patch <- aux[dtmin_clust:dtmax_clust, rtmin_clust:rtmax_clust]
        value2imput <- round(compute_integral2(patch), digits = 0)
        peaktable[j, i+1] <- value2imput
      }
    }
  }
  peak_table<- peaktable
  return(peak_table)
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


