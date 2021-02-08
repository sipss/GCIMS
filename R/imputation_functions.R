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
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_peak_imputation <- function(dir_in, dir_out, prop_samples){
  roi_cluster <- volume <- NULL

  # 0) Set working directory
  setwd(dir_in)
  print("Performing Peak Imputation")

  # 1) Read the roi table in long format
  roi_table_long <- readRDS("all_roi_df.rds")

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
  invisible(capture.output(roi_table_wide_imp <- missForest(roi_table_wide_mis)$ximp))

  # 7) Save results
  setwd(dir_out)
  saveRDS(roi_table_wide_imp, file = "roi_table.rds")

  # 8) Go to the input directory
  setwd(dir_in)
}
