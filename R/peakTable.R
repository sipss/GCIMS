#' Build a peak table
#'
#' @description Extract the volume of each ROI across samples to create a peak table.
#'
#' @param peak_list_clustered The output of [gcims_compute_fom()]. Also, you can create your own peak table
#' and use it as input value for `peak_list_clustered` (see first example below)
#' @param aggregate_conflicting_peaks `NULL` or a function. What to do, in case two peaks from the same sample
#' have been assigned to the same cluster. If `NULL`, throw an error. If `mean`, `max` or any other function,
#' we will summarize all the conflicting volumes into that number (e.g. "take the maximum of the peaks")
#'
#' @return A list with three fields: `peak_table`, `peak_table_matrix`, and `peak_table_duplicity`.
#' `peak_table`, and `peak_table_matrix`, provide information of the peak table. `peak_table` is a dataframe
#' containing cluster volumes, whose columns represent samples and rows clusters. `peak_table_matrix` presents
#' the same information content as `peak_table` but in matrix form. Note that in `peak_table` columns represent
#' clusters and rows samples. Finally, `peak_table_duplicity` is a dataframe that shows ROI duplicity information
#' among clusters. Ideally, only one peak per sample should belong to a cluster.
#'
#' @export
#' @examples
#' \donttest{
#' # Create your peak table from scratch:
#' pl <- data.frame(
#'   SampleID = c("S1", "S1", "S2", "S2"),
#'   cluster = c(1, 2, 1, 2),
#'   Volume = c(10, 20, 8, 18)
#' )
#' peakTable(pl)
#'
#' # Create a peak table from the output of the function gcims_compute_fom()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' peak_list_fom <- readRDS(file.path(dir_in, "peak_list_fom.rds"))
#' peak_table <- peakTable(peak_list_fom, aggregate_conflicting_peaks = max)
#'
#' peak_table$peak_table_matrix
#' }
#'
peakTable <- function(peak_list_clustered, aggregate_conflicting_peaks = NULL) {
  if (!"Volume" %in% colnames(peak_list_clustered)) {
    abort("Please compute a 'Volume' column in peak_list_clustered")
  }

  peak_table_duplicity <- peak_list_clustered %>%
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "Volume"))) %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of("SampleID"),
      values_from = dplyr::all_of("Volume"),
      values_fn = length
    )

  peak_table <- peak_list_clustered %>%
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "Volume"))) %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of("SampleID"),
      values_from = dplyr::all_of("Volume"),
      values_fn = aggregate_conflicting_peaks
    )

  peak_table_mat <- peak_table %>%
    tidyr::pivot_longer(cols = -1, names_to = "SampleID", values_to = "Volume") %>%
    tidyr::pivot_wider(names_from = "cluster", values_from = "Volume") %>%
    tibble::column_to_rownames("SampleID") %>%
    as.matrix()


  # Missing values still need to be filled
  list(
    peak_table = peak_table,
    peak_table_matrix = peak_table_mat,
    peak_table_duplicity = peak_table_duplicity
  )
}
