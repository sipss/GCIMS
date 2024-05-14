#' Build a peak table
#'
#' @description Extract the volume of each ROI across samples to create a peak table.
#'
#' @param peak_list_clustered A peak list with clusters assigned. Also, you can create your own peak table
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
#' # Create your peak table from scratch:
#' pl <- data.frame(
#'   SampleID = c("S1", "S1", "S2", "S2"),
#'   cluster = c("Cluster1", "Cluster2", "Cluster1", "Cluster2"),
#'   Volume = c(10, 20, 8, 18)
#' )
#' peak_table <- peakTable(pl)
#'
#' peak_table$peak_table_matrix
#'
#' # You can use imputePeakTable() to fill in the missing values
#'
#' # If the clustering doesn't work great, you may end up with two peaks
#' # from the same sample on the same cluster. This does not make sense
#' # empirically, because it's either one or the other. In case of such
#' # ambiguity, peakTable() will give an error.
#' #
#' # If you want, you can override the error by taking the average volume
#' # of those ambiguous peaks, or the maximum, using,
#' # e.g. `aggregate_conflicting_peaks = max`.
#' #
#' # In any case, you will get information on how many peaks were aggregated
#' # in the `peak_table_duplicity` field (ideally should be full of `1`):
#' peak_table$peak_table_duplicity
#'
#'
#'
peakTable <- function(peak_list_clustered, aggregate_conflicting_peaks = NULL) {
  if (!"Volume" %in% colnames(peak_list_clustered)) {
    cli_abort("Please compute a 'Volume' column in peak_list_clustered")
  }

  peak_list_clustered <- peak_list_clustered |>
    dplyr::filter(!is.na("cluster"))

  peak_table_duplicity <- peak_list_clustered |>
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "Volume"))) |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of("SampleID"),
      values_from = dplyr::all_of("Volume"),
      values_fn = length
    )

  peak_table <- peak_list_clustered |>
    dplyr::select(dplyr::all_of(c("SampleID", "cluster", "Volume"))) |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of("SampleID"),
      values_from = dplyr::all_of("Volume"),
      values_fn = aggregate_conflicting_peaks
    )

  peak_table_mat <- peak_table |>
    tidyr::pivot_longer(cols = -1, names_to = "SampleID", values_to = "Volume") |>
    tidyr::pivot_wider(names_from = "cluster", values_from = "Volume") |>
    tibble::column_to_rownames("SampleID") |>
    as.matrix()


  # Missing values still need to be filled
  list(
    peak_table = peak_table,
    peak_table_matrix = peak_table_mat,
    peak_table_duplicity = peak_table_duplicity
  )
}


#' Omit ROIs present in certain retention and drift times
#'
#' @description Extract the volume of each ROI across samples to create a peak table.
#'
#' @param peak_list The output of [peaks()]. Also, you can create your own peak table
#' and use it as input value for `peak_list`
#' @param rt_time_2_omit A vector including a set of retention times where ROIs detected
#' should not be considered. As default is is set as `NULL`
#' @param dt_time_2_omit A vector including a set of drift times where ROIs detected
#' should not be considered. As default is is set as `NULL`
#'
#' @return A `peak_list` without the ROIs present in the retention and drift times
#' not desired.
#'
#' @export
#' @examples
#'
#' peak_list <- data.frame(
#'   rt_apex_s =  c(1, 2, 3, 3, 4, 4,  5, 5, 6, 6),
#'   dt_apex_ms = c(2, 4, 6, 4, 8, 4, 10, 4, 4, 12)
#' )
#' peak_list_filt <- omit_times(peak_list, dt_time_2_omit = 4)
#'
#'
omit_times <- function(peak_list, rt_time_2_omit = NULL, dt_time_2_omit = NULL){
  keep_roi <- rep(TRUE, nrow(peak_list))
  if(is.null(rt_time_2_omit) == FALSE){
    for (i in seq_along(rt_time_2_omit)) {
      rt_omit <- rt_time_2_omit[i]
      roi_2_omit <- which(peak_list$rt_apex_s == rt_omit)
      if(length(roi_2_omit) >= 1){
        keep_roi[roi_2_omit] <- FALSE
      }
    }
  }
  if(is.null(dt_time_2_omit) == FALSE){
    for (j in seq_along(dt_time_2_omit)) {
      dt_omit <- dt_time_2_omit[j]
      roi_2_omit <- which(peak_list$dt_apex_ms == dt_omit)
      if(length(roi_2_omit) >= 1){
        keep_roi[roi_2_omit] <- FALSE
      }
    }
  }
  peak_list[keep_roi,,drop = FALSE]
}
