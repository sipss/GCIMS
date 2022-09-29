#' Impute a Peak table
#'
#' @param peak_table A matrix, with samples in rows and clusters in columns. It must have row names and column names.
#' @param dataset The dataset object to extract samples from
#' @param cluster_stats A data frame with the `dt_min_ms` and family columns
#'
#' @return The imputed peak_table
#' @export
imputePeakTable <- function(peak_table, dataset, cluster_stats) {
  num_samples <- nrow(peak_table)
  for (i in seq_len(num_samples)) {
    sample_name <- rownames(peak_table)[i]
    sample <- getSample(dataset, sample = sample_name)
    intmat <- intensity(sample)
    dt <- dtime(sample)
    rt <- rtime(sample)
    dt_step_ms <- dt[2L] - dt[1L]
    rt_step_s <- rt[2L] - rt[1L]
    for (j in which(is.na(peak_table[i,,drop = TRUE]))) {
      cluster_name <- colnames(peak_table)[j]
      cluster_info_row <- which(cluster_stats$cluster == cluster_name)
      dtmin_ms <- cluster_stats$dt_min_ms[cluster_info_row]
      dtmax_ms <- cluster_stats$dt_max_ms[cluster_info_row]
      rtmin_s <- cluster_stats$rt_min_s[cluster_info_row]
      rtmax_s <- cluster_stats$rt_max_s[cluster_info_row]
      patch <- intmat[
        dt >= dtmin_ms & dt <= dtmax_ms,
        rt >= rtmin_s & rt <= rtmax_s,
        drop = FALSE
      ]
      peak_table[i, j] <- sum(patch)*rt_step_s*dt_step_ms

    }
  }
  peak_table
}
