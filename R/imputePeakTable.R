#' Impute a Peak table
#'
#' @param peak_table A matrix, with samples in rows and clusters in columns. It must have row names and column names.
#' @param dataset The dataset object to extract samples from
#' @param cluster_stats A data frame with the `[dt|rt]_[min|max]_[ms|s]` columns
#'
#' @return The imputed peak_table
#' @examples
#'
#' # We are going to create a peak table matrix, typically resulting from [peakTable()]
#' # The peak table may have some missing values
#' # Since the missing values correspond to peaks that have not been detected
#' # in those particular samples, we can integrate the region where they should appear
#' # to get a value different than zero that reflects the noise level.
#' #
#' # Ingredients:
#' # - The GCIMSSample objects, so we can integrate the regions of interest (given as a GCIMSDataset)
#' # - The peak table matrix we want to impute
#' # - The definition of the regions corresponding to each cluster (cluster_stats)
#' #
#' # We will prepare here a synthetic example, check the vignette for a real use case
#' #
#' # Imagine we have information on the location of Cluster1 and Cluster2
#' cluster_stats <- data.frame(
#'   cluster = c("Cluster1", "Cluster2"),
#'   dt_min_ms = c(8, 10),
#'   dt_max_ms = c(9, 12),
#'   rt_min_s = c(120, 300),
#'   rt_max_s = c(128, 320)
#' )
#'
#' # We have a peak table for two samples and two peaks
#' peak_table <- matrix(NA_real_, nrow = 2, ncol = 2)
#' rownames(peak_table) <- c("Sample1", "Sample2")
#' colnames(peak_table) <- c("Cluster1", "Cluster2")
#'
#' # where we previously integrated Cluster2 in Sample 1 and Cluster1 in Sample2:
#' peak_table["Sample1", "Cluster2"] <- 9.5
#' peak_table["Sample2", "Cluster1"] <- 3.6
#' # The table has missing values, because some peaks were not detected.
#' # Maybe they are close to the noise level, or maybe they do not exist
#' peak_table
#'
#' # We will fill the missing values by integrating whatever we find
#' # (typically noise or small peaks) in the cluster regions of each sample. So we
#' # need the sample matrices.
#' #
#' # Let's build dummy Sample1 and Sample2:
#' ## Create drift time and retention time vectors:
#' dt <- seq(from = 0, to = 13, by = 0.1)  # ms
#' rt <- seq(from = 0, to = 350, by = 1)   # s
#'
#' ## Create matrices with random gaussian noise.
#' set.seed(42)
#' s1_intensity <- matrix(
#'   rnorm(length(dt)*length(rt), sd = 0.1),
#'   nrow = length(dt),
#'   ncol = length(rt)
#' )
#' s2_intensity <- matrix(
#'   rnorm(length(dt)*length(rt), sd = 0.1),
#'   nrow = length(dt),
#'   ncol = length(rt)
#' )
#'
#' # The matrix will have a pleateau in a region where the peak is supposed
#' # to be, so when we impute the region corresponding to Sample1-Cluster1 we see a
#' # higher value:
#' s1_intensity[dt > 8.25 & dt < 8.75, rt > 122 & rt < 126] <- 1
#'
#' ## Create GCIMSSample objects
#' s1 <- GCIMSSample(
#'   drift_time = dt,
#'   retention_time = rt,
#'   data = s1_intensity
#' )
#' s2 <- GCIMSSample(
#'   drift_time = dt,
#'   retention_time = rt,
#'   data = s2_intensity
#' )
#' ## And a dataset with the samples:
#' dataset <- GCIMSDataset_fromList(list(Sample1 = s1, Sample2 = s2))
#'
#' # Now we can impute the table
#' peak_table_imp <- imputePeakTable(
#'   peak_table = peak_table,
#'   dataset = dataset,
#'   cluster_stats = cluster_stats
#' )
#' peak_table_imp
#' @export
imputePeakTable <- function(peak_table, dataset, cluster_stats) {
  num_samples <- nrow(peak_table)
  for (i in seq_len(num_samples)) {
    sample_name <- rownames(peak_table)[i]
    sample <- dataset$getSample(sample = sample_name)
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
