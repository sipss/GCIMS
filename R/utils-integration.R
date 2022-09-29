find_regions_rip_saturated <- function(aux, rip_saturation_threshold, verbose = FALSE, retention_time = NULL, drift_time = NULL) {
  the_rip <- find_rip(aux, verbose = verbose, retention_time = retention_time, drift_time = drift_time)
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
  if (!is.null(retention_time)) {
    rt_saturated_regions_s <- matrix(0.0, nrow = length(saturation_list), ncol = 2L)
    rt_saturated_regions_s[,1L] <- retention_time[rt_saturated_regions[,1L]]
    rt_saturated_regions_s[,2L] <- retention_time[rt_saturated_regions[,2L]]
    return(rt_saturated_regions_s)
  } else {
    return(rt_saturated_regions)
  }
}

