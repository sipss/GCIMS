#' Peak integration for a GCIMSSample
#' @param object A [GCIMSSample] object
#' @param peak_list A data frame with the peak list
#' @param integration_size_method If "fixed_size", the ROI integration limits
#' are the same for all the peaks that belong to the same cluster. If "free_size",
#' each ROI has its own integration limits, regardless of the cluster it is assigned to.
#' @param rip_saturation_threshold Used to compute the "Saturation" column. If the ratio of the RIP intensity at the ROI apex
#' with respect to the maximum RIP is below this threshold, the RIP is considered almost depleted, and it's more likely that
#' the ROI suffers from non-linearities.
#' @param verbose If `TRUE`, debug information will be printed
#' @return The modified [GCIMSSample], with an updated peak list
#' @export
setMethod(
  "integratePeaks",
  "GCIMSSample",
  function(object, peak_list, integration_size_method = c("fixed_size", "free_size"), rip_saturation_threshold = 0.1, verbose = FALSE) {
    integration_size_method <- match.arg(integration_size_method)
    # Ensure peak_list only has peaks from object
    sample_name <- description(object)
    peak_list <- dplyr::filter(peak_list, .data$SampleID == sample_name)
    peak_list$Area <- NA_real_
    peak_list$Volume <- NA_real_
    peak_list$Asymmetry <- NA_real_
    peak_list$Saturation <- FALSE

    intmat <- intensity(object)
    rt <- rtime(object)
    dt <- dtime(object)
    dt_step_ms <- dt[2L] - dt[1L]
    rt_step_s <- rt[2L] - rt[1L]
    rt_saturated_regions_s <- find_regions_rip_saturated(
      intmat,
      rip_saturation_threshold = rip_saturation_threshold,
      verbose = verbose,
      retention_time = rt,
      drift_time = dt
    )

    # Area:
    len_dt_ms <- peak_list$dt_max_ms - peak_list$dt_min_ms
    len_rt_s <- peak_list$rt_max_s - peak_list$rt_min_s
    peak_list$Area <- len_dt_ms * len_rt_s

    # ROI asymmetries
    rt_rising_length <- peak_list$rt_apex_s - peak_list$rt_min_s
    rt_falling_length <- peak_list$rt_max_s - peak_list$rt_apex_s
    peak_list$Asymmetry <- round(rt_falling_length/rt_rising_length - 1, digits = 2)

    for (i in seq_len(nrow(peak_list))) {
      # ROI volume
      if (integration_size_method == "fixed_size") {
        dt_range <- c(peak_list$fixedsize_dt_min_ms[i], peak_list$fixedsize_dt_max_ms[i])
        rt_range <- c(peak_list$fixedsize_rt_min_s[i], peak_list$fixedsize_rt_max_s[i])
      } else if (integration_size_method == "free_size") {
        dt_range <- c(peak_list$dt_min_ms[i], peak_list$dt_max_ms[i])
        rt_range <- c(peak_list$rt_min_s[i], peak_list$rt_max_s[i])
      } else {
        stop("Invalid integration_limits parameter")
      }
      patch <- intensity(object, dt_range = dt_range, rt_range = rt_range)
      peak_list$Volume[i] <- sum(patch)*rt_step_s*dt_step_ms

      # ROI saturation
      # Rather than checking if the center of mass is saturated, shouldn't we check if "any region within "some thin peak boundaries" is saturated"?
      for (l in seq_len(nrow(rt_saturated_regions_s))) {
        if (rt_saturated_regions_s[l, 1L] < peak_list$rt_cm_s[i] && rt_saturated_regions_s[l, 2L] > peak_list$rt_cm_s[i]) {
          peak_list$Saturation[i] <- TRUE
          break
        }
      }
    }
    peaks(object) <- peak_list
    object
  }
)
