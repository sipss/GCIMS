get_inv_k0 <- function(drift_time,
                       drift_tube_length_cm,
                       pressure_ims,
                       voltage_ims,
                       ims_temp_k){
  k0 <- (0.0027315 * (drift_tube_length_cm ^ 2) * pressure_ims) /
    (voltage_ims * (drift_time / 1000) * ims_temp_k)
  inv_k0 <- 1 / k0
}


.extract_inv_k0 <- function(object) {
  inv_k0 <- object@inverse_reduced_mobility
  data.frame(
    inv_k0_min = inv_k0[1],
    inv_k0_max = inv_k0[length(inv_k0)],
    inv_k0_step = stats::median(diff(inv_k0))
  )
}

.extract_inv_k0_fun_aggregate <- function(object, extracted_objects) {
  inv_k0_samples <- dplyr::bind_rows(extracted_objects, .id = "SampleID")

  max_inv_k0_min <- max(inv_k0_samples$inv_k0_min)
  min_inv_k0_max <- min(inv_k0_samples$inv_k0_max)
  min_inv_k0_step <- min(inv_k0_samples$inv_k0_step)
  max_inv_k0_length <- floor(round((min_inv_k0_max - max_inv_k0_min)/min_inv_k0_step, digits = 8)) + 1L
  inv_k0_ref <- seq(from = max_inv_k0_min, to = min_inv_k0_max, length.out = max_inv_k0_length)

  object$inv_k0_ref <- inv_k0_ref
}
