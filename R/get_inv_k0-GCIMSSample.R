methods::setMethod(
  "get_inv_k0",
  "GCIMSSample",
  function(object)  {
    if (is.null(object@params$`nom Drift Potential Difference`$value) |
        is.null(object@params$`Temp 1 setpoint`$value)|
        is.null(object@params$`Pressure Ambient`$value) |
        length(object@drift_tube_length) == 0
        ){
      stop("IMS Voltage, IMS Tube length, IMS Tempreature and pressure are
      needed to calculate the inverse reduced mobility ")
    }
    drift_time <- object@drift_time / 1000 # in ms to s
    Voltage_ims <- object@params$`nom Drift Potential Difference`$value # in V
    tube_length_cm <- object@drift_tube_length / 10 # in mm to cm
    ims_temp_k <- object@params$`Temp 1 setpoint`$value + 273.15 # in °C
    pressure_ims <- mean(as.numeric(strsplit(
      object@params$`Pressure Ambient`$value, " ")[[1]]), na.rm = TRUE) #in Pa
    K0 <- (0.0027315 * (tube_length_cm ^ 2) * pressure_ims) /
          (Voltage_ims* drift_time * ims_temp_k)
    inv_k0 <- 1 / K0
    object@inverse_reduced_mobility <- inv_k0
    invisible(object)
  }
)
