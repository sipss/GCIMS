#' Figures Of Merit Calculation

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where files containing the Figures of Meir of each ROI are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to which their Figures of Merit (FOM) have to be caulated.
#' @details `gcims_figures_of_merit` calculates a set of figures of merit for each
#' ROI of the sample. The FOMs calculated are: the area, the volume, the assymetry
#' and the saturation.
#' @return A Set of S3 objects and a ".csv" file with a table that contains the figures
#' of merit for each ROI and the ROI information.
#' @family Utility functions
#' @export
#' @examples
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example of Calculating the Figures of Merit
#' gcims_figures_of_merit(dir_in = roispp, dir_out = fom, samples = (1:15))
#' invisible(file.remove(files))
gcims_figures_of_merit <- function(dir_in, dir_out, samples){
  print(" ")
  print("  ////////////////////////////////////////")
  print(" /    Calculating the Figures of Merit  /")
  print("/////////////////////////////////////////")
  print(" ")

  setwd(dir_in)
  s = 0
  for (i in samples){
    s = s + 1
    print(paste0("Sample ", s, " of ", length(samples)))

    # 1. Data load

    aux_string <- paste0("M", i, ".rds") # Generate file name
    aux_list <- readRDS(aux_string) # Load RDS file
    aux <- (as.matrix(aux_list$data$data_df))
    ROIs <- aux_list$data$ROIs

    # 2. Search of RIP position

    total_ion_spectrum <- rowSums(aux) # Sum per rows
    rip_position <- which.max(total_ion_spectrum) # Find maximum for every column
    rt_idx_with_max_rip <- which.max(aux[rip_position,])
    minima <- pracma::findpeaks(-total_ion_spectrum)[, 2] # Find local minima
    rip_end_index <- minima[min(which((minima - rip_position) > 0))] # Find ending index of RIP
    rip_start_index <- minima[max(which((rip_position - minima) > 0))] # Find starting index of RIP

    labels <- dim(ROIs)[1]
    AsF <- NULL
    volume <- NULL
    area <- NULL
    saturation <- rep(0, length(labels))

    # Search saturation regions
    rip_chrom <- rowSums(aux[rip_start_index: rip_end_index, ]) / length(rip_start_index: rip_end_index)
    max_rip_chrom <- max(rip_chrom)
    saturation_threshold <- 0.1 * max_rip_chrom

    for (n in seq_len(labels)){
      R1 <- ROIs[n, ]

      len_dt <- R1[5] - R1[4]
      len_rt <- R1[7] - R1[6]

      # roi area
      area_roi <- len_rt * len_dt
      area <- c(area, area_roi)

      # roi volume
      volume <- c(volume, sum(aux[R1[4]:R1[5], R1[6]:R1[7]]))

      # roi center of mass
      ind <- arrayInd(which.max(aux[R1[4]:R1[5], R1[6]:R1[7]]),
                      dim(aux[R1[4]:R1[5], R1[6]:R1[7]]))
      x_cm <- R1[4] + ind[1,1] - 1L
      y_cm <- R1[6] + ind[1,2] - 1L
      x_cm <- round(compute_integral2(aux[R1[4]:R1[5], R1[6]:R1[7]] * (R1[5] - R1[4])) / volume[n])
      y_cm <- round(compute_integral2(aux[R1[4]:R1[5], R1[6]:R1[7]] * (R1[7] - R1[6])) / volume[n])

      # roi asymmetries
      half_down_area  <- length(1:y_cm) * len_rt
      half_up_area    <- length(y_cm:len_dt) * len_rt
      asymetry <- round(((half_down_area - half_up_area) / (area_roi)), 2)
      AsF <- c(AsF, asymetry)

      saturation_regions <- which(rip_chrom <= saturation_threshold)
      if (length(saturation_regions) == 0){
        saturation_minima <- NULL
      } else {
        saturation_list <- split(saturation_regions, cumsum(c(1, diff(saturation_regions)) != 1))
        saturation_minima <- matrix(0, length(saturation_list), 2)
        for (k in seq_along(saturation_list)){
          saturation_minima[k, ] <- c(min(saturation_list[[k]]), max(saturation_list[[k]]))
        }
      }

      # roi saturation
      if (length(saturation_minima) == 0){
        # No saturation. Do nothing
        saturation[n] <- 0
      } else {
        for (l in (1:dim(saturation_minima)[1])){
          if ((saturation_minima[l, 1] < y_cm)
              & (saturation_minima[l, 2] > y_cm)) {
            saturation[n] <- 1
            break
          }
        }
      }
    }

    saturation[is.na(saturation)] <- 0
    peaktable <- cbind(ROIs, area, volume, AsF, saturation)
    colnames(peaktable) <- c("ID", "ApexDT", "ApexRT", "minDT", "maxDT", "minRT", "maxRT", "CenterMassDT", "CenterMassRT", "Area", "Volume", "AsF", "Saturation")
    aux_list$data$FOM <- rbind(AsF, saturation, volume)
    aux_list$data$Peaktable <- peaktable
    setwd(dir_out)
    write.csv(peaktable, file = paste0("PeakTable", i, ".csv"))
    M <- aux_list
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}



