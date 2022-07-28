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
#' roi_selection <- tempfile("dir")
#' fom <- tempfile("dir")
#'
#' # Example of Calculating the Figures of Merit
#' gcims_rois_selection(dir_in, roi_selection, samples = c(3, 7), noise_level = 3)
#' peak_list_fom <- gcims_figures_of_merit(dir_in = roi_selection, dir_out = fom, samples = 3)
#' head(peak_list_fom)
gcims_figures_of_merit <- function(dir_in, dir_out, samples, peak_list){
  print(" ")
  print("  ////////////////////////////////////////")
  print(" /    Calculating the Figures of Merit  /")
  print("/////////////////////////////////////////")
  print(" ")

  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  peak_list$Area <- rep(0, dim(peak_list)[1])
  peak_list$Volume <- rep(0, dim(peak_list)[1])
  peak_list$AsF <- rep(0, dim(peak_list)[1])
  peak_list$Saturation <- rep(0, dim(peak_list)[1])


  s = 0

  peaktables <- list()
  for (i in samples){
    s = s + 1
    print(paste0("Sample ", s, " of ", length(samples)))

    # 1. Data load

    aux_string <- paste0("M", i, ".rds") # Generate file name
    aux_list <- readRDS(file.path(dir_in, aux_string)) # Load RDS file
    aux <- (as.matrix(aux_list$data$data_df))
    ROIs <- peak_list[(peak_list$SampleID == aux_string), ]
    if (!"ROIs" %in% names(aux_list$data)) {
      rlang::abort("Please run gcims_rois_selection step first. Check the vignette")
    }


    # 2. Search of RIP position

    total_ion_spectrum <- rowSums(aux) # Sum per rows
    rip_position <- which.max(total_ion_spectrum) # Find maximum for every column
    rt_idx_with_max_rip <- which.max(aux[rip_position,])
    minima <- pracma::findpeaks(-total_ion_spectrum)[, 2] # Find local minima
    rip_end_index <- minima[min(which((minima - rip_position) > 0))] # Find ending index of RIP
    rip_start_index <- minima[max(which((rip_position - minima) > 0))] # Find starting index of RIP

    labels <- nrow(ROIs)
    AsF <- numeric(nrow(ROIs))
    volume <- numeric(nrow(ROIs))
    area <- numeric(nrow(ROIs))
    saturation <- logical(nrow(ROIs))

    # Search saturation regions
    rip_chrom <- rowSums(aux[rip_start_index: rip_end_index, ]) / length(rip_start_index: rip_end_index)
    max_rip_chrom <- max(rip_chrom)
    saturation_threshold <- 0.1 * max_rip_chrom

    for (n in seq_len(labels)){
      R1 <- ROIs[n, ]

      len_dt <- as.numeric(R1["dt_max_idx"] - R1["dt_min_idx"])
      len_rt <- as.numeric(R1["rt_max_idx"] - R1["rt_min_idx"])

      # roi area
      area_roi <- len_rt * len_dt
      area[n] <- area_roi

      # roi volume
      patch <- aux[as.numeric(R1["dt_min_idx"]):as.numeric(R1["dt_max_idx"]),
                   as.numeric(R1["rt_min_idx"]):as.numeric(R1["rt_max_idx"])]

      volume[n] <- round(compute_integral2(patch), digits = 0)


      # roi asymmetries
      half_down_area  <- as.numeric(R1["rt_cm_s"] - R1["rt_min_s"])
      half_up_area    <- as.numeric(R1["rt_max_s"] - R1["rt_cm_s"])
      asymetry <- round(((half_down_area - half_up_area) - 1 ), 2)
      AsF[n] <- asymetry

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
        saturation[n] <- FALSE
      } else {
        for (l in (1:nrow(saturation_minima))){
          if ((saturation_minima[l, 1] < y_cm)
              & (saturation_minima[l, 2] > y_cm)) {
            saturation[n] <- TRUE
            break
          }
        }
      }
    }

    peaktable <- cbind(ROIs, Area = area, Volume = volume, AsF = AsF, Saturation = saturation)
    aux_list$data$Peaktable <- peaktable
    peaktables <- c(peaktables, list(peaktable))
    utils::write.csv(peaktable, file = file.path(dir_out, paste0("PeakTable", i, ".csv")))
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
    peak_list[(peak_list$SampleID == aux_string), "Area"] <- area
    peak_list[(peak_list$SampleID == aux_string), "Volume"] <- volume
    peak_list[(peak_list$SampleID == aux_string), "AsF"] <- AsF
    peak_list[(peak_list$SampleID == aux_string), "Saturation"] <- saturation

  }
  return(peak_list)
  names(peaktables) <- samples
  dplyr::bind_rows(peaktables, .id = "SampleID")
}


#----------------------#
#   compute_integral2  #
#----------------------#

compute_integral2 <- function(data){

  # Set up dimensions and integral limits
  n <- dim(data)[1]
  m <- dim(data)[2]
  xa <- 1
  xb <- n
  ya <- 1
  yb <- m

  # Set up Gauss-Legendre Method
  cx <- pracma::gaussLegendre(n, xa, xb)
  x <- cx$x
  wx <- cx$w
  cy <- pracma::gaussLegendre(m, ya, yb)
  y <- cy$x
  wy <- cy$w

  # Compute the integral
  I <- 0
  for (i in 1:n) {
    for (j in 1:m) {
      I <- I + wx[i] * wy[j] * data[x[i], y[j]]
    }
  }
  return(I)
}




