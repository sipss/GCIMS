#' ROIs Selection

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where RIP removed data files are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to which their RIP has to be removed.
#' @details \code{gcims_remove_rip} substitutes the RIP by its corresponding
#'   linear approximation to the RIP baseline, for every spectrum in a sample.
#'   This process is repeated for all samples in \code{samples}. Use this
#'   function if you are interested in enhancing the contrast of peaks of sample
#'   images / chromatograms / spectra to be obtained from
#'   \code{gcims_view_sample} / \code{gcims_plot_chrom} /
#'   \code{gcims_plot_spec}.
#' @return A Set of S3 objects.
#' @family Utility functions
#' @export
#' @importFrom signal sgolayfilt
#' @importFrom pracma findpeaks
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example of Reactant Ion Peak removal
#' # Before:
#' gcims_view_sample(dir_in, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' # After:
#' gcims_remove_rip(dir_in, dir_out, samples)
#' gcims_view_sample(dir_out, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)



gcims_rois_selection <- function(dir_in, dir_out, samples, noise_level){
  print(" ")
  print("  //////////////////////////")
  print(" /    Selecting the ROIs  /")
  print("//////////////////////////")
  print(" ")

  setwd(dir_in)
  s = 0
  for (i in samples){
    s = s + 1
    print(paste0("Sample ", s, " of ", length(samples)))

    # 1. Data load

    aux_string <- paste0("M", i, ".rds") # Generate file name
    aux_list <- readRDS(aux_string) # Load RDS file
    aux <- (as.matrix(aux_list$data$data_df)) # The data is in data_df

    drift_time <- aux_list$data$drift_time # Extract drift_time from file
    ret_time <- aux_list$data$retention_time # Extract ret_time from file
    fs = c(1/(drift_time[2]-drift_time[1]),1/(ret_time[2]-ret_time[1])); # Calculate sampling rate

    # 2. Search of RIP position

    total_ion_spectrum <- rowSums(aux) # Sum per rows
    rip_position <- which.max(total_ion_spectrum) # Find maximum for every column
    minima <- as.vector(findpeaks(total_ion_spectrum)[, 2]) # Find local minima
    rip_end_index <- minima[min(which((minima - rip_position) > 0))] # Find ending index of RIP
    rip_start_index <- minima[max(which((rip_position - minima) > 0))] # Find starting index of RIP

    # Compute the 2nd derivative for both axes
    drt <- apply(aux, 1, function(x) -computeDerivative(x, p = 2, n = 21, m = 2, dt = 1/fs))
    ddt <- apply(aux, 2, function(x) -computeDerivative(x, p = 2, n = 11, m = 2, dt = 1/fs))

    daux <- ddt + t(drt)

    patch <- daux[100:200,2000:2100]
    sigmaNoise <- sd(patch)

    # Curve fitting of RIP

    signal <- aux[rip_start_index:rip_end_index, 110]
    template <- -computeDerivative(signal, p = 2, n = tail(c(1:length(signal))[c(T,F)], n=1), m = 2)
    tgauss <- drift_time[rip_start_index:rip_end_index]
    f <- fit_gaussian_density(x = tgauss, y = abs(template)) # Fit RIP into Gaussian
    #gaussianDistr = f.a1*exp(-((tgauss-f.b1)/f.c1).^2) + f.a2*exp(-((tgauss-f.b2)/f.c2).^2) # Fitted Gaussian

    # 5. Peaks and Zero-crossings
    ## 5.a. Retention time

    nNoise <- noise_level # Peaks have to be 3 or 10 times above the noise level (according to IUPAC)
    peaksrt <- vector(mode = "list", length = dim(daux)[2]) # Initialization of vector for peaks
    zeros_rt <- vector(mode = "list", length = dim(daux)[2]) # Initialization of vector for zero crossings

    # For loop that iterates through all the rows
    for(j in (1:dim(daux)[2])){
      # Find the max (peaks)
      #locs <- findpeaks(daux[,j], minpeakheight  = nNoise*sigmaNoise,'WidthReference','halfheight', minpeakdistance  = 4*f.c1*fs)
      if (4*sd(f$y)*fs < 1) {
        locs <- findpeaks(daux[,j], minpeakheight = nNoise*sigmaNoise, minpeakdistance = 2)[ ,2]
      } else {
        locs <- findpeaks(daux[,j], minpeakheight = nNoise*sigmaNoise, minpeakdistance = 4*sd(f$y)*fs)[ ,2]
      }

      # Find the zero-crossing points
      posrt <- findZeroCrossings(daux[,j])
      tmp <- NULL
      locs_tmp <- NULL
      for (k in 1:length(locs)){
        dist <- locs[k] - posrt
        peakaround <- findZeroCrossings(dist)
        if (length(peakaround) >= 1){
          idx1 <- posrt[peakaround]
          idx2 <- posrt[peakaround+1]
          tmp <- cbind(tmp, rbind(idx1, idx2))
          locs_tmp <- rbind(locs_tmp, locs[k])
        }
      }
      zeros_rt[[j]] <- tmp
      peaksrt[[j]] <- locs_tmp
    }

    ## 5.b. Peaks and Zero-Crossing for Drift Time
    peaksdt <- vector(mode = "list", length = dim(daux)[1]) # Initialization of vector for peaks
    zeros_dt <- vector(mode = "list", length = dim(daux)[1]) # Initialization of vector for zero crossings

    # For loop that iterates through all the columns
    for(j in (1:dim(daux)[1])){
      # Find the max (peaks)
      locs <- findpeaks(daux[j, ], minpeakheight = nNoise*sigmaNoise)[ ,2]

      # Find the zero-crossing points
      posdt <- findZeroCrossings(daux[j,])
      tmp <- NULL
      locs_tmp <- NULL
      for (k in 1:length(locs)){
        dist <- locs[k] - posdt
        peakaround <- findZeroCrossings(dist)
        if (length(peakaround) >= 1){
          idx1 <- posdt[peakaround]
          idx2 <- posdt[peakaround+1]
          tmp <- cbind(tmp, rbind(idx1, idx2))
          locs_tmp <- rbind(locs_tmp, locs[k])
        }
      }
      zeros_dt[[j]] <- tmp
      peaksdt[[j]] <- locs_tmp
    }


    # Compute intersection

    peaks <- NULL
    ROIs <- NULL
    for (row in (1:length(peaksrt))){
      columns <- peaksrt[[row]]
      zeros_columns <- zeros_rt[[row]]
      for (col in columns){
        c <- peaksdt[[col]]
        zeros_c <- zeros_dt[[col]]

        if(length(intersect(row,c)) >= 1 & (col > rip_start_index)) {
          peaks <- rbind(peaks, c(row, col))
          minY <- zeros_c[1, c == row]
          maxY <- zeros_c[2, c == row]

          minX <- zeros_columns[columns == col][1]
          maxX <- zeros_columns[columns == col][2]

          width <- abs(maxX[1] - minX[1])
          height <- abs(maxY[1] - minY[1])

          if (minX[1] != maxX[1] & minY[1] != maxY[1] & minX[1] < maxX[1] & minY[1] < maxY[1]){
            ROIs <- rbind(ROIs, c(minX[1], maxX[1], minY[1], maxY[1]))
          }
        }
      }
    }


    # Merging algorithm

    thrOverlap <- 0.2
    aff <- c(1:dim(ROIs)[1])

    done <- NULL
    for (j in (1:dim(ROIs)[1])){
      done <- c(done, j)
      R1 <- ROIs[j, ]
      for (k in c((1:dim(ROIs)[1])[- done])){
        R2 <- ROIs[k, ]
        if (aff[k] != j){
          if (overlapPercentage(R1,R2) > thrOverlap){
            aff[aff == k] <- j
          }
        }
      }
    }


    ROIs_overlap <- NULL
    peaks_overlap <- NULL
    uniqueness <- NULL
    labels <- unique(aff)
    AsF <- NULL
    volume <- NULL
    saturation <- rep(0, length(labels))

    # Search saturation regions
    rip_chrom <- rowSums(aux[rip_start_index: rip_end_index, ]) / length(rip_start_index: rip_end_index)
    max_rip_chrom <- max(rip_chrom)
    saturation_threshold <- 0.1 * max_rip_chrom

    for (n in (1:length(labels))){
      idx <- which(aff == labels[n])
      uniqueness <- c(uniqueness, length(idx))
      R1 <- ROIs[idx[1], ]
      if (length(idx) > 1) {
        for (m in (2:length(idx))){
          R2 <- ROIs[idx[m], ]
          R1 <- c(min(R1[1], R2[1]), max(R1[2], R2[2]), min(R1[3], R2[3]), max(R1[4], R2[4]))
        }
      }

      ROIs_overlap <- rbind(ROIs_overlap, R1)

      patch <- aux[R1[1]:R1[2], R1[3]:R1[4]]
      idx <- which.max(patch)
      r <- ((idx-1) %% nrow(aux)) + 1 # retention = y
      c <- floor((idx-1) / nrow(aux)) + 1 # drift = x
      x <- R1[1] + c
      y <- R1[3] + r
      peaks_overlap <- rbind(peaks_overlap, c(x, y)) # Maximo del ROI

      len_dt <- ROIs_overlap[n, 2] - ROIs_overlap[n, 1]
      len_rt <- ROIs_overlap[n, 4] - ROIs_overlap[n, 3]

      # roi area
      area <- len_rt * len_dt

      # roi volume
      volume <- c(volume, sum(aux[R1[1]:R1[2], R1[3]:R1[4]]))

      # roi center of mass
      x_cm <- round(compute_integral2(aux[R1[1]:R1[2], R1[3]:R1[4]] * (ROIs_overlap[n, 2] - ROIs_overlap[n, 1])) / volume[n])
      y_cm <- round(compute_integral2(aux[R1[1]:R1[2], R1[3]:R1[4]] * (ROIs_overlap[n, 4] - ROIs_overlap[n, 3])) / volume[n])
      dt_mc <- ROIs_overlap[n, 1] + x_cm - 1 #min_dt
      rt_mc <- ROIs_overlap[n, 3] + y_cm - 1 #min_rt


      # roi asymmetries
      half_down_area  <- length(1:y_cm) * len_rt
      half_up_area    <- length(y_cm:len_dt) * len_rt
      half_left_area  <- len_dt * length(1:x_cm)
      half_right_area <- len_dt * length(x_cm:len_rt)
      asymetry <- round(((half_down_area - half_up_area) / (area)), 2)
      AsF <- c(AsF, asymetry)

      saturation_regions <- which(rip_chrom <= saturation_threshold)
      if (length(saturation_regions) == 0){
        saturation_minima <- NULL
      } else {
        saturation_list <- split(saturation_regions, cumsum(c(1, diff(saturation_regions)) != 1))
        saturation_minima <- matrix(0, length(saturation_list), 2)
        for (k in 1:length(saturation_list)){
          saturation_minima[k, ] <- c(min(saturation_list[[k]]), max(saturation_list[[k]]))
        }
      }


      # roi saturation
      if (length(saturation_minima) == 0){
        # No saturation. Do nothing
      } else {
        for (l in (1:dim(saturation_minima)[1])){
          if ((saturation_minima[l, 1] < rt_mc)
              & (saturation_minima[l, 2] > rt_mc)) {
            saturation[n] <- 1
            break
          }
        }
      }
    }

    colnames(ROIs_overlap) <- c("minRT", "maxRT", "minDT", "maxDT")

    aux_list$data$data_df <- round(aux)
    aux_list$data$ROIs <- ROIs_overlap
    aux_list$data$Peaks <- peaks_overlap
    aux_list$data$Parameters <- rbind(AsF, saturation, uniqueness, volume)
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}


#---------------#
#   FUNCTIONS   #
#---------------#

#--------------------------#
#   fit_gaussian_density  #
#------------------------#

# https://stats.stackexchange.com/a/83029/62083
fit_gaussian_density <- function(x, y) {
  # Guess for  mu, k and sigma
  mu0 <- which.max(y)
  k0 <- y[mu0]
  x_in_peak <- x[(y - k0/2) > 0]
  if (length(x_in_peak) == 0)  {
    fwhm <- stats::IQR(x)
  } else {
    fwhm <- max(x_in_peak) - min(x_in_peak)
    if (fwhm == 0) {
      fwhm <- stats::IQR(x)
    }
  }
  sigma0 <- fwhm/2.35482004503
  # Fit with nls:
  fit <- stats::nls(
    y ~ k*exp(-1/2*(x-mu)^2/sigma^2),
    start=c(mu=mu0,sigma=sigma0,k=k0),
    data = data.frame(x = x, y = y)
  )
  # Get coefficients and fitted values:
  out <- as.list(stats::coef(fit)) # list(mu = ***, sigma = ***, k = ***)
  out$x <- x
  out$y <- stats::predict(fit)
  # Return the numeric vector
  out
}

#----------------------#
#   computeDerivative  #
#----------------------#


computeDerivative <- function(x, p, n, m, dt){
  dx = sgolayfilt(x, p, n, m) # From package: signal
  return(dx)
}

#----------------------#
#   findZeroCrossings  #
#----------------------#

findZeroCrossings <- function(x){
  signs <- sign(x)
  pos_plus <- which(diff(signs) >= 2)
  pos_minus <- which(diff(signs) == -2)
  pos <- sort(as.vector(union(pos_plus + 1, pos_minus)))
  return(pos)
}

#----------------------#
#   overlapPercentage  #
#----------------------#

overlapPercentage <- function(ROI1, ROI2){
  x_left <- max(ROI1[3], ROI2[3])
  y_top <- max(ROI1[1], ROI2[1])
  x_right <- min(ROI1[4], ROI2[4])
  y_bottom <- min(ROI1[2], ROI2[2])

  if (!(x_right <= x_left | y_bottom <= y_top)){
    area1 <- (ROI1[2] - ROI1[1])*(ROI1[4] - ROI1[3])
    area2 <- (ROI2[2] - ROI2[1])*(ROI2[4] - ROI2[3])
    overlapping_area <- (x_right - x_left)*(y_bottom - y_top)
    p <- overlapping_area / (area1 + area2 - overlapping_area)
  } else {
    p <- 0
  }
  return(p)
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
  cx <- gaussLegendre(n, xa, xb)
  x <- cx$x
  wx <- cx$w
  cy <- gaussLegendre(m, ya, yb)
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


gcims_view_ROIs <- function(dir_in, sample_num, rt_range = NULL, dt_range = NULL){

  Retention_Time <- Drift_Time <- Value <- NULL

  print(" ")
  print("  ///////////////////////////////////")
  print(" /    Sample ROIs Visualization    /")
  print("///////////////////////////////////")
  print(" ")

  setwd(dir_in)
  print(paste0("Visualizing sample ", sample_num))
  aux_string <- paste0("M", sample_num, ".rds")
  aux_list <- readRDS(aux_string) #new
  aux <- (as.matrix(aux_list$data$data_df)) #new
  ROIs <- aux_list$data$ROIs
  colnames(ROIs) <- c("dt1", "dt2", "rt1", "rt2")

  #SOME CHECKS
  retention_time <- aux_list$data$retention_time
  drift_time <- aux_list$data$drift_time
  cond_1_rt <- (rt_range[1] - retention_time[1]) < 0
  cond_2_rt <- (rt_range[2] - retention_time[length(retention_time)]) > 0
  cond_1_dt <-(dt_range[1] - drift_time[1]) < 0
  cond_2_dt <-(dt_range[2] - drift_time[length(drift_time)]) > 0


  if(is.null(rt_range)){# old
    rt_ind <- c(1, dim(aux)[2]) #New

  } else{
    if(cond_1_rt | cond_2_rt){
      stop("Retention time range out of bounds.")
    }
    rt_ind  <- c(which.min(abs(retention_time - rt_range[1])), which.min(abs(retention_time - rt_range[2])))
    if( rt_ind[1] == rt_ind[2]){
      stop("Initial and Final retention time values can't be equal in the variable rt_range.")
    }##New
  }



  if(is.null(dt_range)){# old
    dt_ind <- c(1, dim(aux)[1]) #New
  } else{
    if(cond_1_dt | cond_2_dt){
      stop("Drift time range out of bounds.")
    }
    dt_ind  <- c(which.min(abs(drift_time - dt_range[1])), which.min(abs(drift_time - dt_range[2]))) #New
    if( dt_ind[1] == dt_ind[2]){
      stop("Initial and Final drift time values can't be equal in the variable dt_range.")
    }#
  }

  sel_index_rt <- rt_ind[1]: rt_ind[2]
  sel_index_dt <- dt_ind[1]: dt_ind[2]

  if(is.null(rt_range)){

  } else if((class(sel_index_rt) == "integer") & (sel_index_rt[2] > sel_index_rt[1])){
  } else {
    stop("Possible errors: 1) The selected vector of indexes corresponding to the provided retention time range is not an integer vector, 2) or rt_range[2] <= rt_range[1])")
  }

  if(is.null(dt_range)){

  } else if((class(sel_index_dt) == "integer") & (sel_index_dt[2] > sel_index_dt[1])){
  } else {
    stop("Possible errors: 1) The selected vector of indexes corresponding to the provided drift time range is not an integer vector, 2) or dt_range[2] <= dt_range[1])")
  }

  retention_time <- retention_time[sel_index_rt]
  drift_time <- drift_time[sel_index_dt]

  aux <- aux[sel_index_dt, sel_index_rt]#old
  rownames(aux) <- round(drift_time, digits = 2) #old
  colnames(aux) <- retention_time

  moltaux <- melt((aux))
  colnames(moltaux) <- c("Drift_Time", "Retention_Time", "Value")

  ROIsTime <- matrix(NA, nrow = dim(ROIs)[1], ncol = 4)
  moltaux$ROIsmoltxmin <- rep(NA, dim(moltaux)[1])
  moltaux$ROIsmoltxmax <- rep(NA, dim(moltaux)[1])
  moltaux$ROIsmoltymin <- rep(NA, dim(moltaux)[1])
  moltaux$ROIsmoltymax <- rep(NA, dim(moltaux)[1])

  for (j in (1:dim(ROIsTime)[1])) {
    dt <- drift_time[2] - drift_time [1]
    rt <- retention_time[2] - retention_time[1]
    ROIsTime[,1:2] <- round(ROIs[,1:2] * rt, digits = 2)
    ROIsTime[,3:4] <- round((ROIs[,3:4] * dt) + 6, digits = 2)

    moltaux[which(moltaux[,2] == ROIsTime[j, 1]), "ROIsmoltxmin"] <- ROIsTime[j, 1]
    moltaux[which(moltaux[,2] == ROIsTime[j, 2]), "ROIsmoltxmax"] <- ROIsTime[j, 2]
    moltaux[which(moltaux[,1] == ROIsTime[j, 3]), "ROIsmoltymin"] <- ROIsTime[j, 3]
    moltaux[which(moltaux[,1] == ROIsTime[j, 4]), "ROIsmoltymax"] <- ROIsTime[j, 4]
  }

  #We do this in order to plot the data using geom_raster that is faster than geom_tile
  #perhaps a previous interpolation is needed to avoid this patch:
  rep_dt_index <- rep(seq(from = 1, to = dim(aux)[1], by = 1), times = dim(aux)[2])
  # drift_time_period <- mean(diff(drift_time))
  # corr_drift_time <- seq(from = drift_time[1], by = drift_time_period, length.out = length(drift_time))
  # moltaux$Drift_Time <- corr_drift_time[rep_dt_index]
  moltaux$Drift_Time <- drift_time[rep_dt_index]
  tt <- na.omit(moltaux)
  tt <- tt[,-c(1:3)]

  rm(aux, aux_string)

  p <- ggplot(moltaux, aes(x = Drift_Time, y = Retention_Time, fill = Value)) +
    geom_raster(interpolate = FALSE) +
    scale_fill_viridis(discrete = FALSE, option = "A", direction = -1) +
    theme_minimal()
  # labs(x="Drift Time (ms)",
  #      y="Retention Time (s)",
  #      title = "Sample Matrix Image",
  #      fill = "Intensity") +

  p +
    geom_rect(data=d, mapping=aes(xmin=dt1, xmax=dt2, ymin=rt1, ymax=rt2), color="black", alpha=0.5)

  print(p)
}
