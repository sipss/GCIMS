
#' Two dimensional peak picking

#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         Samples to which perform the peak picking.
#' @param rem_baseline    Perform a baseline removal in retention / drift time
#'                        axes. TRUE by default.
#' @param min_length_tr   Minimum peak length in retention time (now indexes).
#' @param min_length_td   Minimum peak length in drift time (now indexes).
#' @param cor_threshold   Correlation threshold for selecting ROIs.
#' @return A peak table per sample.
#' @family Peak Peaking functions
#' @export
#' @importFrom stats sd
#' @importFrom pracma meshgrid
#' @importFrom chemometrics sd_trim
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }



gcims_peak_picking <- function(dir_in, dir_out, samples, rem_baseline = TRUE,
                               min_length_tr = 10, min_length_td = 10, cor_threshold = 0.6) {


  # Set of functions needed to run gcims_peak_picking:

  #--------------#
  # find_peaks2d #
  #--------------#

  find_peaks2d <- function(z, min_length_tr, min_length_td, cor_threshold){
    gauss2d <- function(x, sigma1, y,  sigma2){
      my_grid <- meshgrid(x, y)
      gaussian <-  exp(-(1/2) * ((my_grid$X) / sigma1) ^ 2) * exp(-(1/2) * ((my_grid$Y) / sigma2) ^ 2)
    }
    x <- seq(from = -round(min_length_tr / 2), to = round(min_length_tr / 2), by = 1)
    y <- seq(from = -round(min_length_td / 2), to = round(min_length_td / 2), by = 1)
    pattern <- gauss2d(x, sigma1 = (min_length_tr / 4), y, sigma2 = (min_length_td / 4))

    N <- dim(z)[1]
    M <- dim(z)[2]
    R <- dim(pattern)[1]
    S <- dim(pattern)[2]

    cross_corr <- matrix(0, N, M)
    for (i in 1: (N-R+1)){
      for(j in 1: (M-S+1)){
        segment <- z[i: (i + R - 1), j: (j + S - 1)]
        cross_corr[i + round(R/2) - 1, j + round(S/2) - 1] = (sum((pattern-mean(pattern))*(segment - mean(segment)))/((R-1)*(S-1)))/(sd(pattern)*sd(segment))
      }
    }
    cross_corr[cross_corr < cor_threshold] <- 0
    return(cross_corr)
  }

  #----------------#
  # find.contigous #
  #----------------#

  find.contiguous <- function(img, x, bg) {
    ## we need to deal with a single (row,col) matrix index
    ## versus a collection of them in a two column matrix separately.
    if (length(x) > 2) {
      lbl <- img[x][1]
      img[x] <- bg
      xc <- x[,1]
      yc <- x[,2]
    } else {
      lbl <- img[x[1],x[2]]
      img[x[1],x[2]] <- bg
      xc <- x[1]
      yc <- x[2]
    }
    ## find all neighbors of x
    xmin <- ifelse((xc-1) < 1, 1, (xc-1))
    xmax <- ifelse((xc+1) > nrow(img), nrow(img), (xc+1))
    ymin <- ifelse((yc-1) < 1, 1, (yc-1))
    ymax <- ifelse((yc+1) > ncol(img), ncol(img), (yc+1))
    ## find all neighbors of x
    x <- rbind(cbind(xmin, ymin),
               cbind(xc  , ymin),
               cbind(xmax, ymin),
               cbind(xmin, yc),
               cbind(xmax, yc),
               cbind(xmin, ymax),
               cbind(xc  , ymax),
               cbind(xmax, ymax))
    ## that have the same label as the original x
    x <- x[img[x] == lbl,]
    ## if there is none, we stop and return the updated image
    if (length(x)==0) return(img);
    ## otherwise, we call this function recursively
    find.contiguous(img,x,bg)
  }

  #-------------------------#
  # extract_contiguous_ones #
  #-------------------------#

  extract_contiguous_ones <- function(img){
    ## as matrices of 1's are extracted by the process
    img <- corr_2d
    ## get all pixel coordinates that are objects
    x <- which(img==1, arr.ind=TRUE)
    ## loop until there are no more pixels that are objects
    ##the output is in the list out
    count <- 0
    out <- list()
    while (length(x) > 0) {
      ## choose a single (e.g., first) pixel location. This belongs to the current
      ## object that we will grow and remove from the image using find.contiguous
      if (length(x) > 2) {
        x1 <- x[1,]
      }
      ## make the call to remove the object from img
      img <- find.contiguous(img, x1, 0)
      ## find the remaining pixel locations belonging to objects
      xnew <- which(img==1, arr.ind=TRUE)
      count <- count + 1
      ## extract the indices for the 1's found by diffing new with x
      out.ind <- x[!(x[,1] %in% xnew[,1] & x[,2] %in% xnew[,2]),]
      ## set it as a matrix in the output
      out[[count]]  <- list(X = min(out.ind[, 1]):max(out.ind[, 1]), Y = min(out.ind[, 2]):max(out.ind[, 2]))
      x <- xnew
    }
    return(out)

  }

  #---------#
  # ind2sub #
  #---------#

  ind2sub <- function(n,ind){
    r <- ((ind-1) %% n) + 1
    c <- floor((ind-1) / n) + 1
    sub <- list(X = r, Y = c)
    return(sub)
  }

  #--------------------#
  # estimate_threshold #
  #--------------------#

  estimate_threshold <- function(x){
    # Compute robust estimation
    # of noise present in a sample
    mu = mean(x, trim = 0.2)
    sigma = sd_trim(x, const = TRUE)
    threshold <- mu + (3 * sigma)
    return(threshold)
  }


  #---------------#
  # compute_power #
  #---------------#

  compute_power <- function (x){
    pow <- sum(x * x) / length(x)
    pow
  }

  #---------------#
  # is_overlapped #
  #---------------#

  is_overlapped <- function(l1x, l1y, r1x, r1y, l2x, l2y,  r2x, r2y){

    # If one rectangle is on left side of other
    if((l1x >= r2x) | (l2x >= r1x)){
      return(FALSE)
    }

    # If one rectangle is above other
    if((l1y <= r2y) | (l2y <= r1y)){
      return(FALSE)
    }
    return(TRUE)
  }

  #------------------#
  # merge_rectangles #
  #------------------#

  merge_rectangles <-function(x, subset){

    # Set of rectangles' diagonal vertex coordinates L,R

    lx <- x[subset,]$roi_lx
    ly <- x[subset,]$roi_ly
    rx <- x[subset,]$roi_rx
    ry <- x[subset,]$roi_ry

    # Diagonal vertex coordinates of the merged rectangle
    # Replace them in the original data

    x[subset,]$roi_lx <- min(lx)
    x[subset,]$roi_ly <- max(ly)
    x[subset,]$roi_rx <- max(rx)
    x[subset,]$roi_ry <- min(ry)

    return(x)
  }



setwd(dir_in)

#m <- 0
#for (i in samples){
# m <- m + 1

m <-1
print(paste0("Performing Peak Picking in sample ", samples[m]))

# Read current data

aux_string <- paste0("M", samples[m], ".rds")
aux_list <- readRDS(aux_string) #new
aux <- t(as.matrix(aux_list$data$data_df))

# Remove baseline in drift time

if (rem_baseline == TRUE){
  aux <- psalsa(data = aux, lambda = 1E7, p = 0.001, k = -1, maxit = 25)$corrected
  # Remove baseline in retention time
  aux <- t(psalsa(data = t(aux), lambda = 1E7, p = 0.001, k = -1, maxit = 25)$corrected)
}

# Find Peaks in 2D

corr_2d <- find_peaks2d(aux, min_length_tr, min_length_td, cor_threshold)

# Select regions with the highest
# correlations and set them to one

corr_2d[corr_2d > 0] <- 1

# Create submatrices within a list
# with the contiguous ones (rois)

out <- extract_contiguous_ones(img = corr_2d) #aichao overflow

# Unfold aux

aux_vector <- sort(as.vector(aux), decreasing = TRUE) #new

# Compute robust estimation
# of noise present in a sample

threshold <- estimate_threshold(aux_vector)

# Compute noise power

noise_power <- compute_power(aux_vector[aux_vector <= threshold])
rm(aux_vector)

#########################
# PEAK TABLE GENERATION #
#########################

# Set to zero matrix values that are below the noise threshold.
# This is done to better locate the limits of the roi
aux[aux <= threshold] <- 0

# Compute a list of tentative peaks.
peak_position_list <- vector(mode = "list", length = length(out))

# Loop for searching peak each peak position
# Data is stored in the intermediate variable peak_position_list

for (k in 1: length(out)){
  max_region <- as.matrix(aux[out[[k]]$X, out[[k]]$Y])
  max_ind <- which.max(max_region)
  relative_peak_position_list <- ind2sub(n = dim(max_region)[1], ind = max_ind)
  peak_position_list[[k]] <- list(X = out[[k]]$X[relative_peak_position_list$X],
                                  Y = out[[k]]$Y[relative_peak_position_list$Y])
  }
rm(relative_peak_position_list, max_region, max_ind)

# Generate an empty peak table
# with 16 colnames (sample, roi and
# peak parameters and figures of merit )

peak_table <- matrix(0, ncol = 16 , nrow = length(out))
colnames(peak_table) <- c("sample_id", "roi_id",
                          "roi_ly", "roi_lx",
                          "roi_ry", "roi_rx",
                          "peak_id",
                          "peak_max", "peak_rt_ind","peak_dt_ind",
                          "peak_sat", "peak_snrdB",
                          "peak_rt_length","peak_dt_length",
                           "peak_rt_assym", "peak_dt_assym")

# Check if the tentative peak comes from the sequence :
# (- -> +) and  (+ -> -) in terms of slope
# If TRUE it's a peak
# else do not selec as a peak
# Apart from that create a first version of the peak_table:

for (k in 1: length(out)){

  # Is this a peak?
  zeros_Y <- which(aux[, peak_position_list[[k]]$Y] == 0)
  left_Y_zeros <- zeros_Y[zeros_Y < peak_position_list[[k]]$X]
  right_Y_zeros <- zeros_Y[zeros_Y > peak_position_list[[k]]$X]

  zeros_X <- which(aux[peak_position_list[[k]]$X, ] == 0)
  left_X_zeros <- zeros_X[zeros_X < peak_position_list[[k]]$Y]
  right_X_zeros <- zeros_X[zeros_X > peak_position_list[[k]]$Y]

  #If there is not any zero the left or right sides of the maximum this is not a peak
  if((length(left_Y_zeros) == 0) | (length(right_Y_zeros) == 0) | (length(left_X_zeros) == 0) | (length(right_X_zeros) == 0)){
    # do nothing
  } else {
    # retention time params
    left_Y_border <- max(left_Y_zeros)
    right_Y_border <- min(right_Y_zeros)

    # drift time params
    left_X_border <- max(left_X_zeros)
    right_X_border <- min(right_X_zeros)

    # rectangle params
    lx <-left_X_border
    ly <- right_Y_border
    rx <- right_X_border
    ry <- left_Y_border

    # create peak table
    peak_table[k, 1] <- samples[m]
    peak_table[k, 3] <- ly
    peak_table[k, 4] <- lx
    peak_table[k, 5] <- ry
    peak_table[k, 6] <- rx
    peak_table[k, 8] <- aux[peak_position_list[[k]]$X, peak_position_list[[k]]$Y]
    peak_table[k, 9] <- peak_position_list[[k]]$X
    peak_table[k, 10] <- peak_position_list[[k]]$Y

  }
}

# Remove false peaks from the table
peak_table <- peak_table[rowSums(peak_table) != 0, ]
peak_table <- as.data.frame(peak_table)

# Set roi and peak ID (defaults)
peak_table$roi_id <- 1:dim(peak_table)[1]
peak_table$peak_id <- 1:dim(peak_table)[1]


#Create and empty list for all possible peaks
overlapped_peaks_list <- vector(mode = "list", length = dim(peak_table)[1])

#check which rectangles are overlapped and group them
peak_id <- 1:dim(peak_table)[1]

remaining_peaks <- peak_id
k <- 0
while(length(remaining_peaks) > 0){
  k <- k + 1
  current_peak_id <- remaining_peaks[1]
  test_result <- as.vector(matrix(FALSE, nrow = 1, ncol = length(remaining_peaks)))
  h <- 0
  for (tested_peak_id in remaining_peaks){
    h <- h + 1
    #params for the rectangle of the current peak
    l1x <- peak_table[current_peak_id,]$roi_lx
    l1y <- peak_table[current_peak_id,]$roi_ly
    r1x <- peak_table[current_peak_id,]$roi_rx
    r1y <- peak_table[current_peak_id,]$roi_ry

    #params for the rectangle of the peak that has to be
    #tested against the current peak
    l2x <- peak_table[tested_peak_id,]$roi_lx
    l2y <- peak_table[tested_peak_id,]$roi_ly
    r2x <- peak_table[tested_peak_id,]$roi_rx
    r2y <- peak_table[tested_peak_id,]$roi_ry

    test_result[h] <- is_overlapped(l1x, l1y, r1x, r1y, l2x, l2y,  r2x, r2y)
  }

  overlapped_peaks_list[[k]] <- remaining_peaks[test_result]
  remaining_peaks  <- remaining_peaks [!test_result]
}

overlapped_peaks_list <-overlapped_peaks_list[lengths(overlapped_peaks_list) != 0]

# Merge the interesting regions and
# upgrade roi coordinates and roi_id

l <- 0
for (h in (1:length(overlapped_peaks_list))){
  l <- l + 1
  peak_table <- merge_rectangles(peak_table, overlapped_peaks_list[[h]])
  peak_table[overlapped_peaks_list[[h]], ]$roi_id <- l

}

# Compute figures of merit for the isolated peaks,
# otherwise set them to NA

for (k in unique(peak_table$roi_id)){
   peaks_per_roi <- length(which(peak_table$roi_id == k))
   if (peaks_per_roi == 1){

      # retention time params
       peak_table[k, ]$peak_rt_length <- peak_table[k, ]$roi_ly - peak_table[k, ]$roi_ry
       srx <-  peak_table[k, ]$peak_rt_ind - peak_table[k, ]$roi_ry
       slx <-  peak_table[k, ]$roi_ly - peak_table[k, ]$peak_rt_ind
       peak_table[k, ]$peak_rt_assym <- srx / slx

      # drift time params
       peak_table[k, ]$peak_dt_length <- peak_table[k, ]$roi_rx - peak_table[k, ]$roi_lx
       sry <-   peak_table[k, ]$roi_rx - peak_table[k, ]$peak_dt_ind
       sly <-  peak_table[k, ]$peak_dt_ind - peak_table[k, ]$roi_lx
       peak_table[k, ]$peak_dt_assym <- sry / sly

      # signal to noise ration in dB
       peak_table[k, ]$peak_snrdB <- 10 * log10((compute_power(aux[peak_table[k, ]$roi_ly: peak_table[k, ]$roi_ry,
                                peak_table[k, ]$roi_rx: peak_table[k, ]$roi_lx]) - noise_power)/ noise_power)

   } else {
     peak_table[peak_table$roi_id == k, 12:16] <- NA

   }
   # saturation
   peak_table[k, ]$peak_sat <- as.logical(peak_table[k, ]$peak_max >= 0.9 * (max(aux)))
}
return(peak_table)
# image(1:200, 1:60, t(aux), xlab = "Drift time index", ylab ="Retention time index", main = "Peak Picking Example")
# points(peak_table$dt_ind, peak_table$rt_ind, pch = 4, lw = 2, col= "blue")
# points(peak_table$roi_lx, peak_table$roi_ly, pch = 4, lw = 2, col= "red")
# points(peak_table$roi_lx, peak_table$roi_ry, pch = 4, lw = 2, col= "red")
# points(peak_table$roi_rx, peak_table$roi_ry,pch = 4, lw = 2, col= "red")
# points(peak_table$roi_rx, peak_table$roi_ly,pch = 4, lw = 2, col= "red")



# setwd(dir_out)
# saveRDS(aux_list, file = paste0("M", i, ".rds"))
# setwd(dir_in)

}


