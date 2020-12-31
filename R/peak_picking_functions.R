
#' Two dimensional peak picking

#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         Samples to which perform the peak picking.
#' @param rem_baseline    Perform a baseline removal in retention / drift time
#'                        axes. TRUE by default.
#' @param min_length_tr   Minimum peak length in retention time (now indexes).
#' @param min_length_td   Minimum peak length in drift time (now indexes).
#' @param preprocess      Preprocess the samples?. TRUE or FALSE
#' @return A peak table per sample.
#' @family Peak Peaking functions
#' @export
#' @importFrom stats sd
#' @importFrom dbscan dbscan kNNdist
#' @importFrom pracma meshgrid findpeaks gaussLegendre
#' @importFrom purrr transpose
#'
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
                               min_length_tr = 50, min_length_td = 10, preprocess = TRUE) {


  #---------------#
  #   FUNCTIONS   #
  #---------------#

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
        condition <- sum(segment)
        if (condition == 0){
          cross_corr[i + round(R/2) - 1, j + round(S/2) - 1] <- 0
        } else{
        cross_corr[i + round(R/2) - 1, j + round(S/2) - 1] = (sum((pattern-mean(pattern))*(segment - mean(segment)))/((R-1)*(S-1)))/(sd(pattern)*sd(segment))
        }
      }
    }
    cross_corr[cross_corr < cor_threshold] <- 0
    return(cross_corr)
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


  #---------#
  # sub2ind #
  #---------#

  sub2ind <- function(sub, n){
    r <- sub[1]
    c <- sub[2]
    ind = ((c-1) * n) + r
    return(ind)
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
  #    cmenger    #
  #---------------#

  cmenger <- function(x1,y1,x2,y2,x3,y3){
    sqrt(abs(4*((x1-x2)^2+(y1-y2)^2)*((x2-x3)^2+(y2-y3)^2)-
               ((x1-x2)^2+(y1-y2)^2+(x2-x3)^2+(y2-y3)^2-(-x3+x1)^2-
                  (-y3+y1)^2)^2))/(sqrt((x1-x2)^2+(y1-y2)^2)*
                                     sqrt((x2-x3)^2+(y2-y3)^2)*sqrt((x3-x1)^2+(y3-y1)^2));
  }

  #---------------#
  #    findknee   #
  #---------------#

  findknee <- function(x,y){
    n <- length(x)
    c <- rep(0,n-2)
    for (j in (1:(n-2)))
    {
      c[j] <- cmenger(x[j],y[j],x[j+1],y[j+1],x[j+2],y[j+2])
    }
    jm  <-which.max(c);cm=max(c)
    ikn <- jm+1
    xkn <-x[ikn]
    c(cm, ikn, x[ikn])
  }



  #-------------------------#
  #   compute_inner_border  #
  #-------------------------#

  compute_inner_border <- function(x, indexes){
    border <- c(max(x[indexes, 2]), max(x[indexes, 1]), min(x[indexes, 2]), min(x[indexes, 1]))
    return(border)

  }


  #-------------------------#
  #   compute_outer_border  #
  #-------------------------#

  compute_outer_border <- function(x, rectangle){
    all_col <- 1:dim(x)[2]
    all_row <- 1:dim(x)[1]
    left <- c(rectangle[2], round((rectangle[3] + rectangle[5])/2))
    bottom <- c(round((rectangle[2] + rectangle[4])/2), rectangle[3])
    right <- c(rectangle[4], round((rectangle[3] + rectangle[5])/2))
    top <- c(round((rectangle[2] + rectangle[4])/2), rectangle[5])

    if(any(x[left[2], (all_col > left[1])] == 0)){
      valleys <- findpeaks(- x[left[2], (all_col > left[1])])[, 2] + left[1]
      zero <- min(which(x[left[2], (all_col > left[1])] == 0)) + left[1]
      if(length(valleys) == 0){
        a <- zero
      } else {
        a <- min(c(zero, valleys[1]))
      }
    } else {
      a <- all_col[length(all_col)]
    }
    if(any(x[(all_row > bottom[2]), bottom[1]] == 0)){
      valleys <- findpeaks(- x[(all_row > bottom[2]), bottom[1]])[, 2] + bottom[2]
      zero <- min(which(x[(all_row > bottom[2]), bottom[1]] == 0)) + bottom[2]
      if(length(valleys) == 0){
        b <- zero
      } else {
        b <- min(c(zero, valleys[1]))
      }
    } else {
      b <- all_row[length(all_row)]
    }

    if(any(x[right[2], (all_col < right[1])] == 0)){
      valleys <- findpeaks(- x[right[2], (all_col < right[1])])[, 2]
      zero <- max(which(x[right[2], (all_col < right[1])] == 0))
      if(length(valleys) == 0){
        c <- zero
      } else {
        c <- max(c(zero, valleys[length(valleys)]))
      }
    } else {
      c <- 1
    }

    if(any(x[(all_row < top[2]), top[1]] == 0)){
      valleys <- findpeaks(-x[(all_row < top[2]), top[1]])[, 2]
      zero <- max(which(x[(all_row < top[2]), top[1]] == 0))
      if(length(valleys) == 0){
        d <- zero
      } else {
        d <- max(c(zero, valleys[length(valleys)]))
      }
    }
    else {
      d <- 1
    }
    border <- c(rectangle[1], a, b, c, d)
    return(border)
  }

  #----------------#
  #  find_indexes  #
  #----------------#

  find_indexes <- function(z){
    indy <- z[4]:z[2]
    indx <- z[5]:z[3]
    my_grid <- meshgrid(indx , indy )
    ind_list <- list(X = my_grid$X, Y = my_grid$Y)
    return(ind_list)
  }


  #-------------------------#
  #     compute_integral2   #
  #-------------------------#

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

 #-----------------------------------------------------------#

 #---------------#
 #   CONSTANSTS  #
 #---------------#

  # digital smothing:

  filter_length <- 19
  polynomial_order <- 2


  # baseline removal:

  lambda <- 1E7
  p <-  0.001
  k <- -1

  # peak picking:

  cor_threshold <- 0.5

  # peak clustering:

  neighbors <- 4


  #-----------------------------------------------------------#


  #-------------#
  #     MAIN    #
  #-------------#

  # 1)   read data
  # 2)   a. search RIP position
  #      b. search saturation regions
  # 3)   remove baseline
  # 4)   a. compute intensity threshold
  #      b. compute noise power
  # 5)   digital smoothing
  # 6)   a. remove data below threshold
  #      b. remove data before the RIP
  # 7)   find peaks in 2D (ROIs if convoluted)
  # 8)   cluster data in ROIs
  # 9)   remove data out of the ROIs
  # 10)  generate ROI table
  # 11)  save results

  setwd(dir_in)
  m <- 0
  for (i in samples){
    m <- m + 1
    print(paste0("Performing Peak Picking in sample ", samples[m]))

    # 1)   read data
    aux_string <- paste0("M", samples[m], ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- t(as.matrix(aux_list$data$data_df))

    # 2)   a. search RIP position
    if (preprocess  == TRUE){
      total_ion_spectrum <- colSums(aux)
      rip_position <- which.max(total_ion_spectrum)
      minima <- as.vector(findpeaks(-total_ion_spectrum)[, 2])
      rip_end_index <- minima[min(which((minima - rip_position) > 0))]
      rip_start_index <- minima[max(which((rip_position - minima) > 0))]
    }
    # 2)   b. search saturation regions
    rip_chrom <- rowSums(aux[, rip_start_index: rip_end_index]) / length(rip_start_index: rip_end_index)
    max_rip_chrom <- max(rip_chrom)
    saturation_threshold <- 0.1 * max_rip_chrom
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

    # 3)   remove baseline
    if (preprocess  == TRUE){
      aux <- psalsa(data = aux, lambda = lambda, p = p, k = k, maxit = 25)$corrected
      # Remove baseline in retention time
      aux <- t(psalsa(data = t(aux), lambda = lambda, p = p, k = k, maxit = 25)$corrected)
    }
    aux_vector <- sort(as.vector(aux), decreasing = TRUE) #new

    # 4)   a. compute intensity threshold
    threshold <- estimate_threshold(aux_vector)
    # 4)  b. compute noise power
    noise_power <- compute_power(aux_vector[aux_vector <= threshold])
    rm(aux_vector)

    # 5)   digital smoothing
    if (preprocess == TRUE){
      aux <- t(aux)
      n <- dim(aux)[1]
      for (j in (1:n)){
        aux[j, ] <- sgolayfilt(aux[j, ], p = polynomial_order, n = filter_length)
      }
      aux <- t(aux)
      n <- dim(aux)[1]
      for (j in (1:n)){
        aux[j, ] <- sgolayfilt(aux[j, ], p = polynomial_order, n = filter_length)
      }
    }

    # 6)   a. remove data below threshold
    aux[aux <= threshold] <- 0
    # 6)   b. remove data before the RIP
    if (preprocess  == TRUE){
    aux[, 1:rip_end_index] <- 0
    }

    # 7)   find peaks in 2D (ROIs if convoluted)
    corr_2d <- find_peaks2d(aux, min_length_tr, min_length_td, cor_threshold)
    coord_linear <- which(corr_2d > 0)
    coord_mat <- matrix(0, nrow = length(coord_linear), ncol = 2)
    coord_row_max <- dim(coord_mat)[1]
    aux_row_max <- dim(aux)[1]
    for (t in (1:coord_row_max)){
      coord_mat[t, ]  <- unlist(ind2sub(aux_row_max, coord_linear[t]))
    }
    rm(corr_2d)

    # 8)   cluster data in ROIs
    sorted_distances <-sort(kNNdist(coord_mat, k = neighbors))
    eps <- sorted_distances[findknee(1:length(sorted_distances), sorted_distances)[3]]
    out <- dbscan(coord_mat, eps, neighbors)$cluster
    clust_ind <- unique(out)
    inner_rectangles <- matrix(0, nrow = length(clust_ind), ncol = 4)
    h <- 0
    for (k in clust_ind){
      h <- h + 1
      indexes <- which(out == k)
      inner_rectangles[h, ] <- compute_inner_border(coord_mat[indexes,])
    }
    inner_rectangles <- cbind(clust_ind, inner_rectangles)
    inner_rectangles <- inner_rectangles[inner_rectangles[,1] != 0, ]
    outer_rectangles <- matrix(0, nrow = dim(inner_rectangles)[1], ncol = dim(inner_rectangles)[2])
    h <- 0
    for (k in inner_rectangles[, 1]){
      h <- h + 1
      indexes <- which(inner_rectangles[, 1] == k)
      outer_rectangles[h, ] <- compute_outer_border(aux, inner_rectangles[indexes,])
    }
    roi_coord_list <- transpose(lapply(split(outer_rectangles, row(outer_rectangles)), find_indexes))
    roi_coord_sub <- cbind(unlist(roi_coord_list$X), unlist(roi_coord_list$Y))
    n <- dim(aux)[1]
    roi_coord_ind <- apply(roi_coord_sub, sub2ind, n, MARGIN = 1)

    # 9)   remove data out of the ROIs
    aux[!roi_coord_ind] <- 0

    # 10)  generate ROI table
    roi_table <- matrix(0, ncol = 16 , nrow = dim(outer_rectangles)[1])
    colnames(roi_table) <- c("sample_id", "roi_id",
                             "max_rt", "min_dt",
                             "min_rt", "max_dt",
                             "rt_length","dt_length",
                             "area", "volume",
                             "rt_mc","dt_mc",
                             "rt_asym", "dt_asym",
                             "saturation", "snr")

    for (k in outer_rectangles[, 1]){
      roi_table[k, 1] <- samples[m]
      roi_table[k, 2] <- k
      roi_table[k, 3] <- outer_rectangles[k, 3]
      roi_table[k, 4] <- outer_rectangles[k, 4]
      roi_table[k, 5] <- outer_rectangles[k, 5]
      roi_table[k, 6] <-outer_rectangles[k, 2]
      # current roi
      roi <- aux[outer_rectangles[k, 5]:outer_rectangles[k, 3],
                 outer_rectangles[k, 4]:outer_rectangles[k, 2]]
      # roi lengths
      len_x <-  outer_rectangles[k, 3] - outer_rectangles[k, 5]
      len_y <- outer_rectangles[k, 2] - outer_rectangles[k, 4]
      roi_table[k, 7] <- len_x
      roi_table[k, 8] <- len_y
      # roi area
      area <- len_x * len_y
      roi_table[k, 9] <- area
      # roi volume
      volume <- compute_integral2(roi)
      roi_table[k, 10] <- volume
      # roi center of mass
      x_cm <- round(compute_integral2(roi * (1:dim(roi)[1])) / volume)
      y_cm <- round(compute_integral2(t(roi) * (1:dim(roi)[2])) / volume)
      roi_table[k, 11] <- outer_rectangles[k, 5] + x_cm - 1
      roi_table[k, 12] <- outer_rectangles[k, 4] + y_cm - 1
      # roi asymmetries
      half_down_area  <- length(1:x_cm) * len_y
      half_up_area    <- length(x_cm:len_x) * len_y
      half_left_area  <- len_x * length(1:y_cm)
      half_right_area <- len_x * length(y_cm:len_y)
      roi_table[k, 13]<- round(((half_up_area - half_down_area) / (area)), 2)
      roi_table[k, 14] <- round(((half_right_area - half_left_area) / (area)), 2)
      # roi saturation
      if (length(saturation_minima) == 0){
        # No saturation. Do nothing
      } else {
        for (l in (1:dim(saturation_minima)[1])){
          if ((saturation_minima[l, 1] < roi_table[k, 11])
              & (saturation_minima[l, 2] > roi_table[k, 11])) {
            roi_table[k, 15] <- 1
            break
          }
        }
      }
      # roi signal to noise ratio
      roi_table[k, 16]  <- compute_power(roi)/ noise_power
    }

    roi_table <- as.data.frame(roi_table)
    aux_list$data$data_df <- aux
    aux_list$data$roi_df <- roi_table
    M <- aux_list

    # 11)  save results
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}


