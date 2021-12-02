#' ROIs Clustering

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where peak table data file is
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples that are going to be included in the peak table.
#' @details \code{gcims_peaks_clustering}  does a clustering along the ROIs
#'   for a peak table creation. The Figures of merits of each ROI are also
#'   reported. In this table are included all samples in \code{samples}. Use this
#'   function if you are interested in obtaining a final peak table for future
#'   classification techniques.
#' @return A Set of S3 objects.
#' @family Peaks CLustering function
#' @export
#' @importFrom cluster pam
#' @importFrom plyr count
#' @importFrom Hotelling hotelling.test
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 1:3
#'
#' # Example of ROIs Clustering for Peak Table Creation
#' gcims_peaks_clustering <- function(dir_in, dir_out, samples)



gcims_peaks_clustering <- function(dir_in, dir_out, samples){
  print(" ")
  print("  ///////////////////////////")
  print(" /    Clustering the ROIs  /")
  print("///////////////////////////")
  print(" ")

  setwd(dir_in)

  imgs <- NULL
  AsF <- NULL
  volume <- NULL
  saturation <- NULL
  peak <- NULL
  nrois <- NULL

  m = 0
  for (i in samples){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))

    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- (as.matrix(aux_list$data$data_df))

    imgs[[i]] <- aux
    AsF[[i]] <- aux_list$data$Parameters["AsF", ]
    volume[[i]] <- aux_list$data$Parameters["volume", ]
    saturation[[i]] <- aux_list$data$Parameters["saturation", ]
    peak[[i]] <- aux_list$data$Peaks
    nrois[[i]] <- aux_list$data$ROIs
  }

  AsF <- unlist(AsF)
  volume <- unlist(volume)
  saturation <- unlist(saturation)
  datapoints <- do.call(rbind.data.frame, peak)
  ROIs <- do.call(rbind.data.frame, nrois)


  lengthroissamples <- NULL
  for(j in 1:length(samples)){
    lengthrois <- dim(nrois[[j]])[1]
    lengthroissamples <- c(lengthroissamples, lengthrois)
  }

  N_clusters <- max(lengthroissamples)

  D <- kmed::distNumeric(as.matrix(datapoints), as.matrix(datapoints), method = "sev") ##Revisar
  ## cluster <- clue::kmedoids(D, N_clusters) Muy pesado
  cluster <- cluster::pam(datapoints,N_clusters, metric = "manhattan") #Programme from scratch
  idx <- cluster$clustering
  idx_post <- matrix(0, length(idx), 3)
  k <- length(unique(idx))

  for (n in 1:k){
    idx_tmp <- which(idx == n)
    acc = 1
    for (j in 1:length(samples)){
      boolean <- idx_tmp >= acc & idx_tmp <= (acc-1) + length(nrois[[j]])
      if (sum(boolean) >= 1){
        dist <- D[idx_tmp[boolean], n]
        mindist <- which.min(dist)
        pos <- which(boolean == TRUE)
        idx_post[idx_tmp[pos[mindist]], 1] <- n
        idx_post[idx_tmp[pos], 2] <- j
      }
      acc <- acc + lengthroissamples[j]
    }
  }

  ## Filtering

  # Filter widths
  widths <- ROIs[, 2] - ROIs[, 1]
  quartiles <- quantile(widths)
  idx_post[c(which(widths > quartiles[2]), which(widths > quartiles[4])), 3] <- 1 #1s son outliers, seran descartados

  # Filter heights
  heights <- ROIs[, 4] - ROIs[, 3]
  medianheight <- median(heights)
  quartiles <- quantile(heights)
  idx_post[c(which(heights > medianheight + quartiles[2]), which(heights > medianheight + quartiles[4])), 3] <- 1 #1s son outliers, seran descartados

  # Filter clusters with 2 or fewer ROIs (cambiar a tanto por cierto)
  N_min_clusters <- 2
  clustfreqs <- plyr::count(idx_post[idx_post[, 3] != 1, 2])

  for (j in clustfreqs[clustfreqs[,2] <= N_min_clusters, 1]){
    idx_post[idx_post[,1] == j, 3] <- 1
  }

  # Repescamiento
  pos_zeros <- which(idx_post[,1] == 0 & idx_post[,3] != 1) #Positions of ROIs where there are zeros and are not discarded

  for (j in pos_zeros){
    possible_candidates <- setdiff(1:dim(idx_post)[1], idx_post[idx_post[, 2] == idx_post[j, 2] & idx_post[, 3] != 1, 1])
    k <- intersect(setdiff(possible_candidates, 0), clustfreqs[clustfreqs[,2] > N_min_clusters, 1]) #transponer
    p_vals <- NULL
    for (n in (1:length(k))){
      clustern <- datapoints[idx_post[,1] == k[n] & idx_post[,3] != 1, ]
      if (dim(clustern)[1] > 2){
        p_val <- Hotelling::hotelling.test(datapoints[j, ], clustern)
        p_vals <- c(p_vals, p_val$pval)
      }
    }
    maxpval <- which.max(p_vals)
    if(p_vals[maxpval] > 0.05){ # If greater than p-value (0.05)
      idx_post[j, 1] <- maxpval
    }
  }

  ## Imputation and statistics

  pos <- idx_post[, 3] != 1 & idx_post[, 1] != 0
  indices_clusters <- unique(idx_post[pos, 1])
  K <- length(indices_clusters)

  stats <- vector(mode = "list", length = K)
  table <- matrix(0, length(samples), K)

  for (k in (1:K)){
    n <- indices_clusters[k]
    pos_n <- idx_post[, 1] == n & idx_post[, 3] != 1
    samples_n <- idx_post[pos_n, 2]
    rest_samples_n <- setdiff((1:length(samples)), samples_n)
    ROI_median <- floor(apply(ROIs[pos_n, ], 2, median))

    AsF_tmp <- rep(0, length(samples))
    AsF_tmp[samples_n] <- AsF[pos_n]
    AsF_tmp[rest_samples_n] <- median(AsF[pos_n])

    volume_tmp <- rep(0, length(samples))
    volume_tmp[samples_n] <- volume[pos_n]

    I <- NULL
    if (length(rest_samples_n) >= 1){
      for (l in rest_samples_n){
        I[[l]] <- imgs[[l]][ROI_median[1]:ROI_median[2], ROI_median[3]:ROI_median[4]]
        volume_tmp[rest_samples_n] <- sum(apply(I[[l]], 1, sum)) #[1, dim(I[[l]], 3)]) # size que?
      }
    }

    ROI_tmp <- matrix(0, length(samples), 4)
    ROI_tmp[samples_n, ] <- as.matrix(ROIs[pos_n, ])
    ROI_tmp[rest_samples_n, ] <- matrix(rep(ROI_median, each = length(rest_samples_n)), nrow = length(rest_samples_n))

    peaks_tmp <- matrix(0, length(samples), 2)
    peaks_tmp[samples_n, ] <- as.matrix(datapoints[pos_n, ])
    peaks_tmp[rest_samples_n, ]  <- matrix(rep(floor(apply(datapoints[pos_n, ], 2, median)), each = length(rest_samples_n)), nrow = length(rest_samples_n))

    stats[[k]]$Name <- paste("Cluster", k)
    stats[[k]]$AsF <- AsF_tmp
    stats[[k]]$volume <- volume_tmp
    table[,k] <- volume_tmp

    stats[[k]]$AsF.mean  <- mean(AsF_tmp)
    stats[[k]]$volume_mean <- mean(volume_tmp)

    stats[[k]]$AsF_std <- sd(AsF_tmp)
    stats[[k]]$volume_std <- sd(volume_tmp)

    stats[[k]]$ROIs <- ROI_tmp
    stats[[k]]$peaks <- peaks_tmp

  }


  ## Make table

  drift_time <- aux_list$data$drift_time
  fs = 1/(drift_time[2]-drift_time[1])
  peak_table <- matrix(0, (length(samples) * K), 7)

  acc <- 1
  for (i in (1:K)){
    c <- 1
    for (j in (acc : (acc + length(samples) - 1))){
      peak_table[j, 1] <- c
      peak_table[j, 2] <- i
      peak_table[j, 3] <- stats[[i]]$peaks[c, 1]
      peak_table[j, 4] <- stats[[i]]$peaks[c, 2]
      peak_table[j, 5] <- stats[[i]]$peaks[c, 1]/fs
      peak_table[j, 6] <- stats[[i]]$peaks[c, 2]/fs
      peak_table[j, 7] <- stats[[i]]$volume[c]

      c <- c + 1
    }
    acc <- acc + length(samples)
  }

  colnames(peak_table) <- c("Sample", "Cluster", "PicoX (indice)", "PicoY (indice)", "PicoX (tiempo)", "PicoY (tiempo)", "Intensity")
  setwd(dir_out)
  write.csv(peak_table, "Peaktable.csv")
  saveRDS(stats, "ROIsParameters.rds")
}


