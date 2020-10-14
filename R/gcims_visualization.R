#' Sample matrix visualization


#' @param dir_in          The input directory.
#' @param sample_num      The number corresponding to the sample
#'                        to be visualized.
#' @param rt_range        Min a Max retention time values. NULL by
#'                        default.
#' @param dt_range        Min a Max drift time values. NULL by
#'                        default.
#' @return An image of the sample matrix.

#' @family Visualization functions
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme_minimal geom_raster
#' @importFrom viridis scale_fill_viridis
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_view_sample <- function(dir_in, sample_num, rt_range = NULL, dt_range = NULL){

  Retention_Time <- Drift_Time <- Value <- NULL

  print(" ")
  print("  /////////////////////////////////////")
  print(" /    Sample Matrix Visualization    /")
  print("/////////////////////////////////////")
  print(" ")

  setwd(dir_in)
  print(paste0("Visualizing sample ", sample_num))
  aux_string <- paste0("M", sample_num, ".rds")
  #aux <- readRDS(aux_string) #OLD
  aux_list <- readRDS(aux_string) #new
  aux <- t(as.matrix(aux_list$data$data_df)) #new

#SOME CHECKS
  retention_time <- aux_list$data$retention_time
  drift_time <- aux_list$data$drift_time

  if(is.null(rt_range)){# old
    rt_ind <- c(1, dim(aux)[1]) #New

  } else{
    rt_ind  <- c(which.min(abs(retention_time - rt_range[1])), which.min(abs(retention_time - rt_range[2]))) #New
  }



  if(is.null(dt_range)){# old
    dt_ind <- c(1, dim(aux)[2]) #New
  } else{
    dt_ind  <- c(which.min(abs(drift_time - dt_range[1])), which.min(abs(drift_time - dt_range[2]))) #New
  }

  sel_index_rt <- rt_ind[1]: rt_ind[2]#old
  sel_index_dt <- dt_ind[1]: dt_ind[2]#old

  retention_time <- retention_time[sel_index_rt]
  drift_time <- drift_time[sel_index_dt]

  aux <- aux[sel_index_rt, sel_index_dt]#old
  rownames(aux) <- retention_time #old

  moltaux <- melt(t(aux))
  colnames(moltaux) <- c("Drift_Time", "Retention_Time", "Value")


#We do this in order to plot the data using geom_raster that is faster than geom_tile
#perhaps a previous interpolation is needed to avoid this patch:
  rep_dt_index <- rep(seq(from = 1, to = dim(aux)[2], by = 1), times = dim(aux)[1])
  #rep_dt_index <- rep(sel_index_dt, times = dim(aux)[1])
  drift_time_period <- mean(diff(drift_time))
  #corr_drift_time <- seq(from = sel_index_dt[1], by = drift_time_period, length.out = length(drift_time))
  corr_drift_time <- seq(from = drift_time[1], by = drift_time_period, length.out = length(drift_time))
  moltaux$Drift_Time <- corr_drift_time[rep_dt_index]

  rm(aux, aux_string)
  p <- ggplot(moltaux, aes(x = Drift_Time, y = Retention_Time, fill = Value)) +
    geom_raster(interpolate = FALSE) +
    scale_fill_viridis(discrete = FALSE, option = "A", direction = -1) +
    labs(x="Drift Time (ms)",
         y="Retention Time (s)",
         title = "Sample Matrix Image",
         fill = "Intensity") +
    theme_minimal()
  print(p)


  # else if((class(rt_ind[1]: rt_ind[2]) == "integer") & (rt_ind[2] > rt_ind[1])){
  #   #else if((class(rt_range[1]: rt_range[2]) == "integer") & (rt_range[2] > rt_range[1])){ # old
  #   if((rt_ind[1] >= 1) & (rt_ind[2] <= dim(aux)[1])){ # NEw
  #   #if((rt_range[1] >= 1) & (rt_range[2] <= dim(aux)[1])){ OLD
  #   } else {
  #     stop("Index out of bounds")
  #   }
  # } else {
  #   stop("Possible errors: 1) rt_range is not a valid (numeric) vector, 2) or rt_range[2] <= rt_range[1])")
  # }
  #
  # if(is.null(dt_range)){ #old
  # } else if((class(dt_ind[1]: dt_ind[2]) == "integer") & (dt_ind[2] > dt_ind[1])){
  #   #else if((class(dt_range[1]: dt_range[2]) == "integer") & (dt_range[2] > dt_range[1])){ #OLD
  #   if((dt_ind[1] >= 1) & (dt_ind[2] <= dim(aux)[2])){ #MODIFICAR ESTO PARA TIEMPOS NO PARA INDICES (O ANTES)
  #   #if((dt_range[1] >= 1) & (dt_range[2] <= dim(aux)[2])){ OLD
  #   } else {
  #     stop("Index out of bounds")
  #   }
  # } else {
  #   stop("Possible errors: 1) dt_range is not a valid (numeric) vector, 2) or dt_range[2] <= dt_range[1])")
  # }

  #END CHECKS

  # if (is.null(rt_range)){#old
  #   sel_index_rt <- 1:dim(aux)[1]#old
  #   #retention_time <- c(0:(dim(aux)[1] - 1))#old
  # } else{#New
  #   sel_index_rt <- rt_ind[1]: rt_ind[2]
  #   #sel_index_rt <- rt_range[1]: rt_range[2]#old
  #   #retention_time <- c((rt_range[1] - 1):(rt_range[2] - 1)) # old#New (OJOOOO: que despues seran tiempos, no indices)
  #   retention_time <- retention_time[sel_index_rt]   #New
  # }#New
  #
  #
  # if (is.null(dt_range)){#old
  #   sel_index_dt <- 1:dim(aux)[2]#old
  # } else{#New
  #   sel_index_dt <- dt_ind[1]: dt_ind[2] #New
  #   drift_time <- drift_time[sel_index_dt] #New
  #   #New
  #   #sel_index_dt <- dt_range[1]: dt_range[2]#old
  # }#New

  # aux <- aux[sel_index_rt, sel_index_dt]#old
  # rownames(aux) <- retention_time #old
  #
  # moltaux <- melt(t(aux))
  # colnames(moltaux) <- c("Drift_Time", "Retention_Time", "Value")
  #
  # rm(aux, aux_string)
  # p <- ggplot(moltaux, aes(x = Drift_Time, y = Retention_Time, fill = Value)) +
  #   geom_raster() +
  #   scale_fill_viridis(discrete = FALSE, option = "A", direction = -1) +
  #   labs(x="Drift Time (ms)",
  #        y="Retention Time (s)",
  #        title = "Sample Matrix Image",
  #        fill = "Intensity") +
  #   theme_minimal()
  # print(p)
}



#' Plots TICs or EICs.

#' @param dir_in          The input directory.
#' @param samples         The set of samples to be
#'                        to be visualized.
#' @param rt_range        Min a Max retention time values. NULL by
#'                        default.
#' @param td_ind          Selected drift time value. NULL by default.
#' @return A set of Chromatograms.
#' @family Visualization functions
#' @export
#' @importFrom dplyr mutate
#' @importFrom magrittr '%>%'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme_minimal geom_line
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_plot_chrom <- function(dir_in, samples, td_ind = NULL, rt_range = NULL){

  Retention_Time <- Index <- Value <- Sample <- NULL

  print(" ")
  print("  ///////////////////////////////")
  print(" /     Plottig Chromatograms   /")
  print("///////////////////////////////")
  print(" ")

  if(is.null(td_ind)){
    #it's OK if it's null
  } else if((is.numeric(td_ind)) & (length(td_ind) == 1)){ #MODIFICAR ESTO PARA TIEMPOS NO PARA INDICES (O ANTES)
    if(td_ind > 0){
    #numeric, positive and scalar it's: OK
    } else{
      stop ("The variable td_ind can't be negative")
    }
  } else{
    stop ("Incorrect input values for td_ind: It must be either NULL or a positive scalar integer")
  }

  setwd(dir_in)
  aux_string <- paste0("M", samples[1], ".rds")
  aux <- readRDS(aux_string)

  if(is.null(rt_range)){
    num_of_rows <- dim(aux)[1]
    sel_index <- 1:num_of_rows
  } else if((class(rt_range[1]: rt_range[2]) == "integer") & (rt_range[2] > rt_range[1])){
    if((rt_range[1] >= 1) & (rt_range[2] <= dim(aux)[1])){ #MODIFICAR ESTO PARA TIEMPOS NO PARA INDICES (O ANTES)
      num_of_rows <- length(rt_range[1]: rt_range[2])
      sel_index <- rt_range[1]: rt_range[2]
      } else {
        stop("Index out of bounds")
        }
  } else {
        stop("Possible errors: 1) rt_range is not an integer vector, 2) or rt_range[2] <= rt_range[1])")

  }

  rm(aux_string, aux)

  m <- 0
  chroms <- matrix(0, nrow = num_of_rows, ncol = length(samples))
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    if (is.null(td_ind)){
      chroms[, m] <- rowSums(aux[sel_index, ])
    } else {
      chroms[, m] <- aux[sel_index, td_ind]
    }
    rm(aux_string, aux)
  }
  rm(m)

  retentiontime <- c(0:(dim(chroms)[1]-1))
  rownames(chroms) <- retentiontime
  moltchroms <- melt(chroms)
  colnames(moltchroms) <- c( "Retention_Time", "Index", "Value")

  if (is.null(td_ind)){
    plot_title <- "Total Ion Chromatogram"
  } else{
    plot_title <- "Extracted Ion Chromatogram"
  }

  moltchroms <- moltchroms %>%
                mutate (Sample = as.factor(samples[Index]))


  p <- ggplot(moltchroms, aes(x = Retention_Time, y = Value, color = Sample)) +
    geom_line() +
    labs(x="Retention Time (a.u.)",
         y="Intensity (a.u.)",
         color = "Sample",
         title = plot_title) +
    theme_minimal()
  print(p)
}



#' Plots TISa or spectra

#' @param dir_in          The input directory.
#' @param samples         The set of samples to be
#'                        to be visualized.
#' @param dt_range        Min a Max drift time values. NULL by
#'                        default.
#' @param tr_ind          Selected retention time value. NULL by default.
#' @return A set of spectra.
#' @family Visualization functions
#' @export
#' @importFrom dplyr mutate
#' @importFrom magrittr '%>%'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme_minimal geom_line
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_plot_spec <- function(dir_in, samples, tr_ind = NULL, dt_range = NULL){

  Drift_Time <- Index <- Value <- Sample <- NULL

  print(" ")
  print("  ///////////////////////////")
  print(" /     Plotting Spectra    /")
  print("///////////////////////////")
  print(" ")

  if(is.null(tr_ind)){
    #it's OK if it's null
  } else if((is.numeric(tr_ind)) & (length(tr_ind) == 1)){ #MODIFICAR ESTO PARA TIEMPOS NO PARA INDICES (O ANTES)
    if(tr_ind > 0){
      #numeric, positive and scalar it's: OK
    } else{
      stop ("The variable tr_ind can't be negative")
    }
  } else{
    stop ("Incorrect input values for tr_ind: It must be either NULL or a positive scalar integer")
  }

  setwd(dir_in)
  aux_string <- paste0("M", samples[1], ".rds")
  aux <- readRDS(aux_string)

  if(is.null(dt_range)){
    num_of_rows <- dim(aux)[2]
    sel_index <- 1:num_of_rows
  } else if((class(dt_range[1]: dt_range[2]) == "integer") & (dt_range[2] > dt_range[1])){
    if((dt_range[1] >= 1) & (dt_range[2] <= dim(aux)[1])){ #MODIFICAR ESTO PARA TIEMPOS NO PARA INDICES (O ANTES)
      num_of_rows <- length(dt_range[1]: dt_range[2])
      sel_index <- dt_range[1]: dt_range[2]
    } else {
      stop("Index out of bounds")
    }
  } else {
    stop("Possible errors: 1) dt_range is not an integer vector, 2) or dt_range[2] <= dt_range[1])")

  }

  rm(aux_string, aux)

  m <- 0
  specs <- matrix(0, nrow = num_of_rows, ncol = length(samples))
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    if (is.null(tr_ind)){
      specs[, m] <- colSums(aux[, sel_index])
    } else {
      specs[, m] <- aux[tr_ind, sel_index]
    }
    rm(aux_string, aux)
  }
  rm(m)

  driftime <- c(0:(dim(specs)[1]-1))
  rownames(specs) <- driftime
  moltspecs <- melt(specs)
  colnames(moltspecs) <- c("Drift_Time", "Index", "Value")

  if (is.null(tr_ind)){
    plot_title <- "Total Ion Spectrum"
  } else{
    plot_title <- "Spectrum"
  }

  moltspecs <- moltspecs %>%
    mutate (Sample = as.factor(samples[Index]))


  p <- ggplot(moltspecs, aes(x = Drift_Time, y = Value, color = Sample)) +
    geom_line() +
    labs(x="Retention Time (a.u.)",
         y="Intensity (a.u.)",
         color = "Sample",
         title = plot_title) +
    theme_minimal()
  print(p)
}

