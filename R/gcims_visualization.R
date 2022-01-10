#' Sample matrix visualization


#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param sample_num      Numeric. Identifies the sample to be visualized from
#'   the dataset.
#' @param rt_range        Min and Max retention time values. If NULL the
#'   complete retention time range is used. NULL by default.
#' @param dt_range        Min and Max drift time values. If NULL the complete
#'   drift time range is used. NULL by default.
#' @return An image of selected data sample.
#' @details `gcims_view_sample` uses represents gcims data as a raster. In
#'   this plot, retention time increases along the y-axis from down to up, while
#'   drift time does it along the x-axis left to right. `gcims_view_sample`
#'   provides to the user qualitative information about what are regions of
#'   interest to be analyzed in a sample.
#'
#' @note `gcims_view_sample` can't provide a reliable visualization of a
#'   data sample if its sampling frequencies  along drift and/or retention time
#'   are not constants. To overcome this problem consider using the function
#'   `gcims_interpolate` (check also the example of use of this function).
#' @family Visualization functions
#' @references {
#' Wickham, Hadley. "ggplot2." Wiley Interdisciplinary Reviews: Computational Statistics 3.2 (2011): 180-185.
#'  }
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme_minimal geom_raster
#' @importFrom viridis scale_fill_viridis
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' sample_num <- 3
#' gcims_view_sample(dir_in, sample_num, rt_range = NULL, dt_range = NULL)
#' setwd(current_dir)
#'
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
  aux_list <- readRDS(aux_string) #new
  aux <- (as.matrix(aux_list$data$data_df)) #new

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

  } else if(methods::is(sel_index_rt, "integer") & (sel_index_rt[2] > sel_index_rt[1])){
  } else {
    stop("Possible errors: 1) The selected vector of indexes corresponding to the provided retention time range is not an integer vector, 2) or rt_range[2] <= rt_range[1])")
  }

  if(is.null(dt_range)){

  } else if(methods::is(sel_index_dt, "integer") & (sel_index_dt[2] > sel_index_dt[1])){
  } else {
    stop("Possible errors: 1) The selected vector of indexes corresponding to the provided drift time range is not an integer vector, 2) or dt_range[2] <= dt_range[1])")
  }

  retention_time <- retention_time[sel_index_rt]
  drift_time <- drift_time[sel_index_dt]

  aux <- aux[sel_index_dt, sel_index_rt]#old
  rownames(aux) <- drift_time #old
  colnames(aux) <- retention_time

  moltaux <- melt((aux))
  colnames(moltaux) <- c("Drift_Time", "Retention_Time", "Value")


  #We do this in order to plot the data using geom_raster that is faster than geom_tile
  #perhaps a previous interpolation is needed to avoid this patch:
   rep_dt_index <- rep(seq(from = 1, to = dim(aux)[1], by = 1), times = dim(aux)[2])
  # drift_time_period <- mean(diff(drift_time))
  # corr_drift_time <- seq(from = drift_time[1], by = drift_time_period, length.out = length(drift_time))
  # moltaux$Drift_Time <- corr_drift_time[rep_dt_index]
  moltaux$Drift_Time <- drift_time[rep_dt_index]

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
}



#' Plots TICs or EICs.
#'
#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   sample to be visualized from the dataset.
#' @param rt_range        Min and Max retention time values. If NULL the
#'   complete retention time range is used. NULL by default.
#' @param dt_value        Numeric. Selects the drift time to show the
#'   corresponding Extracted Ion Chromatograms (EIC) of samples. If NULL, the
#'   Total Ion Chromatogram (TIC) of samples is shown. NULL by default.
#' @param colorby         A variable included in the metadata used to color
#'   sample chromatograms.
#' @return A plot of sample chromatograms colored according to some variable
#'   included in the metadata.
#' @details `gcims_plot_chrom` uses the function plots a set of sample
#'   chromatograms at given drift time (or their TIC) and for a specific range
#'   of retention times. Use `gcims_plot_chrom` to visualize the effects of
#'   digital smoothing, baseline correction and signal alignment algorithms
#'   along the retention time axis.
#' @family Visualization functions
#' @references {
#' Wickham, Hadley. "ggplot2." Wiley Interdisciplinary Reviews: Computational Statistics 3.2 (2011): 180-185.
#'  }
#' @family Visualization functions
#' @export
#' @importFrom dplyr mutate
#' @importFrom magrittr '%>%'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme_minimal geom_line
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' samples <- c(3, 7 ,8, 14, 20, 21)
#' gcims_plot_chrom(dir_in, samples, dt_value = NULL,  rt_range = NULL, colorby = "Class")
#' setwd(current_dir)
#'

gcims_plot_chrom <- function(dir_in, samples, dt_value = NULL, rt_range = NULL, colorby = "Name"){

  Retention_Time <- Index <- Value <- Sample <- Class <- NULL

  print(" ")
  print("  ///////////////////////////////")
  print(" /     Plotting Chromatograms  /")
  print("///////////////////////////////")
  print(" ")


  if(is.null(dt_value)){
    #it's OK if it's null
  } else if((is.numeric(dt_value)) & (length(dt_value) == 1)){ #MODIFICAR ESTO PARA TIEMPOS NO PARA INDICES (O ANTES)
    if(dt_value >= 0){
      #numeric, positive or zero and scalar it's: OK
    } else{
      stop ("The variable dt_value can't be negative")
    }
  } else{
    stop ("Incorrect input value for dt_value: It must be either NULL or a non-negative scalar numeric")
  }

  setwd(dir_in)
  aux_string <- paste0("M", samples[1], ".rds")
  aux_list <- readRDS(aux_string) #new
  aux <- (as.matrix(aux_list$data$data_df)) #new

  retention_time <- aux_list$data$retention_time
  drift_time <- aux_list$data$drift_time
  cond_1_rt <- (rt_range[1] - retention_time[1]) < 0
  cond_2_rt <- (rt_range[2] - retention_time[length(retention_time)]) > 0
  cond_1_dt <-(dt_value - drift_time[1]) < 0
  cond_2_dt <-(dt_value - drift_time[length(drift_time)]) > 0



  if(is.null(rt_range)){# old
    rt_ind <- c(1, dim(aux)[2]) #New

  } else{
    if(cond_1_rt | cond_2_rt){
      stop("Retention time range out of bounds.")
    }
    rt_ind  <- c(which.min(abs(retention_time - rt_range[1])), which.min(abs(retention_time - rt_range[2])))
    if( rt_ind[1] == rt_ind[2]){
      stop("Initial and Final retention time values can't be equal in the variable rt_range.")
    }#New
  }

  sel_index_rt <- rt_ind[1]: rt_ind[2]#old

  if(is.null(dt_value)){# old
    dt_ind <- c(1, dim(aux)[1])
    sel_index_dt <- dt_ind[1]: dt_ind[2]#New
  } else{
    if(cond_1_dt | cond_2_dt){
      stop("Drift time range out of bounds.")
    }
    dt_ind  <- which.min(abs(drift_time - dt_value)) #New
    sel_index_dt <- dt_ind
  }

  if(is.null(rt_range)){

  } else if(methods::is(sel_index_rt, "integer") & (sel_index_rt[2] > sel_index_rt[1])){
  } else {
    stop("Possible errors: 1) The selected vector of indexes corresponding to the provided retention time range is not an integer vector, 2) or rt_range[2] <= rt_range[1])")
  }

  rm(aux_string, aux_list, aux)

  m <- 0
  colorp <- NULL
  chroms <- matrix(0, nrow = length(sel_index_rt), ncol = length(samples))
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- (as.matrix(aux_list$data$data_df)) #new
    if (is.null(dt_value)){
      chroms[, m] <- colSums(aux[sel_index_dt, sel_index_rt])#new
    } else {
      chroms[, m] <- aux[sel_index_dt,  sel_index_rt]
    }
    colorp <- c(colorp, as.character(aux_list$metadata[,colorby]))
    rm(aux_string, aux)
  }
  rm(m)


  rownames(chroms) <- retention_time[sel_index_rt]
  moltchroms <- melt(chroms)
  colnames(moltchroms) <- c( "Retention_Time", "Index", "Value")

  if (is.null(dt_value)){
    plot_title <- "Total Ion Chromatogram"
  } else{
    plot_title <- "Extracted Ion Chromatogram"
  }

  colorp <- as.vector(as.factor(colorp))
  moltchroms <- moltchroms %>%
    mutate (Sample = as.factor(samples[Index])) %>%
    mutate (Class = as.factor(colorp[Index]))

  p <- ggplot(moltchroms, aes(x = Retention_Time, y = Value, color = Class)) +
    geom_line() +
    labs(x="Retention Time (s)",
         y="Intensity (a.u.)",
         color = "Class",
         title = plot_title) +
    theme_minimal()
  print(p)
}



#' Plots TIS or spectra


#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   sample to be visualized from the dataset.
#' @param dt_range        Min and Max drift time values. If NULL the complete
#'   drift time range is used. NULL by default.
#' @param rt_value        Numeric. Selects the retention time to show the
#'   corresponding sample spectra . If NULL, the Total Ion Spectra (TIS) of
#'   samples is shown. NULL by default.
#' @param colorby         A variable included in the metadata used to color
#'   sample spectra.
#' @return A plot of spectra colored according to some variable included in the
#'   metadata.
#' @details `gcims_plot_spec` plots a set of sample spectra at given
#'   retention time (or their TIS) and for a specific range of drift times. Use
#'   `gcims_plot_spec` to visualize the effects of digital smoothing,
#'   baseline correction and signal alignment algorithms along the drift time
#'   axis.
#' @family Visualization functions
#' @references {
#' Wickham, Hadley. "ggplot2." Wiley Interdisciplinary Reviews: Computational Statistics 3.2 (2011): 180-185.
#'  }
#' @export
#' @importFrom dplyr mutate
#' @importFrom magrittr '%>%'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes labs theme_minimal geom_line
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' samples <- c(3, 7 ,8, 14, 20, 21)
#' gcims_plot_spec(dir_in, samples, rt_value = NULL, dt_range = NULL, colorby = "Class")
#' setwd(current_dir)
#'
gcims_plot_spec <- function(dir_in, samples, rt_value = NULL, dt_range = NULL, colorby = "Name"){

  Drift_Time <- Index <- Value <- Sample <- Class <- NULL

  print(" ")
  print("  ///////////////////////////")
  print(" /     Plotting Spectra    /")
  print("///////////////////////////")
  print(" ")


  if(is.null(rt_value)){
    #it's OK if it's null
  } else if((is.numeric(rt_value)) & (length(rt_value) == 1)){ #MODIFICAR ESTO PARA TIEMPOS NO PARA INDICES (O ANTES)
    if(rt_value >= 0){
      #numeric, positive or zero and scalar it's: OK
    } else{
      stop ("The variable rt_value can't be negative")
    }
  } else{
    stop ("Incorrect input value for rt_value: It must be either NULL or a non-negative scalar numeric")
  }

  setwd(dir_in)
  aux_string <- paste0("M", samples[1], ".rds")
  aux_list <- readRDS(aux_string) #new
  aux <- (as.matrix(aux_list$data$data_df)) #new


  retention_time <- aux_list$data$retention_time
  drift_time <- aux_list$data$drift_time
  cond_1_rt <- (rt_value - retention_time[1]) < 0
  cond_2_rt <- (rt_value - retention_time[length(retention_time)]) > 0
  cond_1_dt <- (dt_range[1] - drift_time[1]) < 0
  cond_2_dt <- (dt_range[2] - drift_time[length(drift_time)]) > 0


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

  sel_index_dt <- dt_ind[1]: dt_ind[2]#old

  if(is.null(rt_value)){# old
    rt_ind <- c(1, dim(aux)[2])
    sel_index_rt <- rt_ind[1]: rt_ind[2]#New
  } else{
    if(cond_1_rt | cond_2_rt){
      stop("Retention time range out of bounds.")
    }
    rt_ind  <- which.min(abs(retention_time - rt_value)) #New
    sel_index_rt <- rt_ind
  }

  if(is.null(dt_range)){

  } else if (methods::is(sel_index_dt, "integer") & (sel_index_dt[2] > sel_index_dt[1])){
  } else {
    stop("Possible errors: 1) The selected vector of indexes corresponding to the provided drift time range is not an integer vector, 2) or dt_range[2] <= dt_range[1])")
  }

  rm(aux_string, aux_list, aux)


  m <- 0
  colorp <- NULL
  specs <- matrix(0, nrow = length(sel_index_dt), ncol = length(samples))
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- (as.matrix(aux_list$data$data_df)) #new

    if (is.null(rt_value)){
      specs[, m] <- rowSums(aux[sel_index_dt, sel_index_rt])
    } else {
      specs[, m] <- aux[sel_index_dt, sel_index_rt]
    }
    colorp <- c(colorp, as.character(aux_list$metadata[,colorby]))
    rm(aux_string, aux)
  }
  rm(m)

  rownames(specs) <- drift_time[sel_index_dt]
  moltspecs <- melt(specs)
  colnames(moltspecs) <- c("Drift_Time", "Index", "Value")

  if (is.null(rt_value)){
    plot_title <- "Total Ion Spectrum"
  } else{
    plot_title <- "Spectrum"
  }

  colorp <- as.vector(as.factor(colorp))
  moltspecs <- moltspecs %>%
    mutate (Sample = as.factor(samples[Index])) %>%
    mutate (Class = as.factor(colorp[Index]))

  p <- ggplot(moltspecs, aes(x = Drift_Time, y = Value, color = Class)) +
    geom_line() +
    labs(x="Drift Time (ms)",
         y="Intensity (a.u.)",
         color = "Class",
         title = plot_title) +
    theme_minimal()
  print(p)
}

