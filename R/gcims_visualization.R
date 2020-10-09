#' Sample matrix visualization


#' @param dir_in          The input directory.
#' @param sample_num      The number corresponding to the sample
#'                        to be visualized.
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

gcims_visualization <- function(dir_in, sample_num){

  Retention_Time <- Drift_Time <- Value <- NULL

  print(" ")
  print("  /////////////////////////////////////")
  print(" /    Sample Matrix Visualization    /")
  print("/////////////////////////////////////")
  print(" ")

  setwd(dir_in)
    print(paste0("Visualizing sample ", sample_num))
    aux_string <- paste0("M", sample_num, ".rds")
    aux <- readRDS(aux_string)
    retentiontime <- c(0:(dim(aux)[1]-1))
    rownames(aux) <- retentiontime
    moltaux <- melt(t(aux))
    colnames(moltaux) <- c("Drift_Time", "Retention_Time", "Value")

    rm(aux, aux_string)
    p <- ggplot(moltaux, aes(x = Drift_Time, y = Retention_Time, fill = Value)) +
          geom_raster() +
          scale_fill_viridis(discrete = FALSE, option = "A", direction = -1) +
          labs(x="Drift Time (a.u.)",
               y="Retention Time (a.u.)",
               title = "Sample Matrix Image",
               fill = "Intensity") +
          theme_minimal()
    print(p)
}



#' Plots TICs or EICs.

#' @param dir_in          The input directory.
#' @param samples         The set of samples to be
#'                        to be visualized.
#' @param rt_range        Min a Max retention time values. NULL by
#'                        default.
#' @param td_ind          Selected drift time value. NULL by default.
#' @return A plot of the Total Ion Chromatograms.
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
      sel_index <- 1:num_of_rows
      } else {
        stop("Index out of bounds")
        }
  } else {
        stop("Possible errors: 1) rt_range is not an integer vector, 2) or rt_range[2] <= rt_range[1])")

  }

  rm(aux_string, aux)

  m <- 0
  tics <- matrix(0, nrow = num_of_rows, ncol = length(samples))
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    if (is.null(td_ind)){
      tics[, m] <- rowSums(aux[sel_index, ])
    } else {
      tics[, m] <- aux[sel_index, td_ind]
    }
    #tics[, m] <- rowSums(aux[sel_index, ])
    rm(aux_string, aux)
  }
  rm(m)

  retentiontime <- c(0:(dim(tics)[1]-1))
  rownames(tics) <- retentiontime
  moltics <- melt(tics)
  colnames(moltics) <- c( "Retention_Time", "Index", "Value")

  if (is.null(td_ind)){
    plot_title <- "Total Ion Chromatogram"
  } else{
    plot_title <- "Extracted Ion Chromatogram"
  }

  moltics <- moltics %>%
                mutate (Sample = as.factor(samples[Index]))


  p <- ggplot(moltics, aes(x = Retention_Time, y = Value, color = Sample)) +
    geom_line() +
    labs(x="Retention Time (a.u.)",
         y="Intensity (a.u.)",
         color = "Sample",
         title = plot_title) +
    theme_minimal()
  print(p)
}


#' Total Ion Spectrum visualization


#' @param dir_in          The input directory.
#' @param samples         The set of samples to be
#'                        to be visualized.

#' @return A plot of the Total Ion Spectra.
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

gcims_plot_tisa <- function(dir_in, samples){

  Drift_Time <- Index <- Value <- Sample <- NULL

  print(" ")
  print("  //////////////////////")
  print(" /     Plottig TISa   /")
  print("//////////////////////")
  print(" ")

  setwd(dir_in)
  aux_string <- paste0("M", samples[1], ".rds")
  aux <- readRDS(aux_string)
  num_of_rows <- dim(aux)[2]
  rm(aux_string, aux)

  m <- 0
  tisa <- matrix(0, nrow = num_of_rows, ncol = length(samples))
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    tisa[, m] <- colSums(aux)
    rm(aux_string, aux)
  }
  rm(m)

  driftime <- c(0:(dim(tisa)[1]-1))
  rownames(tisa) <-  driftime
  moltisa <- melt(tisa)
  colnames(moltisa) <- c( "Drift_Time", "Index", "Value")

  moltisa <- moltisa %>%
    mutate (Sample = as.factor(samples[Index]))


  p <- ggplot(moltisa, aes(x = Drift_Time, y = Value, color = Sample)) +
    geom_line() +
    labs(x="Drift Time (a.u.)",
         y="Intensity (a.u.)",
         color = "Samples",
         title = "Total Ion Spectra") +
    theme_minimal()
  print(p)
}


#' Plot Sample Spectra


#' @param dir_in          The input directory.
#' @param samples         The set of samples to be
#'                        to be visualized.
#' @param tr_ind          Retention time Index (to be modified in future)

#' @return A plot of the Sample Spectra.
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

gcims_plot_spec <- function(dir_in, samples, tr_ind){

  Drift_Time <- Index <- Value <- Sample <- NULL

  print(" ")
  print("  /////////////////////////")
  print(" /     Plottig Spectra   /")
  print("/////////////////////////")
  print(" ")

  setwd(dir_in)
  aux_string <- paste0("M", samples[1], ".rds")
  aux <- readRDS(aux_string)
  num_of_rows <- dim(aux)[2]
  rm(aux_string, aux)

  m <- 0
  spec <- matrix(0, nrow = num_of_rows, ncol = length(samples))
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    spec[, m] <- aux[tr_ind,]
    rm(aux_string, aux)
  }
  rm(m)

  driftime <- c(0:(dim(spec)[1]-1))
  rownames(spec) <-  driftime
  moltspec <- melt(spec)
  colnames(moltspec) <- c( "Drift_Time", "Index", "Value")

  moltspec <- moltspec %>%
    mutate (Sample = as.factor(samples[Index]))


  p <- ggplot(moltspec, aes(x = Drift_Time, y = Value, color = Sample)) +
    geom_line() +
    labs(x="Drift Time (a.u.)",
         y="Intensity (a.u.)",
         color = "Samples",
         title = "Spectra") +
    theme_minimal()
  print(p)
}




