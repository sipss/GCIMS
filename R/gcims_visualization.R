#' SAmples Visualization


#' @param dir_in          The input directory.
#' @param samples         The set of samples to be visualized.

#' @return A set of plots of the samples: full espectra, Intensity VS Retention Time
#' and Intensity VS Drift Time.
#' @family Visualization functions
#' @export
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_visualization <- function(dir_in, samples){

  print(" ")
  print("  /////////////////////////////////////")
  print(" /   Samples Spectra Visualization   /")
  print("/////////////////////////////////////")
  print(" ")

  library(ggplot2)
  library(metR)
  library(reshape2)

  setwd(dir_in)
  m = 0;
  for (i in (samples)){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    retentiontime <- c(0:(dim(aux)[1]-1))
    rownames(aux) <- retentiontime
    moltaux <- melt(t(aux))
    colnames(moltaux) <- c("Dirft Time", "Retention Time", "Value")
    ggplot(moltaux, aes(x = `Dirft Time`, y = `Retention Time`, z = Value)) +
      geom_contour_fill()
    ggplot(moltaux, aes(x = `Dirft Time`, y = Value)) +
      geom_line()
    ggplot(moltaux, aes(x = `Retention Time`, y = Value)) +
      geom_line()
  }
}
