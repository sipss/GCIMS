#' SAmples Visualization


#' @param dir_in          The input directory.
#' @param sample_num      The number corresponding to the sample
#'                        to be visualized.

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

gcims_visualization <- function(dir_in, sample_num){

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
    colnames(moltaux) <- c("Drift Time", "Retention Time", "Value")

    rm(aux, aux_string)
    p <- ggplot(moltaux, aes(x = `Drift Time`, y = `Retention Time`, fill = Value)) +
          geom_raster() +
          scale_fill_viridis(discrete = FALSE, option = "A", direction = -1) +
          labs(x="Drift Time (a.u.)",
               y="Retention Time (a.u.)",
               title = "Sample Image",
               fill = "Intensity") +
          theme_minimal()
    print(p)
}
