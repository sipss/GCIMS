
#' Two dimensional peak picking

#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         Samples to which perform the peak picking.
#' @return A peak table per sample.
#' @family Peak Grouping functions
#' @export
#' @importFrom cluster pam
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }


gcims_peak_grouping <- function(dir_in, dir_out, samples){


  #-------------#
  #     MAIN    #
  #-------------#


  setwd(dir_in)
  print("Performing Peak Clustering among samples")
  roi_list <- vector(mode = "list", length = length(samples))
  m <- 0
  for (i in samples){
    m <- m + 1

    # 1)   read data
    aux_string <- paste0("M", samples[m], ".rds")
    aux_list <- readRDS(aux_string) #new
    roi_list[[m]]   <- aux_list$data$roi_df
  }



  all_roi_df <- do.call(rbind, roi_list)
  all_roi_df <- all_roi_df[all_roi_df$snr >= 1, ]
  to_be_clustered <- (all_roi_df[c("rt_mc", "dt_mc")])

   #necesito el número de rois de la muestra promedio
  cl <- pam(to_be_clustered, k = 209, medoids = 1:209, cluster.only = TRUE)
  plot(all_roi_df$dt_mc, all_roi_df$rt_mc, col = cl, pch =  cl, xlim=c(0,4000))

  #mediods (podrian ser los picos de la muestra promedio!!!)
 # hola <- factoextra::fviz_nbclust(lista[c("rt_mc", "dt_mc")], cluster::pam, method = "silhouette", k.max = 220)
  # 11)  save results
  setwd(dir_out)
  saveRDS(all_roi_df, file = paste0("all_roi_df", i, ".rds"))
  setwd(dir_in)

  return(all_roi_df)
}

#
# MS7 <- ms(mc_df)
# out2 <- MS7$cluster.center
# plot(mc_df$dt_mc, mc_df$rt_mc, col = out2)
# table(out2)
