
#' Two dimensional peak picking

#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param distance        A measure of dissimilarity of Samples. Either
#' @aliases               "mahalanobis" or "overlapping"
#' @param samples         Samples to which perform the peak picking.
#' @return A peak table per sample.
#' @family Peak Grouping functions
#' @export
#' @importFrom pracma meshgrid
#' @importFrom cluster pam
#' @importFrom StatMatch mahalanobis.dist
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }


gcims_peak_grouping <- function(dir_in, dir_out, samples, distance = "overlapping"){
#pam

  #---------#
  # sub2ind #
  #---------#

  sub2ind <- function(sub, n){
    r <- sub[1]
    c <- sub[2]
    ind = ((c-1) * n) + r
    return(ind)
  }


  #-------------#
  #     MAIN    #
  #-------------#


  setwd(dir_in)
  print("Performing Peak Clustering among samples")


  # 0) Look for the number of ROIS in the reference sample
    aux_string <- paste0("M", 0, ".rds")
    aux_list <- readRDS(aux_string)
    aux_df <- aux_list$data$roi_df
    ncluster <-length(aux_df$roi_id)
    rm(aux_string, aux_list, aux_df)

  # 1) Initilize roi_list (length equal the number of samples)
  #    and row_list (number of rows or each sample)
    roi_list <- vector(mode = "list", length = length(samples))
    row_list <- vector(mode = "list", length = length(samples))

  # 2) Loop for generating the roi_list
    m <- 0
    for (i in samples){
      m <- m + 1
      # 2.a) read data
      aux_string <- paste0("M", samples[m], ".rds")
      aux_list <- readRDS(aux_string)
      # 2.b) update roi_list
      roi_list[[m]]   <- aux_list$data$roi_df
      row_list[[m]]   <-nrow(aux_list$data$data_df)
    }

  # 3) Convert the roi_list to a dataframe (all_roi_df)
    all_roi_df <- do.call(rbind, roi_list)
    rows_df <- do.call(rbind, row_list)


  #4) Select only the center of mass coordinates for clustering
    to_be_clustered <- all_roi_df[c("rt_mc", "dt_mc", "max_rt", "min_dt", "min_rt", "max_dt")]

  if(distance == "mahalanobis"){
    to_be_clustered_dist <- mahalanobis.dist(to_be_clustered[c("rt_mc", "dt_mc")])
  } else if(distance == "overlapping"){
    hola <-  to_be_clustered[c("max_rt", "min_dt", "min_rt", "max_dt")]
  } else{
    stop("Distance must be either mahalanobis or overlapping")
  }

  #overlapping.dist
  n_dim <-  nrow(to_be_clustered) #  number of rows and columns (the matrix is square and symmetric)
  print(n_dim)
  overlap_dist <- matrix(0, nrow = n_dim, ncol = n_dim)
  for (i in (1:n_dim)){
    print(i)
    seq_x_i <- hola$min_rt[i]:hola$max_rt[i]
    seq_y_i <- hola$min_dt[i]:hola$max_dt[i]
     grid_i_sub <- meshgrid(x = seq_x_i , y= seq_y_i)
     dim(grid_i_sub$X) <- NULL
     dim(grid_i_sub$Y) <- NULL
     grid_i_sub <- cbind(grid_i_sub$X, grid_i_sub$Y)
     print(grid_i_sub)

     n <- rows_df[i]
     print(n)
     grid_i_ind <- apply(grid_i_sub, sub2ind, n, MARGIN = 1)
     print(grid_i_ind)
     #print(grid_i_ind)
     #print(grid_i_ind)
     for(j in (1: (n_dim - i + 1))){
       seq_x_j <- hola$min_rt[j]:hola$max_rt[j]
       seq_y_j <- hola$min_dt[j]:hola$max_dt[j]
       grid_j_sub <- meshgrid(x = seq_x_j , y= seq_y_j)
       dim(grid_j_sub$X) <- NULL
       dim(grid_j_sub$Y) <- NULL
       grid_j_sub <- cbind(grid_j_sub$X, grid_j_sub$Y)

       n <- rows_df[j]

       grid_j_ind <- apply(grid_j_sub, sub2ind, n, MARGIN = 1)
       #print(j)
       #print(grid_j_ind)

       grid_ij_union <- union(grid_i_ind, grid_j_ind)
       #print(grid_ij_union)
       grid_ij_intersect <- intersect(grid_i_ind, grid_j_ind)
       #print(grid_ij_intersect)
       overlap_dist[i, j] <-  1 - (length(grid_ij_intersect) / length(grid_ij_union))
       #print(1 - (length(grid_ij_intersect) / length(grid_ij_union)))
       #overlap_dist[j, i] <- overlap_dist[i, j]
     }



    }

  to_be_clustered_dist <- overlap_dist


  #my_grid <- meshgrid(x, y)

  # 5) cluster the samples according the initial 'seeds'
  #   or centers obtained from the reference samples

  #cl <- kmeans(to_be_clustered[(ncluster + 1): nrow(to_be_clustered), ], ncluster, centers = to_be_clustered[1:ncluster, ])
  # plot(all_roi_df$dt_mc, all_roi_df$rt_mc, col = cl$cluster, pch =  cl$cluster)
  #points(cl$centers[,2], cl$centers[, 1], col = "black", pch = 4)
  color_sampled <- sample(colors(distinct = TRUE), length(colors(distinct = TRUE)))

  cl <- pam(to_be_clustered_dist, k =  ncluster, medoids = 1:ncluster, do.swap = TRUE, diss = TRUE, cluster.only = TRUE, pamonce = 5)
  plot(all_roi_df$dt_mc, all_roi_df$rt_mc, col = color_sampled[cl], pch =  cl)
  # points(to_be_clustered[1:ncluster,2], to_be_clustered[1:ncluster, 1], col = "black", pch = 4)



  # 6) Include results in the output dataframe
  df_cl <- as.data.frame(cl)
  names(df_cl) <- "roi_cluster"
  all_roi_df <- cbind(df_cl, all_roi_df)

  points(all_roi_df$dt_mc[1:ncluster], all_roi_df$rt_mc[1:ncluster], pch = 4, col ="black")

  # 7)  save results
  setwd(dir_out)
  saveRDS(all_roi_df, file = paste0("all_roi_df", i, ".rds"))
  setwd(dir_in)

  return(all_roi_df)
}


