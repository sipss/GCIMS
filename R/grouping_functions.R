
#' Peak/ROI grouping

#' @param dir_in          Input directory. Where input peak picked data files are loaded
#'   from.
#' @param dir_out         Output directory. Where the grouped data (a peak/ROI
#'   table) is stored.
#' @param samples        A vector. Set of samples whose peaks/ROIs need to be
#'   grouped (e.g.: c(1, 2, 3)).
#' @return A peak/ROI table.
#' @details `gcims_peak_grouping` looks for a correspondence among ROI
#'   regions in samples and the ones present in the reference sample.  This
#'   grouping process is performed in two steps: In first step, the algorithm
#'   checks which of the ROI mass centers in a sample are enclosed within any of
#'   the ROI regions of the reference. In the second step, the overlapping
#'   between ROIs of reference and sample pre-selected in step one are computed.
#'   For each ROI of the reference, the ROI of the sample with maximum
#'   overlapping is selected as a representative of this reference ROI. This
#'   way, only one ROI (or none) per sample can be selected as a representative
#'   of one roi of the reference. This process is repeated for all samples in
#'   `samples`.
#' @family Peak Grouping functions
#' @note The grouping process may produce missing values in the final peak/ROI
#'   table. In such a case use `gcims_peak_imputation` to impute them.
#' @export
#' @importFrom pracma meshgrid eye
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


  #---------------#
  #   FUNCTIONS   #
  #---------------#


  #----------------------#
  # check_coarse_overlap #
  #----------------------#

  check_coarse_overlap <- function(x, i, j){
    cond_1 <- (x$rt_mc[j] >= x$min_rt[i]) & (x$rt_mc[j] <= x$max_rt[i])
    cond_2 <- (x$dt_mc[j] >= x$min_dt[i]) & (x$dt_mc[j] <= x$max_dt[i])
    return(cond_1 & cond_2)
  }

  #------------------------#
  # measure_coarse_overlap #
  #------------------------#

  measure_coarse_overlap <- function(x, ncluster_ref, n){
    overlap_coarse <- cbind(eye(ncluster_ref), matrix(0, nrow = ncluster_ref, ncol = n - ncluster_ref + 1))
    for (i in (1:ncluster_ref)){
      for(j in ((ncluster_ref + 1): n)) {
        cond_coarse_overlap <- check_coarse_overlap(x, i, j)
        if(cond_coarse_overlap){
          overlap_coarse[i, j] <- 1
        }
      }
    }
    return(overlap_coarse)
  }

  #-------------------#
  # coarse_clustering #
  #------------------ #

  coarse_clustering <- function(x, y){
    cl <- apply(x, MARGIN = 1, FUN = function(z) which(z == 1))
    cl <- cl[sapply(cl, FUN = function(z) length(z) > 1)] #1
    cl_index <- vector(mode = "list", length = length(cl))
    for (i in seq_along(cl)){
      cl_index[[i]] <- rep(i, length(cl[[i]]))
    }
    cl_index <- unlist(cl_index)
    reliable_rois <- unlist(cl)
    y <- y[reliable_rois, ]
    cl_index <- as.data.frame(cl_index)
    names(cl_index) <- "roi_cluster"
    y <- cbind(cl_index, y)
    return(y)
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

  #----------#
  # grid2ind #
  #----------#

  grid2ind <- function(x, index, rows_df){
    seq_x_index <- x$min_rt[index]:x$max_rt[index]
    seq_y_index <- x$min_dt[index]:x$max_dt[index]
    grid_index_sub <- meshgrid(x = seq_x_index , y = seq_y_index)
    dim(grid_index_sub$X) <- NULL
    dim(grid_index_sub$Y) <- NULL
    grid_index_sub <- cbind(grid_index_sub$X, grid_index_sub$Y)
    grid_index_ind <- apply(grid_index_sub, sub2ind, rows_df, MARGIN = 1)
    return(grid_index_ind)
  }


  #----------------------#
  # compute_grid_overlap #
  #----------------------#

  compute_grid_overlap <- function(grid_x, grid_y){
    grid_xy_union <- union(grid_x, grid_y)
    grid_xy_intersect <- intersect(grid_x, grid_y)
    return(length(grid_xy_intersect) / length(grid_xy_union))
  }


  #----------------------#
  # measure_fine_overlap #
  #----------------------#

  measure_fine_overlap <- function(x, rows_df){
    ncluster  <- length(unique(x$roi_cluster))
    sample_0_indexes <- which(x$sample_id == 0)
    overlap_fine <- vector(mode = "list", length = ncluster)
    for (i in (1:ncluster)){
      current_ind <- sample_0_indexes[i]
      grid_i_ind <- grid2ind(x, current_ind, rows_df)
      tentative_rois <- which(x$roi_cluster == i)
      sample_id <- x$sample_id[tentative_rois]
      h <- 0
      overlap_value <- vector(mode = "list", length = 1)
      for (j in (tentative_rois)){
        h <- h + 1
        grid_j_ind <- grid2ind(x, j, rows_df)
        overlap_value[[1]][[h]] <- compute_grid_overlap(grid_i_ind, grid_j_ind)
      }
      overlap_fine[[i]] <- as.data.frame(t(rbind(tentative_rois, sample_id, overlapping = unlist(overlap_value))))
    }
    return(overlap_fine)
  }

  #--------------------#
  #  which_max_overlap #
  #--------------------#

  which_max_overlap <- function(x) {
    samples_roi <- unique(x$sample_id)
    y <- vector(mode = "list", length = 1)
    h <- 0
    for (k in samples_roi){
      h <- h + 1
      cond <- x$sample_id == samples_roi[h]
      selected_rois <- x$tentative_rois[cond]
      best_hit <- which.max(x$overlapping[cond])
      y[[1]][[h]] <- selected_rois[best_hit]
    }
    x <- x[x$tentative_rois %in% unlist(y), ]
    return(x)
  }

  #-----------------#
  # fine_clustering #
  #-----------------#

  fine_clustering <- function(x, y) {
    z <- do.call(rbind, lapply(y, which_max_overlap))
    x <- x[z$tentative_rois, ]
    return(x)
  }


  #-------------#
  #     MAIN    #
  #-------------#

  # 0) Look for the number of ROIS in the reference sample
  # 1) Initilize roi_list
  # 2) Loop for generating the roi_list
  #    a. read data
  #    b. update roi_list
  # 3) Convert the roi_list to a dataframe
  # 4) Check (coarse) overlapping
  # 5) Perform a coarse clustering
  # 6) Check (fine) overlapping
  # 7) Perform a fine clustering
  # 8) Save results


  setwd(dir_in)
  print("Performing Peak Clustering among samples")


  # 0) Look for the number of ROIS in the reference sample
    aux_string <- paste0("M", 0, ".rds")
    aux_list <- readRDS(aux_string)
    aux_df <- aux_list$data$roi_df
    ncluster_0 <-length(aux_df$roi_id)
    rows_df   <- nrow(aux_list$data$data_df)
    rm(aux_string, aux_list, aux_df)

  # 1) Initilize roi_list
    roi_list <- vector(mode = "list", length = length(samples))

  # 2) Loop for generating the roi_list
    m <- 0
    for (i in samples){
      m <- m + 1
      # 2.a) read data
      aux_string <- paste0("M", samples[m], ".rds")
      aux_list <- readRDS(aux_string)
      # 2.b) update roi_list
      roi_list[[m]]   <- aux_list$data$roi_df
    }

  # 3) Convert the roi_list to a dataframe (all_roi_df)
    all_roi_df <- do.call(rbind, roi_list)
    n_dim  <- nrow(all_roi_df)

  # 4) Check (coarse) overlapping between each roi
  #    in all samples and the rois of the reference sample
    overlap_coarse <- measure_coarse_overlap(all_roi_df, ncluster_0, n_dim)

  # 5) Perform a coarse clustering
    all_roi_df <- coarse_clustering(overlap_coarse, all_roi_df)

  # 6) Check (fine) overlapping between each roi
  #    in all samples and the rois of the reference sample
    overlap_fine <- measure_fine_overlap(all_roi_df, rows_df)

  # 7) Perform a fine clustering
    all_roi_df <- fine_clustering(all_roi_df, overlap_fine)

  # 8) Save results
    setwd(dir_out)
    saveRDS(all_roi_df, file = "all_roi_df.rds")
    setwd(dir_in)
}


