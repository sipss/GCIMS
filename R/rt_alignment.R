library(R.matlab)
rt_alignment <- function(data, aligntr, aux, NumOfSamples, NumOfRows, inicio, fin, max_tr){ 
  print(" ")
  print("  ////////////////////////////")
  print(" / Retention Time Alignment /")
  print("////////////////////////////")
  print(" ")
  
  # The maximes of the previous variable are looked for in the intervals previous defined (inicio and fin), 
  max_row <- NULL
  for (i in (1:NumOfSamples)){
    max_col <- NULL
    for (j in (1:length(inicio))){
      m_x <- which.max(aux[inicio[j]:fin[j],i])
      max_col <- c(max_col, m_x+inicio[j]-1)
    }
    max_row <- rbind(max_row, max_col)
  }
  max_row <- t(max_row)
  new_sum <- aux[1:max_tr,1]
  
  # Once the positions of the reference peaks are found for ech sample, they are aligned taking as patron the position of the reference peaks of the first sample (M1).
  setwd(data)
  M1 <- readMat("M1.mat")[[1]]
  Mr <- M1[1:max_tr, ]
  setwd(aligntr)
  saveRDS(Mr, file = "Mr1.rds")
  tr_sp <- c(1:max_tr)
  
  # The alignment is made interpoling using cubic splines. For each sample and columns of the sample, first the positions of the peaks are adjusted and after, a second interpolation is made in order to make constant the time passes of the sample.
  
  print(paste0("Sample 1 of ", NumOfSamples))
  for (i in (2:NumOfSamples)){
    print(paste0("Sample ", i, " of ", NumOfSamples))
    setwd(data)
    library(ptw)
    tr_sp <- rbind(tr_sp, interp1(max_row[ , 1], max_row[ , i], (1:max_tr), "pchip")) 
    new_sum <- cbind(new_sum, interp1((1:max_tr), aux[1:max_tr, i], tr_sp[i, ]))
    M_string <- paste0("M", i, ".mat")
    M <- readMat(M_string)[[1]]
    Mr <- NULL
    for (j in (1:NumOfColumns)){
      Mr <- cbind(Mr, interp1((1:max_tr), M[1:max_tr, j], tr_sp[i, ]))
    }
    Mr[(which(is.na(Mr) == TRUE))] <- 0 
    setwd(aligntr)
    saveRDS(Mr, file = paste0("Mr", i, ".rds")) 
  }
  saveRDS(new_sum, file = "new_sum.rds")
}
