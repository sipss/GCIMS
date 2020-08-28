
dt_alignment <- function(data,bslnd,alignd,NumOfSamples,NumOfRows,ind_tr,max_tr,dI,td){
  print(" ")
  print("  ////////////////////////")
  print(" / Drift Time Alignment /")
  print("////////////////////////")
  print(" ")

  setwd(data)
  M1 <- readMat("M1.mat")[[1]]
  acc_ind_max <- NULL
  aux <- M1
  for (k in (ind_tr:NumOfRows)){
    max_x <- which.max(aux[k, ])
    acc_ind_max <- t(c(acc_ind_max, max_x))
  }
  max_ref <- round(median(acc_ind_max))
  library(pracma)
  setwd(bslnd)
  for (i in (1:NumOfSamples)){
    print(paste0("Muestra ", i, " de ", NumOfSamples))
    aux_string <- paste0("Mbd", i, ".rds")
    aux <- readRDS(aux_string)
    aux2 <- t(apply(aux, 1, function(x)approx(td, y=NULL, xout = x, rule = 1)$y))
    for (j in (1:dim(aux2)[1])){
      aux2[j, (which(is.na(aux2[j, ]) == TRUE))] <- 0
    }
    setwd(aligntd)
    saveRDS(aux2, file = paste0("Mad", i, ".rds"))
    setwd(bslnd)
  }
}
