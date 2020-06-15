
rt_smoothing <- function(data, p, n){
  print(" ")
  print("  ////////////////////////////")
  print(" / Retention Time Smoothing /")
  print("////////////////////////////")
  print(" ")
  sum_td <- NULL
  for (i in (1:NumOfSamples)){
    setwd(data)
    aux_string <- paste0("M",i,".mat")
    aux <- readMat(aux_string)[[1]]
    sum_td <- cbind(sum_td, apply(aux, 1, sum))
  }
  matplot(sum_td, type="l")
  library(signal)
  aux <- apply(sum_td, 2, function(x) sgolayfilt(x, p = 2, n = 11))
  setwd(smrt)
  saveRDS(aux, file = paste0("Ms", i, ".rds")) 
}