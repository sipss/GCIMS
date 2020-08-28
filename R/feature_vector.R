
feature_vector <- function(aligntr,featvec,NumOfSamples,tr_interval,td_interval,CalRows,CalColumns){
  print(" ")
  print("  //////////////////////////////////")
  print("  / Creation of the Features Vector/")
  print("  //////////////////////////////////")
  print(" ")

  setwd(aligntd)
  for (i in (1:NumOfSamples)){
    print(paste0("Muestra ", i, " de ", NumOfSamples))
    aux_string <- paste0("Md", i, ".rds")
    aux <- readRDS(aux_string)
    Mv <- as.vector(t(aux[tr_interval, td_interval]))
    setwd(featvec)
    saveRDS(Mv, file = paste0("Mv", i, ".rds"))
    setwd(aligntd)
  }
}
