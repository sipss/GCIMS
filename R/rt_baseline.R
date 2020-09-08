
rt_baseline <- function(filtr,bslnr,NumOfSamples,lambda,p,k,order,niter){
  # For each sample, the baseline will be removed in the drift time axis using psalsa (by Sergio Oller). The parametres are the define ones at the beginning of the programm:

  print(" ")
  print("  //////////////////////////////////////////////////////// /")
  print(" / Eliminación de la línea de base en tiempo de retención /")
  print("//////////////////////////////////////////////////////////")
  print(" ")

  setwd(filtr)
  for (i in (1:NumOfSamples)){
    print(paste0("Muestra ", i, " de ", NumOfSamples))
    aux_string <- paste0("Mfr", i, ".rds")
    aux <- readRDS(aux_string)
    Mb <- NULL
    setwd(wd)
    Mb <- t(apply(aux, 2, function(x) psalsa(x, lambda, p, k)))
    Mb <- aux - Mb
    setwd(bslnr)
    saveRDS(Mb, file = paste0("Mbr", i, ".rds"))
    setwd(filtr)
  }
}
