
dt_baseline <- function(filtd,bslnd,NumOfSamples,lambda,p,k,order,niter){
  # For each sample, the baseline will be removed in the drift time axis using psalsa (by Sergio Oller). The parametres are the define ones at the beginning of the programm:

  print(" ")
  print("  ////////////////////////////////////////////////////////")
  print(" / Eliminaci???n de la l???nea de base en tiempo de deriva /")
  print("////////////////////////////////////////////////////////")
  print(" ")

  setwd(filtd)
  for (i in (1:NumOfSamples)){
    print(paste0("Muestra ", i, " de ", NumOfSamples))
    aux_string <- paste0("Mf", i, ".rds")
    aux <- readRDS(aux_string)
    Mb <- NULL
    setwd(wd)
    Mb <- t(apply(aux, 1, function(x) psalsa(x, lambda, p, k)))
    Mb <- aux - Mb
    setwd(bslnd)
    saveRDS(Mb, file = paste0("Mbd", i, ".rds"))
    setwd(filtd)
  }
}

