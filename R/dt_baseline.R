
dt_baseline <- function(filt,psals,NumOfSamples,lambda,p,k,order,niter){
  # For each sample, the baseline will be removed in the drift time axis using psalsa (by Sergio Oller). The parametres are the define ones at the beginning of the programm:

  print(" ")
  print("  ////////////////////////////////")
  print(" / Baseline Removal in Drift Time/")
  print("//////////////////////////////////")
  print(" ")

  setwd(filt)
  for (i in (1:NumOfSamples)){
    print(paste0("Muestra ", i, " de ", NumOfSamples))
    aux_string <- paste0("Mf", i, ".rds")
    aux <- readRDS(aux_string)
    Mb <- NULL
    setwd(wd)
    Mb <- t(apply(aux, 1, function(x) asysm(x, lambda, p)))


    setwd(psals)
    saveRDS(Mb, file = paste0("Mb", i, ".rds"))
    setwd(filt)
  }
}
