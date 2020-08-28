
dt_smoothing <- function(data,filtd,NumOfSamples,max_tr){
  # For each sample, the data will be filtered in the drift time. The filtration is made with Savitzky-Golay, using the folowing parametres:
  # Size: 19 points
  # Polynomio order: 2
  # Derivative order: 0

  print(" ")
  print("  /////////////////////////")
  print(" / Drift Time filtration /")
  print("/////////////////////////")
  print(" ")

  setwd(data)
  for (i in (1:NumOfSamples)){
    print(paste0("Sample ", i, " of ", NumOfSamples))
    aux_string <- paste0("M", i, ".mat")
    aux <- readMat(aux_string)[[1]]
    library(signal)
    Mf <- t(apply(aux, 1, function(x) sgolayfilt(x, p = 2, n = 19)))
    setwd(filtd)
    saveRDS(Mf, file = paste0("Mf", i, ".rds"))
    setwd(data)
  }
}
