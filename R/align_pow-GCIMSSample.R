align_pow_sample <- function(x,w){
  int <- GCIMS::intensity(x)
  inter <- t(apply(int, 1, Align::interpolation, w = w, return = FALSE))
  sel <- !is.na(inter[1,])
  x@retention_time <- x@retention_time[sel]
  x@data <- inter[,sel]
  return(x)
}
