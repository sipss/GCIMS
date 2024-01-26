#' Warp a GCIMSSample
#'
#' Warps a GCIMSSample with the provided warp
#'
#' @param object A [GCIMSSample] object, modified in-place
#' @param w warp to be applied to the GCIMSSample
#' @return The modified [GCIMSSample]
#' @export

methods::setMethod(
  "align_pow",
  "GCIMSSample",
  function(object,w){
    int <- GCIMS::intensity(object)
    inter <- t(apply(int, 1, pow::interpolation, w = w, return = FALSE))
    sel <- !is.na(inter[1,])
    object@retention_time <- object@retention_time[sel]
    object@data <- inter[,sel]
    return(object)
  }
)
