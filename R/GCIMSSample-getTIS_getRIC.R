setMethod("getTIS", "GCIMSSample", function(object) {
  intmat <- intensity(object)
  tis <- rowSums(intmat)
  tis
})

setMethod("getRIC", "GCIMSSample", function(object) {
  intmat <- intensity(object)
  tis <- rowSums(intmat)
  ric_pos <- which.max(tis)
  ric <- intmat[ric_pos, ]
  ric <- max(ric) - ric
  ric <- ric/sum(ric)
  ric
})


