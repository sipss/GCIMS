rand_index <- function(clustering1, clustering2) {
  n <- length(clustering1)
  if (length(clustering2) != n) {
    stop("clustering1 and clustering2 must be of the same length")
  }
  # pairs in the same group for clustering1:
  clu1_pairs <- vapply(clustering1, function(x) x == clustering1, FUN.VALUE = logical(n))
  # Ignore all in diagonal and below:
  clu1_pairs <- upper.tri(clu1_pairs) & clu1_pairs

  # for clustering2:
  clu2_pairs <- vapply(clustering2, function(x) x == clustering2, FUN.VALUE = logical(n))
  clu2_pairs <- upper.tri(clu2_pairs) & clu2_pairs

  # a: number of pairs in the same cluster in clu1 that are also in the same cluster for clu2
  a <- sum(clu1_pairs & clu2_pairs)
  # b: number of pairs in the different cluster in clu1 that are also in different cluster for clu2
  clu1_pairs <- upper.tri(clu1_pairs) & (!clu1_pairs)
  clu2_pairs <- upper.tri(clu2_pairs) & (!clu2_pairs)
  b <- sum(clu1_pairs & clu2_pairs)

  pairs <- n*(n-1)/2
  (a+b)/pairs
}
