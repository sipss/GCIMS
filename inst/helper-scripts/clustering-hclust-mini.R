x <- matrix(c(1,2,3,4.5,6,7,8), ncol = 1)
rownames(x) <- sprintf("P%d", seq_len(nrow(x)))
peak2peak_dist <- dist(x)

# M1:
# P1, P2, P4
peak2peak_dist[1] <- 20
peak2peak_dist[3] <- 20
peak2peak_dist[8] <- 20

# M2:
# P3, P5, P6
# 12
peak2peak_dist[13:14] <- 20
peak2peak_dist[19] <- 20

# M3, P7
peak2peak_dist


clust <- hclust(peak2peak_dist, method = "complete")
plot(clust, ylim = c(0, 10))
clust$merge
as.data.frame(clust$merge)




clust2 <- mdendro::linkage(
  peak2peak_dist,
  method = "complete"
)
plot(clust2, ylim = c(0, 10))


clust2$merger




