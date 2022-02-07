#s2 <- read_mea("inst/extdata/sample_formats/small.mea.gz")

# Test my peak detection function
M3 <- readRDS("inst/extdata/M3.rds")

s2 <- GCIMSSample(
  drift_time = M3$data$drift_time,
  retention_time = M3$data$retention_time,
  data = as.matrix(M3$data$data_df)
)


intmat <- intensity(s2)

peaks <- GCIMS:::detect_peaks(
  intmat = intmat,
  noise_level = 3,
  area_overlap_frac = 0.2,
  rip_saturation_threshold = 0.1
)
message(nrow(peaks))
#rgl::par3d(mouseMode = c(none = "none", left = "trackball", right = "zoom", middle = "fov", wheel = "pull"))
#rgl::par3d(mouseMode = c(none = "none", left = "trackball", right = "fov", middle = "fov", wheel = "pull"))

# 3D plot of the intensity matrix:
rgl::persp3d(
  x = dtime(s2),
  y = rtime(s2),
  z = intmat,
  col = grDevices::heat.colors(10)[cut(intmat, breaks = 10)],
  xlim = c(8, 12)
)
# Adding to the previous plot the detected peaks:
for (i in seq_len(nrow(peaks))) {
  rgl::polygon3d(
    x = dtime(s2)[c(peaks$dt_min[i], peaks$dt_min[i], peaks$dt_max[i], peaks$dt_max[i])],
    y = rtime(s2)[c(peaks$rt_min[i], peaks$rt_max[i], peaks$rt_max[i], peaks$rt_min[i])],
    z = rep(1.1*peaks$intensity_max[i], 4),
    coords = 1:2,
    fill = FALSE,
    col = "blue"
  )
}



# Function being written...
dt <- dtime(s2)
rt <- rtime(s2)


peaks <- GCIMS:::detect_peaks(
  intmat = intmat,
  noise_level = 3,
  area_overlap_frac = 0.2,
  rip_saturation_threshold = 0.1
)

plot3d_peaks <- function(intmax, peaks) {
  dt <- seq_len(nrow(intmat))
  rt <- seq_len(ncol(intmat))
  rgl::persp3d(
    x = dt,
    y = rt,
    z = intmat,
    col = grDevices::heat.colors(10)[cut(intmat, breaks = 10)]
  )
  for (i in seq_len(nrow(peaks))) {
    rgl::polygon3d(
      x = dt[c(peaks$dt_min[i], peaks$dt_min[i], peaks$dt_max[i], peaks$dt_max[i])],
      y = rt[c(peaks$rt_min[i], peaks$rt_max[i], peaks$rt_max[i], peaks$rt_min[i])],
      z = rep(1.1*peaks$peak_max_intensity[i], 4),
      coords = 1:2,
      fill = FALSE,
      col = "blue"
    )
  }
}

plot3d_peaks(intmax, peaks)
