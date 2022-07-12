#s2 <- read_mea("inst/extdata/sample_formats/small.mea.gz")

# Test my peak detection function
M3 <- readRDS("inst/extdata/M3.rds")

s2 <- GCIMSSample(
  drift_time = M3$data$drift_time,
  retention_time = M3$data$retention_time,
  data = as.matrix(M3$data$data_df)
)

orina <- read_mea("../GCIMS-Samples/2021-4orinas/211214_092633.mea.gz")

dt <- dtime(orina)
rt <- rtime(orina)
intmat <- intensity(orina)

intmatsmooth <- apply(intmat, 2, function(x) signal::sgolayfilt(x, p = 2, n = 21))


rt_idx <- (1+50):2500
dt_idx <- (1175+30+30):2100
zmat <- intmatsmooth[dt_idx,rt_idx]

zmax <- max(zmat)
zmat[zmat > zmax/12] <- zmax/12
zmat <- apply(zmat, 2, function(x) signal::sgolayfilt(x, p = 2, n = 21))
zmat <- t(apply(zmat, 1, function(x) signal::sgolayfilt(x, p = 2, n = 21)))

cosa <- log(zmat-min(zmat) +1)
cosa <- cosa - median(cosa)
rgl::persp3d(
  x = dt[dt_idx],
  y = rt[rt_idx],
  z = cosa,
  col = viridis::viridis(100)[cut(cosa, breaks = 100)],
  xlab = "", ylab = "", zlab = ""
#  xlim = c(8, 12)
)
rgl::rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = c('black'), alpha=1, axes=TRUE)

zmax <- max(zmat)
zmat[zmat > zmax/12] <- zmax/12
zmat <- apply(zmat, 2, function(x) signal::sgolayfilt(x, p = 2, n = 21))
zmat <- t(apply(zmat, 1, function(x) signal::sgolayfilt(x, p = 2, n = 21)))


rt_idx <- (1+400):1250
dt_idx <- (1175+30+30):2100
zmat <- intmat[dt_idx,rt_idx]
#zmat[zmat < 55] <- 55
zmat <- baseline::baseline.als(zmat, lambda = 6, p = 0.001)$corrected

cosa <- (zmat - min(zmat) + 1)^(1/5)
#thres <- median(cosa[800:866,2100:2200])
cosa[cosa < 1.5] <- 1.5

rgl::par3d(windowRect = c(20, 20, 800, 800))
rgl::rgl.viewpoint(theta = 0, phi = 0, fov = 45, zoom = 0.75, type = "userviewpoint")
rgl::persp3d(
  x = dt[dt_idx],
  y = rt[rt_idx],
  z = cosa,
  col = viridis::viridis(100)[cut(log(cosa), breaks = 100)],
  xlab = "", ylab = "", zlab = ""
  #  xlim = c(8, 12)
)
rgl::rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = c('white'), alpha=0, axes=TRUE)
#rgl::play3d(rgl::spin3d(axis = c(0, 0, 1), rpm = 10), duration = 6)
dir.create("/tmp/test4/", recursive = TRUE)
rgl::movie3d(
  rgl::spin3d(axis = c(0, 0, 1), rpm = 5),
  duration = 12,
  fps = 60,
  webshot = FALSE,
  type = "mp4",
  convert = TRUE,
  dir = "/tmp/test4/",
  clean = FALSE
)
# ffmpeg -r 60 -i 'movie%03d.png'  -c:v libx264 -r 30 out_slower.mp4


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
