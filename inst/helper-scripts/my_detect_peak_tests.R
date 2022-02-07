# Plots of the laplacian and the partial derivatives to compare peak detection capabilities


aux_list <- readRDS("inst/extdata/M3.rds")
aux <- as.matrix(aux_list$data$data_df)
drt <- apply(aux, 1, function(x) -signal::sgolayfilt(x, p = 2, n = 21, m = 2))
ddt <- apply(aux, 2, function(x) -signal::sgolayfilt(x, p = 2, n = 11, m = 2))
laplacian <- ddt + t(drt)
noise_level <- 3

patch <- laplacian[100:200,2000:2100]
sigmaNoise <- stats::sd(patch)

total_ion_spectrum <- rowSums(aux) # Sum per rows
rip_position <- which.max(total_ion_spectrum) # Find maximum for every column
rt_idx_with_max_rip <- which.max(aux[rip_position,])
minima <- pracma::findpeaks(-total_ion_spectrum)[, 2] # Find local minima
rip_end_index <- minima[min(which((minima - rip_position) > 0))] # Find ending index of RIP
rip_start_index <- minima[max(which((rip_position - minima) > 0))] # Find starting index of RIP


signal <- aux[rip_start_index:rip_end_index, rt_idx_with_max_rip] # Take RIP
if (length(signal) > 21) {
  filter_length <- 21
} else {
  filter_length <- length(signal) - (length(signal) %% 2 == 0)
}
template <- -signal::sgolayfilt(signal, p = 2, n = filter_length, m = 2) # Compute 2nd derivative of RIP

sigma0 <- GCIMS:::estimate_stddev_peak_width(template)

peakloc_signal <- pracma::findpeaks(aux[,150], minpeakheight = 500, minpeakdistance = 4*sqrt(2)*sigma0)[,2]
peakloc_lapl <- pracma::findpeaks(laplacian[,150], minpeakheight = noise_level*sigmaNoise, minpeakdistance = 4*sqrt(2)*sigma0)[,2]
peakloc_ddt <- pracma::findpeaks(ddt[,150], minpeakheight = 5)[,2]
zeroloc_ddt <- GCIMS:::findZeroCrossings(ddt[,150])
zeroloc_lapl <- GCIMS:::findZeroCrossings(laplacian[,150])

layout(matrix(1:4, ncol=1))
plot(aux_list$data$drift_time, aux[,150], type = "l", xlab = "Drift time (ms)", ylab = "Spectrum at idx 150", main = "inst/extdata/M3.rds")
abline(v = aux_list$data$drift_time[peakloc_signal], col="red")

plot(aux_list$data$drift_time, signal::sgolayfilt(aux[,150], p=2, n = 11, m=2), type = "l", xlab = "Drift time (ms)", ylab = "First part. deriv")


plot(aux_list$data$drift_time, laplacian[,150], type = "l", xlab = "Drift time (ms)", ylab = "Laplacian")
abline(v = aux_list$data$drift_time[peakloc_lapl], col="blue")
abline(v = aux_list$data$drift_time[zeroloc_lapl], col="darkgreen")


plot(aux_list$data$drift_time, ddt[,150], type = "l", xlab = "Drift time (ms)", ylab = "Second part. deriv")
abline(v = aux_list$data$drift_time[peakloc_ddt], col="blue")
abline(v = aux_list$data$drift_time[zeroloc_ddt], col="darkgreen")

