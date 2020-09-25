library (signal)
library (R.matlab)
setwd ('C:/Users/Biosignal/Desktop/LCMS_2019/GCIMS/R')
dir_in <- 'C:/Users/Biosignal/Desktop/Files/HAM/JAMONES/data'
dir_out <- 'C:/Users/Biosignal/Desktop/Files/HAM/JAMONES/filtd'
filter_length <- 19
polynomial_order <- 2
by_rows <- TRUE
samples <- c(2)

gcims_smoothing(dir_in, dir_out, filter_length, polynomial_order, by_rows = FALSE, samples)

11
