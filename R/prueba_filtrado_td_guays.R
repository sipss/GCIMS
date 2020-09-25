library (signal)
library (R.matlab)
wd <- setwd ('C:/Users/Biosignal/Desktop/LCMS_2019/GCIMS/R')
dir_in <- 'C:/Users/Biosignal/Desktop/Files/HAM/JAMONES/data'
dir_out <- 'C:/Users/Biosignal/Desktop/Files/HAM/JAMONES/filtd'
filter_length <- 19
polynomial_order <- 2
by_rows <- TRUE
samples <- c(2)

gcims_smoothing(dir_in, dir_out, samples, by_rows, filter_length, polynomial_order)

setwd ('C:/Users/Biosignal/Desktop/LCMS_2019/GCIMS/R')
dir_in <- 'C:/Users/Biosignal/Desktop/Files/HAM/JAMONES/filtd'
dir_out <- 'C:/Users/Biosignal/Desktop/Files/HAM/JAMONES/bsld'

filter_length <- 19

lambda <- 1E4;
p <- 5E-3;
k <-200;
by_rows <- TRUE
samples <- c(2)


wd <- setwd ('C:/Users/Biosignal/Desktop/LCMS_2019/GCIMS/R')
gcims_baseline_removal(dir_in, dir_out, samples, by_rows,lambda, p, k)
