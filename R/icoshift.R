# #' interval Correlation Optimized shifting
# #'
# #' [xCS,ints,ind,target] = icoshift(xT,xP,inter[,n[,options[,Scal]]])
# #'
# #' Splits a spectral database into "inter" intervals and coshift each vector
# #' left-right to get the maximum correlation toward a reference or toward an
# #' average spectrum in that interval. Missing parts on the edges after
# #' shifting are filled with "closest" value or with "NaNs".
# #'
# #' @param xT        Either a numeric array with the reference spectrum or a character
# #'                  vector as follows:
# #'                    - 'average' if you want to use the average spectrum as a reference
# #'                    - 'median' if you want to use the median spectrum as a reference
# #'                    - 'max' if you want to use for each segment the corresponding actual
# #'                      spectrum having max features as a reference
# #'                    - 'average2' for using the average of the average multiplied for a
# #'                      requested number (default=3) as a reference
# #'
# #' @param xP        Matrix of sample vectors to be aligned as a sample-set
# #'                  towards common reference
# #' @param inter     definition of alignment mode
# #'                  'whole'         : it works on the whole spectra (no intervals).
# #'                  nint            : (numeric) number of many intervals.
# #'                  'ndata'         : (string) length of regular intervals
# #'                                    (remainders attached to the last).
# #'                  [I1s I1e,I2s...]: interval definition. ('I(n)s' interval
# #'                                    n start, 'I(n)e' interval n end).
# #'                  (refs:refe)     : shift the whole spectra according to a
# #'                                    reference signal(s) in the region
# #'                                    refs:refe (in sampling points)
# #'                  'refs-refe'     : shift the whole spectra according to a
# #'                                    reference signal(s) in the region
# #'                                    refs-refe (in Scal units)
# #' @param n         (optional)
# #'                  n = integer n.: maximum shift correction in data
# #'                                  points/Scal units (cf. options(5))
# #'                                  in x/rows. It must be >0
# #'                  n = 'b' (best): the algorithm search for the best n
# #'                                  for each interval (it can be time consuming!)
# #'                  n = 'f' (fast): fast search for the best n for each interval (default)
# #'                  A warning is displayed for each interval if "n" appears too small
# #' @param plot      Do plots
# #' @param fill      'previous' or NA or any number
# #' @param coshift   FALSE: NO Co-shift preprocessing (default)
# #'                  TRUE: Execute a co-shift before icoshift
# #' @param maxshift  max allowed shift for the Co-shift preprocessing (default = equal to n if not specified)
# #' @param units     FALSE: datapoints
# #'                  TRUE: ppm
# #'                  (1) triggers plots & warnings:
# #'                      0 : no on-screen output
# #'                      1 : only warnings (default)
# #'                      2 : warnings and plots
# #'                  (2) selects filling mode
# #'                      0 : using not a number
# #'                      1 : using previous point (default)
# #'                  (3) turns on Co-shift preprocessing
# #'                      0 : no Co-shift preprocessing (default)
# #'                      1 : Executes a Co-shift step before carrying out iCOshift
# #'                  (4) max allowed shift for the Co-shift preprocessing (default = equal to n if not specified)
# #'                      it has to be given in Scal units if option(5)=1
# #'                  (5) 0 : intervals are given in No. of datapoints  (deafult)
# #'                      1 : intervals are given in ppm --> use Scal for inter and n
# #' Scal           : vector of scalars used as axis for plot (optional)
# #' avg_power      : (optional) Multiplier for the 2nd average method
# #'
# #' OUTPUT
# #' xCS  (nP � mT): shift corrected vector or matrix
# #' ints (nI � 4) : defined intervals (Int. No., starting point, ending point, size)
# #' ind  (nP � nI): matrix of indexes reporting how many points each spectrum
# #'                 has been shifted for each interval (+ left, - right)
# #' target (1 x mP): actual target used for the final alignment
# #'
# #' Authors:
# #' Francesco Savorani - Department of Food Science
# #'                      Quality & Technology - Spectroscopy and Chemometrics group
# #'                      Faculty of Sciences
# #'                      University of Copenhagen - Denmark
# #' email: frsa@@life.ku.dk - www.models.life.ku.dk
# #'
# #' Giorgio Tomasi -     Department of Basic Science and Environment
# #'                      Soil and Environmental Chemistry group
# #'                      Faculty of Life Sciences
# #'                      University of Copenhagen - Denmark
# #' email: giorgio.tomasi@@ec.europa.eu - www.igm.life.ku.dk
# #'
# #' @export
# # 170508 (FrSa) first working code
# # 211008 (FrSa) improvements and bugs correction
# # 111108 (Frsa) Splitting into regular intervals (number of intervals or wideness in datapoints) implemented
# # 141108 (GT)   FFT alignment function implemented
# # 171108 (FrSa) options implemented
# # 241108 (FrSa) Automatic search for the best or the fastest n for each interval implemented
# # 261108 (FrSa) Plots improved
# # 021208 (FrSa) 'whole' case added & robustness improved
# # 050309 (GT)   Implentation of interpolation modes (NaN); Cosmetics; Graphics
# # 240309 (GT)   Fixed bug in handling missing values
# # 060709 (FrSa) 'max' target and output 'target' added. Some speed, plot and robustness improvements
# # 241109 (GT)   Interval and band definition in units (added options(5))
# # 021209 (GT)   Minor debugging for the handling of options(5)
# # 151209 (FrSa) Cosmetics and minor debugging for the handling of options(5)
# # 151110 (FrSa) Option 'Max' works now also for alignment towards a reference signal
# # 310311 (FrSa) Bugfix for the 'whole' case when mP < 101
# # 030712 (FrSa) Introducing the 'average2' xT (target) for a better automatic target definition. Plots updated to include also this case
# # 050914 (Sara Rica & Sergio Oller) Modify average2 to allow for non-interactive use of the function
# # 150615 (Sergio Oller) Start port to R
# icoshift <- function(xT, xP, inter = 'whole', n = 'fast', plot=FALSE,
#                      fill='previous', coshift=FALSE, maxshift=n, units=FALSE,
#                      xaxis=FALSE) {
#   if (missing(xT)) {
#     stop("No reference given")
#   }
#   if (missing(xP)) {
#     stop("No spectra to be aligned")
#   }
#   output <- list()
#   return(output)
# }
