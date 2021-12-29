#' Peak alignment using Correlation Optimized Warping (COW)


#' @param dir_in          Input directory. Where input data files are
#'  loaded from.
#' @param dir_out         Output directory. Where data files after alignment
#'   correction are stored.
#' @param samples         A vector. Set of samples to be aligned(e.g.: c(1, 2,
#'   3)).
#' @param time            Sets the dimension to be corrected: drift time or
#'   retention time. Introduce "Retention" for retention time; or "Drift" for
#'   drift time.
#' @param seg_vector      Vector of segment lengths.
#' @param slack_vector    Vector of slacks.
#' @details `gcims_baseline_removal` aligns sample chromatograms (in
#'   retention time) or spectra (in drift time).
#'   The alignment in drift time is achieved by transforming the selected temporal
#'   axis of each sample according to the reference using the Correlation Optimized
#'   Warping (COW).
#'
#'   In COW, the reference signal of length *N* is divided
#'   in *L* segments of length *I* and so it is done for the signal to
#'   be aligned. When the correction is applied along retention time direction,
#'   the signals used to estimate COW model are Total Ion Chromatograms (TIC).
#'   In case of correcting drift time misalignments, COW model is estimated
#'   using Total Ion Spectra (TIS).
#'
#'   The alignment between reference and sample
#'   representatives is performed segment by segment, allowing slight
#'   modifications on the segment lengths of a sample. Segment lengths are
#'   controlled by the so-called slack parameter, *t*. For each of the
#'   segments of a sample, segment length is varied so that the correlation
#'   between sample and reference is maximized. At this point a linear
#'   interpolation is realized to recover sample segments of the original length
#'   and sample period. Once the COW model is known, it is applied to the rest
#'   of chromatrograms (if the correction is retention time) or spectra (if the
#'   correction is in drift time) of a sample. This is done for all samples
#'   selected by the argument `samples`.
#'
#'   Cow method is optimized by selecting the parameters that maximize
#'   the average correlation between the reference and the rest of the samples.
#'   In particular, the arguments `seg_vector` and `slack_vector`
#'   provide, respectively, the sets of possible values of *L* and *t*
#'   for COW parameter optimization.
#'
#'   The alignment in retention time is achieved applying parametric time warping.
#'   The arguments `seg_vector` and `slack_vector` are not necessary for
#'   the retention time alignment.
#' @return A set of S3 objects.
#' @export
#'
#' @references {
#'  Niels-Peter Vest Nielsen, Jens Micheal Carstensen and Jørn Smedegaard 'Aligning of singel and multiple
#'           wavelength chromatographic profiles for chemometric data analysis using correlation optimised warping'
#'           J. Chrom. A 805(1998)17-35
#'
#'  Correlation optimized warping and dynamic time warping as preprocessing methods for chromatographic Data
#'            Giorgio Tomasi, Frans van den Berg and Claus Andersson, Journal of Chemometrics 18(2004)231-241
#'  }
#'
#' @author Sergio Oller-Moreno,  \email{soller@@ibecbarcelona.eu}
#'
#' Inspired on the work by:
#'
#' @author Giorgio Tomasi, \email{gt@@kvl.dk}
#' @author Frans van den Berg 070821 (GT) \email{fb@@kvl.dk}
#'
#' Royal Agricultural and Veterinary University - Department of Food Science
#' Quality and Technology - Spectroscopy and Chemometrics group - Denmark
#' www.models.kvl.dk
#' @importFrom signal interp1
#'
#' @examples
#' \donttest{
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- c(3, 7)
#'
#' # Example of Drift time Alignment:
#' # Before:
#' gcims_plot_spec(dir_in, samples, rt_value = NULL, dt_range = NULL, colorby = "Class")
#'
#' # After:
#' time <- "Drift"
#' seg_vector <- seq(from = 10, to = 30, by = 5)
#' slack_vector <- seq(from = 1, to = (seg_vector[length(seg_vector)] - 4), by = 3)
#' gcims_alignment(dir_in, dir_out, samples, time, seg_vector, slack_vector)
#' gcims_plot_spec(dir_out, samples, rt_value = NULL, dt_range = NULL, colorby = "Class")
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#'}
#'
gcims_alignment <- function(dir_in, dir_out, samples, time, seg_vector, slack_vector){


  print(" ")
  print("  /////////////////////////")
  print(" /      Peak Alignment   /")
  print("/////////////////////////")
  print(" ")


  setwd(dir_in)
  if (time == "Drift"){
    Warping <- optimize_cow(dir_in, dir_out, samples, time, seg_vector, slack_vector)
  }

  m <- -1
  for (i in c(0,samples)){
    m <- m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }

    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new


    aux <- as.matrix(aux_list$data$data_df)

    if (i != 0){
      if (time == "Drift"){
        aux <- t(aux)
        aux <- apply_cow(aux, Warping[m, , ])
      } else if (time == "Retention"){
        aux_string <- paste0("M", 0, ".rds")
        aux_list <- readRDS(aux_string) #new
        aux <- as.matrix(aux_list$data$data_df)
        logindx <- log(c(1:dim(aux)[2]))
        reference <- apply(aux, 1, function(x) round(interp1(logindx, x, xi = logindx)))
        aux_string <- paste0("M", i, ".rds")
        aux_list <- readRDS(aux_string) #new
        aux <- as.matrix(aux_list$data$data_df)
        logindx <- log(c(1:dim(aux)[2]))
        aux <- apply(aux, 1, function(x) round(interp1(logindx, x, xi = logindx)))
        num_dt <- nrow(aux)
        num_ret_time <- ncol(aux)
        timepad <- c(rep(NA, 500), aux_list$data$retention_time, rep(NA, 500))
        aux_with_zero_pad <- matrix(0, nrow = 1000 + num_dt, ncol =  num_ret_time)
        aux_with_zero_pad[501:(nrow(aux_with_zero_pad) - 500), ] <- aux
        aux_to_align <- aux_with_zero_pad
        global_warp <- ptw::ptw(ref = t(reference), samp = t(aux_to_align), warp.type = "global", init.coef = c(0, 1), mode = "backward")
        aux_global_warp <- global_warp$warped.sample
        aux_global_warp[is.na(aux_global_warp)] <- 0
        indiv_warp <- ptw::ptw(ref = t(reference), samp = aux_global_warp, warp.type = "individual", init.coef = c(0,1,0), mode = "backward")
        aux_indiv_warp <- indiv_warp$warped.sample
        aux_indiv_warp[is.na(aux_indiv_warp)] <- 0
        aux <- aux_indiv_warp
      }
    } else {
    }

    if (i != 0){
      if (time == "Drift"){
        aux <- t(aux)
      } else if (time == "Retention"){
      }
    }
    aux_list$data$data_df <- round(aux)
    M <- aux_list


    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}


optimize_cow <- function(dir_in, dir_out, samples, time, seg_vector, slack_vector){
  setwd(dir_in)

  # Load a sample to know its length in the
  # direction of alignment (transpose if needed).
  # Create a matrix for storing the samples to
  # be aligned (curves).
  print(" ")
  print("Optimizing COW parameters, please wait")
  print(" ")


  aux_string <- paste0("M", samples[1], ".rds")
  aux_list <- readRDS(aux_string) #new
  aux <- as.matrix(aux_list$data$data_df)


  if (time == "Drift"){
    aux <- t(aux)
  } else if (time == "Retention"){
  }

  curves <- matrix(0, dim(aux)[2], length(samples) + 1)


  # Load the samples (transpose if needed) and
  # compress them (just an addition) in the direction
  # of alignment. Store them in the variable curves

  m <- 0
  for (i in c(0, samples)){
    m <- m + 1
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- as.matrix(aux_list$data$data_df)


    if (time == "Drift"){
      aux <- t(aux)
    } else if (time == "Retention"){
    }
    curves[,m] <- colSums(aux)
  }
  rm(m, aux, aux_string)

  # Separate the reference curve from the rest.
  # The reference curce must be the first one.

  ref_curve <- curves[, 1]



  curves <- curves[ , 2:dim(curves)[2]]


  # Define the correlation matrix.
  corr_data <- matrix(0,length(seg_vector),
                      length(slack_vector)) #ESTO ESTÁ MAL: ARREGLAR

  # Fill the correlation Matrix. For each seg and slack compute
  # the correlation between the reference curve and the rest of
  # of the curves. Them compute the mean value of correlation.

  n <- 0
  diff_seg <- diff(seg_vector)  #ESTO ES MEJORABLE: PENSARLO MEJOR
  step_y <- diff_seg[1]

  for (Z in (seg_vector)){
    n <- n + 1
    m <- 0
    for (Y in (seq(from = 1, to = (Z - 4), by = step_y))){
      m <- m + 1
      cow_results <- cow(ref_curve, t(curves), Z, Y)


      align_curves <- t(cow_results$XWarped)

      sum_corr <- 0
      for (X in (1:dim(align_curves)[2])){
        sum_corr <- sum_corr + stats::cor(ref_curve, align_curves[, X])
      }
      corr_data[n, m] <- sum_corr / dim(align_curves)[2]
    }

  }
  rm(cow_results, diff_seg, step_y, n, m, Z, X, Y)

  # Find the indexes for the maximal value of the average correlation.
  # The find the best seg and slack.

  opt_indexes <- which(corr_data == max(corr_data), arr.ind = TRUE)
  seg <- seg_vector[opt_indexes[1]]
  slack <-slack_vector[opt_indexes[2]]

  # Use cow to obtain the best Warping

  cow_results <- cow(ref_curve, t(curves), seg, slack)
  rm(seg_vector, slack_vector, opt_indexes, seg, slack)
  Warping <- cow_results$Warping


  return(Warping)
}


cow <- function(T, X, Seg, Slack,
                plot=FALSE, corrpow=1, equal_lengths=FALSE,
                max_correction=NA, debug=FALSE) {
  plot = as.logical(plot)
  if (corrpow < 1 || corrpow > 4) {
    stop("corrpow must be in the range 1:4")
  }
  if (any(is.na(T)) | any(is.na(X))) {
    stop('cow can not handle missing values')
  }
  if (!is.matrix(X)) {
    X <- matrix(X, nrow=1)
  }
  nX <- nrow(X)
  pX <- ncol(X)
  pT <- length(T)
  XWarped <- matrix(0, nrow=nX, ncol=pT)
  Time = 0
  Seg = round(Seg)
  Pred_Bound <- length(Seg) > 1
  LenSeg <- NA
  nSeg <- NA
  if (Pred_Bound) {
    assertthat::assert_that(all(Seg[,1] == 1)) #, msg="Segments must start at 1")
    assertthat::assert_that(all(Seg[,ncol(Seg)] == c(pT, pX)))#, msg="Segments must end at the length of the pattern/target")

    LenSeg = t(diff(t(Seg))) # Length of the segments in the - 1

    assertthat::assert_that(all(LenSeg >= 2))#, msg = "Segments must contain at least two points")

    nSeg <- ncol(LenSeg)
  } else {
    assertthat::assert_that(Seg <= min(pX, pT))#, msg='Segment length is larger than length of the signal')

    if (equal_lengths) { # Segments in the signals can have different length from those in the target
      nSeg             = floor((pT - 1)/Seg)
      LenSeg <- matrix(NA, nrow=2, ncol=nSeg)
      LenSeg[1,1:nSeg] = floor((pT - 1)/nSeg)
      LenSeg[2,1:nSeg] = floor((pX - 1)/nSeg)
      message('Segment length adjusted to best cover the remainders')
    } else {
      nSeg = floor((pT - 1) / (Seg - 1))
      LenSeg <- matrix(NA, nrow=2, ncol=nSeg)
      LenSeg[1:2,1:nSeg] = Seg - 1;
      if (floor((pX - 1) / (Seg - 1)) != nSeg) {
        stop('For non-fixed segment lengths the target and the signal do not have the same number of segments (try equal_lengths)')
      }
    }
    temp <- (pT-1) %% LenSeg[1,1] # The remainders are attached to the last segment in the target and in the reference
    if (temp > 0) {
      LenSeg[1,nSeg] = LenSeg[1,nSeg] + temp
      #message(sprintf('Segments: %d points x %d segments + %d (target)',LenSeg[1,1] + 1, nSeg - 1, LenSeg[1,ncol(LenSeg)] + 1))
    } else {
      #message(sprintf('Segments: %d points x %d segments (target)',LenSeg[2,1] + 1, nSeg))
    }
    temp <- (pX-1) %% LenSeg[2,1]
    if (temp > 0) {
      LenSeg[2, nSeg] <- LenSeg[2, nSeg] + temp
      #message(sprintf('       : %d points x %d segments + %d (signals)',LenSeg[2,1] + 1, nSeg - 1, LenSeg[2,ncol(LenSeg)] + 1))
    } else {
      #message(sprintf('       : %d points x %d segments (signals)',LenSeg[2,1] + 1, nSeg))
    }
  }

  if (any(LenSeg <= Slack +2)) {
    stop('The slack cannot be larger than the length of the segments')
  }

  bT      = cumsum(c(1,LenSeg[1,]))
  bP      = cumsum(c(1,LenSeg[2,]))
  Warping = array(data=0, dim=c(nX, nSeg+1,2))
  #
  # Check slack
  if (length(Slack) > 1) {  # Different slacks for the segment boundaries will be implemented
    if (ncol(Slack) <= nSeg) {
      stop('The number of slack parameters is not equal to the number of optimised segments')
    }
    stop('Multiple slacks have not been implemented yet')
  }
  Slacks_vec  = seq(from=-Slack, to=Slack)  # All possible slacks for a segment boundary
  # Set feasible points for boundaries
  Bounds      = matrix(1, nrow=2, ncol=nSeg+1)
  # Slope Constraints
  offs <- matrix(c(-Slack*seq(0, nSeg),
                   +Slack*seq(0, nSeg)),
                 byrow=TRUE, nrow=2)
  Bounds_a <- matrix(rep(bP[1:(nSeg + 1)], each=2), nrow=2)  + offs
  Bounds_b <- matrix(rep(bP[1:(nSeg + 1)], each=2), nrow=2) + offs[,seq(from=nSeg+1, to=1, by=-1)]
  Bounds[1,] <- pmax(Bounds_a[1,], Bounds_b[1,])
  Bounds[2,] <- pmin(Bounds_a[2,], Bounds_b[2,])

  # % Band Constraints
  if (!is.na(max_correction)) {
    if (abs(pT-pX) > max_correction) {
      stop('The band is too narrow and proper correction is not possible')
    }
    Bounds[1,] <- pmax(Bounds[1,], pmax(0, pX/pT*bT-max_correction))
    Bounds[2,] <- pmin(Bounds[2,], pmin(pX, pX/pT*bT + max_correction))
    if (any(diff(Bounds < 0))) {
      stop('The band is incompatible with the fixed boundaries')
    }
  }

  # Calculate first derivatives for interpolation
  Xdiff = t(diff(t(X)))
  #
  # %% Calculate coefficients and indexes for interpolation
  Int_Coeff <- vector(mode = "list", length = nSeg)
  Int_Index <- vector(mode = "list", length = nSeg)
  if (!Pred_Bound) {
    tmp <- InterpCoeff(LenSeg[1,1] + 1, LenSeg[2,1] + Slacks_vec + 1, Slacks_vec)
    A <- tmp$Coeff
    B <- tmp$Index
    Int_Coeff[1:(nSeg-1)] <- lapply(as.list(1:(nSeg-1)), function(x) return(A))
    Int_Index[1:(nSeg-1)] <- lapply(as.list(1:(nSeg-1)), function(x) return(B))

    tmp <- InterpCoeff(LenSeg[1,nSeg] + 1, LenSeg[2,nSeg] + Slacks_vec + 1, Slacks_vec)
    A <- tmp$Coeff
    B <- tmp$Index
    Int_Coeff[[nSeg]] <- A
    Int_Index[[nSeg]] <- B
  } else {
    for (i_seg in 1:nSeg) {
      tmp <- InterpCoeff(LenSeg[1, i_seg] + 1, LenSeg[2,i_seg] + Slacks_vec + 1, Slacks_vec)
      A <- tmp$Coeff
      B <- tmp$Index
      Int_Coeff[[i_seg]] <- A
      Int_Index[[i_seg]] <- B
    }
  }
  #
  # Dynamic Programming Section
  Table_Index    = cumsum(c(0,diff(Bounds) + 1)) # Indexes for the first node (boundary point) of each segment in Table
  Table <- array(data = 0, dim = c(3, Table_Index[nSeg + 2], nX))
  # Table: each column refer to a node
  # (1,i) position of the boundary point in the signal
  # (2,i) optimal value of the loss function up to node (i)
  # (3,i) pointer to optimal preceding node (in Table)

  Table[2, 2:ncol(Table),] <- -Inf # All loss function values apart from node (1) are set to -Inf
  for (i_seg in 1:(nSeg +1)) { # Initialise Table
    v <-  matrix(Bounds[1,i_seg]:Bounds[2,i_seg], ncol=1)
    Table[1, (Table_Index[i_seg] + 1):Table_Index[i_seg + 1],] = v[,rep(1,nX)]
  }
  # Forward phase
  for (i_seg in 1:nSeg) { # Loop over segments
    a = Slacks_vec + LenSeg[2,i_seg]   # a,b,c: auxiliary values that depend only on segment number and not node
    b = Table_Index[i_seg] + 1 - Bounds[1,i_seg]
    c = LenSeg[1,i_seg] + 1
    Count = 1 # Counter for local table for segment i_seg
    Node_Z = Table_Index[i_seg + 2]   # Last node for segment i_seg
    Node_A = Table_Index[i_seg + 1] + 1 #First node for segment i_seg
    Bound_k_Table = array(data = 0, dim=c(2,Node_Z - Node_A + 1,nX)) # Initialise local table for boundary

    Int_Index_Seg = t(Int_Index[[i_seg]]) - (LenSeg[2,i_seg] + 1)  # Indexes for interpolation of segment i_seg
    Int_Coeff_Seg = t(Int_Coeff[[i_seg]])                          # Coefficients for interpolation of segment i_seg

    TSeg          = T[bT[i_seg]:bT[i_seg + 1]] # Segment i_seg of target T
    TSeg_centred  = TSeg - sum(TSeg)/length(TSeg) # Centred TSeg (for correlation coefficients)
    Norm_TSeg_cen = pracma::Norm(TSeg_centred) # (n - 1) * standard deviation of TSeg
    #
    for (i_node in Node_A:Node_Z) { # Loop over nodes (i.e. possible boundary positions) for segment i_seg
      Prec_Nodes <- Table[1, i_node, 1] - a
      Allowed_Arcs <- Prec_Nodes >= Bounds[1,i_seg] & Prec_Nodes <= Bounds[2,i_seg] # Arcs allowed by local and global constraints
      Nodes_TablePointer <- b + Prec_Nodes[Allowed_Arcs] # Pointer to predecessors in Table
      N_AA               = sum(Allowed_Arcs) # Number of allowed arcs
      if (N_AA > 0) { # Sometimes boundaries are ineffective and few nodes are allowed that cannot be reached
        # It has to be further investigated
        Index_Node = Table[1,i_node, 1] + Int_Index_Seg[,Allowed_Arcs]  # Interpolation signal indexes for all the allowed arcs for node i_node
        Coeff_b = Int_Coeff_Seg[,Allowed_Arcs] # Interpolation coefficients for all the allowed arcs for node i_node
        Coeff_b = matrix(Coeff_b, nrow=1)
        Coeff_b   = Coeff_b[rep(1, nX),]
        Xi_Seg = X[,Index_Node]
        Xi_diff = Xdiff[,Index_Node]
        Xi_Seg  = pracma::Reshape(t(Xi_Seg + Coeff_b * Xi_diff),c, N_AA * nX);  # Interpolate for all allowed predecessors
        Xi_Seg_mean = matrix(apply(Xi_Seg, 2, sum)/nrow(Xi_Seg), nrow=1) # Means of the interpolated segments
        Norm_Xi_Seg_cen <- sqrt(apply(Xi_Seg^2, 2, sum) - nrow(Xi_Seg) * Xi_Seg_mean^2)       # Fast method for calculating the covariance of T and X (no centering of X is needed)
        CCs_Node <- (TSeg_centred %*% Xi_Seg)/(Norm_TSeg_cen %*% Norm_Xi_Seg_cen)     # Correlation coefficients relative to all possible predecessors
        CCs_Node[!is.finite(CCs_Node)] = 0                                                              # If standard deviation is zero, update is not chosen
        CCs_Node = pracma::Reshape(CCs_Node,N_AA,nX)
        if (corrpow == 1) {
          Cost_Fun = pracma::Reshape(Table[2,Nodes_TablePointer,],N_AA,nX) + CCs_Node   # Optimal value of loss function from all predecessors
        } else {
          Cost_Fun = pracma::Reshape(Table[2,Nodes_TablePointer,],N_AA,nX) + CCs_Node^corrpow
        }
        pos = apply(Cost_Fun, 2, which.max)
        ind = apply(Cost_Fun, 2, max)   # Optimal value of loss function from all predecessors
        Bound_k_Table[1, Count,] <- ind
        Bound_k_Table[2, Count,] <- Nodes_TablePointer[pos]
        Count <- Count + 1
      }
    }
    Table[2:3,Node_A:Node_Z,] = Bound_k_Table # Update general table (it turned out to be faster than using Table directly in the loop over nodes
  }
  #
  for (i_sam in 1:nX) {  # Loop over samples/signals
  # Backward phase
    Pointer = ncol(Table) # Backtrace optimal boundaries using the pointers in Table
    Warping[i_sam, nSeg + 1, 1] = pX
    for (i_bound in seq(from=nSeg, to=1, by=-1)) {
      Pointer <- Table[3, Pointer, i_sam]
      Warping[i_sam,i_bound, 1] = Table[1, Pointer, i_sam]
    }
  }

  Warping[,,2] = matrix(bT, nrow=1)[rep(1, nX),]
  #
  output <- list()
  # Output
  # Reconstruct aligned signals
  for (i_seg in 1:nSeg) {
    indT = bT[i_seg]:bT[i_seg + 1]
    lenT = bT[i_seg + 1] - bT[i_seg]
    for (i_sam in 1:nX) {
      indX = Warping[i_sam,i_seg,1]:Warping[i_sam,i_seg + 1,1]
      lenX = Warping[i_sam,i_seg + 1,1] - Warping[i_sam,i_seg,1]
      # NB the right handside expression must be transposed to fit MATLAB version 6.5
      x = indX - Warping[i_sam, i_seg,1] + 1
      y = X[i_sam,indX]
      xi <- (0:lenT)/lenT * lenX + 1
      #print(xi)
      #print(x)
      #print(y)
      XWarped[i_sam, indT] = t(signal::interp1(x, y, xi))
    }
  }
  output$Warping <- Warping
  output$XWarped <- XWarped
  if (debug) {
    output$Diagnos <- list(indexP=bP,
                           indexT=bT,
                           Nsegments=nSeg,
                           options=list(),
                           rangeP=t(Bounds),
                           segment_length=LenSeg,
                           slack=Slack,
                           table=Table)
  }

  #
  # %% Plot
  # if Options(1)
  #    figure
  #    minmaxaxis = [1 max([pT pX]) min([T X(nX,:)]) max([T X(nX,:)])] ;
  #    subplot(2,1,1);
  #    plot(1:pT,T,'b',bT,T(bT),'.b',1:pX,X(nX,:),'g',bP,X(nX,bP),'.g');
  #    hold on
  #    for a = 2:length(Warping(nX,:,1))
  #       plot([bT(a) Warping(nX,a,1)],[T(Warping(nX,a,2)) T(Warping(nX,a,2))],'r');
  #       if (Warping(nX,a,2) > Warping(nX,a,1))
  #          plot(Warping(nX,a,2),T(Warping(nX,a,2)),'>r');
  #       else
  #          plot(Warping(nX,a,2),T(Warping(nX,a,2)),'<r');
  #       end
  #    end
  #    hold off
  #    axis(minmaxaxis)
  #    grid
  #    title(['COW reference = blue, Sample ' num2str(nX) '(/' num2str(nX) ') = green, Segment-boundary movement = red']);
  #    subplot(2,1,2);
  #    plot(1:pT,T,'b',1:pT,XWarped(nX,:),'g');
  #    grid;
  #    axis(minmaxaxis);
  #    title('Warped sample')
  # end
  return(output)
}




apply_cow <- function(X, Wrp) {

  nX <- dim(X)[1]
  pX <- dim(X)[2]

  Xw <- matrix(0, nrow=nX, ncol=pX)
  pW <-dim(Wrp)[1] # dim(Warping)[2]

  for (i_seg in 1:(pW-1)) {

    indT <- Wrp[i_seg, 2]: Wrp[i_seg + 1, 2]
    lenT <- Wrp[i_seg + 1, 2]- Wrp[i_seg, 2]

    for (i_sam in 1:nX) {

      indX = Wrp[i_seg, 1]: Wrp[i_seg + 1, 1]
      lenX = Wrp[i_seg + 1, 1] - Wrp[i_seg, 1]
      x = indX - Wrp[i_seg,1 ] + 1
      y = X[i_sam, indX]
      xi <- (0:lenT)/lenT * lenX + 1
      Xw[i_sam, indT] = t(signal::interp1(x, y, xi))

    }

  }
  return(Xw)

}


InterpCoeff <- function(n,nprime,offs) {
  p     = length(nprime)
  q     = n - 1
  Coeff = matrix(0, nrow = p, ncol = n)
  Index = matrix(0, nrow = p, ncol = n)
  for (i_p in 1:p) {
    pp = 1:nprime[i_p]
    p = (0:q) * (nprime[i_p] - 1)/q + 1;
    k <- pracma::histc(p,pp)$bin
    k[p < 1] = 1
    k[p >= nprime[i_p]] = nprime[i_p] - 1
    Coeff[i_p,] <- p - pp[k]
    Index[i_p,] <- k - offs[i_p]
  }
  return(list(Coeff = Coeff, Index = Index))
}
