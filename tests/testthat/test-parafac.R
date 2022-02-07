test_that("3-way PARAFAC with random data", {

  skip("Toy example, not meant to run.")
  mydim <- c(40,15,6) # Dimensions for A, B and C (rows)
  N <- 3 # Dimensions for columns (3-way PARAFAC)

  A <- matrix(rnorm(mydim[1]*N), mydim[1], N) # Generate (mydim, N) normally distributed RVs for A
  B <- matrix(runif(mydim[2]*N), mydim[2], N) # Generate (mydim, N) uniformly distributed RVs for B
  C <- matrix(runif(mydim[3]*N), mydim[3], N) # Generate (mydim, N) uniformly distributed RVs for C

  Xmat <- array(Matrix::tcrossprod(A, multiway::krprod(C, B)), dim = mydim) # X = sum_f a_f conv b_f conv c_f, where f are the columns of A,B,C (formula in page 152 on the second column here: https://www.cs.cmu.edu/~pmuthuku/mlsp_page/lectures/Parafac.pdf)

  E <- array(rnorm(prod(mydim)), dim = mydim) # Generate normally distributed RVs with the dimensions of X
  E <- multiway::nscale(E, 0, multiway::sumsq(Xmat)) # SNR = 1 (0 -> SNR = 1, 1 -> SNR = 10, 2 -> SNR = 10^2, etc.)
  X <- Xmat + E # X is Xmat plus some error E

  pfac <- multiway::parafac(X, nfac = N) # Fit PARAFAC model

  expect_length(pfac, 13) # Expected list of 13 elements

  Xhat <- stats::fitted(pfac) # Fit solution

  expect_equal(dim(Xhat),dim(X)) # Expected dimensions of Xhat

})

test_that("3-way PARAFAC2 for concentration prediction",{
  skip("Toy example, not meant to run.")

  #### IMPORT ####

  samples = 1:15 # Samples
  f_table <- readMat('table.mat') # Import table
  table <- f_table$table.all # Table in file
  concentrations <- list('2','2','2','5','5','5','10','10','10','15','15','15','20','20','20') # Concentrations (list)
  concentrations_num <- c(rep(2,3),rep(5,3),rep(10,3),rep(15,3),rep(20,3)) # Concentrations (integers)

  # Compute the median of the sizes
  median_width <- round(median(table[,6] - table[,5]))
  median_height <- round(median(table[,8] - table[,7]))

  # Initialize X with zeros
  nROIs <- length(unique(table[,2]))
  X <- array(rep(0,median_width*median_height*length(samples)*nROIs), dim=c(median_width,median_height,length(samples),length(unique(table[,2]))))

  # For loop that iterates through all the samples
  for (s in samples){
    f <- readMat(paste0(s, ".mat")) # Import file
    data <- f$data.df # Extract data from file
    for (r in 1:nROIs){
      idx <- which(table[,1] == s & table[,2] == r) # Position where ROI r is in sample s
      center <- c(table[idx,3],table[idx,4]) # Pick center
      patch <- data[seq(center[1] - ceiling(median_width/2-1), center[1] + floor(median_width/2)),seq(center[2] - ceiling(median_height/2-1),center[2] + floor(median_height/2))] # Extract patch with intensities
      X[,,s,r] <- patch # Put patch inside X
    }
  }

  ### CONCENTRATION PREDICTION FOR ROI 1 ###

  colors_lgnd <- list("red","blue","green","purple","brown") # Colors for plotting
  ROI <- 1 # Concentration prediction for ROI1

  conditions <- c() # Initialize variable to store all the conditions
  nComponents <- 1:4 # Number of components
  layout(matrix(nComponents, ncol=2, byrow = TRUE)) # Layout of plots

  RMSE <- c() # Initialize array to save the RMSE

  # Loop that iterates through all the components
  for (n in nComponents){

    # Run PARAFAC2
    pfac2<- multiway::parafac2(X[,,,1], nfac = n, const=c("ortnon","ortnon","ortnon")) # Fit PARAFAC model for ROI 2 (monomer)

    lgnd <- c() # Initialize legend

    for (i in 1:dim(pfac2$C)[2]){
      x <- as.numeric(unlist(pfac2$C[,i])) # x-axis (C from PARAFAC2)

      fitting <- lm(concentrations_num ~ x + I(x^2)) # Parabola fitting (second-order quadratic regression)

      y_axis <- fitting$coefficients[1] + fitting$coefficients[2]*x + fitting$coefficients[3]*x^2 # Predicted values based on regression

      conditions <- c(conditions,all(abs(y_axis - concentrations_num) < 3)) # Check if the absolute error is less than 3

      # Subplot of predicted values (datapoints)
      if (i == 1){
        plot(concentrations_num, y_axis, type = "p", pch = 19, col = unlist(colors_lgnd[i]), ylab = "predicted concentration (ppb)", xlab = "real concentration (ppb)", ylim=c(0,20), xlim=c(0,20))
      }else{
        lines(concentrations_num, y_axis, type = "p", pch = 19, col = unlist(colors_lgnd[i]), ylim=c(0,20), xlim=c(0,20)) #
      }
      abline(a=0,b=1, col="black", lty=2) # Plot ideal prediction (dotted line)

      RMSE <- c(RMSE,round(sqrt(sum((y_axis - concentrations_num)^2) /length(y_axis)),2))
      lgnd <- c(lgnd,paste0("RMSE = ",RMSE[i])) # Add legend
    }
    title(paste0('Number of components: ', n)) # Title of subplot
    legend("bottomright", cex=0.65, legend=lgnd, pch=15, col=unlist(colors_lgnd[1:dim(pfac2$C)[2]])) # Legend
  }

  expect_equal(any(conditions), TRUE)

})
