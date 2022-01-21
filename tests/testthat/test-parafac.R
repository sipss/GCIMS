test_that("3-way PARAFAC with random data", {

  skip("toy example, not meant to run")
  mydim <- c(40,15,6) # Dimensions for A, B and C (rows)
  N <- 3 # Dimensions for columns (3-way PARAFAC)

  A <- matrix(stats::rnorm(mydim[1]*N), mydim[1], N) # Generate (mydim, N) normally distributed RVs for A
  B <- matrix(stats::runif(mydim[2]*N), mydim[2], N) # Generate (mydim, N) uniformly distributed RVs for B
  C <- matrix(stats::runif(mydim[3]*N), mydim[3], N) # Generate (mydim, N) uniformly distributed RVs for C

  Xmat <- array(tcrossprod(A, multiway::krprod(C, B)), dim = mydim) # X = sum_f a_f conv b_f conv c_f, where f are the columns of A,B,C (formula in page 152 on the second column here: https://www.cs.cmu.edu/~pmuthuku/mlsp_page/lectures/Parafac.pdf)

  E <- array(rnorm(prod(mydim)), dim = mydim) # Generate normally distributed RVs with the dimensions of X
  E <- multiway::nscale(E, 0, multiway::sumsq(Xmat)) # SNR = 1 (0 -> SNR = 1, 1 -> SNR = 10, 2 -> SNR = 10^2, etc.)
  X <- Xmat + E # X is Xmat plus some error E

  pfac <- multiway::parafac(X, nfac = N) # Fit PARAFAC model

  expect_length(pfac, 13) # Expected list of 13 elements

  Xhat <- stats::fitted(pfac) # Fit solution

  expect_equal(dim(Xhat),dim(X)) # Expected dimensions of Xhat

})
