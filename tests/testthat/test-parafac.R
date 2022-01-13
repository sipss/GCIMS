
mydim <- c(40,15,6) # Dimensions for A, B and C (rows)
N <- 3 # Dimensions for columns (3-way PARAFAC)

A <- matrix(rnorm(mydim[1]*N), mydim[1], N) # Generate (mydim, N) normally distributed RVs for A
B <- matrix(runif(mydim[2]*N), mydim[2], N) # Generate (mydim, N) uniformly distributed RVs for B
C <- matrix(runif(mydim[3]*N), mydim[3], N) # Generate (mydim, N) uniformly distributed RVs for C

Xmat <- array(tcrossprod(A, krprod(C, B)), dim = mydim) # X = sum_f a_f conv b_f conv c_f, where f are the columns of A,B,C (formula in page 152 on the second column here: https://www.cs.cmu.edu/~pmuthuku/mlsp_page/lectures/Parafac.pdf)

E <- array(rnorm(prod(mydim)), dim = mydim) # Generate normally distributed RVs with the dimensions of X
E <- nscale(E, 0, sumsq(Xmat)) # SNR = 1 (0 -> SNR = 1, 1 -> SNR = 10, 2 -> SNR = 10^2, etc.)
X <- Xmat + E # X is Xmat plus some error E

pfac <- parafac(X, nfac = N) # Fit PARAFAC model
Xhat <- fitted(pfac) # Fit solution
SSE <- sum((Xmat-Xhat)^2)/prod(mydim) # SSE
SSE
