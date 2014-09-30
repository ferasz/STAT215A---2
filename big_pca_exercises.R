
library("Matrix")
library("irlba")
library("microbenchmark")

GetCov <- function(p, m, max.corr = .5, sparse = T) {
  # Generate a covariance matrix with limited off-diagonal elements.
  #
  # Args:
  #   p: dimensionality of the cov mat
  #   m: number non-zero elements in each direction off the diagonal
  #   max.corr: maximum correlation between variables
  #   sparse: whether to use sparse data structure (Matrix vs matrix)
  #
  # Returns:
  #   A matrix with nonzeros close to the diagonal and zeros everywhere else
  #   Each row will look like 
  #       0 0 0 0 0 .1 .2 ... .9 1 .9  ... .2 .1 0 0 0 0 0

  r <- seq(max.corr, 0, length.out=m + 1)
  r <- r[ -length(r)]
  if (sparse) {
    mat <- Matrix(0, nrow = p, ncol = p, sparse = T)
  } else {
    mat <- matrix(0, nrow = p,ncol = p)
  }
  
  for (i in 1:length(r)) {
    mat[seq(from = i+1, by = p+1, length.out = p-i )] <- r[i]
  }
  
  mat <- mat + t(mat)
  diag(mat) <- 1
  return(mat)
}

# Model parameters.
n <- 100000  # The number of observations.
p <- 100     # The dimension of each observation.
m <- 20      # The number of off-diagonal covariance terms in each direction.

# Get the covariance matrix, set number of non zero off diagonals at 40
# First, use sparse matrices and check the size
m.sparse <- GetCov(p, m, .9, T)
object.size(m.sparse)

# Now use dense matrices and check the size
m.dense <- GetCov(p, m, .5, F)
object.size(m.dense)
