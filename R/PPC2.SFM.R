#' @name PPC2.SFM
#' @title Apply the PPC method to the Skew factor model
#' @description This function performs Projected Principal Component Analysis (PPC) on a given data set to reduce dimensionality. It calculates the estimated values for the loadings, specific variances, and the covariance matrix.
#' @param data The total data set to be analyzed.
#' @param m The number of principal components.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @return A list containing:
#' \item{Ap2}{Estimated factor loadings.}
#' \item{Dp2}{Estimated uniquenesses.}
#' \item{MSESigmaA}{Mean squared error for factor loadings.}
#' \item{MSESigmaD}{Mean squared error for uniquenesses.}
#' \item{LSigmaA}{Loss metric for factor loadings.}
#' \item{LSigmaD}{Loss metric for uniquenesses.}
#' @examples
#' library(SOPC)
#' library(matrixcalc)
#' library(MASS)
#' library(psych)
#' library(sn)
#' n=1000
#' p=10
#' m=5
#' mu=t(matrix(rep(runif(p,0,1000),n),p,n))
#' mu0=as.matrix(runif(m,0))
#' sigma0=diag(runif(m,1))
#' F=matrix(mvrnorm(n,mu0,sigma0),nrow=n)
#' A=matrix(runif(p*m,-1,1),nrow=p)
#' r <- rsn(n*p,0,1)
#' epsilon=matrix(r,nrow=n)
#' D=diag(t(epsilon)%*%epsilon)
#' data=mu+F%*%t(A)+epsilon
#' results <- PPC2.SFM(data, m, A, D)
#' print(results)
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
PPC2.SFM <- function(data, m,A, D) {
  # Standardize the data to have zero mean and unit variance
  X <- scale(data)

  # Get the number of observations and variables
  n <- nrow(X)
  p <- ncol(X)

  # Create a projection matrix P to eliminate noise or specific factors
  P <- as.matrix(diag(c(0, 1), n, n))

  # Project the data using the projection matrix P
  Xpro <- scale(P %*% X)

  # Calculate the covariance matrix of the projected data
  Sigmahatpro <- cov(Xpro)

  # Perform eigenvalue decomposition on the covariance matrix
  eig <- eigen(Sigmahatpro)
  lambdahat <- eig$values[1:m]  # Extract the first 'm' eigenvalues
  ind <- order(lambdahat, decreasing = TRUE)  # Sort the eigenvalues in descending order
  lambdahat <- lambdahat[ind]
  Q <- eig$vectors  # Extract the eigenvectors
  Q <- Q[, ind]
  Qhat <- Q[, 1:m]  # Select the first 'm' eigenvectors

  # Calculate the estimated loading matrix
  Apro <- matrix(0, nrow = p, ncol = m)
  for (j in 1:m) {
    Apro[, j] <- sqrt(lambdahat[j]) * Qhat[, j]
  }

  # Calculate the hat matrix (projection matrix) and the specific variance matrix
  hpro <- diag(Apro %*% t(Apro))
  Dpro <- diag(Sigmahatpro - hpro)

  MSESigmaA = frobenius.norm(Apro - A)^2 / (p^2)
  MSESigmaD = frobenius.norm(Dpro - D)^2 / (p^2)
  LSigmaA = frobenius.norm(Apro - A)^2 / frobenius.norm(A)^2
  LSigmaD = frobenius.norm(Dpro - D)^2 / frobenius.norm(D)^2
  # Return the results as a list
  return(list(Ap2 = Apro,
              Dp2 = Dpro,
              MSESigmaA = MSESigmaA,
              MSESigmaD = MSESigmaD,
              LSigmaA = LSigmaA,
              LSigmaD = LSigmaD))
}
