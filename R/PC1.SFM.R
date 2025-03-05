#' @name PC1.SFM
#' @title Apply the PC method to the Laplace factor model
#' @description This function performs Principal Component Analysis (PCA) on a given data set to reduce dimensionality. It calculates the estimated values for the loadings, specific variances, and the covariance matrix.
#' @param data The total data set to be analyzed.
#' @param m The number of principal components to retain in the analysis.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @return A list containing:
#' \item{A1}{Estimated factor loadings.}
#' \item{D1}{Estimated uniquenesses.}
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
#' results <- PC1.SFM(data, m, A, D)
#' print(results)
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
PC1.SFM <- function(data, m, A, D) {
  X <- scale(data, scale = FALSE)
  p <- ncol(X)
  n <- nrow(X)
  Sigmahat <- cov(X)
  eig <- eigen(Sigmahat)
  lambdahat <- eig$values[1:m]
  ind <- order(lambdahat, decreasing = TRUE)
  lambdahat <- lambdahat[ind]
  Q <- eig$vectors
  Q <- Q[, ind]
  Qhat <- Q[, 1:m]
  Ahat <- matrix(0, nrow = p, ncol = m)
  for (j in 1:m) {
    Ahat[, j] <- sqrt(lambdahat[j]) * Qhat[, j]
  }
  h0 <- diag(Ahat %*% t(Ahat))
  Dhat <- diag(Sigmahat - h0)
  S2 <- Ahat %*% t(Ahat) + Dhat

  MSESigmaA = frobenius.norm(Ahat - A)^2 / (p^2)
  MSESigmaD = frobenius.norm(Dhat - D)^2 / (p^2)
  LSigmaA = frobenius.norm(Ahat - A)^2 / frobenius.norm(A)^2
  LSigmaD = frobenius.norm(Dhat - D)^2 / frobenius.norm(D)^2
  return(list(A1 = Ahat,
              D1 = Dhat,
              MSESigmaA = MSESigmaA,
              MSESigmaD = MSESigmaD,
              LSigmaA = LSigmaA,
              LSigmaD = LSigmaD))
}
