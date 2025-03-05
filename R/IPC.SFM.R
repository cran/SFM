#' @name IPC.SFM
#' @title Apply the IPC method to the Skew factor model
#' @description This function performs Incremental Principal Component Analysis (IPC) on the provided data. It updates the estimated factor loadings and uniquenesses as new data points are processed, calculating mean squared errors and loss metrics for comparison with true values.
#' @param x The data used in the IPC analysis.
#' @param m The number of common factors.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @param p The number of variables.
#' @return A list of metrics including:
#' \item{Ai}{Estimated factor loadings updated during the IPC analysis, a matrix of estimated factor loadings.}
#' \item{Di}{Estimated uniquenesses updated during the IPC analysis, a vector of estimated uniquenesses corresponding to each variable.}
#' \item{MSESigmaA}{Mean squared error of the estimated factor loadings (Ai) compared to the true loadings (A).}
#' \item{MSESigmaD}{Mean squared error of the estimated uniquenesses (Di) compared to the true uniquenesses (D).}
#' \item{LSigmaA}{Loss metric for the estimated factor loadings (Ai), indicating the relative error compared to the true loadings (A).}
#' \item{LSigmaD}{Loss metric for the estimated uniquenesses (Di), indicating the relative error compared to the true uniquenesses (D).}
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
#' result <- IPC.SFM(data, m = m, A = A, D = D, p = p)
#' print(result)
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
IPC.SFM <- function(x, m, A, D, p) {
  frobenius.norm <- function(matrix) {
    matrix <- as.matrix(matrix)
    return(norm(matrix, type = "F"))
  }
  Ai = IPC(data = x, m = m, eta = 0.8)$Ai
  Di = IPC(data = x, m = m, eta = 0.8)$Di
  MSESigmaA = frobenius.norm(Ai - A)^2 / (p^2)
  MSESigmaD = frobenius.norm(Di - D)^2 / (p^2)
  LSigmaA = frobenius.norm(Ai - A)^2 / frobenius.norm(A)^2
  LSigmaD = frobenius.norm(Di - D)^2 / frobenius.norm(D)^2
  return(c('Ai' = Ai,
           'Di' = Di,
           'MSESigmaA' = MSESigmaA,
           'MSESigmaD' = MSESigmaD,
           'LSigmaA' = LSigmaA,
           'LSigmaD' = LSigmaD))
}
