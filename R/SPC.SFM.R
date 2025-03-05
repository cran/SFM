#' @name SPC.SFM
#' @title Apply the SPC method to the Skew factor model
#' @description This function performs Sparse Principal Component Analysis (SPC) on the input data. It estimates factor loadings and uniquenesses while calculating mean squared errors and loss metrics for comparison with true values.
#'
#' @param data The data used in the SPC analysis.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @param m The number of common factors.
#' @param p The number of variables.
#' @return A list containing:
#' \item{As}{Estimated factor loadings, a matrix of estimated factor loadings from the SPC analysis.}
#' \item{Ds}{Estimated uniquenesses, a vector of estimated uniquenesses corresponding to each variable.}
#' \item{MSESigmaA}{Mean squared error of the estimated factor loadings (As) compared to the true loadings (A).}
#' \item{MSESigmaD}{Mean squared error of the estimated uniquenesses (Ds) compared to the true uniquenesses (D).}
#' \item{LSigmaA}{Loss metric for the estimated factor loadings (As), indicating the relative error compared to the true loadings (A).}
#' \item{LSigmaD}{Loss metric for the estimated uniquenesses (Ds), indicating the relative error compared to the true uniquenesses (D).}
#' \item{tau}{Proportion of zero factor loadings in the estimated loadings matrix (As).}
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
#' results <- SPC.SFM(data, A, D, m, p)
#' print(results)
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
SPC.SFM <- function(data, A, D, m, p) {
  As <- SPC(data, m = m, gamma = 0.1)$As
  Ds <- SPC(data, m = m, gamma = 0.1)$Ds
  MSESigmaA <- frobenius.norm(As - A)^2 / (p^2)
  MSESigmaD <- frobenius.norm(Ds - D)^2 / (p^2)
  LSigmaA <- frobenius.norm(As - A)^2 / frobenius.norm(A)^2
  LSigmaD <- frobenius.norm(Ds - D)^2 / frobenius.norm(D)^2
  tau <- as.vector(table(As == 0) / (p * m))[2]

  return(list('As' = As,
              'Ds' = Ds,
              'MSESigmaA' = MSESigmaA,
              'MSESigmaD' = MSESigmaD,
              'LSigmaA' = LSigmaA,
              'LSigmaD' = LSigmaD,
              'tau' = tau))
}
