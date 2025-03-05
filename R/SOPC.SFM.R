#' @name SOPC.SFM
#' @title SOPC Estimation Function
#' @description This function processes Skew Factor Model (SFM) data using the Sparse Online Principal Component (SOPC) method.
#' @param data A numeric matrix containing the data used in the SOPC analysis.
#' @param m An integer specifying the number of subsets or common factors.
#' @param p An integer specifying the number of variables in the data.
#' @param A A numeric matrix representing the true factor loadings.
#' @param D A numeric matrix representing the true uniquenesses.
#' @usage SOPC.SFM(data, m, p, A, D)
#' @return A list containing the following metrics:
#' \item{Aso}{Estimated factor loadings matrix.}
#' \item{Dso}{Estimated uniquenesses matrix.}
#' \item{MSEA}{Mean squared error of the estimated factor loadings (Aso) compared to the true loadings (A).}
#' \item{MSED}{Mean squared error of the estimated uniquenesses (Dso) compared to the true uniquenesses (D).}
#' \item{LSA}{Loss metric for the estimated factor loadings (Aso), indicating the relative error compared to the true loadings (A).}
#' \item{LSD}{Loss metric for the estimated uniquenesses (Dso), indicating the relative error compared to the true uniquenesses (D).}
#' \item{tauA}{Proportion of zero factor loadings in the estimated loadings matrix (Aso), representing the sparsity.}
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
#' results <- SOPC.SFM(data, m, p, A, D)
#' print(results)
#' @export
#' @importFrom SOPC SOPC
#' @importFrom matrixcalc frobenius.norm

SOPC.SFM <- function(data, m, p, A, D){
  Aso=SOPC(data,m=m,gamma=0.1,eta=0.8)$Aso
  Dso=SOPC(data,m=m,gamma=0.1,eta=0.8)$Dso
  MSEA=frobenius.norm(Aso-A)^2/(p^2)
  MSED=frobenius.norm(Dso-D)^2/(p^2)
  LSA=frobenius.norm(Aso-A)^2/frobenius.norm(A)^2
  LSD=frobenius.norm(Dso-D)^2/frobenius.norm(D)^2
  tauA=as.vector(table(Aso==0)/(p*m))[2]
  return(c('Aso'=Aso,
           'Dso'=Dso,
           'MSEA'=MSEA,
           'MSED'=MSED,
           'LSA'=LSA,
           'LSD'=LSD,
           'tauA'=tauA))
}
