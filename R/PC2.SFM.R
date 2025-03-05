#' @name PC2.SFM
#' @title Apply the PC method to the Laplace factor model
#' @description This function performs Principal Component Analysis (PCA) on a given data set to reduce dimensionality. It calculates the estimated values for the loadings, specific variances, and the covariance matrix.
#' @param data The total data set to be analyzed.
#' @param m The number of principal components to retain in the analysis.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @return A list containing:
#' \item{A2}{Estimated factor loadings.}
#' \item{D2}{Estimated uniquenesses.}
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
#' results <- PC2.SFM(data, m, A, D)
#' print(results)
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
PC2.SFM <- function(data, m, A, D) {
  X<-scale(data)
  R<-cor(X)
  S<-R
  eig<-eigen(S)
  p<-nrow(S)
  diag_S<-diag(S)
  sum_rank<-sum(diag_S)
  rowname<-paste("X",1:p,sep="")
  colname<-paste("Factor",1:m,sep="")
  Ahat<-matrix(0,nrow=p,ncol=m, dimnames=list(rowname,colname))
  rowname<-c("SS loadings","Proportion Var","Cumulative Var")
  for (i in 1:m){
    Ahat[,i]<-sqrt(eig$values[i])*eig$vectors[,i]
  }
  h2<-diag(Ahat%*%t(Ahat))
  Dhat<-diag(S-h2)

  MSESigmaA = frobenius.norm(Ahat - A)^2 / (p^2)
  MSESigmaD = frobenius.norm(Dhat - D)^2 / (p^2)
  LSigmaA = frobenius.norm(Ahat - A)^2 / frobenius.norm(A)^2
  LSigmaD = frobenius.norm(Dhat - D)^2 / frobenius.norm(D)^2
  return(list(A2 = Ahat,
              D2 = Dhat,
              MSESigmaA = MSESigmaA,
              MSESigmaD = MSESigmaD,
              LSigmaA = LSigmaA,
              LSigmaD = LSigmaD))
}
