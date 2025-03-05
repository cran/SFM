#' @name OPC.SFM
#' @title Apply the OPC method to the Skew factor model
#' @description This function computes Online Principal Component Analysis (OPC) for the provided input data, estimating factor loadings and uniquenesses. It calculates mean squared errors and sparsity for the estimated values compared to true values.
#' @param data A matrix of input data.
#' @param m The number of principal components.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @param p The number of variables.
#' @return A list containing:
#' \item{Ao}{Estimated factor loadings.}
#' \item{Do}{Estimated uniquenesses.}
#' \item{MSEA}{Mean squared error for factor loadings.}
#' \item{MSED}{Mean squared error for uniquenesses.}
#' \item{tau}{The sparsity.}
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
#' results <- OPC.SFM(data, m, A, D, p)
#' print(results)
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
OPC.SFM <- function(data,m=m, A, D, p){
Ao=OPC(data,m=m,eta=0.5)$Ao
Do=OPC(data,m=m,eta=0.5)$Do
MSEA=frobenius.norm(Ao-A)^2/(p^2)
MSED=frobenius.norm(Do-diag(D))^2/(p^2)
tau=as.vector(table(Ao==0)/(p*m))[2]

  return(list(Ao = Ao,
              Do = Do,
              MSEA = MSEA,
              MSED = MSED,
              tau = tau))
}

