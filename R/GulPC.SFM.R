#' @name GulPC.SFM
#' @title Apply the GulPC method to the Skew factor model
#' @description This function performs General Unilateral Loading Principal Component (GulPC) analysis on a given data set. It calculates the estimated values for the first layer and second layer loadings, specific variances, and the mean squared errors.
#' @param data A matrix of input data.
#' @param m The number of principal components.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @return A list containing:
#' \item{AU1}{The first layer loading matrix.}
#' \item{AU2}{The second layer loading matrix.}
#' \item{DU3}{The estimated specific variance matrix.}
#' \item{MSESigmaD}{Mean squared error for uniquenesses.}
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
#' results <- GulPC.SFM(data, m, A, D)
#' print(results)
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
GulPC.SFM <- function(data, m, A, D) {
    X=scale(data)
    n=nrow(X)
    p<-ncol(X)
    SigmaU1hat=cor(X)
    eig5<-eigen(SigmaU1hat)
    lambda1hat = eig5$values[1:m]
    ind<-order(lambda1hat,decreasing=T)
    lambda1hat<-lambda1hat[ind]
    Q1<-eig5$vectors
    Q1=Q1[,ind]
    Q1hat<- Q1[, 1:m]
    AU1 <- matrix(0, nrow = p, ncol = m)
    for (j in 1:m) {AU1[, j] <- sqrt(lambda1hat[j]) * Q1hat[, j]}; AU1
    hU1 <- diag(AU1 %*% t(AU1))
    DU1 <- diag(SigmaU1hat - hU1)
    pc=2
    F1hat=X%*%AU1
    F1star<-F1hat/sqrt(n)
    SigmaU2hat=cov(F1star)
    eig6<-eigen(SigmaU2hat)
    lambda2hat =eig6$values[1:pc]
    ind<-order(lambda2hat,decreasing=T)
    lambda2hat<-lambda2hat[ind]
    Q2<-eig6$vectors
    Q2=Q2[,ind]
    Q2hat<- Q2[, 1:pc]
    AU2<- matrix(0, nrow = m, ncol = pc)
    for (j in 1:pc) {AU2[, j] <- sqrt(lambda2hat[j]) * Q2hat[, j]}; AU2
    hU2 <- diag(AU2%*% t(AU2))
    DU2 <- diag(SigmaU2hat - hU2)
    Fhat=F1star%*%AU2
    Xhat=Fhat%*%t(AU2)%*%t(AU1)
    S1hat=cov(Xhat)
    hU3 <- diag(t(t(AU2)%*%t(AU1))%*%(t(AU2)%*%t(AU1)))
    DU3 <- diag(S1hat - hU3)

  MSESigmaD = frobenius.norm(DU3 - D)^2 / (p^2)
  LSigmaD = frobenius.norm(DU3 - D)^2 / frobenius.norm(D)^2
  # Return the results as a list
  return(list(AU1 = AU1,
              AU2 = AU2,
              DU3 = DU3,
              MSESigmaD = MSESigmaD,
              LSigmaD = LSigmaD))
}
