#' The SFM function is to generate Skew Factor Models data.
#'
#' The function supports various distribution types for generating the data,
#' including: Skew-Normal Distribution, Skew-Cauchy Distribution, Skew-t Distribution.
#'
#' @name SFM
#'
#' @usage SFM(n, p, m, xi, omega, alpha, distribution_type)
#' @param n Sample size.
#' @param p Sample dimensionality.
#' @param m Number of factors.
#' @param xi A numerical parameter used exclusively in the "Skew-t" distribution, representing the distribution's xi parameter.
#' @param omega A numerical parameter representing the omega parameter of the distribution, which affects the degree of skewness in the distribution.
#' @param alpha A numerical parameter representing the alpha parameter of the distribution, which influences the shape of the distribution.
#' @param distribution_type The type of distribution.
#'
#' @return A list containing:
#' \item{data}{A matrix of generated data.}
#' \item{A}{A matrix representing the factor loadings.}
#' \item{D}{A diagonal matrix representing the unique variances.}
#'
#' @examples
#' library(MASS)
#' library(SOPC)
#' library(sn)
#' library(matrixcalc)
#' library(psych)
#' n <- 100
#' p <- 10
#' m <- 5
#' xi <- 5
#' omega <- 2
#' alpha <- 5
#' distribution_type <- "Skew-Normal Distribution"
#' X <- SFM(n, p, m, xi, omega, alpha, distribution_type)
#' 
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats cov
SFM <- function(n, p, m, xi, omega, alpha, distribution_type) {
  mu <- t(matrix(rep(runif(m, 0, 100), n), m, n))
  mu0 <- as.matrix(runif(p, 0))
  sigma0 <- diag(runif(p, 1))
  F <- matrix(MASS::mvrnorm(n, mu0, sigma0), nrow = n)
  A <- matrix(runif(p * m, -1, 1), nrow = m)
  
  if (distribution_type == "Skew-Normal Distribution") {
    epsilon <- matrix(sn::rsn(n * m, omega = omega, alpha = alpha), nrow = n)
  } else if (distribution_type == "Skew-Cauchy Distribution") {
    epsilon <- matrix(sn::rsc(n * m, omega = omega, alpha = alpha), nrow = n)
  } else if (distribution_type == "Skew-t Distribution") {
    epsilon <- matrix(sn::rst(n * m, xi = xi, omega = omega, alpha = alpha), nrow = n)
  } else {
    stop("Unsupported distribution_type")
  }
  
  D0 <- omega * diag(p)
  D <- diag(D0)
  data <- mu + F %*% t(A) + epsilon
  
  return(list(data = data, A = A, D = D))
}

