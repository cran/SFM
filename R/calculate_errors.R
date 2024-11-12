#' @name calculate_errors
#' @title calculate_errors Function
#' @description This function calculates the Mean Squared Error (MSE) and relative error for factor loadings and uniqueness estimates obtained from factor analysis.
#'
#' @usage calculate_errors(data, A, D)
#' @param data Matrix of SFM data.
#' @param A Matrix of true factor loadings.
#' @param D Matrix of true uniquenesses.
#' @return A named vector containing:
#' \item{MSEA}{Mean Squared Error for factor loadings.}
#' \item{MSED}{Mean Squared Error for uniqueness estimates.}
#' \item{LSA}{Relative error for factor loadings.}
#' \item{LSD}{Relative error for uniqueness estimates.}
#' @export
#' @examples
#' set.seed(123) # For reproducibility
#' # Define dimensions
#' n <- 10  # Number of samples
#' p <- 5   # Number of factors
#' 
#' # Generate matrices with compatible dimensions
#' A <- matrix(runif(p * p, -1, 1), nrow = p)  # Factor loadings matrix (p x p)
#' D <- diag(runif(p, 1, 2))  # Uniquenesses matrix (p x p)
#' data <- matrix(runif(n * p), nrow = n)  # Data matrix (n x p)
#' 
#' # Calculate errors
#' errors <- calculate_errors(data, A, D)
#' print(errors)


calculate_errors <- function(data, A, D) {
  # Estimate factor loadings and uniquenesses using SOPC
  estimation_results <- SOPC(data, m = ncol(data), gamma = 0.1, eta = 0.8)
  Aso <- estimation_results$Aso
  Dso <- estimation_results$Dso
  
  frobenius.norm <- function(mat) {
    return(sqrt(sum(mat^2)))
  }
  
  # Calculate Mean Squared Error (MSE) and relative error
  MSEA <- frobenius.norm(Aso - A)^2 / (ncol(A)^2)
  MSED <- frobenius.norm(Dso - D)^2 / (ncol(D)^2)
  LSA <- frobenius.norm(Aso - A)^2 / frobenius.norm(A)^2
  LSD <- frobenius.norm(Dso - D)^2 / frobenius.norm(D)^2
  
  return(c(MSEA = MSEA, MSED = MSED, LSA = LSA, LSD = LSD))
}