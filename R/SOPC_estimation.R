#' @name SOPC_estimation
#' @title SOPC Estimation Function
#' @description This function processes Skew Factor Model (SFM) data using the Sparse Online Principal Component (SOPC) method. 
#' 
#' @param data Matrix of SFM data.
#' @param gamma Tuning parameter for the sparseness of the loadings matrix.
#' @param eta Tuning parameter for the sparseness of the common factors  matrix.
#' @usage SOPC_estimation(data, gamma, eta)
#'
#' @return A list containing:
#' \item{Aso}{Estimated factor loadings.}
#' \item{Dso}{Estimated common factors .}
#' \item{tauA}{Sparsity of the loadings matrix, calculated as the proportion of zeros.}
#' @export
#' @examples
#' set.seed(123) # For reproducibility
#' data <- matrix(runif(200), nrow = 20) # Skew Factor Model data
#' sopc_results <- SOPC_estimation(data, 0.1, 0.8)
#' print(sopc_results)
#' 
SOPC_estimation <- function(data, gamma, eta) {
  # Estimate factor loadings and uniquenesses using SOPC
  estimation_results <- SOPC(data, m = ncol(data), gamma = gamma, eta = eta)
  Aso <- estimation_results$Aso
  Dso <- estimation_results$Dso
  
  # Calculate the sparsity of the loadings matrix
  tauA <- as.vector(table(Aso == 0) / (ncol(data) * nrow(data)))[2]
  
  return(list(Aso = Aso, Dso = Dso, tauA = tauA))
}

