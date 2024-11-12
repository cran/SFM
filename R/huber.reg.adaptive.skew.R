#' Adaptive Huber Regression for Skew Factor Models
#'
#' Performs adaptive Huber regression tailored for skew factor models, 
#' and returns the estimated regression coefficients in a matrix (loading matrix) format.
#'
#' @param X A matrix of predictor variables.
#' @param Y A vector of response variables.
#' @param tau Initial robustification parameter (default is 1.35).
#' @param max_iterations Maximum number of iterations (default is 100).
#' @param tolerance Convergence tolerance (default is 1e-6).
#' @param n_factors The number of factors (columns) for the loading matrix (default is 1).
#' @return A matrix of estimated regression coefficients with dimensions `p x n_factors`.
#' @examples
#' # Generate some example data for skew factor models
#' set.seed(123)
#' n <- 200
#' d <- 10
#' beta <- rep(1, d)
#' skew_factor <- rnorm(n)  # Adding a skew factor
#' X <- matrix(rnorm(n * d), n, d)
#' err <- rnorm(n)
#' Y <- 1 + skew_factor + X %*% beta + err
#'
#' # Perform adaptive Huber regression for skew factor model
#' loading_matrix <- huber.reg.adaptive.skew(X, Y, n_factors = 3)
#' print(loading_matrix)
#'
#' @export
huber.reg.adaptive.skew <- function(X, Y, tau = 1.35, max_iterations = 100, tolerance = 1e-6, n_factors = 1) {
  # Check for dimension compatibility
  if (nrow(X) != length(Y)) {
    stop("The number of rows in X must match the length of Y.")
  }
  
  # Calculate the initial residuals using OLS
  ols_fit <- lm(Y ~ X - 1)
  residuals <- residuals(ols_fit)
  
  # Adaptive Huber regression algorithm
  iteration <- 0
  
  # Huber loss function
  huber_loss <- function(r, tau) {
    sapply(r, function(value) {
      if (abs(value) <= tau) {
        return(value^2)
      } else {
        return(tau * abs(value) - 0.5 * tau^2)
      }
    })
  }
  
  while (iteration < max_iterations) {
    iteration <- iteration + 1
    
    # Apply the Huber loss to each residual
    huber_residuals <- huber_loss(residuals, tau)
    
    # Update the weights based on the Huber loss
    adaptive_weights <- 1 / pmax(huber_residuals, 1e-6)  # Avoid division by zero
    
    # Perform weighted least squares regression
    wls_fit <- lm(Y ~ X - 1, weights = adaptive_weights)
    
    # Calculate new residuals
    new_residuals <- residuals(wls_fit)
    
    # Check for convergence
    if (max(abs(new_residuals - residuals)) < tolerance) {
      break
    }
    
    # Update residuals for the next iteration
    residuals <- new_residuals
  }
  
  # Get the coefficients and reshape as a loading matrix
  coefficients <- coef(wls_fit)
  
  # Reshape coefficients into a p x n_factors matrix
  loading_matrix <- matrix(coefficients, nrow = ncol(X), ncol = n_factors, byrow = TRUE)
  
  return(loading_matrix)
}





