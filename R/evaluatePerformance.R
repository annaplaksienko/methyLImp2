#' Evaluate root mean square error (RMSE) and mean absolute error (MAE)
#'
#' This function computes root mean square error and mean absolute error as an element-wise difference between two matrices (apart from NA elements): RMSE = sqrt(sum_i (beta_i - hat_beta_i)^2 / #artificial_NAs), MAE = sum_i |beta_i - hat_beta_i| / #artificial_NAs.
#'
#' @param beta_true first numeric data matrix.
#' @param beta_imputed second numeric data matrix
#' 
#' @return A number, root mean square error.
#'
#' @export

evaluatePerformance <- function(beta_true, beta_imputed, na_positions) {

  nna <- length(na_positions)
  total_nna <- sum(sapply(na_positions, function(x) length(x$na_rows)))
  
  #initiate the sum that in the end we'll divide by total_nna        
  RMSE_sum <- MAE_sum <- 0
  
  for (l in 1:nna) {
    col <- na_positions[[l]]$na_col
    rows <- na_positions[[l]]$na_rows 
    
    True <- beta_true[rows, col]
    Pred <- beta_imputed[rows, col]
    
    RMSE_sum <- RMSE_sum + sum((Pred - True)^2)
    MAE_sum <- MAE_sum + sum(abs(Pred - True))
  }
  
  RMSE <- sqrt(RMSE_sum / total_nna)
  MAE <- MAE_sum / total_nna
  
  performance <- c("RMSE" = RMSE, "MAE" = MAE)
  
  return(performance)
}