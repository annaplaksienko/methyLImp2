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
  total_nna <- 0
  
  #initiate the sum that in the end we'll divide by total_nna        
  RMSE_sum <- MAE_sum <- 0
  
  for (l in 1:nna) {
    col <- na_positions[[l]]$na_col
    rows <- na_positions[[l]]$na_rows 
    
    True <- beta_true[rows, col]
    Pred <- beta_imputed[rows, col]
    
    #note that if some values were note imputed 
    #because there were too many NAs in one column and we didn't consider it
    #then we don't evaluate the performance for that value here
    total_nna <- total_nna + sum(!is.na(Pred))
    
    RMSE_sum <- RMSE_sum + sum((Pred - True)^2, na.rm = TRUE)
    MAE_sum <- MAE_sum + sum(abs(Pred - True), na.rm = TRUE)
  }
  
  RMSE <- sqrt(RMSE_sum / total_nna)
  MAE <- MAE_sum / total_nna
  
  performance <- c("RMSE" = RMSE, "MAE" = MAE)
  
  return(performance)
}