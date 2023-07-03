#' Evaluate root mean square error (RMSE) and mean absolute error (MAE)
#'
#' This function computes root mean square error and mean absolute error as 
#' an element-wise difference between two matrices (apart from NA elements): 
#' \eqn{RMSE = \sqrt{\sum_i (true_i - est_i)^2 / \#NAs)}}, 
#' \eqn{MAE = \sum_i |true_i - est_i| / \#NAs}.
#'
#' @param beta_true first numeric data matrix.
#' @param beta_imputed second numeric data matrix
#' @param na_positions a list where each element is a list of two elements: 
#' column id and ids of rows with NAs in that column. 
#' We need this because some NAs in the dataset are from real data and 
#' not artificial, so we can't evaluate the performance of the method 
#' on them since we do not know real value. 
#' Therefore, we need to know the positions of artificial NAs. 
#' 
#' @return A numerical vector of two numbers, root mean square error 
#' and mean absolute error.
#'
#' @export
#' 
#' @examples
#' {
#' data(beta)
#' with_missing_data <- generateMissingData(beta, lambda = 3.5)
#' beta_with_nas <- with_missing_data$beta_with_nas
#' na_positions <- with_missing_data$na_positions
#' beta_imputed <- methyLImp2(input = beta_with_nas, type = "EPIC", 
#'                           minibatch_frac = 0.5, ncores = 1)
#' evaluatePerformance(beta, beta_imputed, na_positions)
#' }

evaluatePerformance <- function(beta_true, beta_imputed, na_positions) {

    nna <- length(na_positions)
    total_nna <- 0
  
    #initiate the sum that in the end we'll divide by total_nna        
    RMSE_sum <- MAE_sum <- 0
  
    for (l in seq_len(nna)) {
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