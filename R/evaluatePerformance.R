#' Extract a vector of beta values of interest out of a matrix, 
#' given the positions of artificial NAs
#'
#' @param beta numeric data matrix
#' @param na_positions a list where each element is a list of two elements: 
#' column id and ids of rows with NAs in that column (structure matches the 
#' output of generateMissingData function). We need this because some NAs in 
#' the dataset are from real data and not artificial, so we can't evaluate the 
#' performance of the method on them since we do not know real value. 
#' Therefore, we need to know the positions of artificial NAs. 
#' 
#' @return A numerical vector of beta values
#'
#' @keywords internal

extract_values <- function(beta, na_positions) {
    
    nna <- length(na_positions)
    
    beta_vector <- c()
    
    for (l in seq_len(nna)) {
        col <- na_positions[[l]]$na_col
        rows <- na_positions[[l]]$na_rows 
        
        curr_values <- beta[rows, col]
        beta_vector <- c(beta_vector, curr_values)
    }
    
    names(beta_vector) <- NULL
    
    return(beta_vector)
}

#' Evaluate performance metrics: root mean square error (RMSE), 
#' mean absolute error (MAE), mean absolute percentage error (MAPE) and
#' Pearson correlation coefficient (PCC).
#'
#' This function computes performance metrics as an element-wise difference 
#' between two matrices (skipping NA elements that were not imputed): 
#' * \eqn{RMSE = \sqrt{\sum_i (true_i - est_i)^2 / \#NAs)}};
#' * \eqn{MAE = \sum_i |true_i - est_i| / \#NAs};
#' * \eqn{MAPE = \dfrac{100}{n} \sum_i |true_i - est_i / true_i|} 
#' (here we omit the true beta-values equal to 0 and their  predicted values
#' to  avoid an indeterminate measure);
#' * \eqn{PCC = \dfrac{\sum_i (true_i - \bar{true_i}) \sum_i(est_i - \bar{est}_i)}{\sqrt{\sum_i (true_i - \bar{true}_i)^2}} \sqrt{\sum_i(est_i - \bar{est_i})^2}}.
#'
#' @param beta_true first numeric data matrix.
#' @param beta_imputed second numeric data matrix
#' @param na_positions a list where each element is a list of two elements: 
#' column id and ids of rows with NAs in that column (structure matches the 
#' output of generateMissingData function). We need this because some NAs in 
#' the dataset are from real data and not artificial, so we can't evaluate the 
#' performance of the method on them since we do not know real value. 
#' Therefore, we need to know the positions of artificial NAs. 
#' 
#' @return A numerical vector of four numbers: root mean square error (RMSE), 
#' mean absolute error (MAE), mean absolute percentage error (MAPE) and
#' Pearson correlation coefficient (PCC).
#'
#' @importFrom stats cor
#' 
#' @export
#' 
#' @examples
#' data(beta)
#' with_missing_data <- generateMissingData(beta, lambda = 3.5)
#' beta_with_nas <- with_missing_data$beta_with_nas
#' na_positions <- with_missing_data$na_positions
#' beta_imputed <- methyLImp2(input = beta_with_nas, type = "EPIC", 
#'                           minibatch_frac = 0.5, ncores = 1)
#' evaluatePerformance(beta, beta_imputed, na_positions)

evaluatePerformance <- function(beta_true, beta_imputed, na_positions) {

    beta_true_vec <- extract_values(beta_true, na_positions)
    beta_imputed_vec <- extract_values(beta_imputed, na_positions)
    
    RMSE <- sqrt(mean((beta_true_vec - beta_imputed_vec)^2, na.rm = TRUE))
    MAE <- mean(abs(beta_true_vec - beta_imputed_vec), na.rm = TRUE)
    PCC <- cor(beta_true_vec, beta_imputed_vec, use = "pairwise.complete.obs")
    
    Filt <- abs(beta_true_vec) != 0
    MAPE <- mean(abs((beta_true_vec[Filt] - beta_imputed_vec[Filt]) / 
                         beta_true_vec[Filt]), na.rm = TRUE)
    
    performance <- c("RMSE" = RMSE, "MAE" = MAE, 
                     "PCC" = PCC, "MAPE" = MAPE)
    return(performance)
}