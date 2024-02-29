test_that("Extraction of values from the matrix is correct", {
    expect_equal({
    A <- matrix(c(11, 12, 13,
                  21, 22, 23,
                  31, 32, 33), byrow = TRUE, nrow = 3)
    na_positions <- list("1" = list("na_col" = c(1), "na_rows" = c(1)),
                         "3" = list("na_col" = c(3), "na_rows" = c(2, 3)))
    A_vec <- methyLImp2:::extract_values(A, na_positions)
    
}, c(11, 23, 33))
    })

test_that("Performance evaluation is correct", {
          expect_equal({
              A <- matrix(c(11, 12, 13,
                            21, 22, 23,
                            31, 32, 33), byrow = TRUE, nrow = 3)
              na_positions <- list("1" = list("na_col" = c(1), "na_rows" = c(1)),
                                   "3" = list("na_col" = c(3), "na_rows" = c(2, 3)))
              
              A_imputed <- matrix(c(11, 12.2, 13,
                                    21, 22, 23,
                                    31, 32.2, 33.2), byrow = TRUE, nrow = 3)
              evaluatePerformance(A, A_imputed, na_positions)
}, c(RMSE = 0.115470053837927, MAE = 0.0666666666666676, PCC = 0.999983980905238, 
     MAPE = 0.00202020202020205))
    })
