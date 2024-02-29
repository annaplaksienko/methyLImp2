test_that("Missing value generation is correct", {
    expect_equal({
        set.seed(1994)
        A <- matrix(rep(1, 1000), nrow = 10, ncol = 100)
        unlist(generateMissingData(A, lambda = 1)$na_positions)
    }, c(`88.na_col` = 88L, `88.na_rows1` = 4L, `88.na_rows2` = 6L, 
         `88.na_rows3` = 5L, `24.na_col` = 24L, `24.na_rows` = 5L, `11.na_col` = 11L, 
         `11.na_rows1` = 2L, `11.na_rows2` = 10L, `11.na_rows3` = 9L))
})