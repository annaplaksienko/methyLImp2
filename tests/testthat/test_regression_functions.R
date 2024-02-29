test_that("Pseudo-logit function is correct", {
    expect_equal({
        A <- matrix(c(0.11, 0.12, 
                      0.21, 0.22), byrow = TRUE, nrow = 2)
        round(plogit(A), digits = 2)
    },
    matrix(c(-2.09, -1.99,
             -1.32, -1.27), byrow = TRUE, nrow = 2))
})

test_that("Inverse of pseudo-logit function is correct", {
    expect_equal({
        B <- matrix(c(-2.09, -1.99,
                      -1.32, -1.27), byrow = TRUE, nrow = 2)
        round(inv.plogit(B), digits = 2)
    },
    matrix(c(0.11, 0.12, 
             0.21, 0.22), byrow = TRUE, nrow = 2))
})

test_that("Moore-Penrose generalized inverse of a matrix is correct", {
    expect_equal({
        A <- matrix(c(0.11, 0.12, 
                      0.21, 0.22), byrow = TRUE, nrow = 2)
        pinvr(A)
    },
    matrix(c(-220, 120,
             210, -110), byrow = TRUE, nrow = 2))
})