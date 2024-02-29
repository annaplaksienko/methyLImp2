test_that("Chromosome split is correct", {
    expect_equal({
        A <- matrix(rep(1, 1000), nrow = 10, ncol = 100)
        colnames(A) <- paste("col", 1:100, sep = "")
        anno <- data.frame(cpg = paste("col", 1:100, sep = ""), 
                           chr = c(rep("A", 50), rep("B", 50)))
        A_chr <- methyLImp2:::split_by_chromosomes(A, type = "user", 
                                                   annotation = anno)
        length(A_chr)
    },
    2)
    expect_equal({
        A <- matrix(rep(1, 1000), nrow = 10, ncol = 100)
        colnames(A) <- paste("col", 1:100, sep = "")
        anno <- data.frame(cpg = paste("col", 1:100, sep = ""), 
                           chr = c(rep("A", 50), rep("B", 50)))
        A_chr <- methyLImp2:::split_by_chromosomes(A, type = "user", 
                                                   annotation = anno)
        sapply(A_chr, function(x) class(x)[1])
    },
    c(A = "matrix", B = "matrix"))
})

test_that("Imputation is correct", {
    expect_equal({
        set.seed(1994)
        A <- matrix(rep(1, 1000), nrow = 10, ncol = 100)
        rownames(A) <- paste("row", 1:10, sep = "")
        colnames(A) <- paste("col", 1:100, sep = "")
        anno <- data.frame(cpg = paste("col", 1:100, sep = ""), 
                           chr = c(rep("A", 50), rep("B", 50)))
        A_with_nas <- generateMissingData(A, lambda = 1)
        A_imp <- methyLImp2(input = A_with_nas$beta_with_nas,
                            type = "user",
                            annotation = anno)
        methyLImp2:::extract_values(A_imp, A_with_nas$na_positions)
    },
    rep(1, 7))
})

