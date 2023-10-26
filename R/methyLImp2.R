#' Impute missing values in methylation dataset
#'
#' This function performs missing value imputation specific for DNA methylation 
#' data. The method is based on linear regression since methylation levels 
#' show a high degree of inter-sample correlation. Implementation is 
#' parallelised over chromosomes to improve the running time.
#'
#' @param input either a numeric data matrix with missing values, 
#' with samples in rows and variables (probes) in named columns, or a 
#' SummarizedExperiment object, from which the first assays slot (with 
#' variables in rows and samples in columns, as standard) will be imputed.
#' @param type a type of data, 450K or EPIC. Type is used to split CpGs across 
#' chromosomes. Match of CpGs to chromosomes is taken from Illumina website. 
#' However, it is a part of the package and is not dynamically 
#' updated. Therefore, if you wish to provide your own match, specify "user" in 
#' this argument and provide a data frame in the next argument. 
#' @param annotation a data frame, user provided match between CpG sites and 
#' chromosomes. Must contain two columns: cpg and chr. Choose "user" in the 
#' previous argument to be able to provide user annotation.
#' @param range a vector of two numbers, \eqn{min} and \eqn{max}, 
#' specifying the range of values in the data. 
#' Since we assume the beta-value representation of the methylation data, 
#' the default range is \eqn{[0, 1]}. 
#' However, if a user wishes to apply the method to other kind of data, 
#' they can change the range in this argument.
#' @param groups a vector of the same length as the number of samples that 
#' identifies what groups does each sample correspond. Unique elements of the 
#' vector will be identified as groups and data will be split according to 
#' them. Imputation will be done for each group separately. The default is NULL, 
#' so all is considered as one group. 
#' @param skip_imputation_ids a numeric vector of ids of the columns with NAs for which 
#' \emph{not} to perform the imputation. If \code{NULL}, all columns are considered.
#' @param ncores number of cores to use in parallel computation. 
#' If \code{NULL}, set to \eqn{\#physical cores - 1}.
#' @param minibatch_frac a number, what percentage of samples to use for 
#' mini-batch computation. The default is 1 (i.e., 100\% of samples are used, 
#' no mini-batch). Remember that if your data has several groups, mini-batch 
#' will be applied to each group separately but with the same \%, 
#' so choose it accordingly. However, if your chosen \% will be smaller than 
#' matrix dimension for some groups, mini-batch will be just ignored.
#' @param minibatch_reps a number, how many times repeat computations with 
#' a fraction of samples (more times -> better performance). The default is 1 
#' (as a companion to default fraction of 100\%. i.e. no mini-batch).
#'
#' @importFrom BiocParallel SnowParam bplapply bpstart bpstop
#' @importFrom parallel detectCores
#' @importFrom SummarizedExperiment assays
#' @importFrom methods is
#'
#' @return Either a numeric matrix with imputed data or a SummarizedExperiment 
#' object where the assay slot is replaced (!) with imputed data. 
#'
#' @export
#' 
#' @examples
#' data(beta)
#' beta_with_nas <- generateMissingData(beta, lambda = 3.5)$beta_with_nas
#' beta_imputed <- methyLImp2(input = beta_with_nas, type = "EPIC", 
#'                           minibatch_frac = 0.5, ncores = 1)

methyLImp2 <- function(input,
                      type = c("450K", "EPIC", "user"),
                      annotation = NULL,
                      range = NULL,
                      groups = NULL,
                      skip_imputation_ids = NULL,
                      ncores = NULL,
                      minibatch_frac = 1, minibatch_reps = 1) {
  
    #check the class of input
    if (is.matrix(input)) {
        if (is.numeric(input)) {
            data_full <- input
            flag <- "matrix"
        } else {
            stop("The input matrix is not numeric.") 
        }
    } else if (is(input, "SummarizedExperiment")) {
        data_full <- t(assays(input)[[1]])
        flag <- "SE"
    } else {
        stop("Input is not a matrix or a SummarizedExperiment object.")
    }
    
    #load annotation
    type <- match.arg(type)
    if (type == "450K") {
        data("anno450K")
        anno <- anno450K
    } else if (type == "EPIC") {
        data("annoEPIC")
        anno <- annoEPIC
    } else if (type == "user") {
        if (!is.null(annotation)) {
            anno <- annotation
        } else stop("Please provide annotation: a data frame with cpg and chr columns.")
    } else stop("Choose type of data correctly to select annotation: 450K or EPIC - or provide your annotation.")
    
    #check the provided data range
    if (is.null(range)) {
        min <- 0
        max <- 1
    } else {
        if (length(range) != 2) {
            stop("Range must be a vector of two numbers.")
        } else {
            if (range[1] >= range[2]) {
                stop("Left entry of the interval must be less than the right one.")
            } else {
                min <- range[1]
                max <- range[2]
            }
        }
    }
    
    #check groups
    if (is.null(groups)) {
        groups <- rep(1, dim(data_full)[1])
    } else if (is.vector(groups)) {
        if (length(groups) != dim(data_full)[1]) {
            stop("Length of the groups vector does not match number of samples in the input argument.")
        } 
    } else {
        stop("Object provided for the groups argument is not a vector.")
    }
  
    unique_groups <- unique(groups)
    ngroups <- length(unique_groups)
    samples_order <- rownames(data_full)
    out_full <- vector(mode = "list", length = ngroups)
    
    for (i in seq_len(ngroups)) {
        curr_group <- unique_groups[i]
        data <- data_full[groups == curr_group, ]
        
        #split the data by chromosomes
        chromosomes <- unique(anno$chr)
        nchr <- length(chromosomes)
        data_chr <- vector(mode = "list", length = nchr)
        names(data_chr) <- chromosomes
        probes <- colnames(data)
        for (c in seq_len(nchr)) {
            curr_cpgs <- anno[anno$chr == chromosomes[c], ]$cpg
            data_chr[[c]] <- data[, colnames(data) %in% curr_cpgs]
            probes <- probes[! probes %in% curr_cpgs]
            print(length(probes))
        }
        if (length(probes) != 0) {
            #answer <- readline("Some probes are not present in our annotation")
            warning("Some probes are not present in oput annotation. They will not be imputed or used for imputation. Consider using a custom annotation.")
        }
        #drop chromosomes that were not present in the the input dataset 
        #(if such exist)
        data_chr <- data_chr[names(data_chr)[vapply(data_chr, 
                                                    function(x) 
                                                dim(x)[2], numeric(1)) != 0]]
        nchr <- length(data_chr)
        #reorder chromosomes from biggest to smallest
        chr_order <- names(sort(vapply(data_chr, function(x) dim(x)[2],
                                       numeric(1)), 
                                decreasing = TRUE))
        data_chr <- data_chr[chr_order]
        
        #run methyLImp in parallel for each chromosome
        if (is.null(ncores)) {
            ncores <- detectCores(logical = FALSE) - 1
        }

        #we set seed inside the parameters so that random sampling for mini-batch
        #is reproducible
        parallel_param <- SnowParam(workers = ncores, type = "SOCK",
                                    tasks = nchr,
                                    exportglobals = FALSE, 
                                    exportvariables = FALSE,
                                    RNGseed = 1994)
        bpstart(parallel_param)
        res <- bplapply(data_chr, methyLImp2_internal, 
                         min, max, skip_imputation_ids,
                         minibatch_frac, minibatch_reps,
                         BPPARAM = parallel_param)
        bpstop(parallel_param)
        
        names(res) <- names(data_chr)
        
        #if some chromosomes didn't have NAs or didn't have enough data for imputation,
        #tell this to the user and then later skip them in the replacement stage
        problem_chr <- c()
        if (sum(vapply(res, function(x) is.character(x), logical(1))) > 0) {
            problem_chr <- names(data_chr)[vapply(res, function(x) is.character(x), 
                                                  logical(1))]
            for (i in length(problem_chr)) {
                message("For group " , curr_group, "and chromosome ", 
                        problem_chr[i], ": ", res[problem_chr[i]])
            }
        }
        
        #replace NAs with imputed values
        out <- data
        for (c in seq_len(nchr)) {
            if (names(res)[c] %in% problem_chr) next
            imputed <- res[[c]]
            rows_id <- rownames(data) %in% rownames(imputed)
            cols_id <- colnames(data) %in% colnames(imputed)
            out[rows_id, cols_id] <- imputed
        }
        
        out_full[[i]] <- out
    }
    
    out_full <- do.call(rbind, out_full)
    out_full <- out_full[samples_order, ]
    
    if (flag == "SE") {
        assays(input)[[1]] <- t(out_full)
        return(input)
    } else {
        return(out_full)
    }
    
}
