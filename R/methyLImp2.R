#' Split methylation dataset into a list by chromosomes
#'
#' This function split a given methylation dataset into a list of datasets according to a given annotation.
#'
#' @param data a numeric data matrix with with samples in rows and 
#' variables (probes) in named columns.
#' @param type a type of data, 450K or EPIC. Type is used to split CpGs across 
#' chromosomes. Match of CpGs to chromosomes is taken from ChAMPdata package. 
#' If you wish to provide your own match, specify "user" in 
#' this argument and provide a data frame in the next argument. 
#' @param annotation a data frame, user provided match between CpG sites and 
#' chromosomes. Must contain two columns: cpg and chr. Choose "user" in the 
#' previous argument to be able to provide user annotation.
#'
#' @import ChAMPdata
#' @importFrom utils data
#'
#' @return A list of numeric data matrices where each matrix contains probes 
#' from one chromosome.

split_by_chromosomes <- function(data,
                       type = c("450K", "EPIC", "user"),
                       annotation = NULL) {
    
    #load annotation
    type <- match.arg(type)
    if (type == "450K") {
        utils::data("hm450.manifest.hg19", package = "ChAMPdata",
                    envir = environment())
        anno <- data.frame(cpg = rownames(hm450.manifest.hg19), 
                           chr = hm450.manifest.hg19$CpG_chrm)
        remove(hm450.manifest.hg19)
    } else if (type == "EPIC") {
        utils::data("EPIC.manifest.hg19", package = "ChAMPdata",
                    envir = environment())
        anno <- data.frame(cpg = rownames(EPIC.manifest.hg19), 
                           chr = EPIC.manifest.hg19$CpG_chrm)
        remove(EPIC.manifest.hg19)
    } else if (type == "user") {
        anno <- annotation
    } 
    
    #split the data by chromosomes
    chromosomes <- unique(anno$chr)
    nchr <- length(chromosomes)
    data_chr <- vector(mode = "list", length = nchr)
    names(data_chr) <- chromosomes
    probes <- colnames(data)
    for (c in seq_len(nchr)) {
        curr_cpgs <- anno[anno$chr == chromosomes[c], ]$cpg
        data_chr[[c]] <- data[, colnames(data) %in% curr_cpgs, drop = FALSE]
        probes <- probes[! probes %in% curr_cpgs]
    }
    if (length(probes) != 0) {
        #answer <- readline("Some probes are not present in our annotation")
        warning("Some probes are not present in the annotation. They will not be imputed or used for imputation. Most probably, these are probes that have no match to a chromosome in Illumina manifest or are not present in the manifest at all. Consider using a custom annotation if you want to impute and/or use these columns in the imputation.")
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
    
    return(data_chr)
}

#' Impute missing values in methylation dataset
#'
#' This function performs missing value imputation specific for DNA methylation 
#' data. The method is based on linear regression since methylation levels 
#' show a high degree of inter-sample correlation. Implementation is 
#' parallelised over chromosomes to improve the running time.
#'
#' @param input either a numeric data matrix with missing values to be, 
#' with named samples in rows and variables (probes) in named columns, or a 
#' SummarizedExperiment object, with an assay with variables in rows 
#' and samples in columns, as standard.
#' @param which_assay a character specifying the name of assay of the 
#' SummarizedExperiment object to impute. By default the first one will be imputed.
#' @param type a type of data, 450K or EPIC. Type is used to split CpGs across 
#' chromosomes. Match of CpGs to chromosomes is taken from ChAMPdata package. 
#' If you wish to provide your own match, specify "user" in 
#' this argument and provide a data frame in the next argument. 
#' @param annotation a data frame, user provided match between CpG sites and 
#' chromosomes. Must contain two columns: cpg and chr. Choose "user" in the 
#' previous argument to be able to provide user annotation.
#' @param groups a vector of the same length as the number of samples that 
#' identifies what groups does each sample correspond, e.g. \code{c(1, 1, 2, 3)}
#' or \code{c("group1", "group1", "group2", "group3")}. Unique elements of the 
#' vector will be identified as groups and data will be split accordingly. 
#' Imputation will be done for each group separately consecutively. 
#' The default is NULL, so all samples are considered as one group. 
#' @param range a vector of two numbers, \eqn{min} and \eqn{max}, 
#' specifying the range of values in the data. 
#' Since we assume the beta-value representation of the methylation data, 
#' the default range is \eqn{[0, 1]}. 
#' However, if a user wishes to apply the method to the other kind of data, 
#' they can change the range in this argument.
#' @param skip_imputation_ids a numeric vector of ids of the columns with NAs  
#' for which \emph{not} to perform the imputation. If \code{NULL}, all columns 
#' are considered.
#' @param BPPARAM set of options for parallelization through BiocParallel package.
#' For details we refer to their documentation. The one thing most users probably
#' wish to customize is the number of cores. By default it is set to 
#' \eqn{\#cores - 2}. If you wish to change is, supply
#' \code{BBPARAM = SnowParam(workers = ncores)} with your desired \code{ncores}.
#' If the default or user-specified number of workers is higher than number of 
#' chromosomes, it will be overwritten.
#' We also recommend setting \code{exportglobals = FALSE} since it can help reduce
#' running time.
#' @param minibatch_frac a number between 0 and 1, what fraction of samples 
#' to use for mini-batch computation. Remember that if your data has several groups, 
#' mini-batch will be applied to each group separately but with the same fraction, 
#' so choose it accordingly. However, if your chosen fraction will be smaller 
#' than a matrix dimension for some groups, mini-batch will be just ignored. 
#' We advise to use mini-batch only if you have large number of samples, 
#' order of hundreds. The default is 1 (i.e., 100\% of samples are used, 
#' no mini-batch). 
#' @param minibatch_reps a number, how many times to repeat computations with 
#' a fraction of samples specified above (more times -> better performance but 
#' more runtime). The default is 1 (as a companion to default fraction of 100\%,
#' i.e. no mini-batch).
#' @param overwrite_res a boolean specifying whether to overwrite an imputed slot
#' of the SummarizedExperiment object or to add another slot with 
#' imputed data. The default is \code{TRUE} to reduced the object size.
#'
#' @importFrom BiocParallel SnowParam bplapply bpstart bpstop
#' @importFrom parallel detectCores
#' @importFrom SummarizedExperiment assays assay assayNames
#' @importFrom methods is
#'
#' @return Either a numeric matrix with imputed data or a SummarizedExperiment 
#' object.
#'
#' @export
#' 
#' @examples
#' data(beta)
#' beta_with_nas <- generateMissingData(beta, lambda = 3.5)$beta_with_nas
#' beta_imputed <- methyLImp2(input = beta_with_nas, type = "EPIC", 
#'                           minibatch_frac = 0.5, 
#'                           BPPARAM = BiocParallel::SnowParam(workers = 1))

methyLImp2 <- function(input, which_assay = NULL,
                      type = c("450K", "EPIC", "user"),
                      annotation = NULL,
                      groups = NULL,
                      range = NULL,
                      skip_imputation_ids = NULL,
                      BPPARAM = BiocParallel::bpparam(),
                      minibatch_frac = 1, minibatch_reps = 1,
                      overwrite_res = TRUE) {
  
    #check the class of input
    if (is.matrix(input)) {
        if (is.numeric(input)) {
            data_full <- input
            flag <- "matrix"
            if (is.null(rownames(input))) {
                stop("Please provide rownames for the input matrix.")
            }
            if (is.null(colnames(input))) {
                stop("Please provide colnames for the input matrix.")
            }
        } else {
            stop("The input matrix is not numeric.") 
        }
    } else if (is(input, "SummarizedExperiment")) {
        if (is.null(which_assay)) {
            data_full <- t(assay(input))
        } else {
            if (which_assay %in% assayNames(input)) {
                data_full <- t(assay(input, which_assay))
            } else {
                stop("Please provide correct assay name.")
            }
        }
        flag <- "SE"
    } else {
        stop("Input is not a numeric matrix or a SummarizedExperiment object.")
    }
    
    #check annotation
    type <- match.arg(type)
    if (! type %in% c("450K", "EPIC", "user")) {
        stop("Choose type of data correctly to select annotation: 450K or EPIC - or provide your own annotation.")
    }
    if (type == "user") {
        if (is.null(annotation)) {
            stop("Please provide annotation: a data frame with cpg and chr columns.")
        } else if (sum(colnames(data_full) %in% annotation$cpg) == 0) {
            stop ("CpG names in the annotation and in the input data do not match.")
        }
    } 
    
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
    
    #check the parallelization set-up
    {
        #split the data by chromosomes to see the nchr
        data_chr_full <- split_by_chromosomes(data = data_full, type = type,
                                              annotation = annotation)
        nchr <- length(data_chr_full)
        
        if (BPPARAM$tasks != 0 & BPPARAM$tasks != nchr) {
            BPPARAM$tasks <- nchr
            warning("BPPARAM$tasks parameter must be equal to the number of chromosomes. Your input was overwritten.")
        }
        
        if (nchr < BPPARAM$workers) {
            BPPARAM$workers <- nchr
            warning("Number of chromosomes is less than specified BPPARAM$workers. Your input was overwritten since more cores than number of chromosomes is not necessary.")
        }
        
        if (is.null(BPPARAM$RNGseed)) {
            BPPARAM$RNGseed <- 1994
        }
        
        if (BPPARAM$exportglobals) {
            warning("BPPARAM$exportglobals = TRUE. Setting it to FALSE may reduce the running time.")
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
        data_group <- data_full[groups == curr_group, ]
        
        data_chr <- vector(mode = "list", length = nchr)
        names(data_chr) <- names(data_chr_full)
        for (j in seq_len(nchr)) {
            data_chr[[j]] <- data_chr_full[[j]][groups == curr_group, ]
        }
        
        #run methyLImp in parallel for each chromosome
        bpstart(BPPARAM)
        res <- bplapply(data_chr, methyLImp2_internal, 
                         min, max, skip_imputation_ids,
                         minibatch_frac, minibatch_reps,
                         BPPARAM = BPPARAM)
        bpstop(BPPARAM)
        
        names(res) <- names(data_chr)
        
        #if some chromosomes didn't have NAs or didn't have enough data for imputation,
        #tell this to the user and then later skip them in the replacement stage
        problem_chr <- c()
        if (sum(vapply(res, function(x) is.character(x), logical(1))) > 0) {
            problem_chr <- names(data_chr)[vapply(res, function(x) is.character(x), 
                                                  logical(1))]
            for (i in 1:length(problem_chr)) {
                message("Sample group " , curr_group, ", ", 
                        problem_chr[i], ": ", res[problem_chr[i]])
            }
        }
        
        #replace NAs with imputed values
        out <- data_group
        for (c in seq_len(nchr)) {
            if (names(res)[c] %in% problem_chr) next
            imputed <- res[[c]]
            rows_id <- rownames(data_group) %in% rownames(imputed)
            cols_id <- colnames(data_group) %in% colnames(imputed)
            out[rows_id, cols_id] <- imputed
        }
        
        out_full[[i]] <- out
    }
    
    out_full <- do.call(rbind, out_full)
    out_full <- out_full[samples_order, ]
    
    if (flag == "SE") {
        if (overwrite_res) {
            if (is.null(which_assay)) {
                assay(input) <- t(out_full)
            } else {
                assay(input, which_assay) <- t(out_full)
            }    
            return(input)
        } else {
            assay(input, paste(which_assay, "_imputed", sep = "_")) <- 
                t(out_full)
            return(input)
        }
    } else {
        return(out_full)
    }
    
}

#' Impute missing values in methylation dataset
#'
#' @param dat a numeric data matrix with missing values, 
#' with samples in rows and variables (probes) in columns.
#' @param min a number, minimum value for bounded-range variables. 
#' Default is 0 (we assume beta-value representation of the methylation data). 
#' Can be user provided in case of other types of data.
#' @param max a number, maximum value for bounded-range variables. 
#' Default is 1 (we assume beta-value representation of the methylation data). 
#' Can be user provided in case of other types of data.
#' @param skip_imputation_ids a numeric vector of ids of the columns with NAs for which 
#' \emph{not} to perform the imputation. If \code{NULL}, all columns are considered.
#' @param minibatch_frac a number, what percentage of samples to use for 
#' mini-batch computation. The default is 1 (i.e., 100\% of samples are used, 
#' no mini-batch).
#' @param minibatch_reps a number, how many times repeat computations with a 
#' fraction of samples (more times - better performance). 
#' The default is 1 (as a companion to default fraction of 100\%. i.e. no mini-batch).
#'
#' @return A numeric matrix \eqn{out} with imputed data is returned.
#' 
#' @keywords internal

methyLImp2_internal <- function(dat,
                                min, max,
                                skip_imputation_ids,
                                minibatch_frac, minibatch_reps) {
    
    #let's identify columns with same patterns of NAs
    {
        colnames_dat <- colnames(dat)
        
        #detect NAs
        dat_na <- is.na(dat)
        #exclude columns with no NAs
        dat_na <- dat_na[, colSums(dat_na) > 0, drop = FALSE]
        if (dim(dat_na)[2] == 0) {
            return("No columns with missing values detected.")
        }
        #save all columns with NA
        all_NA_cols <- which(colnames_dat %in% colnames(dat_na))
        # exclude from the imputation columns with all NAs or 
        # a single not NA value: not enough information for imputation
        dat_na <- dat_na[, colSums(dat_na) < (dim(dat_na)[1] - 1), drop = FALSE]
        
        #If all the columns have missing values we cannot do anything
        if (dim(dat_na)[2] == dim(dat)[2]) {
            return("No CpGs without any missing values, not possible to conduct imputation.")
        } 
        
        unique_patterns <- unique(t(dat_na))
        npatterns <- dim(unique_patterns)[1]
        
        ids <- vector(mode = "list", length = npatterns)
        
        for (i in seq_len(npatterns)) {
            curr_pattern <- unique_patterns[i, ]
            
            col_match <- apply(dat_na, 2, function(x) identical(x, curr_pattern))
            col_match <- colnames(dat_na)[col_match]
            NAcols_id <- which(colnames_dat %in% col_match)
            #if some of the chosen columns are restricted by user, exclude them
            if(!is.null(skip_imputation_ids)) {
                NAcols_id <- NAcols_id[!(NAcols_id %in% skip_imputation_ids)]
            }
            
            row_id <- which(curr_pattern == TRUE)
            
            ids[[i]] <- list(row_id = row_id, NAcols_id = NAcols_id)
        }
        
        names(ids) <- paste("group", seq_len(npatterns), sep = "_")
        
    }
    
    out <- dat
    for (j in seq_len(npatterns)) {
        row_id <- ids[[j]]$row_id
        NAcols_id <- ids[[j]]$NAcols_id
        
        C <- dat[row_id, -all_NA_cols]
        
        A_full <- dat[-row_id, -all_NA_cols]
        B_full <- dat[-row_id, NAcols_id]
        
        imputed_list <- vector(mode = "list", length = minibatch_reps)
        for (r in seq_len(minibatch_reps)) {
            sample_size <- ifelse(dim(A_full)[1] > ceiling(dim(A_full)[1] * 
                                                               minibatch_frac),
                                  ceiling(dim(A_full)[1] * minibatch_frac),
                                  dim(A_full)[1])
            chosen_rows <- sort(sample(seq_len(dim(A_full)[1]), 
                                       size = sample_size))
            A <- A_full[chosen_rows, , drop = FALSE]
            if (is.null(dim(B_full))) {
                B <- B_full[chosen_rows]
            } else {
                B <- B_full[chosen_rows, ]
            }
            
            # Updates or computes max.sv from A. Negative or zero value not allowed
            max.sv <- NULL
            max.sv <- max(ifelse(is.null(max.sv), min(dim(A)), max.sv), 1)
            
            if(is.null(min) || is.null(max)) {
                # Unrestricted-range imputation
                # X <- pinvr(A, rank) %*% B (X = A^-1 * B)
                # O <- C %*% X             (O = C*X)
                imputed_list[[r]] <- C %*% (methyLImp2:::pinvr(A, max.sv) %*% B)
            } else {
                # Bounde-range imputation
                # X <- pinvr(A, rank) %*% logit(B, min, max) (X = A^-1 * logit(B))
                # P <- inv.logit(C %*% X, min, max)         (P = logit^-1 (C * X))
                imputed_list[[r]] <- methyLImp2:::inv.plogit(C %*% 
                                                                 (methyLImp2:::pinvr(A, max.sv) %*%
                                                                      methyLImp2:::plogit(B, min, max)), min, max)
            }
        }
        
        imputed <- Reduce('+', imputed_list) / minibatch_reps
        out[row_id, NAcols_id] <- imputed
    }
    
    return(out)
    
}
