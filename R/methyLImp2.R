#' Impute missing values in methylation dataset
#'
#' This function performs missing value imputation specific for DNA methylation data. The method is based on linear regression since methylation levels show a high degree of inter-sample correlation. Implementation is parallelised over chromosomes since probes on different chromosomes are usually independent.
#'
#' @param data a numeric data matrix with missing values, with samples in rows and variables (probes) in columns.
#' @param type a type of data. It is used to provide correct annotation: 450K, EPIC or user-provided. If the latter is chosen, annotation must be provided in the next argument.
#' @param annotation a data frame, user provided annotation. Must contain two columns: cpg and chr, and provide a match between CpG sites and chromosomes. Choose "user" in the previous argument to be able to provide user annotation.
#' @param range a vector of two numbers, \eqn{min} and \eqn{max}, specifying the range of values in the data. Since we assume the beta-value representation of the methylation data, the default range is \eqn{[0, 1]}. However, if a user wishes to apply the method to other kind of data, they can change the range in this argument.
#' @param col.list a numeric vector, restricts the imputation only to the specified columns. If \code{NULL}, all columns are considered.
#' @param ncores number of cores to use in parallel computation. If \code{NULL}, set to \eqn{\#physical cores - 1}.
#' @param minibatch_frac a number, what percentage of samples to use for mini-batch computation. The default is 1 (i.e., 100\% of samples are used, no mini-batch).
#' @param minibatch_reps a number, how many times repeat computations with a fraction of samples (more times -> better performance). The default is 1 (as a companion to default fraction of 100\%. i.e. no mini-batch).
#'
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport clusterApply stopCluster
#' @importFrom dplyr distinct
#'
#' @return A numeric matrix \eqn{out} with imputed data is returned.
#'
#' @export

methyLImp2 <- function(data,
                      type = c("450K", "EPIC", "user"),
                      annotation = NULL,
                      range = NULL,
                      col.list = NULL,
                      ncores = NULL,
                      minibatch_frac = 1, minibatch_reps = 1) {
  
  #check if the input data is of correct format
  if(!is.matrix(data)) {
    stop("The data object is not a matrix.")
  } else if (!is.numeric(data)) {
    stop("The data matrix is not numeric.")
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
        stop("Left interval must be less than the right one.")
      } else {
        min <- range[1]
        max <- range[2]
      }
    }
  }
  
  #split the data by chromosomes
  type <- match.arg(type)
  if (type == "450K") {
    anno <- anno450K
  } else if (type == "EPIC") {
    anno <- annoEPIC
  } else if (type == "user") {
    if (!is.null(annotation)) {
      anno <- annotation
    } else stop("Please provide annotation: a data frame with cpg and chr columns.")
  } else stop("Choose type of data correctly to select annotation: 450K or EPIC - or provide your annotation.")

  chromosomes <- unique(anno$chr)
  nchr <- length(chromosomes)
  data_chr <- vector(mode = "list", length = nchr)
  names(data_chr) <- chromosomes
  for (c in 1:nchr) {
    curr_cpgs <- anno[anno$chr == chromosomes[c], ]$cpg
    data_chr[[c]] <- data[, colnames(data) %in% curr_cpgs]
  }
  #drop chromosomes that were not present in the the input dataset (if such exist)
  data_chr <- data_chr[names(data_chr)[sapply(data_chr, 
                                              function(x) dim(x)[2]) != 0]]
  nchr <- length(data_chr)
  #reorder chromosomes from biggest to smallest
  chr_order <- names(sort(sapply(data_chr, function(x) dim(x)[2]), 
                          decreasing = TRUE))
  data_chr <- data_chr[chr_order]

  #run methyLImp in parallel for each chromosome
  if (is.null(ncores)) {
    ncores <- detectCores(logical = FALSE) - 1
  }
  cl <- makeCluster(ncores)
  clusterEvalQ(cl, library(methyLImp2))
  clusterEvalQ(cl, library(dplyr))
  res <- parallel::clusterApplyLB(cl, data_chr, methyLImp2_internal,
                                  min, max, col.list,
                                  minibatch_frac, minibatch_reps)
  stopCluster(cl)
  names(res) <- names(data_chr)
  
  #if some chromosomes didn't have NAs or didn't have enough data for imputation,
  #tell this to the user and then later skip them in the replacement stage
  problem_chr <- c()
  if (sum(sapply(res, function(x) is.character(x))) > 0) {
    problem_chr <- names(data_chr)[sapply(res, function(x) is.character(x))]
    for (i in length(problem_chr)) {
      message(paste0("For chromosome ", problem_chr[i], ": ", 
                     res[problem_chr[i]]))
    }
  }

  #replace NAs with imputed values
  out <- data
  for (c in 1:nchr) {
    if (names(res)[c] %in% problem_chr) next
    imputed <- res[[c]]
    rows_id <- rownames(data) %in% rownames(imputed)
    cols_id <- colnames(data) %in% colnames(imputed)

    out[rows_id, cols_id] <- imputed
  }

  return(out)

}
