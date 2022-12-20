#' Impute missing values in methylation dataset
#'
#' This function performs missing value imputation specific for DNA methylation data. The method is based on linear regression since methylation levels show a high degree of inter-sample correlation. Implementation is parallelised over chromosomes since probes on different chromosomes are usually independent.
#'
#' @param data a numeric data matrix with missing values, with samples in rows and variables (probes) in columns.
#' @param type a type of data. It is used to provide correct annotation: 450K, EPIC or user-provided. If the latter is chosen, annotation must be provided in the next argument.
#' @param annotation a data frame, user provided annotation. Must contain two columns: cpg and chr, and provide a match between CpG sites and chromosomes. Choose "user" in the previous argument to be able to provide user annotation.
#' @param min a number, minimum value for bounded-range variables. Default: 0 (we assume beta-value representation of the methylation data). Unrestricted range if \code{min} or \code{max = NULL}.
#' @param max a number, maximum value for bounded-range variables. Default: 1 (we assume beta-value representation of methylation data). Unrestricted range if \code{min} or \code{max = NULL}.
#' @param max.sv a number, maximum of a singular value to be used in the pseudoinverse. If \code{NULL}, use all singular values.
#' @param col.list a numeric vector, restricts the imputation only to the specified columns. If \code{NULL}, all columns are considered.
#' @param ncores number of cores to use in parallel computation. If \code{NULL}, set to \code{#physical cores} - 1.
#' @param minibatch_frac a number, what percentage of samples to use for mini-batch computation. The default is 1 (i.e., 100% of samples are used, no mini-batch).
#' @param minibatch_reps a number, how many times repeat computations with a fraction of samples (more times - better performance). The default is 1 (as a companion to default fraction of 100%. i.e. no mini-batch).
#'
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport clusterApply stopCluster
#' @importFrom dplyr distinct
#'
#' @return A numeric matrix \code{out} with imputed data is returned.
#'
#' @export

methyLImp <- function(data,
                      type = c("450K", "EPIC", "user"),
                      annotation = NULL,
                      min = 0, max = 1, max.sv = NULL,
                      col.list = NULL,
                      ncores = NULL,
                      minibatch_frac = 1, minibatch_reps = 1) {

  #first, split NA imputed data by chromosomes
  type <- match.arg(type)
  if (type == "450K") {
    anno <- anno450K
  } else if (type == "EPIC") {
    anno <- annoEPIC
  } else if (type == "user") {
    if (!is.null(annotation)) {
      anno <- annotation
    } else stop("Please provide annotation: two-column data.frame with CpGs matching chromosomes.")
  } else stop("Choose type of data correctly to select annotation: 450K or EPIC - or provide your annotation.")

  chromosomes <- unique(anno$chr)
  nchr <- length(chromosomes)
  data_chr <- vector(mode = "list", length = nchr)
  names(data_chr) <- chromosomes
  for (c in 1:nchr) {
    curr_cpgs <- anno[anno$chr == chromosomes[c], ]$cpg
    data_chr[[c]] <- data[, colnames(data) %in% curr_cpgs]
  }

  #run methyLImp in parallel for each chromosome
  if (is.null(ncores)) {
    ncores <- detectCores(logical = FALSE) - 1
  }
  cl <- makeCluster(ncores)
  clusterEvalQ(cl, library(methyLImp))
  clusterEvalQ(cl, library(dplyr))
  clusterExport(cl, c("min", "max", "max.sv", "col.list",
                      "minibatch_frac", "minibatch_reps"))
  res <- parallel::clusterApplyLB(cl, data_chr, methyLImp_internal)

  out <- data
  for (c in 1:nchr) {
    imputed <- res[[c]]
    rows_id <- rownames(data) %in% rownames(imputed)
    cols_id <- colnames(data) %in% colnames(imputed)

    out[rows_id, cols_id] <- imputed
  }

  stopCluster(cl)

  return(out)

}
