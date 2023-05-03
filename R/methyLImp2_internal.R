#' Impute missing values in methylation dataset
#'
#' @param dat a numeric data matrix with missing values, with samples in rows and variables (probes) in columns.
#' @param min a number, minimum value for bounded-range variables. Default is 0 (we assume beta-value representation of the methylation data). Can be user provided in case of other types of data.
#' @param max a number, maximum value for bounded-range variables. Default is 1 (we assume beta-value representation of the methylation data). Can be user provided in case of other types of data.
#' @param col.list a numeric vector, restricts the imputation only to the specified columns. If \code{NULL}, all columns are considered.
#' @param minibatch_frac a number, what percentage of samples to use for mini-batch computation. The default is 1 (i.e., 100\% of samples are used, no mini-batch).
#' @param minibatch_reps a number, how many times repeat computations with a fraction of samples (more times - better performance). The default is 1 (as a companion to default fraction of 100\%. i.e. no mini-batch).
#'
#' @importFrom dplyr distinct
#'
#' @return A numeric matrix \eqn{out} with imputed data is returned.
#' 
#' @keywords internal

methyLImp2_internal <- function(dat,
                               min, max,
                               col.list,
                               minibatch_frac, minibatch_reps,
                               max.sv = NULL) {

  #let's identify columns with same patterns of NAs
  {
    colnames_dat <- colnames(dat)

    #detect NAs
    dat_na <- is.na(dat)
    #exclude columns with no NAs
    dat_na <- dat_na[, colSums(dat_na) > 0]
    if (dim(dat_na)[2] == 0) {
      return("No columns with missing values detected.")
    }
    #save all columns with NA
    all_NA_cols <- which(colnames_dat %in% colnames(dat_na))
    # exclude from the imputation columns with all NAs or 
    # a single not NA value: not enough information for imputation
    dat_na <- dat_na[, colSums(dat_na) < (dim(dat_na)[1] - 1)]

    #If all the columns have missing values we cannot do anything
    if (dim(dat_na)[2] == dim(dat)[2]) {
      return("Not enough data without missing values to conduct imputation.")
    } else {
      cat("#cols with #NAs < (#samples - 1):", dim(dat_na)[2], "\n")
    }

    #if some columns are restricted by user, exclude them
    #if(!is.null(col.list)) dat_na <- dat_na[, !col.list]

    unique_patterns <- as.matrix(distinct(as.data.frame(t(dat_na))))
    ngroups <- dim(unique_patterns)[1]
    cat("#regression groups:", ngroups, "\n")

    ids <- vector(mode = "list", length = ngroups)

    for (i in 1:ngroups) {
      curr_pattern <- unique_patterns[i, ]

      col_match <- apply(dat_na, 2, function(x) identical(x, curr_pattern))
      col_match <- colnames(dat_na)[col_match]
      NAcols_rowid <- which(colnames_dat %in% col_match)

      row_id <- which(curr_pattern == TRUE)

      ids[[i]] <- list(row_id = row_id, NAcols_rowid = NAcols_rowid)
    }

    names(ids) <- paste("group", c(1:ngroups), sep = "_")

  }

  out <- dat
  for (i in 1:ngroups) {
    print(i)
    row_id <- ids[[i]]$row_id
    NAcols_rowid <- ids[[i]]$NAcols_rowid

    C <- dat[row_id, -all_NA_cols]

    A_full <- dat[-row_id, -all_NA_cols]
    B_full <- dat[-row_id, NAcols_rowid]

    imputed_list <- vector(mode = "list", length = minibatch_reps)
    for (r in 1:minibatch_reps) {
      sample_size <- ifelse(dim(A_full)[1] > ceiling(dim(A_full)[1] * minibatch_frac),
                            ceiling(dim(A_full)[1] * minibatch_frac),
                            dim(A_full)[1])
      chosen_rows <- sort(sample(1:dim(A_full)[1], size = sample_size))
      A <- A_full[chosen_rows, , drop = FALSE]
      if (is.null(dim(B_full))) {
        B <- B_full[chosen_rows]
      } else {
        B <- B_full[chosen_rows, ]
      }

      # Updates or computes max.sv from A. Negative or zero value not allowed
      max.sv <- max(ifelse(is.null(max.sv), min(dim(A)), max.sv), 1)
    
      if(is.null(min) || is.null(max)) {
        # Unrestricted-range imputation
        # X <- pinvr(A, rank) %*% B (X = A^-1 * B)
        # O <- C %*% X             (O = C*X)
        imputed_list[[r]] <- C %*% (pinvr(A, max.sv) %*% B)
      } else {
        # Bounde-range imputation
        # X <- pinvr(A, rank) %*% logit(B, min, max) (X = A^-1 * logit(B))
        # P <- inv.logit(C %*% X, min, max)         (P = logit^-1 (C * X))
        imputed_list[[r]] <- inv.plogit(C %*% (pinvr(A, max.sv) %*%
                                       plogit(B, min, max)), min, max)
      }
    }

    imputed <- Reduce('+', imputed_list) / minibatch_reps
    #test2[[i]] <- imputed
    out[row_id, NAcols_rowid] <- imputed
  }

   return(out)
  
}
