#' Pseudo-logit function with range restricted to \eqn{[-X, X]} where \eqn{X} depends on the double machine precision.
#'
#' @param x a numeric matrix
#' @param min a number, default value is 0
#' @param max a number, default value is 1
#'
#' @return A numeric matrix 
#'
#' @keywords internal

plogit <- function(x, min = 0, max = 1)
{
  p <- (x - min) / (max - min)
  # fix -Inf
  p <- ifelse(p < .Machine$double.neg.eps, .Machine$double.neg.eps, p)
  # fix +Inf
  p <- ifelse(p > 1 - .Machine$double.neg.eps, 1 - .Machine$double.neg.eps, p)
  log(p / (1 - p))
}


#' Inverse of the pseudo-logit function.
#'
#' @param x a numeric matrix
#' @param min a number, default value is 0
#' @param max a number, default value is 1
#'
#' @return A numeric matrix 
#' 
#' @keywords internal

inv.plogit <- function(x, min = 0, max = 1)
{
  p <- exp(x) / (1 + exp(x))
  # fix problems with +Inf
  p <- ifelse(is.na(p) & !is.na(x), 1, p )
  # fix 0 rounding
  p <- ifelse(p <= exp(plogit(0)) / (1 + exp(plogit(0))), 0, p)
  p * (max - min) + min
}

#' Computes the Moore-Penrose generalized inverse of a matrix. Allows rank reduction of the generalized inverse.
# '
#' This function is directly taken from MASS package (code on GPLv3 license) and modified in order to include the rank reduction option. The added code for rank reduction is commented in the implementation.
#'
#' @param X a numeric matrix
#' @param tol a number, default value is \eqn{sqrt(.Machine$double.eps)}.
#'
#' @return A numeric matrix 
#'
#' @keywords internal

pinvr <- function(X, max.sv = min(dim(X)), tol = sqrt(.Machine$double.eps))
{
  # based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
  
  if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if(!is.matrix(X)) X <- as.matrix(X)
  Xsvd <- svd(X)
  if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)

  # Rank reduction extension: START
  max.sv <- min(ifelse(max.sv < 0, 1, max.sv), min(dim(X)))
  L <- logical(length(Positive))
  L[seq_len(max.sv)] <- TRUE
  Positive <- Positive & L
  # Rank reduction extension: END

  if (all(Positive)) Xsvd$v %*% (1 / Xsvd$d * t(Xsvd$u))
  else if(!any(Positive)) array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1 / Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}