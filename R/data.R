#' A subset of GSE199057 dataset for vignette demonstration
#'
#' The GSE199057 Gene Expression Omnibus dataset contains 68 mucosa samples 
#' from non-colon-cancer patients, from which we randomly sampled 24. 
#' Methylation data were measured on EPIC arrays and after removal of 
#' sex chromosomes and SNPs loci, it contains 816 126 probes. 
#' Pre-processing can be found on the _methyLImp2_ simulation github page 
#' https://github.com/annaplaksienko/methyLImp2_simulation_studies.
#' Here we subset only a quarter of probes from two smallest chromosomes (18 and 21) 
#' for the sake of demonstration. 
#'
#' @docType data
#'
#' @usage data(beta)
#'
#' @format A numeric matrix
#' @return A numeric matrix
#'
#' @keywords datasets
#'
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
"beta"

#' Metadata information for GSE199057 dataset for vignette demonstration
#'
#' The GSE199057 Gene Expression Omnibus dataset contains 68 mucosa samples 
#' from non-colon-cancer patients, from which we randomly sampled 24. 
#' Dataset contains metadata for those.
#'
#' @docType data
#'
#' @usage data(beta_meta)
#'
#' @format A data.frame
#' @return A data.frame.
#'
#' @keywords datasets
#'
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
"beta_meta"
