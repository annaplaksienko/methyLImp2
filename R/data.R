#' A subset of GSE199057 dataset for vignette demonstration
#'
#' The GSE199057 Gene Expression Omnibus dataset contains 68 mucosa samples 
#' from non-colon-cancer patients, from which we randomly sampled 24. 
#' Methylation data were measured on EPIC arrays and after removal of 
#' sex chromosomes and SNPs loci, it contains 816 126 probes. 
#' Pre-processing can be found on the _methyLImp2_ simulation github page 
#' https://github.com/annaplaksienko/methyLImp2_simulation_studies.
#' Here we subset only probes from two chromosomes (18 and 21) 
#' to reduce the size for the sake of demonstration. 
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


#' Annotation for 450K
#'
#' Data matches each probe of Illumina 450K array to a chromosome. 
#' Information is from https://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html -> 
#' Infinium HumanMethylation450K v1.2 Product Files -> 
#' HumanMethylation450 v1.2 Manifest File (CSV Format). 
#' Pre-processing can be found on the _methyLImp2_ simulation github page 
#' https://github.com/annaplaksienko/methyLImp2_simulation_studies.
#'
#' @docType data
#'
#' @usage data(anno450K)
#'
#' @format A data frame.
#' @return A data frame.
#'
#' @keywords datasets
#'
#' @source https://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html
"anno450K"


#' Annotation for EPIC array
#'
#' Data matches each probe of Illumina EPIC array to a chromosome. 
#' Information is from https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html -> 
#' Infinium MethylationEPIC v1.0 Product Files (Legacy BeadChip) -> 
#' Infinium MethylationEPIC v1.0 B5 Manifest File (CSV Format). 
#' Pre-processing can be found on the _methyLImp2_ simulation github page 
#' https://github.com/annaplaksienko/methyLImp2_simulation_studies.
#'
#' @docType data
#'
#' @usage data(annoEPIC)
#'
#' @format A data frame
#' @return A data frame
#'
#' @keywords datasets
#'
#' @source https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html
"annoEPIC"
