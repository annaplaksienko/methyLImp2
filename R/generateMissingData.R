#' Generation of artifical missing values
#'
#' This function generates missing values for the simulation purposes (to apply _methyLImp_ method and then compare the imputed values with the true ones that have been replaced by NAs). First, we randomly choose 3\% of all probes. Then for each of the chosen probes, we randomly define the number of NAs from a Poisson distribution with \eqn{\lambda}, appropriate to the sample size of the dataset (unless specified by the user, here we use \eqn{lambda = 0.15 * \#samples + 0.2}). Finally, these amount of NAs is randomly placed among the samples.
#'
#' @param beta a numeric data matrix into which one wants to add some missing values
#' @param lambda a number, parameter of the Poisson distribution that will indicated how many samples will have missing values in each selected probe.
#' 
#' @return A numeric data matrix with generated NAs in some entries.
#' 
#' @importFrom stats rpois
#'
#' @export

generateMissingData <- function(beta, lambda = NULL) {
  
  #how many NAs do we have to start with?
  na_start <- sum(is.na(beta))
  message(paste("The input dataset has", na_start, "missing entries."))
  
  #let's save number of samples and number of probes
  nsamples <- dim(beta)[1]
  nprobes <- dim(beta)[2]
  
  #how many probes will have NAs?
  #according to the literature, around 3%
  nna <- floor(nprobes * 0.03)
  #which probes will have NAs?
  na_cols <- sample(1:nprobes, size = nna)
  
  #make a list where each element corresponds to a probe and it will store which samples are NA for that probe
  na_values <- vector(mode = "list", length = nna)
  names(na_values) <- na_cols
  
  #we assume that number of NAs in the chosen probes follows Poisson distrubution.
  #if lambda is not provided by the user, we choose it by lambda = 0.15 * #samples - 0.2.
  #it's a linear model made to fit the parameters we have chosen in the simulation, which were
  # 1 for 9 samples
  # 2.5 for 17 samples
  # 5 for 34 samples
  # 7.5 for 51 samples
  # 10 for 68 samples
  if (is.null(lambda)) {
    lambda <- 0.15 * nsamples - 0.2
  }
  
  beta_with_nas <- beta
  
  for (i in 1:nna) {
    #save column id
    na_values[[i]]$na_col <- na_cols[i]
    
    #for each probe, choose how many NAs it will have. 
    #If the number is 0, choose again
    curr_nna <- 0
    while (curr_nna == 0) {
      curr_nna <- rpois(n = 1, lambda = lambda)
    }
    #now sample, which samples are NA in the current probe
    na_values[[i]]$na_rows <- sample(1:nsamples, size = curr_nna)
    #if those probes already have NAs (coming from original dataset), sample again
    while (sum(is.na(beta_with_nas[na_values[[i]]$na_rows, na_cols[i]])) != 0) {
      na_values[[i]]$na_rows <- sample(1:nsamples, size = curr_nna)
    }
    #assign artificial NA
    beta_with_nas[na_values[[i]]$na_rows, na_cols[i]] <- NA
  }
  
  #what's the number of NAs now?
  na_imputed <- sum(is.na(beta_with_nas))
  message(paste("After NA generation, the dataset has", na_imputed, "missing entries."))
  
  return(list("beta_with_nas" = beta_with_nas,
              "na_positions" = na_values))
}