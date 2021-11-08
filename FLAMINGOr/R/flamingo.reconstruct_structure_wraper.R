#' flamingo.reconstruct_structure_wraper
#'
#' A wraper of the low-rank matrix completion algorithm to process all the domains
#' @param index Index of a domain
#' @param sw Downsample rates, suggesting the fraction of the observed data to be used in the model.
#' @param lambda The penalty term for the diagnal entries.
#' @param max_dist The maximum allowed distance for the diagnal entries.
#' @param nThread Number of thread to be used in the model, default is 28.
#' @keywords FLAMINGO
#' @return A list of flamingo_prediction object
#' @export

flamingo.reconstruct_structure_wraper<- function(index,sw,lambda,max_dist,nThread=28){
  pd <- read.table(paste0("./Domain_data/PD_domain_",index,'.txt'))
  input_if <- read.table(paste0("./Domain_data/IF_domain_",index,'.txt'))
  file <- dir('./Domain_data/')
  if(check_data_availability(pd)){
    res =  flamingo.reconstruct_structure_worker(input_if,pd,sw,lambda,max_dist,nThread=1)
  }else{
    res = NULL
  }
  return(res)
}