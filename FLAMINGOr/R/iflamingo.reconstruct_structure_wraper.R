#' iflamingo.reconstruct_structure_wraper
#'
#' A wraper of the low-rank matrix completion algorithm to process all the domains using iFLAMINGO
#' @param index Index of a domain
#' @param sw Downsample rates, suggesting the fraction of the observed data to be used in the model.
#' @param lambda The penalty term for the diagnal entries.
#' @param lambda_epi The penalty term for the epigenomic term.
#' @param max_dist The maximum allowed distance for the diagnal entries.
#' @param nThread Number of thread to be used in the model, default is 28.
#' @keywords FLAMINGO
#' @return A list of flamingo_prediction object
#' @export

iflamingo.reconstruct_structure_wraper<- function(index,sw,lambda,lambda_epi,max_dist,nThread=28){
  pd <- read.table(paste0("./Domain_data/PD_domain_",index,'.txt'))
  input_if <- read.table(paste0("./Domain_data/IF_domain_",index,'.txt'))
  epi_feature <- read.table(paste0("./Genomic_loc/3D_Dist_impute_Epi_",index,'.txt'))
  file <- dir('./Domain_data/')
  if(check_data_availability(pd)){
    res =  iflamingo.reconstruct_structure_worker(input_if,pd,epi_feature,sw,lambda,lambda_epi,max_dist,nThread=1)
  }else{
    res = NULL
  }
  return(res)
}
