#' flamingo.reconstruct_backbone_structure
#'
#' A wraper of the low-rank matrix completion algorithm to reconstruct the backbone structure
#' @param flamingo_data_obj A flamingo object containing the IF and PD matrix of the backbone, e.g 1MB resolution.
#' @param sw Downsample rates, suggesting the fraction of the observed data to be used in the model.
#' @param lambda The penalty term for the diagnal entries.
#' @param max_dist The maximum allowed distance for the diagnal entries.
#' @param nThread Number of thread to be used in the model, default is 28.
#' @keywords FLAMINGO
#' @return A flamingo_prediction object containing the 3D structure of the backbone.
#' @export
flamingo.reconstruct_backbone_structure<- function(flamingo_data_obj,sw,lambda,max_dist,nThread=28){
  res = flamingo.reconstruct_structure_worker(flamingo_data_obj@IF,flamingo_data_obj@PD,sw,lambda,max_dist,nThread)
  return(res)
}
