#' iflamingo.reconstruct_structure
#'
#' A wraper of the low-rank matrix completion algorithm to process all the domains using iflamingo
#' @param sw Downsample rates, suggesting the fraction of the observed data to be used in the model.
#' @param lambda The penalty term for the diagnal entries.
#' @param max_dist The maximum allowed distance for the diagnal entries.
#' @param mThread Number of thread to be used in the model, default is 28.
#' @keywords FLAMINGO
#' @return A list of flamingo_prediction object
#' @examples
#' iflamingo.reconstruct_structure(0.75,10,0.01)

iflamingo.reconstruct_structure<- function(sw,lambda,lambda_epi,max_dist,nThread=28){
  library(parallel)
  cl <- makeCluster(nThread)
  clusterCall(cl, function() devtools::load_all('/mnt/home/wangha73/FLAMINGO'))
  file <- dir('./Domain_data/')
  n <- length(grep('PD',file))
  res = list()
  domain_id <- c()
  clusterExport(cl,c("sw","lambda","max_dist"),envir=environment())
  worker_res = parSapply(cl,1:n,function(x){
    iflamingo.reconstruct_structure_wraper(x,sw,lambda,max_dist,nThread)
  })
  for(i in 1:length(worker_res)){
    if(!is.null(worker_res[[i]])){
      res[[i]] <- worker_res[[i]]
      domain_id <- c(domain_id,i)
    }
  }
  res <- res[lengths(res)>0]
  names(res) <- domain_id
  return(res)
}