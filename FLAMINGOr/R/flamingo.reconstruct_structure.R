#' flamingo.reconstruct_structure
#'
#' A wraper of the low-rank matrix completion algorithm to process all the domains
#' @param sw Downsample rates, suggesting the fraction of the observed data to be used in the model.
#' @param lambda The penalty term for the diagnal entries.
#' @param max_dist The maximum allowed distance for the diagnal entries.
#' @param nThread Number of thread to be used in the model, default is 28.
#' @keywords FLAMINGO
#' @return A list of flamingo_prediction object
#' @export
flamingo.reconstruct_structure<- function(sw,lambda,max_dist,nThread=28){
  cl <- parallel::makeCluster(nThread)
  parallel::clusterCall(cl, function() library(FLAMINGOr))
  file <- dir('./Domain_data/')
  n <- length(grep('PD',file))
  res = list()
  domain_id <- c()
  parallel::clusterExport(cl,c("sw","lambda","max_dist"),envir=environment())
  worker_res = parallel::parSapply(cl,1:n,function(x){
    flamingo.reconstruct_structure_wraper(x,sw,lambda,max_dist,nThread=1)
  })
  for(i in 1:length(worker_res)){
    if(!is.null(worker_res[[i]])){
      res[[i]] <- worker_res[[i]]
      domain_id <- c(domain_id,i)
    }
  }
  res <- res[lengths(res)>0]
  names(res) <- domain_id
  parallel::stopCluster(cl)
  return(res)
}
