#' flamingo.reconstruct_structure_worker
#'
#' Reconstruct the 3D genome structure using low-rank tensor completion
#' @param input_if Input interaction frequency matrix.
#' @param pd Input pairwise distance matrix.
#' @param sw Downsample rates, suggesting the fraction of the observed data to be used in the model.
#' @param lambda The penalty term for the diagnal entries.
#' @param max_dist The maximum allowed distance for the diagnal entries.
#' @param nThread Number of thread to be used in the model, default is 28.
#' @keywords FLAMINGO
#' @return A flamingo_prediction object containing the fragment id and 3D coordinates
#' @export
flamingo.reconstruct_structure_worker <- function(input_if,pd,sw,lambda,max_dist,nThread=2){
  library(parallel)
  library(mgcv)
  library(Matrix)
  # process input data
  input_if <- as.matrix(input_if)
  input_if[which(is.na(input_if))] <- NA
  pd <- as.matrix(pd)
  pd <- pd^2
  pd[which(pd==Inf|is.na(pd))] <- 3
  diag(pd) <- 0
  rm_id <- which(apply(input_if,1,max)==0)
  n <- dim(pd)[1]
  # identify the gram matrix of the input data
  H = Diagonal(x=rep(1,n))-1/n*rep(1,n)%*%t(rep(1,n))
  M = -1/2*H%*%pd%*%H
  M <- as.matrix(M)
  # define measurement set omega
  cl <- makeCluster(nThread)
  input_if <- as.matrix(input_if)
  input_if <- as(input_if,'sparseMatrix')
  if_th <- 0
  col_j = findInterval(seq(input_if@x)-1,input_if@p[-1])+1
  row_i = input_if@i+1
  x_x = input_if@x
  df = data.frame(row_i,col_j,x_x)
  omega <- subset(df[,1:2],df[,3] > if_th & df[,1] != df[,2])
  omega <- as.matrix(omega)
  n_omega <- dim(omega)[1]
  diag_term <- which(omega[,2]-omega[,1]==1)
  omega_diag = omega[diag_term,]
  omega <<- omega[unique(c(diag_term,sample(1:n_omega,sw*n_omega))),]
  n_omega <- dim(omega)[1]
  if(length(diag_term)==1){
    n_omega_diag=1
    omega_diag = matrix(omega_diag,ncol=2)
  }
  else if(length(diag_term)==0){
    return(NULL)
  }else{
    n_omega_diag <- dim(omega_diag)[1]
  }
  # prepare for the optimization
  clusterExport(cl,"omega",envir=environment())
  clusterExport(cl,"omega_diag",envir=environment())
  clusterExport(cl,"convert_index",envir=environment())
  w_list<- parLapply(cl,1:n_omega,function(x){cbind(convert_index(omega[x,]),x)})
  func_list <- list()
  tmp_w_list <- do.call(rbind,w_list)
  all_elements <- unique(data.frame(tmp_w_list[,1:2]))
  loc <- prodlim::row.match(data.frame(tmp_w_list[,1:2]),all_elements)
  func_list <- split(data.frame(tmp_w_list[,c(3,4)]),loc)
  clusterExport(cl,"func_list",envir=environment())
  clusterExport(cl,"omega_diag",envir=environment())
  w_diag_list<- parLapply(cl,1:n_omega_diag,function(x){cbind(convert_index(omega_diag[x,]),x)})
  func_list_diag <- list()
  tmp_w_list <- do.call(rbind,w_diag_list)
  all_elements_diag <- unique(data.frame(tmp_w_list[,1:2]))
  loc <- prodlim::row.match(data.frame(tmp_w_list[,1:2]),all_elements_diag)
  func_list_diag <- split(data.frame(tmp_w_list[,c(3,4)]),loc)
  clusterExport(cl,"func_list_diag",envir=environment())
  b = cal_proj(omega,M)
  d = cal_proj_diag(omega_diag,M)
  for(i in 1:length(d)){
    d[i] = min(d[i],max_dist)
  }
  # initialize parameters
  gamma_p <- rep(0,n_omega)
  q=3
  re_error <- 10
  r = 1
  p_t <- sub_prob(gamma_p,cl,b,d,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,q,lambda,r)
  l <- as.matrix(p_t%*%t(p_t))
  gamma_t <- gamma_p+cal_proj(omega,l)-b
  error = sum(diag(l)) + r/2*norm(cal_proj(omega,l)-b+gamma_t,"2")^2
  re_error <- norm(p_t%*%t(p_t)-M,"2")
  gamma_p <- gamma_t
  
  if(length(rm_id)>0){
    frag_id <- (1:n)[-rm_id]
  }else{
    frag_id <- 1:n
  }
  stopCluster(cl)
  return(new('flamingo_prediction',id=frag_id,coordinates=p_t,input_n=n))
}

