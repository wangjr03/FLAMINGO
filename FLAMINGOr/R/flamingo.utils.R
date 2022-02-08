
cal_proj <- function(omega,x){

  proj <- apply(omega,1,function(y){
    
    y <- as.numeric(y)
    
    p <- x[y[1],y[1]]+x[y[2],y[2]]-x[y[1],y[2]]-x[y[2],y[1]]
    
    return(p)
  })
  
  return(proj)
  
}

cal_proj_diag <- function(omega_diag,x){
  
  proj <- apply(omega_diag,1,function(y){
    
    y <- as.numeric(y)
    
    p <-x[y[1],y[1]]+x[y[2],y[2]]-x[y[1],y[2]]-x[y[2],y[1]]
    
    return(p)
  })
  
  return(proj)
  
}


convert_index <- function(x){
  
  # col_j <- findInterval(seq(x@x)-1,x@p[-1])+1
  # row_i <- x@i+1
  
  df <- expand.grid(as.vector(x),as.vector(x))
  
  df[,3] <- (df[,2] == df[,1])*2-1
  
  return(df)
  
  
  
}



new_calculate_adj <- function(x,func_list,all_elements,n,cl){
  
  clusterExport(cl,'x',envir = environment())
  
  tmp <- parSapply(cl,func_list,function(y){
    
    sum(x[y[,2]] * y[,1])
    
  })
  
  mat <- sparseMatrix(i=all_elements[,1],j=all_elements[,2],x=tmp,dims=c(n,n))
  
  return(mat)
  
}




new_calculate_adj_diag <- function(x,func_list_diag,all_elements_diag,n){
  
  tmp <- sapply(func_list_diag,function(y){
    
    sum(x[y[,2]] * y[,1])
    
  })
  
  mat <- sparseMatrix(i=all_elements_diag[,1],j=all_elements_diag[,2],x=tmp,dims=c(n,n))
  
  return(mat)
  
  
}



grad <- function(P,gamma,b,d,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,cl,lambda,r){
  tmp_prod <- as.matrix(P%*%t(P))
  a <- cal_proj(omega,tmp_prod)
  l <- a-b+gamma
  tmp_b <- new_calculate_adj(l,func_list,all_elements,n,cl)
  diag_proj = cal_proj_diag(omega_diag,tmp_prod)
  diag_adj = new_calculate_adj_diag(diag_proj-d,func_list_diag,all_elements_diag,n)
  pal = 2*lambda*diag_adj%*%P
  g <- 2*P+2*r*tmp_b%*%P+pal
  return(g)
}


sub_prob <- function(gamma,cl,b,d,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,q,lambda,r){
  P_p <- matrix(rnorm(n*q),n,q)
  P_pp <- matrix(rnorm(n*q),n,q)
  error <- 10
  error.2 <- 100
  iter <- 1
  while(error>1e-3){
    #print(paste('iter:',iter))
    s = P_p-P_pp
    if(iter > 1){
      g_2 <- g_1
      g_1 <- grad(P_p,gamma,b,d,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,cl,lambda,r)
    }else{
      g_1 = grad(P_p,gamma,b,d,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,cl,lambda,r)
      g_2 <- grad(P_pp,gamma,b,d,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,cl,lambda,r)
    }
    y <- g_1 - g_2
    s = as.matrix(s)
    y=as.matrix(y)
    t_k = (sum(diag(t(s)%*%s)))/(sum(diag(t(s)%*%y)))
    P_n <- P_p-t_k*g_1
    P_n <- as.matrix(P_n)
    P_pp <- P_p
    P_p <- P_n
    error <- norm(P_p-P_pp,"f")
    d_mat <- P_n%*%t(P_n)
    d_mat <- as.matrix(d_mat)
    error.2 <- sum(diag(d_mat))+r/2*norm(cal_proj(omega,d_mat)-b+gamma,"2")^2
    #print(paste0('Loss: ',error.2))
    #print(paste0('Coordinate change: ',error))
    iter <- iter+1
    if(iter>200) break
  }
  return(P_n)
}


iflamingo_grad <- function(P,gamma,b,d,epi_vec,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,cl,lambda,lambda_2,r){
  tmp_prod <- as.matrix(P%*%t(P))
  a <- cal_proj(omega,tmp_prod)
  l <- a-b+gamma
  tmp_b <- new_calculate_adj(l,func_list,all_elements,n,cl)
  diag_proj = cal_proj_diag(omega_diag,tmp_prod)
  diag_adj = new_calculate_adj_diag(diag_proj-d,func_list_diag,all_elements_diag,n)
  pal = 2*lambda*diag_adj%*%P
  d_term = 2*lambda_2*new_calculate_adj(a-epi_vec,func_list,all_elements,n,cl)%*%P
  g <- 2*P+2*r*tmp_b%*%P+pal+d_term
  return(g)
}


iflamingo_sub_prob <- function(gamma,cl,b,d,epi_vec,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,q,lambda,lambda_2,r){
  P_p <- matrix(rnorm(n*q),n,q)
  P_pp <- matrix(rnorm(n*q),n,q)
  error <- 10
  error.2 <- 100
  iter <- 1
  while(error>1e-3){
    #print(paste('iter:',iter))
    s = P_p-P_pp
    if(iter > 1){
      g_2 <- g_1
      g_1 <- iflamingo_grad(P_p,gamma,b,d,epi_vec,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,cl,lambda,lambda_2,r)
    }else{
      g_1 = iflamingo_grad(P_p,gamma,b,d,epi_vec,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,cl,lambda,lambda_2,r)
      g_2 <- iflamingo_grad(P_pp,gamma,b,d,epi_vec,omega,omega_diag,func_list,func_list_diag,all_elements,all_elements_diag,n,cl,lambda,lambda_2,r)
    }
    y <- g_1 - g_2
    s = as.matrix(s)
    y=as.matrix(y)
    t_k = (sum(diag(t(s)%*%s)))/(sum(diag(t(s)%*%y)))
    P_n <- P_p-t_k*g_1
    P_n <- as.matrix(P_n)
    P_pp <- P_p
    P_p <- P_n
    error <- norm(P_p-P_pp,"f")
    d_mat <- P_n%*%t(P_n)
    d_mat <- as.matrix(d_mat)
    error.2 <- sum(diag(d_mat))+r/2*norm(cal_proj(omega,d_mat)-b+gamma,"2")^2
    #print(paste0('Loss: ',error.2))
    #print(paste0('Coordinate change: ',error))
    iter <- iter+1
    if(iter>500) break
  }
  return(P_n)
  
}






rotation_matrix = function(x,y){

  u=x/sqrt(sum(x^2))
  
  v=y-sum(u*y)*u
  v=v/sqrt(sum(v^2))
  
  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
  
  sint=sqrt(max(1-cost^2,0));
  #print(cost)
  diag(length(x)) - u %*% t(u) - v %*% t(v) + 
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
}

rotate <- function(x,c,r){
  
  tmp_c <- t(apply(x,1,function(y){y-c}))
  z = as.matrix(tmp_c) %*% r
  z = t(apply(z,1,function(y){y+c}))
  return(z)
  
}

get_start_point <- function(all_points){
  n <- length(all_points)
  t(sapply(all_points[2:(n)],function(x){
    as.matrix(head(x,n=1))
  }))
  
}

get_end_point <- function(all_points){
  n <- length(all_points)
  t(sapply(all_points[1:(n-1)],function(x){
    as.matrix(tail(x,n=1))
  }))
  
}

ave_dist <- function(r,pd){
  n = dim(pd)[1]
  res = pd[row(pd)-col(pd)==r]
  a = mean(res,na.rm=T)
  return(mean(a))
}

get_dist_vec <- function(all_points,id_list,pd){
  n <- length(all_points)
  start_id <- sapply(id_list[2:(n)],function(y){y[1]})
  end_id <- sapply(id_list[1:(n-1)],function(y){tail(y,n=1)})
  dist_vec <- c()
  for(i in 1:(n-1)){
    tmp_dist = pd[start_id[i],end_id[i]]
    if(is.na(tmp_dist)){
      tmp_dist = ave_dist(start_id[i]-end_id[i],pd)
    }else if(tmp_dist == Inf){
      tmp_dist = ave_dist(start_id[i]-end_id[i],pd)
    }else if(tmp_dist == 0){
      tmp_dist = 0.001
    }
    dist_vec = c(dist_vec,tmp_dist)
    
  }
  dist_vec[which(dist_vec > 3)] <- max(dist_vec[which(dist_vec<3)])
  dist_vec[which(is.na(dist_vec))] <- mean(dist_vec,na.rm=T)
  #dist_vec <- dist_vec/radius*scaler
  return(dist_vec)
}


#start optimization

evaluate_dist <- function(start_point,end_point){
  
  apply(start_point-end_point,1,function(x){norm(x,'2')})
  
}

smt <- function(o){
  o_smt <- o
  for(i in 2:c(dim(o_smt)[1]-2)){
    
    o_smt[i,] <- apply(o_smt[(i-1):(i+1),],2,mean)
    
    
  }
  return(o_smt)
}

check_data_availability <- function(x){
  x <- as.matrix(x)
  invalid_id <- which(apply(x,1,min)==Inf)
  if(length(invalid_id) == dim(x)[1]){
    return(F)
  }else{
    return(T)
  }
}

convert_huge_mat <- function(sparse_mat){
  print('Contact map is too large, large matrix mod is on')
  n <- dim(sparse_mat)[1]
  res <- matrix(0,n,n)
  n_bin <- ceiling(n/1000)
  bin_size = 1000
  pb <- progress::progress_bar$new(total = n_bin*n_bin)
  for(i in 1:n_bin){
    row_idx = (1+(i-1)*bin_size):min(n,i*bin_size)
    for(j in 1:n_bin){
      col_idx = (1+(j-1)*bin_size):min(n,j*bin_size)
      res[row_idx,col_idx] <- as.matrix(sparse_mat[row_idx,col_idx])
      pb$tick()
    }
  }
  return(res)
}

write_huge_mat <- function(mat,file_path,nThread=10){
  n = dim(mat)[1]
  block_size=1000
  n_block = ceiling(n/block_size)
  pb <- progress::progress_bar$new(total = n_block)
  for(row in 1:n_block){
    start_id = 1+(row-1)*1000
    end_id = min(row*1000,n)
    data.table::fwrite(mat[start_id:end_id,],file_path,col.names=F,row.names=F,sep='\t',quote=F,append=T,nThread=nThread)
    pb$tick()
  }
}
