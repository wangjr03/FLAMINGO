args <- commandArgs(T)

#1. fraction of basis to use
#2. penalty
#3 distance threshold
#4 input pairwise distance matrix
#5 input interaction frequency matrix
#6 output file name
#7 output file dir



library(parallel)
library(mgcv)
library(Matrix)
#given a pairwise distance matrix, calculate gram matrix

pd <- read.table(args[4])
#load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/data/chr21_5kb_rao_dist_KR.Rdata")
#load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/data/chr21_5kb_rao_IF_KR.Rdata")
input_if <- read.table(args[5])
input_if <- as.matrix(input_if)
input_if[which(is.na(input_if))] <- NA
pd <- as.matrix(pd)

pd <- pd^2

pd[which(pd==Inf)] <- 3
pd[which(is.na(pd))] <- 3
rm_id <- which(apply(input_if,1,max)==0)
n <- dim(pd)[1]

diag(pd) <- 0

H = Diagonal(x=rep(1,n))-1/n*rep(1,n)%*%t(rep(1,n))

lambda = as.numeric(args[2])

M = -1/2*H%*%pd%*%H
M <- as.matrix(M)

sw <- as.numeric(args[1])
max_dist = as.numeric(args[3])
output_pth <- args[7]


########define functions

cal_proj <- function(x){
  
  proj <- apply(omega,1,function(y){
    
    y <- as.numeric(y)
    
    p <- x[y[1],y[1]]+x[y[2],y[2]]-x[y[1],y[2]]-x[y[2],y[1]]
    
    return(p)
  })
  
  return(proj)
  
}

cal_proj_diag <- function(x){
  
  proj <- apply(omega_diag,1,function(y){
    
    y <- as.numeric(y)
    
    p <-x[y[1],y[1]]+x[y[2],y[2]]-x[y[1],y[2]]-x[y[2],y[1]]
    
    return(p)
  })
  
  return(proj)
  
}


#calculate adjoint

def_smat <- function(x,y){
  
  sparseMatrix(i=c(x,y,x,y),j=c(x,y,y,x),x=c(1,1,-1,-1),dims=c(n,n))
  #matrix(x=c(x,y,x,y),y=c(x,y,y,x),x=c(1,1,-1,-1),dims=c(n,n))
  
}


convert_index <- function(x){
  
  # col_j <- findInterval(seq(x@x)-1,x@p[-1])+1
  # row_i <- x@i+1
  
  df <- expand.grid(as.vector(x),as.vector(x))
  
  df[,3] <- (df[,2] == df[,1])*2-1
  
  return(df)
  
  
  
}



new_calculate_adj <- function(x){
  
  clusterExport(cl,'x',envir = environment())
  
  tmp <- parSapply(cl,func_list,function(y){
    
    sum(x[y[,2]] * y[,1])
    
  })
  
  mat <- sparseMatrix(i=all_elements[,1],j=all_elements[,2],x=tmp,dims=c(n,n))
  
  return(mat)
  
}




new_calculate_adj_diag <- function(x){
  
  clusterExport(cl,'x',envir = environment())
  
  tmp <- sapply(func_list_diag,function(y){
    
    sum(x[y[,2]] * y[,1])
    
  })
  
  mat <- sparseMatrix(i=all_elements_diag[,1],j=all_elements_diag[,2],x=tmp,dims=c(n,n))
  
  return(mat)
  
  
}



grad <- function(P,gamma){
  t1 <- Sys.time()
  tmp_prod <- as.matrix(P%*%t(P))
  
  a <- cal_proj(tmp_prod)
  t2 <- Sys.time()
  l <- a-b+gamma
  
  #clusterExport(cl,'l',envir = environment())
  
  tmp_b <- new_calculate_adj(l)
  t3 <- Sys.time()
  pal = 2*lambda*new_calculate_adj_diag(cal_proj_diag(tmp_prod)-d)%*%P
  t4 <- Sys.time()
  g <- 2*P+2*r*tmp_b%*%P+pal
  #g <- 2*P+2*r*tmp_b%*%P
  print(norm(pal,'f'))
  print(c(t1,t2,t3,t4))
  return(g)
  
}


sub_prob <- function(gamma,cl){
  
  P_p <- matrix(rnorm(n*q),n,q)
  
  P_pp <- matrix(rnorm(n*q),n,q)
  
  error <- 10
  
  error.2 <- 100
  
  iter <- 1
  
  while(error>1e-2){
    
    s = P_p-P_pp
    
    if(iter > 1){
      
      g_2 <- g_1
      
      g_1 <- grad(P_p,gamma)
      
    }else{
      
      g_1 = grad(P_p,gamma)
      
      g_2 <- grad(P_pp,gamma)
      
    }
    y <- g_1 - g_2
    
    t_k = (sum(diag(t(s)%*%s)))/(sum(diag(t(s)%*%y)))
    
    P_n <- P_p-t_k*g_1
    
    P_pp <- P_p
    
    P_p <- P_n
    
    error <- norm(P_p-P_pp,"f")
    
    d_mat <- P_n%*%t(P_n)
    d_mat <- as.matrix(d_mat)
    error.2 <- sum(diag(d_mat))+r/2*norm(cal_proj(d_mat)-b+gamma,"2")^2
    
    print(error)
    
    #print("########")
    
    print(error.2)
    
    print("########")
    
    iter <- iter+1
    
    if(iter>500) break
    
  }
  return(P_n)
  
}


###############pre-calculate basis and adjoint matching##################

n_cores = 28

cl <- makeCluster(n_cores)


n <- dim(M)[1]

#omega <- t(combn(1:n,2))
#filter out low quality data
#input_if <- read.table("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/data/chr21_5kb_rao_IF_KR.txt")
input_if <- as.matrix(input_if)
input_if <- as(input_if,'sparseMatrix')
if_th <- 0

#omega_if <- apply(omega, 1, function(x){

#  input_if[x[1],x[2]]

#})

col_j = findInterval(seq(input_if@x)-1,input_if@p[-1])+1
row_i = input_if@i+1
x_x = input_if@x
df = data.frame(row_i,col_j,x_x)

omega <- subset(df[,1:2],df[,3] > if_th & df[,1] != df[,2])
omega <- as.matrix(omega)
#set sample size
n_omega <- dim(omega)[1]
diag_term <- which(omega[,2]-omega[,1]==1)

omega_diag = omega[diag_term,]

#sel_term <-  which(omega[,2]-omega[,1] <= iter)
omega <<- omega[unique(c(diag_term,sample(1:n_omega,sw*n_omega))),]
n_omega <- dim(omega)[1]
n_omega_diag <- dim(omega_diag)[1]


clusterExport(cl,"omega")
clusterExport(cl,"omega_diag")
clusterExport(cl,"convert_index")

w_list<- parLapply(cl,1:n_omega,function(x){cbind(convert_index(omega[x,]),x)})

func_list <- list()

tmp_w_list <- do.call(rbind,w_list)


all_elements <- unique(data.frame(tmp_w_list[,1:2]))


loc <- prodlim::row.match(data.frame(tmp_w_list[,1:2]),all_elements)

func_list <- split(data.frame(tmp_w_list[,c(3,4)]),loc)
clusterExport(cl,"func_list")


clusterExport(cl,"omega_diag")

w_diag_list<- parLapply(cl,1:n_omega_diag,function(x){cbind(convert_index(omega_diag[x,]),x)})

func_list_diag <- list()

tmp_w_list <- do.call(rbind,w_diag_list)

all_elements_diag <- unique(data.frame(tmp_w_list[,1:2]))

loc <- prodlim::row.match(data.frame(tmp_w_list[,1:2]),all_elements_diag)

func_list_diag <- split(data.frame(tmp_w_list[,c(3,4)]),loc)

clusterExport(cl,"func_list_diag")

b = cal_proj(M)

d = cal_proj_diag(M)

for(i in 1:length(d)){
  
  d[i] = min(d[i],max_dist)
  
}

gamma_p <- rep(0,n_omega)

q=3

error=10

re_error <- 10

r = 1
iter <- 1

setwd(output_pth)

while(iter < 2){
  
  p_t <- sub_prob(gamma_p,cl)
  l <- as.matrix(p_t%*%t(p_t))
  
  gamma_t <- gamma_p+cal_proj(l)-b
  
  
  error[1+iter] = sum(diag(l)) + r/2*norm(cal_proj(l)-b+gamma_t,"2")^2
  
  re_error[1+iter] <- norm(p_t%*%t(p_t)-M,"2")
  
  gamma_p <- gamma_t
  print(iter)
  
  print(error[1+iter])
  
  P_p = p_t
  iter <- iter + 1
  #p_t <- p_t[-rm_id,]
  save.image(paste0(args[6],'.Rdata'))
  
}

