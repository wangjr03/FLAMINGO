# the whole structure
newargs <- commandArgs(T)
#1. input 1mb interaction frequency
#2. backbone structure
#3. 5kb result output directory
#4. 5kb pairwise distance matrix
#5. output result name

input_if <- read.table(newargs[1])

#load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/output/chr21_5kb_frag/0727_chr21_5kb_fragAll.Rdata")


#load('/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/data/chr21_1mb_backbone_K562.Rdata')

#get which point should be removed
rm_id <- apply(input_if,1,sum)
rm_id <- which(rm_id == 0)
n <- dim(input_if)[1]

#total_frag <- dim(input_if)[1]

#total_struc <- p_t[-rm_id,]

frag_id <- c(1:n)[-rm_id]
curr_dir <- getwd()
#load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/output/mds_total_struc.Rdata")


load(newargs[2])
#load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/data/chr21_1mb_backbone_K562.Rdata")

#load each sub fragment and calculate coordinates
p_t <- as.matrix(p_t)
total_frag <- dim(input_if)[1]

total_struc <- p_t
all_points <- list()
setwd(newargs[3])
rm_id_list <- list()
scaler <- as.matrix(dist(total_struc))

scaler <- mean(scaler[which(row(scaler)+1 == col(scaler))])
frag_id_valid <- c()
acc_counter <- 1
for(counter in 1:length(frag_id)){
  if(!paste0('5kb_frag',frag_id[counter],'.Rdata') %in% dir()){
    next
  }
  
  load(paste0('5kb_frag',frag_id[counter],'.Rdata'))
  tmp_center <- total_struc[counter,]
  
  # 
  # prev_center <- total_struc[i-1,]
  # 
  # next_center <- total_struc[i+1,]
  #
  input_if <- as.matrix(input_if)
  input_if[which(is.na(input_if))] <- 0
  rm_id <- apply(input_if,1,sum)
  rm_id <- which(rm_id == 0)
  rm_id_list[[acc_counter]] <- rm_id
  if(length(rm_id) > 0){
    
    p_t <- p_t[-rm_id,]}
  
  if(dim(p_t)[1] == 0){
    next
  }
  frag_id_valid <- c(frag_id_valid,frag_id[counter])
  old_center <- apply(p_t,2,mean)
  
  radius <- max(apply(p_t,1,function(x){norm(x-old_center,'2')}))
  
  new_loc <- t(apply(p_t,1,function(x){
    
    (x-old_center)/radius*scaler+tmp_center
    #(x-old_center)+tmp_center
    
  }))
  
  new_loc <- as.data.frame(new_loc)
  # new_loc$id <- counter
  # all_points <- rbind(all_points,new_loc)
  all_points[[acc_counter]] <- as.data.frame(new_loc)
  acc_counter = acc_counter+1
}
total_struc <- total_struc[match(frag_id_valid,frag_id),]
frag_id <- frag_id_valid
names(all_points) <- frag_id

setwd(curr_dir)

#prepare input data

#calculate point id
id_list <- list()
for(i in 1:length(frag_id)){
  tmp_id <- 1:200 
  if(length(rm_id_list[[i]])>0){
    tmp_id <- tmp_id[-rm_id_list[[i]]]
  }
  
  tmp_id <- (frag_id[i]-1)*200+tmp_id
  id_list[[i]] <- tmp_id[1:dim(all_points[[i]])[1]]
  
  
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


#allocate the first point
direction <- (total_struc[2,]-total_struc[1,])
tail_direction <- tail(all_points[[1]],n=1)-total_struc[1,]
direction <- direction/norm(direction,'2')
tail_direction <- tail_direction/norm(tail_direction,'2')
direction <- as.vector(direction)
tail_direction <- as.vector(as.matrix(tail_direction))
r_mat <- rotation_matrix(matrix(tail_direction,ncol=1),as.matrix(direction,nrow=1))

all_points[[1]] <- rotate(all_points[[1]],total_struc[1,],r_mat)

#get all point location and pair-wise distance
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

get_dist_vec <- function(all_points,id_list){
  n <- length(all_points)
  start_id <- sapply(id_list[2:(n)],function(y){y[1]})
  end_id <- sapply(id_list[1:(n-1)],function(y){tail(y,n=1)})
  dist_vec <- c()
  for(i in 1:(n-1)){
    
    dist_vec = c(dist_vec,pd[start_id[i],end_id[i]])
    
  }
  dist_vec[which(dist_vec > 3)] <- max(dist_vec[which(dist_vec<3)])
  #dist_vec <- dist_vec/radius*scaler
  return(dist_vec)
}

#start optimization

evaluate_dist <- function(start_point,end_point){
  
  apply(start_point-end_point,1,function(x){norm(x,'2')})
  
}

#load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/data/chr21_Hi_C_PD_5kb.Rdata")


load(newargs[4])

#target function l = \sum (||y_s_i+1 - y_e_i || - d_i,i+1)^2
# gradient 4 (\tild(D) - D)^2(y_s-y_e)
#update total_structure, update y_s,y_e,\tild(D)

error <- 10

d = get_dist_vec(all_points,id_list)

t = 1e-4
#all_points_bp <- all_points
g_p <- matrix(0,length(all_points)-1,3)
y_s_p <-  matrix(0,length(all_points)-1,3)
error_p = error_change <- 1
while(error>1e-4&error_change>1e-5){
  
  y_s <- get_start_point(all_points)
  y_e <- get_end_point(all_points)
  d_tild <- evaluate_dist(y_s,y_e)
  
  g <- 4*(d_tild-d)*(y_s-y_e)
  
  y = g-g_p
  s = y_s-y_s_p
  t_k = (sum(diag(t(s)%*%s)))/(sum(diag(t(s)%*%y)))
  
  y_s_prim <- y_s-t_k*g
  
  g_p <- g
  y_s_p <- y_s
  
  #update structure
  
  for(i in 2:length(all_points)){
    tmp_points <- all_points[[i]]
    tmp_start <- tmp_points[1,]
    tmp_center <- total_struc[i,]
    old_direction <- as.vector(as.matrix(tmp_start-tmp_center))
    new_direction <- y_s_prim[i-1,] - tmp_center
    r_mat <- rotation_matrix(old_direction,new_direction)
    all_points[[i]] <- rotate(all_points[[i]],total_struc[i,],r_mat)
  }
  
  error <- norm(d_tild-d,'2')/norm(d,'2')
  error_change <- abs(error-error_p)
  error_p <- error
  print(error)
  
}

setwd(newargs[3])

all_points <- do.call(rbind,all_points)

save(all_points,id_list,file=paste0(newargs[5],'.Rdata'))

