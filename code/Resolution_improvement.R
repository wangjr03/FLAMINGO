args <- commandArgs(T)

curr_dir <- getwd()
#divide each 10kb fragment into two points
#three arguments: 1. data file; 2. output file 3. DNase file 4. resolution 5 numeric representation of resolution
setwd(args[1])

dir.create(args[2])

res <- (args[4])

res_num = as.numeric(args[5])

chr <- args[6]

raw_count <- read.table(paste0(chr,"_",res,".RAWobserved"))

kr_file <- read.table(paste0(chr,"_",res,".KRnorm"))


n <- floor(max(raw_count[,c(1,2)]))/res_num+1

normalized_mat <- matrix(0,n,n)

normalized_num <- c()

for(i in 1:dim(raw_count)[1]){
  
  tmp_i <- floor(raw_count[i,1]/res_num)+1
  
  tmp_j <- floor(raw_count[i,2]/res_num)+1
  
  normalized_mat[tmp_i,tmp_j] <- raw_count[i,3]/(kr_file[tmp_i,]*kr_file[tmp_j,])
  
  normalized_num[i] <- raw_count[i,3]/(kr_file[tmp_i,]*kr_file[tmp_j,])
  
}

input_if <- normalized_mat
input_if <- input_if + t(input_if)
diag(input_if) <- diag(input_if)/2


#pd <- input_if^(-1/4)

#expand the matrix into 2n
n = dim(input_if)[1]
input_if_expand = matrix(0,res_num/5e3*n,res_num/5e3*n)
diag(input_if) <- 0
for(i in 1:n){
  
  input_if_expand[((i-1)*res_num/5e3+1):(i*res_num/5e3),]<- rep(rep(input_if[i,],each=res_num/5e3),each = res_num/5e3)
  #input_if_expand[(i-1)*2+1,i*2] <- input_if_expand[i*2,(i-1)*2+1] <- 0
  
}

input_if <- input_if_expand
pd <- input_if^(-1/4)

save(input_if,file=paste0(args[2],'/',chr,"_",res,'_expanded_IF.Rdata'))
save(pd,file=paste0(args[2],"/",chr,"_",res,'_expanded_PD.Rdata'))

#generate fragment

setwd(args[2])

n <- ceiling(dim(input_if)/100)[1]

dir.create(paste0(chr,"_",res,'_frag'))
setwd(paste0(chr,"_",res,'_frag'))

for(i in 1:n){
  
  start_id <- (i-1)*100+1
  end_id <- min(dim(pd)[1],i*100)
  
  tmp_pd <- pd[start_id:end_id,start_id:end_id]
  tmp_input_if <- input_if[start_id:end_id,start_id:end_id]
  write.table(as.matrix(tmp_pd),paste0("Dist_frag",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
  write.table(as.matrix(tmp_input_if),paste0("IF_frag",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
  
}

#generate fragment
setwd(args[2])

load(paste0(chr,"_",res,'PD.Rdata'))

options(scipen = 999)
N = dim(pd)[1]/100

dir.create(paste0('./Genomic_fragment_',res,'/'))

for(i in 1:N){
  start_id <- (i-1)*100+1
  end_id <- min(dim(pd)[1],i*100)
  start_loc <- (seq(start_id,end_id)-1)*10000+1
  end_loc <- start_loc+9999
  frag_list <- data.frame(chr=chr,start=start_loc,end=end_loc)
  write.table(frag_list,paste0("./Genomic_fragment_",res,"/Genomic_frag",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
  
}

# overlap DNase
dnase_file <- args[3]
setwd(curr_dir)
command <- paste('./4_generate_DNase_profile.sh',paste0(args[2],"/Genomic_fragment_",res,"/"),dnase_file)
system(command)


#impute DNase
  
load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/output/1001/exponential_linear_model_-0.4.Rdata")

setwd(paste0(args[2],"/Genomic_fragment_expanded_",res,"/"))

file = dir()
n = length(grep('DNase_Genomic_frag',file))

for(id in 1:n){
  
  Dnase_signal <- read.table(paste0("DNase_Genomic_frag",id,".txt"))
  
  N = dim(Dnase_signal)[1]
  
  Dnase_mat = matrix(0,N,N)
  
  od_vec <- abs(row(Dnase_mat)-col(Dnase_mat))*5000
  
  data <- data.frame(genomic_distance=as.vector(od_vec),dnase_1 = Dnase_signal[row(Dnase_mat),4],dnase_2 = Dnase_signal[col(Dnase_mat),4])
  
  pred_val <- predict(l_model,data)
  
  predict_dist <- matrix(pred_val,byrow=F,N)
  
  write.table(predict_dist,paste0("DNase_3D_impute",id,".txt"),col.names = F,row.names = F,sep='\t',quote=F)
}









