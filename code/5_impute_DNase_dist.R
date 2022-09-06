# load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/output/1001/exponential_linear_model_-0.4.Rdata")
load("../predictions/exponential_linear_model_-0.4.Rdata")

args <- commandArgs(T)
#one arguments: 1. directory of DNase-genomic-frag data

setwd(args)

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

