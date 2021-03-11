args <- commandArgs(T)

#two arguments: 1. data file; 2. output file
curr_dir <- getwd()
setwd(args[1])
chr <- args[3]
raw_count <- read.table(paste0(chr,"_1mb.RAWobserved"))

kr_file <- read.table(paste0(chr,"_1mb.KRnorm"))


n <- floor(max(raw_count[,c(1,2)]))/1e6+1

normalized_mat <- matrix(0,n,n)

normalized_num <- c()

for(i in 1:dim(raw_count)[1]){
  
  tmp_i <- floor(raw_count[i,1]/1e6)+1
  
  tmp_j <- floor(raw_count[i,2]/1e6)+1
  
  normalized_mat[tmp_i,tmp_j] <- raw_count[i,3]/(kr_file[tmp_i,]*kr_file[tmp_j,])
  
  normalized_num[i] <- raw_count[i,3]/(kr_file[tmp_i,]*kr_file[tmp_j,])
  
}

input_if <- normalized_mat
input_if <- input_if + t(input_if)
diag(input_if) <- diag(input_if)/2


pd <- input_if^(-1/4)
setwd(curr_dir)
write.table(input_if,paste0(args[2],'/',chr,'_1mb_IF.txt'),col.names = F,row.names = F,sep='\t',quote=F)
write.table(pd,paste0(args[2],'/',chr,'_1mb_PD.txt'),col.names = F,row.names = F,sep='\t',quote=F)





