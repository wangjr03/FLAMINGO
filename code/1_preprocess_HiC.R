args <- commandArgs(T)

curr_dir <- getwd()

#two arguments: 1. data file; 2. output file; 3 chromosome information


dir.create(args[2])

setwd(args[1])

raw_count <- read.table(paste0(args[3],"_5kb.RAWobserved"))

kr_file <- read.table(paste0(args[3],"_5kb.KRnorm"))


n <- floor(max(raw_count[,c(1,2)]))/5e3+1

normalized_mat <- matrix(0,n,n)

normalized_num <- c()



i_ind <- floor(raw_count[,1]/5e3)+1

j_ind <- floor(raw_count[,2]/5e3)+1


for(i in 1:dim(raw_count)[1]){
  tmp_i <- i_ind[i]
  tmp_j <- j_ind[i]
  
 
  normalized_mat[tmp_i,tmp_j] <- raw_count[i,3]/(kr_file[tmp_i,]*kr_file[tmp_j,])
  
  #normalized_num[i] <- raw_count[i,3]/(kr_file[tmp_i,]*kr_file[tmp_j,])
  
}

input_if <- normalized_mat
input_if <- input_if + t(input_if)
diag(input_if) <- diag(input_if)/2


pd <- input_if^(-1/4)
setwd(curr_dir)
save(input_if,file=paste0(args[2],paste0("/",args[3],'_5kb_IF.Rdata')))
save(pd,file=paste0(args[2],paste0("/",args[3],'_5kb_PD.Rdata')))
