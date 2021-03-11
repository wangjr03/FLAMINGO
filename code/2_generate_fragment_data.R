args <- commandArgs(T)

#two arguments: 1. data dir
setwd(args[1])
chr <- args[2]


load(paste0(chr,'_5kb_IF.Rdata'))
load(paste0(chr,'_5kb_PD.Rdata'))
n <- ceiling(dim(input_if)/200)[1]

dir.create(paste0(chr,'_5kb_frag'))
setwd(paste0(chr,'_5kb_frag'))

for(i in 1:n){
  
  start_id <- (i-1)*200+1
  end_id <- min(dim(pd)[1],i*200)
  
  tmp_pd <- pd[start_id:end_id,start_id:end_id]
  tmp_input_if <- input_if[start_id:end_id,start_id:end_id]
  write.table(as.matrix(tmp_pd),paste0("Dist_frag",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
  write.table(as.matrix(tmp_input_if),paste0("IF_frag",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
  
}

