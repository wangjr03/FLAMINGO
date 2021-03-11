args <- commandArgs(T)

#one arguments: path to the data dir
setwd(args[1])
chr <- args[2]
load(paste0(chr,'_5kb_PD.Rdata'))

options(scipen = 999)
N = dim(pd)[1]/200

dir.create('./Genomic_fragment/')

for(i in 1:N){
  start_id <- (i-1)*200+1
  end_id <- min(dim(pd)[1],i*200)
  start_loc <- (seq(start_id,end_id)-1)*5000+1
  end_loc <- start_loc+4999
  frag_list <- data.frame(chr=chr,start=start_loc,end=end_loc)
  write.table(frag_list,paste0("./Genomic_fragment/Genomic_frag",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
  
}
