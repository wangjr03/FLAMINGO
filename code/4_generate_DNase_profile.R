args <- commandArgs(T)
input_file_dir <- args[1]
dnase_file <- args[2]
setwd(input_file_dir)

file <- dir()
file <- file[grep("Genomic_frag",file)]

n <- length(file)

total_file <- c()

for(i in 1:n){
  
  tmp_file <- read.table(paste0('Genomic_frag',i,".txt"))
  total_file <- rbind(total_file,tmp_file)
  
  
}

write.table(total_file,"total_file.txt",col.names = F,sep = '\t',row.names = F,quote=F)

command <- paste("bedtools map -a total_file.txt -b",dnase_file, "-c 4 -o mean > total_Dnase.txt")
system(command)


total_DNase <- read.table('total_Dnase.txt')
dnase_list <- split(total_DNase,floor((0:(dim(total_DNase)[1]-1))/200))

for(i in 1:length(dnase_list)){
  
  write.table(dnase_list[[i]],paste0('DNase_Genomic_frag',i,'.txt'),quote=F,col.names=F,row.names=F,sep='\t')
  
}
