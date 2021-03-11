#compare correlations of two cell lines
args <- commandArgs(T)

c_1 <- args[1]
c_2 <- args[2]
chr <- args[3]
res <- args[4]
frac <- args[5]



#cell line 1 within cell line result
load(paste0('../../data/',c_1,'/',c_1,'_',c_1,'/',chr,"_",res,"_",'result/',frac,'.Rdata'))


#cell line 1 IF and PD
load(paste0('../../data/',c_1,'/',chr,"_",res,"_PD.Rdata"))
load(paste0('../../data/',c_1,'/',chr,"_",res,"_IF.Rdata"))


c1_pd <- pd
id_list_1 <- id_list

#cell line 2 within cell line result
load(paste0('../../data/',c_2,'/',c_2,'_',c_2,'/',chr,"_",res,"_",'result/',frac,'.Rdata'))


#cell line 2 IF and PD
load(paste0('../../data/',c_2,'/',chr,"_",res,"_PD.Rdata"))
load(paste0('../../data/',c_2,'/',chr,"_",res,"_IF.Rdata"))
c2_pd <- pd
id_list_2 <- id_list


#cell line 2 within cell line result
load(paste0('../../data/',c_2,'/',c_1,'_',c_2,'/',chr,"_",res,"_",'result/',frac,'.Rdata'))

res_id_list <- id_list
result <- all_points

#match points
common_points <- intersect(unlist(id_list_2),intersect(unlist(id_list_1),unlist(res_id_list)))

#process two datasets
c1_pd_common <- c1_pd[unlist(common_points),unlist(common_points)]
c2_pd_common <- c2_pd[unlist(common_points),unlist(common_points)]

res_common <- result[na.omit(match(common_points,unlist(res_id_list))),]

res_pd <- as.matrix(dist(res_common))

missing_id <- which(c1_pd_common>=3|is.na(c1_pd_common)|c2_pd_common>=3|is.na(c2_pd_common))
c2_pd_common[which(c2_pd_common > 3)] <- 3
c2_pd_common[which(is.na(c2_pd_common))] <- 3
c1_pd_common[which(c1_pd_common > 3)] <- 3
c1_pd_common[which(is.na(c1_pd_common))] <- 3
all_cor_cross <- cor(as.vector(res_pd),as.vector(c2_pd_common))
valid_cor_cross <- cor(as.vector(res_pd[-missing_id]),as.vector(c2_pd_common[-missing_id]))
all_cor_wi <- cor(as.vector(res_pd),as.vector(c1_pd_common))
valid_cor_wi <- cor(as.vector(res_pd[-missing_id]),as.vector(c1_pd_common[-missing_id]))


intra_tad_dist <- list()
intra_tad_obs_wi <- list()
intra_tad_obs_cross <- list()
inter_tad_obs_cross <- list()
inter_tad_obs_wi <-list()
inter_tad_dist <- list()
common_id_list <- list()

for(i in 1:length(res_id_list)){
  
  common_id_list[[i]] <- intersect(unlist(id_list_2),intersect(res_id_list[[i]],unlist(id_list_1)))
  
}

for(i in 1:length(res_id_list)){
  
  n <- unlist(common_id_list[0:(i-1)])
  n <- length(n)
  
  start <- n+1
  end <- n+length(id_list[[i]])
  #id <- which(c1_pd_common[start:end,start:end]>=3)
  intra_tad_dist[[i]] <- res_pd[start:end,start:end]
  intra_tad_obs_wi[[i]] <- c1_pd_common[start:end,start:end]
  intra_tad_obs_cross[[i]] <- c2_pd_common[start:end,start:end]
  
  inter_tad_dist[[i]] <- res_pd[start:end,-c(start:end)]
  inter_tad_obs_cross[[i]] <-c2_pd_common[start:end,-c(start:end)]
  inter_tad_obs_wi[[i]] <-c1_pd_common[start:end,-c(start:end)]
  inter_tad_obs_wi[[i]] <-c1_pd_common[start:end,-c(start:end)]
}

rm_id_wi <- which(unlist(intra_tad_obs_wi)==3)
rm_id_cross <- which(unlist(intra_tad_obs_cross)==3)

rm_id_inter_cross <- which(unlist(inter_tad_obs_cross)==3)
rm_id_inter_wi <- which(unlist(inter_tad_obs_wi)==3)

intra_cor_wi <- cor(unlist(intra_tad_dist)[-rm_id_wi],unlist(intra_tad_obs_wi)[-rm_id_wi])
intra_cor_cross <- cor(unlist(intra_tad_dist)[-rm_id_cross],unlist(intra_tad_obs_cross)[-rm_id_cross])
inter_cor_cross <- cor(unlist(inter_tad_dist)[-rm_id_inter_cross],unlist(inter_tad_obs_cross)[-rm_id_inter_cross])
inter_cor_wi <- cor(unlist(inter_tad_dist)[-rm_id_inter_wi],unlist(inter_tad_obs_wi)[-rm_id_inter_wi])


#calculate hic similarity
common_valid <- which(c1_pd_common<3&c2_pd_common<3)
hic_cor <- cor(as.vector(c1_pd_common[common_valid]),as.vector(c2_pd_common[common_valid]))

#calculate DNase similarity
#get c1 dnase data
curr_dir <- getwd()
setwd(paste0('../../data/',c_1,'/Genomic_fragment'))
file <- dir()
file <- file[grep('DNase_Genomic_frag',file)]
N <- length(file)
c1_dnase <- list()

for(i in 1:N){
  
  c1_dnase[[i]] <- read.table(paste0('DNase_Genomic_frag',i,'.txt'),sep='\t')
  
}

c1_dnase <- do.call(rbind,c1_dnase)


setwd(curr_dir)

setwd(paste0('../../data/',c_2,'/Genomic_fragment'))
file <- dir()
file <- file[grep('DNase_Genomic_frag',file)]
N <- length(file)
c2_dnase <- list()

for(i in 1:N){
  
  c2_dnase[[i]] <- read.table(paste0('DNase_Genomic_frag',i,'.txt'),sep='\t')
  
}

c2_dnase <- do.call(rbind,c2_dnase)

dnase_cor <- cor(c1_dnase[,4],c2_dnase[,4])



df <- data.frame(c_1,c_2,chr,res,frac,all_cor_wi,all_cor_cross,valid_cor_wi,valid_cor_cross,intra_cor_wi,intra_cor_cross,inter_cor_wi,inter_cor_cross,hic_cor,dnase_cor)

write.table(df,'/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/data/correlation_summary_no_DNase.txt',append = T,col.names = F,row.names = F,sep='\t',quote=F)
