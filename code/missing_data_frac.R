args <- commandArgs(T)
cell <- args[1]
res = args[2]
chr <- args[3]
load(paste0("../../data/",cell,"/",chr,"_5kb_PD.Rdata"))
load(paste0("../../data/",cell,"/",chr,"_5kb_IF.Rdata"))


load(paste0("../../data/",cell,"/result_",res,"/ensemble_structure.Rdata"))

pd <- pd[unlist(id_list),unlist(id_list)]
missing_id <- which(pd>3|is.na(pd))
pd[missing_id] <- 3
input_if <- input_if[unlist(id_list),unlist(id_list)]


our_pd <- as.matrix(dist(all_points))
all_cor <- cor(as.vector(our_pd),as.vector(pd))
valid_cor <- cor(as.vector(our_pd[-missing_id]),as.vector(pd[-missing_id]))

intra_tad_dist <- list()
intra_tad_obs <- list()
inter_tad_obs <- list()
inter_tad_dist <- list()
for(i in 1:length(id_list)){
  
  n <- unlist(id_list[0:(i-1)])
  n <- length(n)
  
  start <- n+1
  end <- n+length(id_list[[i]])
  intra_tad_dist[[i]] <- our_pd[start:end,start:end]
  intra_tad_obs[[i]] <- pd[start:end,start:end]
  inter_tad_dist[[i]] <- our_pd[start:end,-c(start:end)]
  inter_tad_obs[[i]] <-pd[start:end,-c(start:end)]
}

rm_id <- which(unlist(intra_tad_obs)==3)
rm_id_inter <- which(unlist(inter_tad_obs)==3)
intra_cor <- cor(unlist(intra_tad_dist)[-rm_id],unlist(intra_tad_obs)[-rm_id])
inter_cor <- cor(unlist(inter_tad_dist)[-rm_id_inter],unlist(inter_tad_obs)[-rm_id_inter])

df = data.frame(cell,res,all_cor,intra_cor,inter_cor)

write.table(df,"/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj10_recon_3D/data/missing_data_frac.txt",col.names = F,row.names = F,sep='\t',quote=F,append = T)
