args <- commandArgs(T)
#1. 5kb data
#2. cell line
#3. resolution
chr <- args[4]



load(paste0("../../data/",args[1],'/',chr,'_5kb_PD.Rdata'))
high_pd <- pd

load(paste0("../../data/",args[1],"/",chr,"_",args[2],"_expanded_no_DNase/0.75.Rdata"))

low_id_list = id_list
low_coord = all_points

high_id <- 1:dim(high_pd)[1]

common_id <- intersect(high_id,unlist(low_id_list))

high_pd_sub <- high_pd[common_id,common_id]

low_coord_sub <- low_coord[match(common_id,unlist(low_id_list)),]

low_res_pd <- as.matrix(dist(low_coord_sub))

id <- which(high_pd_sub>3 | is.na(high_pd_sub))

all_cor_no_dnase = cor(as.vector(high_pd_sub[-id]),as.vector(low_res_pd[-id]))


intra_tad_dist <- list()
intra_tad_obs<- list()
inter_tad_obs<-list()
inter_tad_dist <- list()
common_id_list <- list()

low_pd <- as.matrix(dist(low_coord))

for(i in 1:length(low_id_list)){
  
  common_id_list[[i]] <- intersect(high_id,low_id_list[[i]])
  
}


for(i in 1:length(common_id_list)){
  tmp_loc = na.omit(match(common_id_list[[i]],unlist(low_id_list)))
  final_id <- unlist(low_id_list)[tmp_loc]
  high_pd_id <- match(final_id,high_id)
  
  intra_tad_obs[[i]] <- high_pd[high_pd_id,high_pd_id]
  
  intra_tad_dist[[i]] <- low_pd[tmp_loc,tmp_loc]
  
  high_pd_sub_id <- which( common_id%in% final_id)
  
  inter_tad_obs[[i]] <- high_pd_sub[high_pd_sub_id,-high_pd_sub_id]
  inter_tad_dist[[i]] <-low_res_pd[tmp_loc,-tmp_loc]
}


rm_id <- which(unlist(intra_tad_obs)>3 | is.na(unlist(intra_tad_obs)))

rm_id_inter <- which(unlist(inter_tad_obs)>3 | is.na(unlist(inter_tad_obs)))

intra_cor_no_dnase <- cor(unlist(intra_tad_dist)[-rm_id],unlist(intra_tad_obs)[-rm_id])
inter_cor_no_dnase <- cor(unlist(inter_tad_dist)[-rm_id_inter],unlist(inter_tad_obs)[-rm_id_inter])



load(paste0("../../data/",args[1],"/",chr,"_",args[2],"_expanded_DNase/0.75.Rdata"))

low_id_list = id_list
low_coord = all_points

high_id <- 1:dim(high_pd)[1]

common_id <- intersect(high_id,unlist(low_id_list))

high_pd_sub <- high_pd[common_id,common_id]

low_coord_sub <- low_coord[match(common_id,unlist(low_id_list)),]

low_res_pd <- as.matrix(dist(low_coord_sub))

id <- which(high_pd_sub>3 | is.na(high_pd_sub))

all_cor_dnase = cor(as.vector(high_pd_sub[-id]),as.vector(low_res_pd[-id]))


#check intra-TAD cor

intra_tad_dist <- list()
intra_tad_obs<- list()
inter_tad_obs<-list()
inter_tad_dist <- list()
common_id_list <- list()

low_pd <- as.matrix(dist(low_coord))

for(i in 1:length(low_id_list)){
  
  common_id_list[[i]] <- intersect(high_id,low_id_list[[i]])
  
}


for(i in 1:length(common_id_list)){
  tmp_loc = na.omit(match(common_id_list[[i]],unlist(low_id_list)))
  final_id <- unlist(low_id_list)[tmp_loc]
  high_pd_id <- match(final_id,high_id)
  
  intra_tad_obs[[i]] <- high_pd[high_pd_id,high_pd_id]
  
  intra_tad_dist[[i]] <- low_pd[tmp_loc,tmp_loc]
  
  high_pd_sub_id <- which( common_id%in% final_id)
  
  inter_tad_obs[[i]] <- high_pd_sub[high_pd_sub_id,-high_pd_sub_id]
  inter_tad_dist[[i]] <-low_res_pd[tmp_loc,-tmp_loc]
}



rm_id <- which(unlist(intra_tad_obs)>3 | is.na(unlist(intra_tad_obs)))

rm_id_inter <- which(unlist(inter_tad_obs)>3 | is.na(unlist(inter_tad_obs)))

intra_cor_dnase <- cor(unlist(intra_tad_dist)[-rm_id],unlist(intra_tad_obs)[-rm_id])
inter_cor_dnase <- cor(unlist(inter_tad_dist)[-rm_id_inter],unlist(inter_tad_obs)[-rm_id_inter])

df = data.frame(cell_line = args[1],resolution = args[2],all_cor_no_dnase,all_cor_dnase,intra_cor_no_dnase,intra_cor_dnase,inter_cor_no_dnase,inter_cor_dnase)

write.table(df,'../../data/Resolution_summary_V2.txt',append = T,col.names = F,row.names = F,sep='\t',quote=F)


