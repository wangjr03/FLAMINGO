#' flamingo.impute_dist_from_epi
#'
#' Impute 3D distance between DNA fragments based on the 1D epigenomic signal
#' @param frag_res Size of the small DNA fragment.
#' @keywords FLAMINGO
#' @return Write the imputed 3D distances between DNA fragments.
#' @export

flamingo.impute_dist_from_epi <- function(frag_res){
  file <- dir('./Genomic_loc')
  file <- file[grep("Epi_sig_domain_",file)]
  n <- length(file)
  pb <- progress::progress_bar$new(total = n)
  for(id in 1:n){
    Dnase_signal <- read.table(paste0("./Genomic_loc/Epi_sig_domain_",id,".txt"))
    N = dim(Dnase_signal)[1]
    Dnase_mat = matrix(0,N,N)
    od_vec <- abs(row(Dnase_mat)-col(Dnase_mat))*frag_res
    data <- data.frame(genomic_distance=as.vector(od_vec),dnase_1 = Dnase_signal[row(Dnase_mat),4],dnase_2 = Dnase_signal[col(Dnase_mat),4])
    pred_val <- predict(l_model,data)
    predict_dist <- matrix(pred_val,byrow=F,N)
    write.table(predict_dist,paste0("./Genomic_loc/3D_Dist_impute_Epi_",id,".txt"),col.names = F,row.names = F,sep='\t',quote=F)
    pb$tick()
  }
  
}
