#' flamingo.calculating_epigenomic_signal
#'
#' Calculate the averaged epigenomic signal for each DNA fragment.
#' @param epigenomic_file The path to the epigenomic file, which should has four columns (chr,start,end,score).
#' @keywords FLAMINGO
#' @return Write out a list of files containing the genomic location and averaged epigenomic signals for DNA fragments.
#' @export
flamingo.calculating_epigenomic_signal <- function(epigenomic_file){
  options(scipen = 999)
  all_score = bedtoolsr::bt.map(a="./Genomic_loc/Genomic_loc_all_domain.txt",b=epigenomic_file,c=4,o='mean')
  file <- dir('./Genomic_loc')
  file <- file[grep("Genomic_loc_domain_",file)]
  n <- length(file)
  pb <- progress::progress_bar$new(total = n)
  print('Calculating Epigenomic signal...')
  for(i in 1:n){
    tmp_file <- read.table(paste0('./Genomic_loc/Genomic_loc_domain_',i,".txt"))
    tmp_row_id <- prodlim::row.match(tmp_file,all_score[,1:3])
    write.table(all_score[tmp_row_id,],paste0('./Genomic_loc/Epi_sig_domain_',i,'.txt'),quote=F,col.names=F,row.names=F,sep='\t')
    pb$tick()
  }
}
