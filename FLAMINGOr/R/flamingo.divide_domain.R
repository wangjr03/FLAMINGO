#' flamingo.divide_domain
#'
#' Divide the whole chromosome into large domains to perform hierarchical reconstruction
#' @param flamingo_obj FLAMINGO object containing the chromosome-wide interaction frequency matrix and pairwise distance matrix.
#' @param domain_res Size of the large domain.
#' @param frag_res Size of the small DNA fragment.
#' @keywords FLAMINGO
#' @return Write out a list of Interaction Frequency matrix (IF) and Pairwise Distance matrix (PD)
#' @export
flamingo.divide_domain <- function(flamingo_obj,domain_res,frag_res){
  options(scipen = 999)
  dir.create("Domain_data")
  dir.create('Genomic_loc')
  n = flamingo_obj@n_frag
  bin_size = domain_res/frag_res
  if(bin_size != round(bin_size)){
    stop('The domain resolution is not an integer multiple of the fragment resolution!')
  }
  res = list()
  chr = flamingo_obj@chr_name
  pd = flamingo_obj@PD
  input_if = flamingo_obj@IF
  print('Processing Fragments...')
  n_domain = ceiling(n/bin_size)
  pb <- progress::progress_bar$new(total = n_domain)
  for(i in 1:n_domain){
    start_id <- (i-1)*bin_size+1
    end_id <- min(dim(pd)[1],i*bin_size)
    
    tmp_pd <- pd[start_id:end_id,start_id:end_id]
    tmp_input_if <- input_if[start_id:end_id,start_id:end_id]
    start_loc = ((start_id:end_id)-1)*frag_res+1
    end_loc = (start_id:end_id)*frag_res
    tmp_frag = data.frame(chr=chr,start=start_loc,end=end_loc)
    write.table(as.matrix(tmp_pd),paste0("./Domain_data/PD_domain_",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
    write.table(as.matrix(tmp_input_if),paste0("./Domain_data/IF_domain_",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
    write.table(tmp_frag,paste0("./Genomic_loc/Genomic_loc_domain_",i,'.txt'),col.names = F,row.names = F,sep="\t",quote=F)
    pb$tick()
  }
}

