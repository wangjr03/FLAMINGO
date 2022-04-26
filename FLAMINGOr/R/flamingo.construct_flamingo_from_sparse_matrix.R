#' flamingo.construct_flamingo_from_sparse_matrix
#'
#' Generate FLAMINGO input data from the sparse format of Hi-C data from Rao et al 2014.
#' @param raw_count_file Path to the RawCount file from Rao et al.
#' @param norm_file Path to the noralization vector file from Rao et al.
#' @param resolution Target resolution.
#' @param chr_name Name of the chromosome, e.g. chr1.
#' @param chr_size Chromosome size
#' @param alpha convertion factor between interaction frequency and distance, default is -0.25
#' @keywords FLAMINGO
#' @return A FLAMINGO object containing the resulted Interaction Frequency (IF) matrix, Pairwise Distance (PD) matrix and number of fragments.
#' @export

flamingo.construct_flamingo_from_sparse_matrix <- function(raw_count_file,norm_file=NULL,resolution,chr_name,chr_size,alpha=-0.25){
  raw_count <- read.table(raw_count_file)
  norm_file <- read.table(norm_file)
  n <- ceiling(chr_size/resolution)
  normalized_mat <- matrix(0,n,n)
  normalized_num <- c()
  i_ind <- 1+floor(raw_count[,1]/resolution)
  j_ind <- 1+floor(raw_count[,2]/resolution)
  row_count = dim(raw_count)[1]
  pb <- progress::progress_bar$new(total = row_count)
  print('Constructing matrix...')
  for(i in 1:row_count){
    tmp_i <- i_ind[i]
    tmp_j <- j_ind[i]
    normalized_mat[tmp_i,tmp_j] <- raw_count[i,3]/(norm_file[tmp_i,]*norm_file[tmp_j,])
    if(is.na(normalized_mat[tmp_i,tmp_j])){
    normalized_mat[tmp_i,tmp_j] <- 0
    }
    pb$tick()
  }
  
  input_if <- normalized_mat                                                                                             
  input_if <- input_if + t(input_if)
  diag(input_if) <- diag(input_if)/2
  pd <- input_if^(alpha)
  res = new('flamingo',IF=input_if,PD=pd,n_frag=n,chr_name=chr_name)
  return(res)
}
