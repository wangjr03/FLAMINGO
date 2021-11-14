#' flamingo.construct_flamingo_from_mcool
#'
#' Generate FLAMINGO input data from the .mcool format data.
#' @param mcool_file Path to the .mcool file.
#' @param normalization Method of normalization.
#' @param resolution Target resolution.
#' @param chr_name chromosome name
#' @param alpha convertion factor between interaction frequency and distance, default is -0.25
#' @param n_row number of rows to use the large matrix mode, default to 45000
#' @keywords FLAMINGO
#' @return A FLAMINGO object containing the resulted Interaction Frequency (IF) matrix, Pairwise Distance (PD) matrix and number of fragments.
#' @export
flamingo.construct_flamingo_from_mcool <- function(mcool_file,normalization,resolution,chr_name,alpha=-0.25,n_row=45000){
  options(scipen = 999)
  all_dir = rhdf5::h5ls(mcool_file)
  parent_dir = all_dir[1,2]
  target_dir = paste(c("",parent_dir,resolution),collapse='/')
  mcool_dat = rhdf5::h5read(mcool_file,target_dir)
  available_normalization = setdiff(names(mcool_dat$bins),c('chrom','start','end','weight'))
  if(!normalization %in% available_normalization){
    stop(
      paste('Normalization method not exist in the .mcool data! Must be one of: ',paste(available_normalization,collapse=', '))
         )
  }
  csr_rawcount = data.frame(bin_1 = mcool_dat$pixels$bin1_id,
                            bin_2 = mcool_dat$pixels$bin2_id,
                            value = mcool_dat$pixels$count)
  normalization_file = mcool_dat$bins[[normalization]]
  chr_id <- which(mcool_dat$bins$chrom == chr_name)
  offset <- mcool_dat$indexes$chrom_offset[which(mcool_dat$chroms$name==chr_name)]
  n <- length(chr_id)
  csr_rawcount <- subset(csr_rawcount,csr_rawcount[,1] %in% chr_id & csr_rawcount[,2] %in% chr_id)
  csr_rawcount[,1] <- csr_rawcount[,1]-offset+1
  csr_rawcount[,2] <- csr_rawcount[,2]-offset+1
  csr_rawcount <- as.matrix(csr_rawcount)
  normalization_file <- normalization_file[chr_id]
  for(i in 1:(dim(csr_rawcount)[1])){
    csr_rawcount[i,3] <- csr_rawcount[i,3]/(normalization_file[csr_rawcount[i,1]]*normalization_file[csr_rawcount[i,2]])
  }
  input_if <- Matrix::sparseMatrix(i=csr_rawcount[,1],j=csr_rawcount[,2],x=csr_rawcount[,3],dims=c(n,n))
  if(n<n_row){
    input_if <- as.matrix(input_if)
  }else{
    input_if <- convert_huge_mat(input_if)
  }
  input_if <- input_if + t(input_if)
  diag(input_if) <- diag(input_if)/2
  pd <- input_if^(alpha)
  res = new('flamingo',IF=input_if,PD=pd,n_frag=n,chr_name=chr_name)
  return(res)
  
}
