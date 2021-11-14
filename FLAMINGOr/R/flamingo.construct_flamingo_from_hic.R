#' flamingo.construct_flamingo_from_hic
#'
#' Generate FLAMINGO input data from the .hic format data.
#' @param hic_file Path to the hic file.
#' @param normalization Method of normalization. Must be one of "NONE", "VC", "VC_SQRT", "KR".
#' @param resolution Target resolution.
#' @param chr_name chromosome name
#' @param chr_size chromsome size
#' @param alpha convertion factor between interaction frequency and distance, default is -0.25
#' @param n_row number of rows to use the large matrix mode, default to 45000
#' @keywords FLAMINGO
#' @return A FLAMINGO object containing the resulted Interaction Frequency (IF) matrix, Pairwise Distance (PD) matrix and number of fragments.
#' @export
flamingo.construct_flamingo_from_hic <- function(hic_file,normalization,resolution,chr_name,chr_size,alpha=-0.25,n_row=45000){
  options(scipen = 999)
  chr_number <- gsub("chr","",chr_name)
  normalized_data = strawr::straw(normalization,hic_file,chr_number,chr_number,unit='BP',binsize=resolution)
  n <- ceiling(chr_size/resolution)
  i_ind <- 1+(normalized_data[,1]/resolution)
  j_ind <- 1+(normalized_data[,2]/resolution)
  input_if = Matrix::sparseMatrix(i=i_ind,j=j_ind,x=normalized_data[,3],dims=c(n,n))
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
