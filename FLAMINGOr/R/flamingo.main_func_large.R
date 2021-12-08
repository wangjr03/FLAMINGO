#' flamingo.main_func_large
#'
#' Main function of FLAMINGO, a wraper for all steps for super large data, i.e. 1kb Hi-C
#' @param hic_data_low Input Hi-C data in .hic format, .mcool format or sparse matrix format.
#' @param file_format foramt of the input hic data. Could be .hic, .mcool or sparse matrix format. 
#' @param domain_res Size of the domains in bps, e.g. 1e6.
#' @param frag_res Size of the fragment in bps, e.g. 5e3.
#' @param chr_size Size of the chromosome in bps.
#' @param chr_name Name of the chromosome, e.g. chr1.
#' @param normalization Normalization method.
#' @param downsampling_rates Fraction of contacts to be used during the reconstruction.
#' @param lambda Lagrigian coefficient.
#' @param max_dist Maximum allowed distance betwee two consecutive points.
#' @param nThread Number of thread avalable for the reconstruction.
#' @param alpha Convertion factor between interaction frequency and pairwise distance. default -0.25.
#' @param max_iter Maximum iteration for the assembling algorithm. default 500.
#' @param hic_data_high Optional. The high resolution HiC data in sparse matrix format. Only required if the file_format is 'sparse matrix'.
#' @param norm_low Optional. The normalization vector for the low resolution Hi-C data. Only required if the file_format is 'sparse matrix'.
#' @param norm_high Optional. The normalization vector for the high resolution Hi-C data. Only required if the file_format is 'sparse matrix'.
#' @param n_row Optional. Number of rows to use the large matrix format. reduce the number if memory is limited. default t0 45000 rows.
#' @keywords FLAMINGO
#' @return A data.frame containing the FLAMINGO predicted 3D structure.
#' @export

flamingo.main_func_large <- function(hic_data_low,file_format,domain_res,frag_res,chr_size,chr_name,normalization,downsampling_rates,lambda,max_dist,nThread,alpha=-0.25,max_iter=500,hic_data_high=NULL,norm_low=NULL,norm_high=NULL,n_row=45000){
  if(!file_format %in% c('sparse matrix','hic','mcool')){
    stop("file format must be one of: 'sparse matrix'','hic' or 'mcool'")
  }
  if(file_format=='sparse matrix'){
    if(is.null(hic_data_high)){
      stop('To use the sparse matrix format input data, please provide both low resolution data and high resolution data by setting parameters hic_data_low and hic_data_high')
    }
    if(is.null(norm_low) | is.null(norm_high)){
      stop('Please provide the normalizing vector file for both high resolution data and low resolution data')
    }
  }
  if(file_format=='hic'){
    flamingo_low_res_data_obj <- flamingo.construct_flamingo_from_hic(hic_file = hic_data_low,
                                                                  resolution = domain_res,
                                                                  chr_size=chr_size,
                                                                  chr_name=chr_name,
                                                                  normalization = normalization,
                                                                  alpha=alpha,
                                                                  n_row=n_row)
    flamingo_high_res_data_obj <- flamingo.construct_flamingo_from_hic(hic_file = hic_data_low,
                                                                  resolution = frag_res,
                                                                  chr_size=chr_size,
                                                                  chr_name=chr_name,
                                                                  normalization = normalization,
                                                                  alpha=alpha,
                                                                  n_row=n_row)
  }
  if(file_format=='mcool'){
    flamingo_low_res_data_obj <- flamingo.construct_flamingo_from_mcool(mcool_file = hic_data_low,
                                                                        resolution = domain_res,
                                                                        chr_name=chr_name,
                                                                        normalization = normalization,
                                                                        alpha=alpha,
                                                                        n_row=n_row)
    flamingo_high_res_data_obj <- flamingo.construct_flamingo_from_mcool(mcool_file = hic_data_low,
                                                                        resolution = frag_res,
                                                                        chr_name=chr_name,
                                                                        normalization = normalization,
                                                                        alpha=alpha,
                                                                        n_row=n_row)
  }
  if(file_format=='sparse matrix'){
    flamingo_low_res_data_obj <- flamingo.construct_flamingo_from_sparse_matrix(raw_count_file = hic_data_low,
                                                                                norm_file = norm_low,
                                                                                resolution = domain_res,
                                                                                chr_name=chr_name,
                                                                                chr_size=chr_size,
                                                                                alpha=alpha)
    flamingo_high_res_data_obj <- flamingo.construct_flamingo_from_sparse_matrix(raw_count_file = hic_data_high,
                                                                                norm_file = norm_high,
                                                                                resolution = frag_res,
                                                                                chr_name=chr_name,
                                                                                chr_size=chr_size,
                                                                                alpha=alpha)
  }
  print('Dividing domains...')
  flamingo.divide_domain(flamingo_obj = flamingo_high_res_data_obj,domain_res=domain_res,frag_res=frag_res)
  print('caching datasets...')
  write_huge_mat(flamingo_high_res_data_obj@IF,"tmp_IF.txt.gz",nThread=nThread)
  write_huge_mat(flamingo_high_res_data_obj@PD,"tmp_PD.txt.gz",nThread=nThread)
  high_res_n_frag = flamingo_high_res_data_obj@n_frag
  high_res_chr_name = flamingo_high_res_data_obj@chr_name
  rm(flamingo_high_res_data_obj)
  print('Reconstructing backbones...')
  flamingo_backbone_prediction = flamingo.reconstruct_backbone_structure(flamingo_data_obj = flamingo_low_res_data_obj,sw=downsampling_rates,lambda=lambda,max_dist = max_dist,nThread=1)
  print('Reconstructing intra-domain structures...')
  flamigo_intra_domain_prediction = flamingo.reconstruct_structure(sw=downsampling_rates,lambda = lambda,max_dist = max_dist,nThread=nThread)
  save(flamigo_intra_domain_prediction,file="intra_domain.Rdata")
  print('Assembling structures...')
  flamingo_high_res_IF = data.table::fread("tmp_IF.txt.gz",sep='\t',header=F,nThread=nThread)
  flamingo_high_res_IF = as.matrix(flamingo_high_res_IF)
  flamingo_high_res_PD = data.table::fread("tmp_PD.txt.gz",sep='\t',header=F,nThread=nThread)
  flamingo_high_res_PD = as.matrix(flamingo_high_res_PD)
  flamingo_high_res_data_obj = new('flamingo',IF=flamingo_high_res_IF,PD=flamingo_high_res_PD,n_frag=high_res_n_frag,chr_name=high_res_chr_name)
  res = flamingo.assemble_structure(flamingo_backbone_prediction_obj=flamingo_backbone_prediction,
                              flamingo_final_res_data_obj=flamingo_high_res_data_obj,
                              list_of_flamingo_domain_prediction_obj=flamigo_intra_domain_prediction,
                              max_iter=max_iter)
  return(res)
  
}