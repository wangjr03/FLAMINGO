#' flamingo.estimate_alpha
#'
#' Estimate an alpha value based on the input Hi-C data.
#' @param flamingo_obj FLAMINGO object containing the 5kb resolution interaction frequency matrix and pairwise distance matrix.
#' @keywords FLAMINGO
#' @return An alpha value.
#' @export
flamingo.estimate_alpha <- function(flamingo_obj){
	input_if = flamingo_obj@IF
	median_diag_if = median(input_if[row(input_if)-col(input_if)==1],na.rm=T)
	alpha = -0.25*log(48)/log(median_diag_if)
	return(alpha)
}

