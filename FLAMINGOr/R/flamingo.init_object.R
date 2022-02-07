#' flamingo.construct_flamingo_from_sparse_matrix
#'
#' Create an object named flamingo for data organization and an object named flamingo_prediction for result organization
#' @keywords FLAMINGO
#' @return Objects named flamingo and flamingo_prediction are ready for use.
#' @export
setClass("flamingo", slots=list(IF="matrix", PD='matrix',n_frag='numeric',chr_name='character'))
setClass("flamingo_prediction", slots=list(id="numeric", coordinates='matrix',input_n='numeric'))
setClass("flamingo_large", slots=list(IF="dgCMatrix", PD='dgCMatrix',n_frag='numeric',chr_name='character'))