#=========================================================================

#' Log2 transform exprs
#' @param object SummarizedExperiment
#' @param verbose TRUE or FALSE
#' @return Updated SummarizedExperiment
#' @noRd
log2transform <- function(object, verbose = FALSE){
    if (verbose) message('\t\tLog2 transform exprs')
    exprs(object) %<>% log2()
    if (!is.null(occupancies(object))){
    if (verbose) message('\t\tLog2 transform occupancies')
        occupancies(object) %<>% log2()
    }
    return(object)
}


