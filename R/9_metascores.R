#' Progeny
#'
#' Find pathway responsive genes for activity inference
#'
#' @param  object     SummarizedExperiment
#' @param  symbolvar  string: gene symbol var
#' @return SummarizedExperiment (sdata columns added)
#' @seealso https://saezlab.github.io/progeny/
#' @examples
#' require(magrittr)
#' file <- autonomics::download_data('billing16.rnacounts.txt')
#' object <- read_rnaseq_counts(file)
#' object %<>% progeny(symbolvar = 'gene_name')
#' object
#' @export
fit_progeny <- function(object, symbolvar){
# Assert
    if (!requireNamespace('progeny', quietly = TRUE)){
        message("\tBiocManager::install('progeny'")
        return(object) }
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset('log2counts', assayNames(object))
    assert_is_subset(symbolvar, fvars(object))
# Run
    uniquesymbols <- data.table(fdata(object))[, feature_id[1], by = symbolvar]$V1
    obj <- object[uniquesymbols, ]
    fnames(obj) <- as.character(fdata(obj)[[symbolvar]])
    model_human_full <- progeny::model_human_full
    results <- progeny::progeny(log2counts(obj))
# Merge
    results %<>% mat2dt('sample_id')
    names(results)[-1] %<>% paste0('progeny:', .)
    object %<>% merge_sdata(results)
    object
}

