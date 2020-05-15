#========================
# AFFYMETRIX MICROARRAYS
#========================

# https://stackoverflow.com/a/4090208
# install_if_required <- function(pkgs){
#     pkgs %<>% extract(!(pkgs %in% installed.packages()[,"Package"]))
#     if(length(pkgs)) BiocManager::install(pkgs)
# }


add_affy_fdata <- function(object){

    # Extract entrez identifiers
    entrezgs <- vapply(
      stri_split_fixed(fnames(object), '_'), extract, character(1), 1)

    # Get annotation db
    pkgname <- paste0(metadata(object)$annotation, '.db')
    #install_if_required(pkgname)
    db <- getFromNamespace(pkgname, pkgname)

    # Map
    rowData(object) <- DataFrame(
        feature_id    = fnames(object),
        feature_name  = suppressMessages(mapIds(
                    db, entrezgs, column = 'SYMBOL',   keytype = 'ENTREZID')),
        feature_descr = suppressMessages(mapIds(
                    db, entrezgs, column = 'GENENAME', keytype = 'ENTREZID')),
        row.names     = fnames(object)
    )

    # Return
    object

}

#' Read affymetrix microarray
#' @param celfiles string vector: CEL file paths
#' @return RangedSummarizedExperiment
#' @examples
#' require(magrittr)
#' url <- paste0('http://www.bioconductor.org/help/publications/2003/',
#'                 'Chiaretti/chiaretti2/T33.tgz')
#' localfile <- file.path('~/importomicscache', basename(url))
#' if (!file.exists(localfile)){
#'     download.file(url, destfile = localfile)
#'     untar(localfile, exdir = path.expand('~/importomicscache'))
#' }
#' localfile %<>% substr(1, nchar(.)-4)
#' read_affymetrix(celfiles = list.files(localfile, full.names = TRUE))
#' @export
read_affymetrix <- function(celfiles){

    # read
    message('Read Affymetrix CEL files: ', basename(celfiles)[1], ', ...')
    suppressWarnings(eset1 <- just.rma(filenames = celfiles))
    object <- makeSummarizedExperimentFromExpressionSet(eset1)

    # sdata
    snames(object) %<>% stri_replace_first_fixed('.CEL', '')
    sdata(object) <- data.frame(sample_id = snames(object),
                                row.names = snames(object),
                                stringsAsFactors = FALSE)
    # fdata
    object %<>% add_affy_fdata()

    # return
    return(object)

}



