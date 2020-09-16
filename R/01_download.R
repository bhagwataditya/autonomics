#' Available autonomics datasets
#' @noRd
IMPORTOMICS_DATASETS <- c(
    'uniprot_hsa_20140515.fasta',
    'stemcells.bam.zip',
    'stemcells.rnacounts.txt',
    'stemcells.proteinGroups.txt',
    'stemcells.somascan.adat',
    'differentiation.rnacounts.txt',
    'differentiation.proteinGroups.txt',
    'differentiation.phosphoSites.txt',
    'hypoglycemia.somascan.adat',
    'hypoglycemia.metabolon.xlsx',
    'glutaminase.metabolon.xlsx')

#' Download autonomics data
#' @param file      name of file to download
#' @param localdir  directory where results will be saved
#' @param unzip     TRUE (default) or FALSE: whether to unzip
#' @return return   localfile invisibly
#' @examples
#' print(download_data(file = 'stemcells.somascan.adat'))
#' @export
download_data <- function(
    file, localdir = '~/importomicscache/datasets', unzip = TRUE
){
    assert_is_subset(file, IMPORTOMICS_DATASETS)
    . <- NULL

    bitbucket <- 'https://bitbucket.org/graumannlabtools/importomics/downloads'
    localdir  %<>% paste(vapply(stri_split_fixed(file, '.'), extract, character(1), 1), sep = '/')
    dir.create(localdir, showWarnings = FALSE, recursive = TRUE)
    localfile <- paste0(localdir,  '/', file)
    if (file.exists(localfile)){
        message('Use already available file: ', localfile)
    } else {
        download.file(paste0(bitbucket, '/', file), localfile, mode = 'wb')
    }

    if (file_ext(file) == 'zip'){
        unzip(localfile, exdir = dirname(localfile))
        localfile %<>% substr(1, nchar(.)-4)
    }

    return(invisible(localfile))
}


