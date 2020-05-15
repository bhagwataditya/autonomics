#' Available autonomics datasets
#' @noRd
AUTONOMICS_DATASETS <-
c( 'stemcells_rna.txt', 'stemcells_proteinGroups.txt', 'stemcells_soma.adat',
    'stemcells.bam.zip',
    'diff_rna.txt', 'diff_proteinGroups.txt', 'diff_phosphoSites.txt',
    'hypo_soma.adat', 'hypo_metab.xlsx',
    'glutaminase_metab.xlsx')

#' Download autonomics data
#' @param file      name of file to download
#' @param localdir  directory where results will be saved
#' @param unzip     TRUE (default) or FALSE: whether to unzip
#' @return return   localfile invisibly
#' @examples
#' print(download_autonomics_data(file = 'stemcells_soma.adat'))
#' @export
download_autonomics_data <- function(
    file, localdir = '~/autonomicscache', unzip = TRUE
){
    assert_is_subset(file, AUTONOMICS_DATASETS)
    dir.create(localdir, showWarnings = FALSE, recursive = TRUE)

    bitbucket <- 'https://bitbucket.org/graumannlabtools/autonomics/downloads'
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


