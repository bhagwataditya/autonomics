#' Available autonomics datasets
#' @noRd
IMPORTOMICS_DATASETS <- c(
    'uniprot_hsa_20140515.fasta',
    'billing16.bam.zip',
    'billing16.rnacounts.txt',
    'billing16.proteingroups.txt',
    'billing16.somascan.adat',
    'billing19.rnacounts.txt',
    'billing19.proteingroups.txt',
    'billing19.phosphosites.txt',
    'fukuda20.proteingroups.txt',
    'atkin18.somascan.adat',
    'atkin18.metabolon.xlsx',
    'halama18.metabolon.xlsx')

#' Download example datasets
#' @param file      name of file to download
#' @param localdir  directory where results will be saved
#' @param unzip     TRUE (default) or FALSE: whether to unzip
#' @return return   localfile invisibly
#' @examples
#' # atkin18 - hypoglycemia - pubmed 30525282
#'     download_data('atkin18.somascan.adat')        # somascan  intensities
#'     download_data('atkin18.metabolon.xlsx')       # metabolon intensities
#'
#' # billing16 - stemcell characterization - pubmed 26857143  pride PXD001856
#'     download_data('billing16.proteingroups.txt')  # proteingroup ratios
#'     download_data('billing16.somascan.adat')      # somascan     intensities
#'     download_data('billing16.rnacounts.txt')      # rnaseq       counts
#'     download_data('billing16.bam.zip')          # rnaseq       alignments
#'
#' # billing19 - stemcell differentiation - pubmed 31332097  pride PXD004652
#'     download_data('billing19.proteingroups.txt')  # proteingroup ratios
#'     download_data('billing19.phosphosites.txt')   # phosphosite  ratios
#'     download_data('billing19.rnacounts.txt')      # rnaseq       counts
#'
#' # fukuda20 - heart regeneration - pubmed PXD016235  pride PXD016235
#'     download_data('fukuda20.proteingroups.txt')   # proteingroup LFQ intens.
#'
#' # halama18 - glutaminase inhibition - pubmed 30525282
#'     download_data('halama18.metabolon.xlsx')      # metabolon intensities
#' @export
download_data <- function(
    file, localdir = '~/importomicscache/datasets', unzip = TRUE
){
    assert_is_subset(file, IMPORTOMICS_DATASETS)
    . <- NULL

    bitbucket <- 'https://bitbucket.org/graumannlabtools/importomics/downloads'
    localdir  %<>% paste(vapply(
        stri_split_fixed(file, '.'), extract, character(1), 1), sep = '/')
    dir.create(localdir, showWarnings = FALSE, recursive = TRUE)
    localfile <- paste0(localdir,  '/', file)
    if (!file.exists(localfile))    download.file(
        paste0(bitbucket, '/', file), localfile, mode = 'wb')

    if (file_ext(file) == 'zip'){
        if (!dir.exists(dirname(localfile))){
            unzip(localfile, exdir = dirname(localfile))
        }
        localfile %<>% substr(1, nchar(.)-4)
    }

    return(invisible(localfile))
}


