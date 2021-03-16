# Data in examples/vignettes/tests
CORE_DATASETS  <- c('atkin18.somascan.adat',
                    'atkin18.metabolon.xlsx',
                    'billing16.bam.zip',
                    'billing16.rnacounts.txt',
                    'billing16.proteingroups.txt',
                    'billing16.somascan.adat',
                    'billing19.rnacounts.txt',
                    'billing19.proteingroups.txt',
                    'billing19.phosphosites.txt',
                    'fukuda20.proteingroups.txt',
                    'halama18.metabolon.xlsx')

#' Data used in examples/vignette/tests/longtests
#' @export
AUTONOMICS_DATASETS <- c(CORE_DATASETS, 
                        'uniprot_hsa_20140515.fasta')

#' Download autonomics example data
#' 
#' @param filename file name
#' \itemize{
#'     \item \href{https://pubmed.ncbi.nlm.nih.gov/26857143/}{
#'                 Billing 2016: stemcell comparison:} E, EM, BM
#'     \itemize{ 
#'         \item \code{'billing16.bam.zip'}
#'         \item \code{'billing16.rnacounts.txt'}
#'         \item \code{'billing16.somascan.adat'}
#'         \item \code{'billing16.proteingroups.txt'}
#'     }
#'     \item \href{https://pubmed.ncbi.nlm.nih.gov/30525282/}{
#'                 Atkin 2018: hypoglycemia:} t0, t1, t2, t3
#'     \itemize{
#'         \item \code{'atkin18.somascan.adat'}
#'         \item \code{'atkin18.metbolon.xlsx'}
#'     }
#'     \item \href{https://pubmed.ncbi.nlm.nih.gov/29777783/}{
#'                 Halama 2018: glutaminase inhibition: } 4 conc, 4 timepoints
#'     \itemize{
#'         \item \code{'halama18.metabolon.xlsx'}
#'     }
#'     \item \href{https://pubmed.ncbi.nlm.nih.gov/31332097/}{
#'            Billing 2019: stemcell differentiation:} 
#'            E00, E01, E02, E05, EM15, EM30, M00
#'     \itemize{
#'         \item \code{'billing19.rnacounts.txt'    }
#'         \item \code{'billing19.proteingroups.txt'}
#'         \item \code{'billing19.phosphosites.txt'} 
#'     }
#'     \item \href{https://pubmed.ncbi.nlm.nih.gov/32648304/}{
#'           Fukuda 2020: zebrafish development:} X30dpt, Adult
#'     \itemize{
#'         \item \code{'fukuda20.proteingroups.txt'}
#'     }
#' }
#' @param verbose TRUE / FALSE
#' @return local file path
#' @examples
#' # atkin18 - hypoglycemia - pubmed 30525282
#'     download_data('atkin18.somascan.adat')          # somascan  intensities
#'     download_data('atkin18.metabolon.xlsx')         # metabolon intensities
#'
#' # billing16 - stemcell characterization - pubmed 26857143
#'     download_data('billing16.proteingroups.txt')  # proteingroup ratios
#'     download_data('billing16.somascan.adat')      # somascan  intensities
#'     download_data('billing16.rnacounts.txt')      # rnaseq    counts
#'     download_data('billing16.bam.zip')            # rnaseq    alignments
#'
#' # billing19 - stemcell differentiation - pubmed 31332097
#'     # download_data('billing19.proteingroups.txt')    # proteingroup ratios
#'     # download_data('billing19.phosphosites.txt')     # phosphosite  ratios
#'     # download_data('billing19.rnacounts.txt')        # rnaseq       counts
#'
#' # fukuda20 - heart regeneration - pubmed PXD016235
#'     download_data('fukuda20.proteingroups.txt')     # proteingroup LFQ
#'
#' # halama18 - glutaminase inhibition - pubmed 30525282
#'     download_data('halama18.metabolon.xlsx')        # metabolon intensities
#' @export
download_data <- function(filename, verbose = TRUE){
    assert_is_subset(filename, AUTONOMICS_DATASETS)
    fileURL <- paste0(
        "https://bitbucket.org/graumannlabtools/autonomics/downloads/", 
        filename)
    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, filename, "rname")$rid
    if (!length(rid)) {
        if(verbose) message( "Downloading ", filename)
        rid <- names(BiocFileCache::bfcadd(bfc, filename, fileURL ))
    }
    #if (!isFALSE(bfcneedsupdate(bfc, rid)))
    #bfcdownload(bfc, rid)
    filepath <- BiocFileCache::bfcrpath(bfc, rids = rid)
    if (tools::file_ext(filename)=='zip'){
        if (verbose)  message('unzip')
        utils::unzip(filepath, exdir = substr(filepath, 1, nchar(filepath)-4))
        filepath %<>% substr(1, nchar(.)-4)
    }
    filepath
}

.get_cache <- function(){
    cache <- rappdirs::user_cache_dir(appname="autonomics")
    BiocFileCache::BiocFileCache(cache)
}

download_data_old <- function(
    file, localdir = '~/autonomicscache/datasets', unzip = TRUE, ntries = 3
){
    assert_is_subset(file, AUTONOMICS_DATASETS)
    . <- NULL
    bitbucket <- 'https://bitbucket.org/graumannlabtools/autonomics/downloads'
    localdir  %<>% paste(vapply(
        stri_split_fixed(file, '.'), extract, character(1), 1), sep = '/')
    dir.create(localdir, showWarnings = FALSE, recursive = TRUE)
    localfile <- paste0(localdir,  '/', file)
    if (!file.exists(localfile)){
        bitbucketfile <- paste0(bitbucket, '/', file)
        while (ntries > 0){  # bioconductor.org/developers/how-to/web-query
            result <- tryCatch(
                        download.file(bitbucketfile, localfile, mode = 'wb'), 
                        error = identity)
            if (!inherits(result, 'error'))  break
            ntries %<>% magrittr::subtract(1) }
        assertive::assert_all_are_not_equal_to(ntries, 0) }
    if (file_ext(file) == 'zip'){
        localdir <- file_path_sans_ext(localfile)
        if (!dir.exists(localdir)) unzip(localfile, exdir = localdir)
        localfile %<>% substr(1, nchar(.)-4) }
    return(invisible(localfile))
}
