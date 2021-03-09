#' Available autonomics datasets
#' @export
AUTONOMICS_DATASETS <- c(
    'atkin18.somascan.adat',
    'atkin18.metabolon.xlsx',
    'billing16.bam.zip',
    'billing16.rnacounts.txt',
    'billing16.proteingroups.txt',
    'billing16.somascan.adat',
    'billing19.rnacounts.txt',
    'billing19.proteingroups.txt',
    'billing19.phosphosites.txt',
    'fukuda20.proteingroups.txt',
    'halama18.metabolon.xlsx',
    'uniprot_hsa_20140515.fasta')

#' Download example dataset
#' 
#' @param file      string in \code{AUTONOMICS_DATASETS}: filename to download
#' @param localdir  directory where results will be saved
#' @param unzip     TRUE (default) or FALSE: whether to unzip
#' @param ntries    no of times to retry
#' @return return localfile invisibly
#' @examples
#' # atkin18 - hypoglycemia - pubmed 30525282
#'     download_data('atkin18.somascan.adat')        # somascan  intensities
#'     download_data('atkin18.metabolon.xlsx')       # metabolon intensities
#'
#' # billing16 - stemcell characterization - pubmed 26857143  pride PXD001856
#'     download_data('billing16.proteingroups.txt')  # proteingroup ratios
#'     download_data('billing16.somascan.adat')      # somascan     intensities
#'     download_data('billing16.rnacounts.txt')      # rnaseq       counts
#'     download_data('billing16.bam.zip')            # rnaseq       alignments
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
download_data <- function(filename, verbose = FALSE ){
    fileURL <- paste0(
        "https://bitbucket.org/graumannlabtools/autonomics/downloads/", 
        filename)
    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, filename, "rname")$rid
    if (!length(rid)) {
     if( verbose )
         message( "Downloading ", filename)
     rid <- names(BiocFileCache::bfcadd(bfc, filename, fileURL ))
    }
    #if (!isFALSE(bfcneedsupdate(bfc, rid)))
    #bfcdownload(bfc, rid)
    filepath <- BiocFileCache::bfcrpath(bfc, rids = rid)
    if (tools::file_ext(filename)=='zip'){
        if (verbose)  cmessage('unzip')
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
