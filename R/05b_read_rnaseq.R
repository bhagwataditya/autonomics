#=========================================================
# RNASEQ
#=========================================================
#download_gtf
release_to_build <- function(release, organism){
    if        (organism == 'Homo sapiens'){
        if (release >= 76)  'GRCh38'   else 'GRCh37'  }

    else if   (organism == 'Mus musculus'){
        if (release >= 68)  'GRCm38'   else 'NCBIM37'  }

    else if   (organism == 'Rattus norvegicus'){
        if (release >= 80)  'Rnor_6.0' else 'Rnor_5.0'  }
}

#-----------------------------------------------------------
# Following GTF functions are soft-deprecated.
# Better to outsource this functionality to biomartr::getGTF
#-----------------------------------------------------------

#' Make link to GTF file
#' @param organism 'Homo sapiens', 'Mus musculus', or 'Rattus norvegicus'
#' @param release   number
#' @examples
#' make_gtf_url(organism = 'Homo sapiens',            release = 95)
#' make_gtf_url(organism = 'Mus musculus',            release = 95)
#' make_gtf_url(organism = 'Sacharomyces cerevisiae', release = 100)
#' @noRd
make_gtf_url <- function(organism, release){
    sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/%s/%s.%s.%s.gtf.gz',
            release,
            stri_replace_first_fixed(tolower(organism), ' ', '_'),
            stri_replace_first_fixed(        organism,  ' ', '_'),
            release_to_build(release, organism),
            release)
}


#' Download GTF file
#'
#' Download GTF file with feature annotations
#' @param organism  'Homo sapiens', 'Mus musculus' or 'Rattus norvegicus'
#' @param release    GTF release. By default release 95 selected
#' @param gtffile    string: path to local GTF file
#' @examples
#' \dontrun{ # requires internet and does not always work:
#'           # https://stackoverflow.com/questions/55532102
#'    download_gtf(organism = 'Homo sapiens')
#'    download_gtf(organism = 'Mus musculus')
#'    download_gtf(organism = 'Rattus norvegicus')
#' }
#' @noRd
download_gtf <- function(
    organism,
    release = 100,
    gtffile = sprintf("~/autonomicscache/gtf/%s",
        basename(make_gtf_url(organism, release) %>% substr(1, nchar(.)-3)))
){
    assert_is_subset(organism,
                    c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'))
    . <- NULL
    remote <- make_gtf_url(organism, release)

    if(file.exists(gtffile)){
        message(sprintf("GTF file already available at %s", gtffile))
    } else {
        message(sprintf("download   %s'" , remote))
        message(sprintf("to         %s", gtffile ))
        dir.create(dirname(gtffile), showWarnings = FALSE, recursive = TRUE)
        tryCatch(expr = { download.file(
                                url = remote, destfile = gtffile, quiet = TRUE)
                            gunzip(gtffile,  remove = TRUE, overwrite = TRUE)},
                error = function(cond){
                    message('failed     repeat manually using browser\n');
                    message(cond)})
    }
    gtffile %<>% stri_replace_last_fixed('.gz', '')
    invisible(gtffile)
}


#' Read GTF file into data.table
#'
#' Read GTF file into a data.table. Filter for particular values of variable.
#'
#' @param gtffile    string: path to gtffile
#' @param var        string: variable on which to filter
#' @param values     filter variable for these values. If NULL, not filtering.
#' @param writefile  string: file to write gtf table to
#' @examples
#' \dontrun{ # requires internet connection
#'    require(magrittr)
#'    gtffile <- download_gtf(organism = 'Homo sapiens', release = 95)
#'    gtfdt <- read_gtf(gtffile, var = 'gene_id', values = 'ENSG00000198947')
#' }
#' @noRd
read_gtf <- function(
    gtffile,
    var       = 'gene_id',
    values    = NULL,
    writefile = NULL
){

    # Assert
    assert_all_are_existing_files(gtffile)
    assert_is_a_string(var)

    # Read
    dt <- rtracklayer::import(gtffile) %>%
          GenomicRanges::as.data.frame() %>% data.table()

    # Filter
    if (!is.null(values)){
        dt %>% setkeyv(var)
        dt %<>% extract(values)
    }

    # Write
    if (!is.null(writefile)){
        message(sprintf("\t\tWrite   %s", writefile))
        dir.create(dirname(writefile), recursive = TRUE, showWarnings = FALSE)
        fwrite(dt, writefile)
    }

    # Return
    dt

}

#' Read BAM files into SummarizedExperiment
#'
#' @param bamdir       string: path to SAM or BAM file directory
#'                    (one SAM or BAM file per sample)
#' @param ispaired     TRUE or FALSE (default): paired end reads?
#' @param gtffile      NULL (use Rsubread's default human/mouse annotations) or
#'                     string (path to GTF file)
#' @param fvars        character vector: GTF variables to include in object.
#' @param nthreads     number of cores to be used by Rsubread::featureCounts()
#' @param ...          passed to Rsubread::featureCounts
#' @return SummarizedExperiment
#' @examples
#' bamdir <- download_autonomics_data("stemcells.bam.zip")
#' read_bam(bamdir, ispaired = TRUE)
#' @export
read_bam <- function(bamdir, ispaired = FALSE, gtffile = NULL,
    fvars = character(0), nthreads   = detectCores(), ...
){
    # Assert
    assert_all_are_existing_files(bamdir)
    assert_is_a_bool(ispaired)
    if (!is.null(gtffile))   assert_all_are_existing_files(gtffile)
    assert_is_a_number(nthreads)

    # Count reads
    files <- list.files(
        bamdir, pattern = ".sam$|.bam$", full.names = TRUE, recursive = TRUE)
    fcounts <-  featureCounts(
                    files               =  files,
                    annot.ext           =  gtffile,
                    isGTFAnnotationFile = !is.null(gtffile),
                    GTF.attrType.extra  =  if(length(fvars)==0) NULL else fvars,
                    isPairedEnd         =  ispaired,
                    nthreads            =  nthreads,
                    ...)

    # Forge SummarizedExperiment
    filenames   <- basename(file_path_sans_ext(files))
    subdirnames <- basename(dirname(files))
    sample_names <- if (has_no_duplicates(filenames)){
                        filenames
                    } else if (has_no_duplicates(subdirnames)){
                        subdirnames
                    } else {
                        paste0(subdirnames, '_', filenames)
                    }
    object <- SummarizedExperiment(assays = list(
        exprs = fcounts$counts %>% set_colnames(sample_names)))

    # Add sdata
    message("\t\tAdd sdata")
    colData(object) <- DataFrame(
                        sample_id = sample_names, row.names = sample_names)
    object$subgroup <- guess_subgroup_values(object, verbose = FALSE)

    # Add fdata
    message("\t\tAdd fdata")
    rowData(object) <- fcounts$annotation[ , c('GeneID', fvars), drop = FALSE]
    rownames(object) <- rowData(object)$GeneID
    fvars(object) %<>% stri_replace_first_fixed('GeneID', 'feature_id')

    # Return
    object

}

#' Read rnaseq counts
#'
#' Read tsv file with rnaseq counts into SummarizedExperiment
#'
#' File format: header row
#'              feature annotations in first few columns
#'              feature counts      in next columns
#'
#' @param file      string: path to rnaseq counts file
#' @param fid_var   string or number: feature id variable
#' @param fname_var string or number: feature name variable
#' @return SummarizedExperiment
#' @examples
#' file <- download_autonomics_data('stemcells_rna.txt')
#' read_counts(file, fid_var = 'gene_id', fname_var = 'gene_name')
#' @seealso merge_sdata, merge_fdata
#' @export
read_counts <- function(
    file,
    fid_var,
    fname_var = character(0)
){
    assert_all_are_existing_files(file)
    dt <- fread(file, integer64='numeric')

    assert_is_subset(fid_var, names(dt))
    fid_col <- which(names(dt)==fid_var)
    expr_cols   <- which(unname(vapply(dt, is.integer, logical(1))))
    fdata_cols  <- c(fid_col,
                     1 + which(unname(!vapply(
                                        dt[, -fid_col, with = FALSE],
                                        is.integer,
                                        logical(1)))))
    object <- read_omics(
                file,
                fid_rows   = 2:nrow(dt),   fid_cols   = fid_col,
                sid_rows   = 1,            sid_cols   = expr_cols,
                expr_rows  = 2:nrow(dt),   expr_cols  = expr_cols,
                fvar_rows  = 1,            fvar_cols  = fdata_cols,
                fdata_rows = 2:nrow(dt),   fdata_cols = fdata_cols,
                transpose  = FALSE,
                verbose    = TRUE)

    sdata(object)$subgroup  <- guess_subgroup_values(object, verbose = TRUE)
    if (length(fname_var)>0){
        assert_is_subset(fname_var, fvars(object))
        fdata(object)$feature_name <- fdata(object)[[fname_var]]
        fdata(object) %<>% pull_columns(c('feature_id', 'feature_name'))
    }

    object

}

