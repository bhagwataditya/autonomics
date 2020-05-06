#' Make vector components unique by appending spaces
#' @param x character or factor vector
#' @param method 'make.unique' or 'make.unique.spaces'
#' @param verbose TRUE (default) or FALSE
#' @return character vector
#' @seealso \code{\link[base]{make.unique}}
#' @examples
#' x <- c('A', 'B', 'C', 'A', 'D')
#' uniquify(x, 'make.unique')
#' uniquify(x, 'make.unique.spaces')
#' @noRd
uniquify <- function(x, method = 'make.unique.spaces', verbose = TRUE){
    idx <- cduplicated(x)
    if (any(idx)){
        uniquefun <- get(method)
        uniquex <- uniquefun(x)
        if (verbose){
            cmessage(   '\t\tUniquify ( %s -> %s ) duplicates of',
                        x[idx][1], uniquefun(x[idx][c(1, 1)])[2])
            cmessage_df("\t\t\t  %s", table(x[idx]) %>% as.list() %>% unlist)
                # unlist(as.list(.)) is to prevent empty line
        }
    } else {
        uniquex <- x
    }
    uniquex
}

#' Convenient (two way) duplicated
#' @param x vector
#' @return logical vector
#' @examples
#' require(magrittr)
#' c(1,2,3,4,5,2) %>% cduplicated()
#' @noRd
cduplicated <- function(x){
  duplicated(x) | duplicated(x, fromLast = TRUE)
}


#=================================================
# GENERIC
#=================================================

is_excel_file <- function(file){
    stringi::stri_detect_fixed(
            tools::file_ext(file), 'xls', case_insensitive = TRUE)
}

nrows <- function(x, sheet=1){
    if (is_excel_file(x)){
        nrow(readxl::read_excel(
                x, sheet = sheet, .name_repair = 'minimal', col_names = FALSE))
    } else {
        length(readLines(x, warn=FALSE))
    }
}


ncols <- function(x, sheet=1){
    if (is_excel_file(x)){
        ncol(readxl::read_excel(
            x, sheet=sheet, .name_repair = 'minimal', col_names = FALSE))
    } else {
        max(utils::count.fields(x, quote = "", sep = '\t'))
    }
}

#' Available autonomics datasets
#' @export
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
    assertive::assert_is_subset(file, AUTONOMICS_DATASETS)
    dir.create(localdir, showWarnings = FALSE, recursive = TRUE)

    bitbucket <- 'https://bitbucket.org/graumannlabtools/autonomics/downloads'
    localfile <- paste0(localdir,  '/', file)
    if (file.exists(localfile)){
        message('Use already available file: ', localfile)
    } else {
        download.file(paste0(bitbucket, '/', file), localfile, mode = 'wb')
    }

    if (tools::file_ext(file) == 'zip'){
        utils::unzip(localfile, exdir = dirname(localfile))
        localfile %<>% substr(1, nchar(.)-4)
    }

    return(invisible(localfile))
}


#' Is this a fixed column file?
#' @param file string: file name
#' @return TRUE or FALSE
#' @examples
#' # SOMASCAN
#' #---------
#' file <- paste0('../../wcmq/atkin.hypo/atkin.hypoglycemia/extdata/soma/',
#'     'WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.original.adat')
#' is_fixed_col_file(file)
#'
#' # METABOLON
#' #----------
#' file <- download_autonomics_data('glutaminase_metab.xlsx')
#' is_fixed_col_file(file)
#' @noRd
is_fixed_col_file <- function(file){
    if (is_excel_file(file)){ TRUE
    } else {
        fields <- utils::count.fields(file, quote = "", sep = '\t')
        all(fields==fields[1])
    }
}

#=================================================
# extract_rectangle
#=================================================

#' Extract rectangle from omics file, data.table, or matrix
#'
#' @param x          omics datafile or datatable
#' @param rows       numeric vector
#' @param cols       numeric vector
#' @param verbose    logical
#' @param transpose  logical
#' @param drop       logical
#' @param sheet      numeric or string
#' @param ...        allow for S3 method dispatch
#' @return matrix
#' @examples
#' # FROM FILE: extract_rectangle.character
#' #=======================================
#' # exprs
#'    require(magrittr)
#'    x <- download_autonomics_data('glutaminase_metab.xlsx')
#'    extract_rectangle(x, rows = 11:401, cols = 15:86, sheet = 2) %>%
#'    extract(1:3, 1:3)
#'
#' # fids
#'    extract_rectangle(x, rows = 11:401, cols = 5, sheet = 2)     %>%
#'    extract(1:3, )
#'
#' # sids
#'    extract_rectangle(x, rows = 2, cols = 15:86, sheet = 2)      %>%
#'    extract(,1:3)
#'
#' # fdata
#'    extract_rectangle(x, rows = 10:401, cols = 1:14,  sheet = 2) %>%
#'    extract(1:3, 1:3)
#'
#' # sdata
#'    extract_rectangle(x, rows = 1:10,   cols = 14:86, sheet = 2,
#'    transpose = TRUE) %>% extract(1:3, 1:3)
#'
#' # FROM MATRIX: extract_rectangle.matrix
#' #======================================
#' # exprs
#'    x <-download_autonomics_data('glutaminase_metab.xlsx') %>%
#'        extract_rectangle(sheet = 2)
#'    extract_rectangle(x, rows = 11:401, cols = 15:86, sheet = 2) %>%
#'    extract(1:3, 1:3)
#'
#' # fids
#'    extract_rectangle(x, rows = 11:401, cols = 5,     sheet = 2) %>%
#'    extract(1:3, )
#'
#' # sids
#'    extract_rectangle(x, rows = 2,      cols = 15:86, sheet = 2) %>%
#'    extract(,1:3)
#'
#' # fdata
#'    extract_rectangle(x, rows = 10:401, cols = 1:14,  sheet = 2) %>%
#'    extract(1:3, 1:3)
#'
#' # sdata
#'    extract_rectangle(x, rows = 1:10,   cols = 14:86, sheet = 2,
#'    transpose = TRUE) %>% extract(1:3, 1:3)
#' @importFrom magrittr %>% %<>%
#' @noRd
extract_rectangle <- function(x, ...){
    UseMethod('extract_rectangle')
}


extract_rectangle.character <- function(
    x, rows = 1:nrows(x, sheet=sheet), cols = 1:ncols(x, sheet=sheet),
    verbose = FALSE, transpose = FALSE, drop = FALSE, sheet = 1, ...
){
    # Assert
    assertive.files::assert_all_are_existing_files(x)
    assertive.types::assert_is_numeric(rows)
    assertive.types::assert_is_numeric(cols)
    row1 <- min(rows); rown <- max(rows)
    col1 <- min(cols); coln <- max(cols)

    # Read file
    if (verbose) message('\t\tRead ',
                if (sheet==1) '' else paste0("sheet '", sheet, "' of"), file)
    dt <- if (is_excel_file(x)){
            data.table::data.table(readxl::read_excel(
                x,
                sheet       = sheet,
                col_names   = FALSE,
                range       = sprintf('R%dC%d:R%dC%d', row1, col1, rown, coln),
                .name_repair = 'minimal'))
        } else {
            data.table::fread(
                x,
                na.strings = "",
                header     = FALSE,
                integer64  = 'numeric',
                skip       = row1-1,
                nrow       = 1+rown-row1,
                select     = col1:coln) }

    # Extract rectangle
    extract_rectangle.data.table(dt, transpose = transpose, drop = drop)

}

# Extract row
extract_dt_row <- function(dt, i) unname(as.matrix(dt[i,])[1,])
extract_dt_col <- function(dt, i) dt[[i]]


extract_rectangle.data.table <- function(
    x,
    rows = 1:nrow(x),
    cols = 1:ncol(x),
    transpose = FALSE,
    drop = FALSE,
    ...
){
    magrittr::extract(x, rows, cols, with = FALSE) %>%
    as.matrix() %>%
    extract_rectangle.matrix(transpose = transpose, drop = drop)
}


extract_rectangle.matrix <- function(
    x,
    rows      = 1:nrow(x),
    cols      = 1:ncol(x),
    transpose = FALSE,
    drop      = FALSE,
    ...
){
    rectangle <- x[rows, cols, drop = FALSE]
    if (transpose) rectangle %<>% t()
    if (drop) if (nrow(rectangle)==1 | ncol(rectangle)==1) rectangle %<>%
                                                        as.vector('character')
    rectangle
}

#=================================================
# extract_(s|f)data
#=================================================

# Leave rownames(fdata1) empty: fids1 may contain non-valid values
# This happens in MaxQuant files, which sometimes contain missing rows
# (probably after opening in excell)
extract_fdata <- function(
    x, sheet, fids, fvar_rows, fvar_cols, fdata_rows, fdata_cols, transpose
){
    fdata1 <- data.frame(feature_id = fids, stringsAsFactors = FALSE)
    fdata_available <- !is.null(fvar_rows) & !is.null(fvar_cols)
    if (fdata_available){
        fvars1  <-  extract_rectangle(x,
                                    rows      = fvar_rows,
                                    cols      = fvar_cols,
                                    transpose = transpose,
                                    drop      = TRUE,
                                    sheet     = sheet)
        fdata1 %<>% cbind(
                extract_rectangle(x,
                                rows       = fdata_rows,
                                cols       = fdata_cols,
                                transpose  = transpose,
                                drop       = FALSE,
                                sheet      = sheet) %>%
                set_colnames(fvars1) %>%
                data.frame(stringsAsFactors = FALSE, check.names = FALSE))
    }
    fdata1
}


# Leave rownames(sdata1) empty: sids may contain non-valid values
# This happens in SOMA files, where CLIENT_IDENTIFIER is not unique for
# calibrator and buffer samples
extract_sdata <- function(
    x, sheet, sids, svar_rows, svar_cols, sdata_rows, sdata_cols, transpose
){
    sdata1 <- data.frame(sample_id = sids, stringsAsFactors = FALSE)
    sdata_available <- !is.null(svar_rows) & !is.null(svar_cols)
    if (sdata_available){
        svars1 <- extract_rectangle(x,
                                    rows      = svar_rows,
                                    cols      = svar_cols,
                                    transpose =  transpose,
                                    drop      = TRUE,
                                    sheet     = sheet)
        sdata1 %<>% cbind(
                extract_rectangle(  x,
                                    rows       = sdata_rows,
                                    cols       = sdata_cols,
                                    transpose  = !transpose,
                                    drop       = FALSE,
                                    sheet      = sheet) %>%
                set_colnames(svars1) %>%
                data.frame(stringsAsFactors = FALSE, check.names = FALSE))
    }
    sdata1
}


#' Read omics data from rectangular file
#' @param file       string: name of text (txt, csv, tsv, adat) or
#'                           excel (xls, xlsx) file
#' @param sheet      integer/string: only relevant for excel files
#' @param fid_rows   numeric vector: featureid rows
#' @param fid_cols   numeric vector: featureid cols
#' @param sid_rows   numeric vector: sampleid rows
#' @param sid_cols   numeric vector: sampleid cols
#' @param expr_rows  numeric vector: expr rows
#' @param expr_cols  numeric vector: expr cols
#' @param fvar_rows  numeric vector: fvar rows
#' @param fvar_cols  numeric vector: fvar cols
#' @param svar_rows  numeric vector: svar rows
#' @param svar_cols  numeric vector: svar cols
#' @param fdata_rows numeric vector: fdata rows
#' @param fdata_cols numeric vector: fdata cols
#' @param sdata_rows numeric vector: sdata rows
#' @param sdata_cols numeric vector: sdata cols
#' @param transpose  TRUE or FALSE (default)
#' @param verbose    TRUE (default) or FALSE
#' @examples
#' # RNASEQ
#'    file <- download_autonomics_data('stemcells_rna.txt')
#'    read_omics(file,fid_rows   = 2:58736,   fid_cols   = 1,
#'                    sid_rows   = 1,         sid_cols   = 4:14,
#'                    expr_rows  = 2:58736,   expr_cols  = 4:14,
#'                    fvar_rows  = 1,         fvar_cols  = 1:3,
#'                    fdata_rows = 2:58736,   fdata_cols = 1:3,
#'                    transpose  = FALSE)
#' # LCMSMS PROTEINGROUPS
#'    file <- download_autonomics_data('diff_proteinGroups.txt')
#'    read_omics(file,fid_rows   = 2:9783,  fid_cols   = 383,
#'                    sid_rows   = 1,       sid_cols   = seq(124, 316, by = 6),
#'                    expr_rows  = 2:9783,  expr_cols  = seq(124, 316, by = 6),
#'                    fvar_rows  = 1,       fvar_cols  = c(2, 6, 7, 383),
#'                    fdata_rows = 2:9783,  fdata_cols = c(2, 6, 7, 383),
#'                    transpose  = FALSE)
#' # SOMASCAN
#'    file <- download_autonomics_data('stemcells_soma.adat')
#'    read_omics(file,fid_rows   = 21,       fid_cols   = 19:1146,
#'                    sid_rows   = 30:40,    sid_cols   = 4,
#'                    expr_rows  = 30:40,    expr_cols  = 19:1146,
#'                    fvar_rows  = 21:28,    fvar_cols  = 18,
#'                    svar_rows  = 29,       svar_cols  = 1:17,
#'                    fdata_rows = 21:28,    fdata_cols = 19:1146,
#'                    sdata_rows = 30:40,    sdata_cols = 1:17,
#'                    transpose  = TRUE)
#' # METABOLON
#'    file <- download_autonomics_data('glutaminase_metab.xlsx')
#'    read_omics(file, sheet = 2,
#'                    fid_rows   = 11:401,    fid_cols   = 5,
#'                    sid_rows   = 3,         sid_cols   = 15:86,
#'                    expr_rows  = 11:401,    expr_cols  = 15:86,
#'                    fvar_rows  = 10,        fvar_cols  = 1:14,
#'                    svar_rows  = 1:10,      svar_cols  = 14,
#'                    fdata_rows = 11:401,    fdata_cols = 1:14,
#'                    sdata_rows = 1:10,      sdata_cols = 15:86,
#'                    transpose  = FALSE)
#' @export
read_omics <- function(file, sheet = 1, fid_rows, fid_cols, sid_rows, sid_cols,
    expr_rows, expr_cols, fvar_rows  = NULL, fvar_cols = NULL, svar_rows = NULL,
    svar_cols  = NULL, fdata_rows = NULL,  fdata_cols = NULL, sdata_rows = NULL,
    sdata_cols = NULL, transpose  = FALSE, verbose    = TRUE
){
# Assert
    assertive.files::assert_all_are_existing_files(file)
    assertive.types::assert_is_a_bool(transpose)
# Read (in one go if fixed col file)
    is_fixed_col <- is_fixed_col_file(file)
    x  <-   if (is_fixed_col){ extract_rectangle.character(file, sheet=sheet)
            } else {           file }
# Extract components
    fids1  <- extract_rectangle(x, rows = fid_rows, cols = fid_cols,
                        transpose = transpose, drop = TRUE, sheet = sheet)
    sids1  <- extract_rectangle(x, rows = sid_rows, cols = sid_cols,
                        transpose = transpose, drop = TRUE,  sheet = sheet)
    exprs1 <- extract_rectangle(x, rows = expr_rows, cols = expr_cols,
                        transpose = transpose, drop = FALSE, sheet = sheet)
    class(exprs1) <- 'numeric'
    fdata1 <- extract_fdata(x,   sheet = sheet, fids = fids1,
        fvar_rows  = fvar_rows,  fvar_cols  = fvar_cols,
        fdata_rows = fdata_rows, fdata_cols = fdata_cols, transpose = transpose)
    sdata1 <- extract_sdata(x,   sheet = sheet, sids = sids1,
        svar_rows  = svar_rows,  svar_cols  = svar_cols,
        sdata_rows = sdata_rows, sdata_cols = sdata_cols, transpose = transpose)
# Rm features/samples with missing ids (eg MaxQuant: empty interspersed lines)
    idx <- !is.na(fids1)
    if (any(!idx)){
        if (verbose) message(
            "\t\tRm ", sum(!idx), " features with missing 'feature_id' values")
        fids1 <- fids1[idx]; fdata1 <- fdata1[idx,]; exprs1 <- exprs1[idx, ] }
    idx <- !is.na(sids1)
    if (any(!idx)){
        if (verbose) message(
            "\t\tRm ", sum(!idx), " samples with missing 'sample_id' values")
        sids1 <- sids1[idx]; sdata1 <- sdata1[idx, ]; exprs1 <- exprs1[, idx] }
# Name features/samples
    fids1 %<>% uniquify('make.unique'); sids1 %<>% uniquify('make.unique')
    rownames(exprs1) <- rownames(fdata1) <- fids1
    colnames(exprs1) <- rownames(sdata1) <- sids1

# Wrap into Sumexp and return
    object <- SummarizedExperiment(assays = list(exprs = exprs1))
    rowData(object) <- as(fdata1, 'DataFrame')
    colData(object) <- as(sdata1, 'DataFrame')
    metadata(object)$analysis  <-list(
        nfeatures = c(all = nrow(exprs1)), nsamples  = c(all = ncol(exprs1)))
    object
}


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
            stringi::stri_replace_first_fixed(tolower(organism), ' ', '_'),
            stringi::stri_replace_first_fixed(        organism,  ' ', '_'),
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
        tryCatch(expr = { utils::download.file(
                                url = remote, destfile = gtffile, quiet = TRUE)
                            R.utils::gunzip(
                                gtffile,  remove = TRUE, overwrite = TRUE)},
                error = function(cond){
                    message('failed     repeat manually using browser\n');
                    message(cond)})
    }
    gtffile %<>% stringi::stri_replace_last_fixed('.gz', '')
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
#' @importFrom magrittr %>%
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
    dt <- rtracklayer::import(gtffile) %>% GenomicRanges::as.data.frame() %>% data.table::data.table()

    # Filter
    if (!is.null(values)){
        dt %>% data.table::setkeyv(var)
        dt %<>% magrittr::extract(values)
    }

    # Write
    if (!is.null(writefile)){
        message(sprintf("\t\tWrite   %s", writefile))
        dir.create(dirname(writefile), recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(dt, writefile)
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
#' @param sumexpfile   string: file which object is saved to (with saveRDS)
#' @param nthreads     number of cores to be used by Rsubread::featureCounts()
#' @param ...          passed to Rsubread::featureCounts
#' @examples
#' bamdir <- download_autonomics_data("stemcells.bam.zip")
#' read_bam(bamdir, ispaired = TRUE)
#' @export
read_bam <- function(bamdir, ispaired = FALSE, gtffile = NULL,
    fvars = character(0), nthreads   = parallel::detectCores(), ...
){
    # Assert
    assert_all_are_existing_files(bamdir)
    assert_is_a_bool(ispaired)
    if (!is.null(gtffile))   assert_all_are_existing_files(gtffile)
    if (!is.null(sumexpfile))    assert_all_are_existing_files(sumexpfile)
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
    filenames   <- basename(tools::file_path_sans_ext(files))
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

#' Pull columns in a dataframe to the front
#' @param df         data.frame
#' @param first_cols character vector: columns to be pulled to the front
#' @param verbose    TRUE (default) or FALSE
#' @return dataframe with re-ordered columns
#' @examples
#' require(magrittr)
#' df <- data.frame(
#'    symbol = c('A1BG', 'A2M'),
#'    id     = c('1',    '2'),
#'    name   = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'),
#'    type   = c('proteincoding', 'proteincoding'))
#' first_cols <- c('id', 'symbol', 'location', 'uniprot')
#' df %>% autonomics.support::pull_columns(first_cols)
#' @noRd
pull_columns <- function(df, first_cols, verbose = TRUE){

    assertive.types::assert_is_data.frame(df)
    assertive.types::assert_is_character(first_cols)
    extract <- magrittr::extract

    idx <- first_cols %in% names(df)
    if (any(!idx)){
        if (verbose) cmessage(
            'pull_columns: ignore absent columns %s',
            paste0(sprintf("'%s'", first_cols[!idx]), collapse = ', '))
        first_cols %<>% magrittr::extract(idx)
    }

  df %>% extract(, c(first_cols, setdiff(names(df), first_cols)), drop = FALSE)
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

#========================
# AFFYMETRIX MICROARRAYS
#========================

# https://stackoverflow.com/a/4090208
install_if_required <- function(pkgs){
    pkgs %<>% extract(!(pkgs %in% installed.packages()[,"Package"]))
    if(length(pkgs)) BiocManager::install(pkgs)
}


add_affy_fdata <- function(object){

    # Extract entrez identifiers
    entrezgs <- vapply(
      stri_split_fixed(fnames(object), '_'), extract, character(1), 1)

    # Get annotation db
    pkgname <- paste0(metadata(object)$annotation, '.db')
    install_if_required(pkgname)
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
#' @return SummarizedExperiment
#' @examples
#' url <- paste0('http://www.bioconductor.org/help/publications/2003/',
#'                 'Chiaretti/chiaretti2/T33.tgz')
#' localfile <- file.path('~/autonomicscache', basename(url))
#' if (!file.exists(localfile)){
#'     download.file(url, destfile = localfile)
#'     untar(localfile, exdir = path.expand('~/autonomicscache'))
#' }
#' localfile %<>% substr(1, nchar(.)-4)
#' read_affy(celfiles = list.files(localfile, full.names = TRUE))
#' @export
read_affy <- function(celfiles){

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



#==========================================================
# GENEX
#==========================================================

#' Read genex file
#' @param file string: path to exiqon genex file
#' @export
read_genex <- function(file){
    assert_all_are_existing_files(file)
    dt <- extract_rectangle(file, sheet=1)
    read_omics(
        file,
        sheet = 1,
        fid_rows   = 1,                       fid_cols   = 2:(ncol(dt)-2),
        sid_rows   = 2:(nrow(dt)-3),          sid_cols   = 1,
        expr_rows  = 2:(nrow(dt)-3),          expr_cols  = 2:(ncol(dt)-2),
        fvar_rows  = (nrow(dt)-2):nrow(dt),   fvar_cols  = 1,
        svar_rows  = 1,                       svar_cols  = (ncol(dt)-1):ncol(dt),
        fdata_rows = (nrow(dt)-2):nrow(dt),   fdata_cols = 2:(ncol(dt)-2),
        sdata_rows = 2:(nrow(dt)-3),          sdata_cols = (ncol(dt)-1):ncol(dt),
        transpose  = TRUE,
        verbose    = TRUE)
}


#=======================================================
# SOMASCAN
#=======================================================

#' Read somascan
#'
#' Read data from somascan adat file
#'
#' @param file         string: path to *.adat file
#' @param fid_var      string: feature_id   variable
#' @param sid_var      string: sample_id    variable
#' @param subgroup_var string: subgroup     variable
#' @param fname_var    string: feature_name variable
#' @param ...          provide backward compatibility to deprecated function load_soma
#' @return Summarizedexperiment
#' @seealso prepare_somascan
#' @examples
#' file <- download_autonomics_data('hypo_soma.adat')
#' read_somascan(file)
#' @export
read_somascan <- function(file, fid_var = 'SeqId', sid_var = 'SampleId',
    subgroup_var = 'SampleGroup', fname_var    = 'EntrezGeneSymbol'
){
    # Assert
    assert_all_are_existing_files(file)
    assert_is_a_string(fid_var)
    assert_is_a_string(sid_var)

    # Understand file structure
    content <- readLines(file)
    n_row <- length(content)
    n_col <- max(count.fields(file, quote='', sep='\t'))
    f_row <- 1 + which(stri_detect_fixed(content, '^TABLE_BEGIN'))
    s_row <- content %>% extract(f_row:length(.)) %>%
        detect_index(function(y) stri_detect_regex(y, '^\t+', negate=TRUE)) %>%
        add(f_row-1)
    f_col <- content %>% extract(f_row) %>%
        stri_extract_first_regex('^\\t+') %>% stri_count_fixed('\t') %>% add(1)
    fid_cols <- (1+f_col):n_col
    fid_rows <- content %>% extract(f_row:(s_row-1)) %>%
        stri_extract_first_words() %>% equals(fid_var) %>% which() %>%
        add(f_row-1)
    sid_rows <- (1+s_row):n_row
    sid_cols <- content %>% extract(s_row) %>% stri_extract_all_words() %>%
                unlist() %>% equals(sid_var) %>% which()
    # Read
    object <- read_omics(
        file,
        fid_rows   =  fid_rows,         fid_cols   =  fid_cols,
        sid_rows   =  sid_rows,         sid_cols   =  sid_cols,
        expr_rows  = (s_row+1):n_row,   expr_cols  = (f_col+1):n_col,
        fvar_rows  =  f_row:(s_row-1),  fvar_cols  =  f_col,
        fdata_rows =  f_row:(s_row-1),  fdata_cols  = (f_col+1):n_col,
        svar_rows  =  s_row,            svar_cols  = 1:(f_col-1),
        sdata_rows = (s_row+1):n_row,   sdata_cols = 1:(f_col-1),
        transpose  = TRUE,
        verbose    = TRUE)

    # Add sdata/fdata
    sdata(object) %<>% (function(y){ y$subgroup <- y[[subgroup_var]]
    y %>% pull_columns(c('sample_id', 'subgroup'))})

    assert_is_subset(fname_var, fvars(object))
    fdata(object) %<>% (function(y){ y$feature_name <- y[[fname_var]]
    y %>% pull_columns(c('feature_id', 'feature_name'))})

    # Return
    object
}


#======================================================
# METABOLON
#======================================================

find_origscale_sheet <- function(file){
    readxl::excel_sheets(file) %>%
    magrittr::extract(stringi::stri_detect_fixed(., 'OrigScale')) %>%
    magrittr::extract2(1)
}


#' Read metabolon
#' @param file          string: path to metabolon xlsx file
#' @param sheet         number/string: xls sheet number or name
#' @param fid_var       string: feature_id variable (ideally transcends dataset)
#' @param sid_var       string: sample_id variable
#' @param subgroup_var  string: subgroup variable (human comprehensible)
#' @param fname_var     string: feature_name variable
#' @param ...           enable backward compatibility to deprecated load_metabolon
#' @examples
#' file <- download_autonomics_data('hypo_metab.xlsx')
#' read_metabolon(file)
#' @export
read_metabolon <- function(file, sheet = find_origscale_sheet(file),
    fid_var      = '(COMP|COMP_ID)', sid_var = '(CLIENT_IDENTIFIER|Client ID)',
    subgroup_var = 'Group', fname_var    = 'BIOCHEMICAL'
){
    assertive.files::assert_all_are_existing_files(file)

    d_f <- read_excel(file, sheet, col_names = FALSE, .name_repair = 'minimal')

    fvar_rows <- which(!is.na(d_f %>% extract_dt_col(1))) %>% extract(1)
    svar_cols <- which(!is.na(d_f %>% extract_dt_row(1))) %>% extract(1)
    fvar_cols <- fdata_cols <- 1:svar_cols
    svar_rows <- sdata_rows <- 1:fvar_rows
    fvar_names <- d_f %>% extract_dt_row(fvar_rows) %>% extract(1:svar_cols)
    svar_names <- d_f %>% extract_dt_col(svar_cols) %>% extract(1:fvar_rows)

    fid_var <- fvar_names %>% extract(stri_detect_regex(., fid_var))
    sid_var <- svar_names %>% extract(stri_detect_regex(., sid_var))
    fid_rows  <- fdata_rows <- expr_rows <- (fvar_rows+1):nrow(d_f)
    sid_cols  <- sdata_cols <- expr_cols <- (svar_cols+1):ncol(d_f)

    fid_cols  <-  fvar_names %>% equals(fid_var) %>% which()
    sid_rows  <-  svar_names %>% is_in(sid_var) %>% which() %>% extract(1)

    object <- read_omics(file,
                        sheet      = sheet,
                        fid_rows   = fid_rows,      fid_cols   = fid_cols,
                        sid_rows   = sid_rows,      sid_cols   = sid_cols,
                        expr_rows  = expr_rows,     expr_cols  = expr_cols,
                        fvar_rows  = fvar_rows,     fvar_cols  = fvar_cols,
                        svar_rows  = svar_rows,     svar_cols  = svar_cols,
                        fdata_rows = fdata_rows,    fdata_cols = fdata_cols,
                        sdata_rows = svar_rows,     sdata_cols = sdata_cols,
                        transpose  = FALSE,
                        verbose    = TRUE)

    # sdata
    is_subgroup_col <- stringi::stri_detect_regex(svars(object), subgroup_var)
    subgroup_var <- if (any(is_subgroup_col)){  svars(object)[is_subgroup_col]
                    } else {                    sid_var }
    sdata(object) %<>% (function(y){
                            y$subgroup <- y[[subgroup_var]]
                            y %>% pull_columns(c('sample_id', 'subgroup'))})
    # fdata
    assertive.sets::assert_is_subset(fname_var, fvars(object))
    fdata(object) %<>% (function(y){ y$feature_name <- y[[fname_var]]
    y %>% pull_columns(c('feature_id', 'feature_name'))})

    # return
    object
}

#====================================
# LCMS PROTEINGROUPS & PHOSPHOSITES
#====================================

#' maxquant patterns
#' @export
maxquant_patterns <- c(
    `Ratio normalized`             =
      '^Ratio ([HM]/[ML]) normalized (.+)$',
    `Ratio`                        =
      '^Ratio ([HM]/[ML]) (?!count|type|variability|iso-count|normalized)(.+)',
    `LFQ intensity`                =
      '^LFQ intensity ([HML])? ?(.+)$',
    `Reporter intensity corrected` =
      '^Reporter intensity corrected ([0-9]+) (.+)$',
    `Reporter intensity`           =
      '^Reporter intensity ([0-9]+) (.+)$',
    `Intensity labeled`            =
      '^Intensity ([HML]) (.+)$',
    `Intensity`                    =
      '^Intensity (.+)$')


#' Guess maxquant quantity from snames
#'
#' character vector, dataframe, or SummarizedExperiment.
#'
#' @param x character vector, dataframe, or SummarizedExperiment
#' @param ... used for proper S3 method dispatch
#' @return  string: value from names(maxquant_patterns)
#' @examples
#' # file
#'     x <- download_autonomics_data('stemcells_proteinGroups.txt')
#'     guess_maxquant_quantity(x)
#'
#' # character vector
#'     x <- "Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "Ratio M/L STD(L)_EM00(M)_EM01(H)_R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "LFQ intensity EM00.R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "Reporter intensity corrected 0 STD(0)EM00(1)EM01(2)_R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "Reporter intensity 0 STD(0)EM00(1)EM01(2)_R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "Intensity H STD(L)_EM00(M)_EM01(H)_R1"
#'     guess_maxquant_quantity(x)
#'
#' # dataframe
#'     file <- download_autonomics_data( 'stemcells_proteinGroups.txt')
#'     x <- data.table::fread(file)
#'     guess_maxquant_quantity(x)
#'
#' # SummarizedExperiment
#'      file <-download_autonomics_data( 'stemcells_proteinGroups.txt'))
#'      x <- read_proteingroups(
#'              x, standardize_snames = FALSE, demultiplex_snames = FALSE)
#'      guess_maxquant_quantity(x)
#' @export
guess_maxquant_quantity <- function(x, ...){
    UseMethod("guess_maxquant_quantity", x)
}

#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.character <- function(x, ...){

    # read if x is filename
    if (is_existing_file(x)){
        x <- names(fread(x, header = TRUE, nrows = 1))
    }

    # guess from character vector
    for (quantity in names(maxquant_patterns)){
        pattern <- maxquant_patterns[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.data.frame <- function(x, ...){
    x <- names(x)
    for (quantity in names(maxquant_patterns)){
        pattern <- maxquant_patterns[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.SummarizedExperiment <- function(x, ...){
    x <- snames(x)
    for (quantity in names(maxquant_patterns)){
        pattern <- maxquant_patterns[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#' proteingroups fvars
#' @export
proteingroups_fvars <- c(
    'id', 'Majority protein IDs', 'Protein names', 'Gene names',
    'Contaminant', 'Potential contaminant', 'Reverse', 'Phospho (STY) site IDs')


#' Standardize maxquant snames
#'
#' Standardize maxquant sample names
#'
#' Drop "Ratio normalized", "LFQ intensity" etc from maxquant sample names
#'
#' @param x        character(.) or SummarizedExperiment
#' @param quantity maxquant quantity
#' @param verbose  logical(1)
#' @param ...      allow for proper S3 method dispatch
#' @examples
#' # character vector
#'    x <- "Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1"
#'    standardize_maxquant_snames(x)
#'
#'    x <- "Ratio M/L STD(L)_EM00(M)_EM01(H)_R1"
#'    standardize_maxquant_snames(x)
#'
#'    x <-'LFQ intensity STD_R1'
#'    standardize_maxquant_snames(x)
#'
#'    x <- 'LFQ intensity L STD(L)_EM00(M)_EM01(H)_R1'
#'    standardize_maxquant_snames(x)
#'
#'    x <-'Reporter intensity 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1'
#'    standardize_maxquant_snames(x)
#'
#'    x <- 'Reporter intensity corrected 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1'
#'    standardize_maxquant_snames(x)
#'
#' # SummarizedExperiment
#'    file <- download_autonomics_data('stemcells_proteinGroups.txt')
#'    x <- read_proteingroups(
#'            file, standardize_snames = FALSE, demultiplex_snames = FALSE)
#'    standardize_maxquant_snames(x)
#' @export
standardize_maxquant_snames <- function (x, ...) {
    UseMethod("standardize_maxquant_snames", x)
}


#' @export
#' @rdname standardize_maxquant_snames
standardize_maxquant_snames.character <- function(
    x,
    quantity = guess_maxquant_quantity(x),
    verbose  = FALSE,
    ...
){
    # x = mix + channel. Return mix if single channel.
    pattern <- maxquant_patterns %>% magrittr::extract2(quantity)

    # Decompose mix and channel
    if (quantity == 'Intensity'){
        mix     <- stri_replace_first_regex(x, pattern, '$1')
        channel <- rep('', length(mix))
    } else {
        mix     <- stri_replace_first_regex(x, pattern, '$2')
        channel <- stri_replace_first_regex(x, pattern, '$1')
    }

    # Standardize
    if (all(channel=='')){
        cleanx <- mix
    } else {
        cleanx <- sprintf('%s{%s}', mix, channel)
    }
    message('\t\tStandardize snames: ', x[1], '  ->  ', cleanx[1])
    return(cleanx)
}

#' @export
#' @rdname standardize_maxquant_snames
standardize_maxquant_snames.SummarizedExperiment <- function(
    x,
    quantity = guess_maxquant_quantity(x),
    verbose  = FALSE,
    ...
){
    newsnames <- standardize_maxquant_snames(
                    snames(x), quantity = quantity, verbose=verbose)
    snames(x) <- sdata(x)$sample_id <- newsnames
    x
}


#' Demultiplex snames
#'
#' Demultiplex sample names
#'
#' @param x        character vector or SummarizedExperiment
#' @param verbose  logical
#' @param ...      allow for S3 dispatch
#' @examples
#' # character vector
#'
#'    # Alternate multiplexing forms supported
#'    demultiplex_snames("STD(L)_EM00(M)_EM01(H)_R1{M/L}") # Label Ratio
#'    demultiplex_snames('A(0)_B(1)_C(2)_D(3)_R1{0}'     ) # Reporter intensity
#'    demultiplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')   # Label Intensity
#'
#'    # Alternate separators supported
#'    demultiplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')   # underscore
#'    demultiplex_snames('STD(L).EM00(M).EM01(H).R1{L}')   # dot
#'    demultiplex_snames('STD(L)EM00(M)EM01(H).R1{L}')     # no separator
#'
#'    # Composite snames supported
#'    demultiplex_snames("WT.t0(L)_WT.t1(M)_WT.t2(H)_R1{H/L}")
#'
#'    # Uniqueness ensured by appending labels when necessary
#'    demultiplex_snames(c("STD(L).BM00(M).BM00(H).R10{M/L}",
#'                        "STD(L).BM00(M).BM00(H).R10{H/L}"))
#'
#'    # Uniplexed snames are returned unchanged
#'    demultiplex_snames(c('STD_R1', 'EM0_R1'))
#'
#' # SummarizedExperiment
#'    file <- download_autonomics_data('stemcells_proteinGroups.txt')
#'    x <- read_proteingroups(
#'            file, standardize_snames = FALSE, demultiplex_snames = FALSE)
#'    x %<>% standardize_maxquant_snames()
#'    demultiplex_snames(x, verbose = TRUE)
#' @export
demultiplex_snames <- function (x, ...) {
    UseMethod("demultiplex_snames", x)
}

#' @rdname demultiplex_snames
#' @export
demultiplex_snames.character <- function(x, verbose = FALSE, ...){

    # Return unchanged if not multiplexed
    # KD(H)WT(L){H/L}
    pattern <- '(.+)\\{(.+)\\}'
    n_open   <- x %>% stringi::stri_count_fixed('(')
    n_closed <- x %>% stringi::stri_count_fixed(')')
    is_multiplexed <- all(stringi::stri_detect_regex(x, pattern) & (n_open==n_closed) & (n_open>0))
    if (!is_multiplexed) return(x)

    # Separate mix and channel
    mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
    channel <- x %>% stringi::stri_replace_first_regex(pattern, '$2')

    # Separate labels and samples
    pattern <- '\\(.+?\\)'
    labels  <- mix %>% stringi::stri_extract_all_regex(pattern) %>% lapply(stringi::stri_replace_first_fixed, '(', '') %>%
        lapply(stringi::stri_replace_first_fixed, ')', '')
    samples <- mix %>% stringi::stri_split_regex(pattern) %>%
        # rm sep from samples (but not from replicate - needed to glue back later!)
        lapply(function(y){y[1:length(labels[[1]])] %<>% stringi::stri_replace_first_regex('^[_. ]', ''); y})

    # Return unchanged if mixes differ in no of labels or samples
    are_all_identical <- function(y) if (length(y)==1) TRUE else all(y[-1] == y[1])
    n_samples <- vapply(samples,  length, integer(1))
    n_labels  <- vapply(labels, length, integer(1))
    if (!are_all_identical(n_samples) | !are_all_identical(n_labels)){
        autonomics.support::cmessage('\t\tCannot demultiplexing snames: mixes differ in number of samples or labels')
        return(x)
    }

    # Extract replicate
    n_samples %<>% unique()
    n_labels  %<>% unique()
    if (n_samples > n_labels){ replicate <- mix %>% stringi::stri_split_regex(pattern) %>% vapply((function(y) y %>% magrittr::extract(length(y))), character(1))
    samples %<>% lapply(magrittr::extract, 1:(n_samples-1))
    } else {                   replicate <- rep('', length(samples))
    }

    # Extract channel samples from mix
    is_ratio <- channel %>% stringi::stri_detect_fixed('/') %>% all()
    samples %<>% mapply(magrittr::set_names, ., labels, SIMPLIFY = FALSE)
    if (is_ratio){
        num_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(magrittr::extract, character(1), 1)
        den_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(magrittr::extract, character(1), 2)
        den_samples <- mapply(magrittr::extract, samples, den_label)
        num_samples <- mapply(magrittr::extract, samples, num_label)
        xdemultiplex <- sprintf('%s_%s%s', num_samples, den_samples, replicate)
    } else {
        samples %<>% mapply(magrittr::extract, ., channel)
        xdemultiplex <- sprintf('%s%s', samples, replicate)
    }
    if (verbose) autonomics.support::cmessage('\t\tDemultiplex snames: %s  ->  %s', x[1], xdemultiplex[1])

    # Ensure uniqueness. Add labels if required.
    idx <- autonomics.support::cduplicated(xdemultiplex) %>% which()
    if (length(idx)>0){
        label_tags <- channel[idx] %>% stringi::stri_replace_first_fixed('/', '')
        if (verbose)   autonomics.support::cmessage('\t\tUniquify snames: %s -> %s%s (for %d/%d snames)',
                                                    xdemultiplex[idx][1], xdemultiplex[idx][1], label_tags[1],
                                                    length(idx), length(xdemultiplex))
        xdemultiplex[idx] %<>% paste0(label_tags)
    }

    # Return
    return(xdemultiplex)
}

#' @rdname demultiplex_snames
#' @importFrom magrittr %>%
#' @export
demultiplex_snames.SummarizedExperiment <- function(x,verbose  = FALSE, ...){
    newsnames <- snames(x) %>% demultiplex_snames(verbose = verbose)
    snames(x) <- sdata(x)$sample_id <- newsnames
    x
}

#' Read proteingroups
#' @param file                character(1). Path to 'proteinGroups.txt'
#' @param quantity            character(1). Expression columns to extract from proteinGroups file: 'Ratio normalized', 'Ratio', 'LFQ intensity',
#'                                          'LFQ intensity', 'Reporter intensity corrected', 'Reporter intensity', 'Intensity labeled', or 'Intensity'.
#' @param fvars               character(n). Annotation columns to extract from proteinGroups file.
#' @param standardize_snames  logical(1).   Standardize maxquant snames: Ratio normalized H/L WT(L).KD(H).R1 -> WT(L).KD(H).R1{H/L} ?
#' @param demultiplex_snames  logical(1).   Demultiplex (standardized) maxquant snames: WT(L).KD(H).R1{H/L} -> KD(H)_WT(L).R1 ?
#' @param verbose             logical(1).   Message progress?
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% read_proteingroups()
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_proteingroups <- function(
    file,
    quantity           = guess_maxquant_quantity(file),
    fvars              = proteingroups_fvars,
    standardize_snames = TRUE,
    demultiplex_snames = TRUE,
    verbose            = TRUE
){

    # Assert
    assertive.files::assert_all_are_existing_files(file)
    assertive.sets::assert_is_subset(quantity, names(maxquant_patterns))
    assertive.types::assert_is_a_bool(verbose)

    # Initial Read
    assertive.files::assert_all_are_existing_files(file)
    dt <- data.table::fread(file, integer64 = 'numeric', header = TRUE)
    fvars %<>% intersect(names(dt))

    # Define components
    fid_rows   <- 2:nrow(dt)
    fid_cols   <- which(names(dt) == 'id')
    sid_rows   <- 1
    sid_cols   <- names(dt) %>% stringi::stri_detect_regex(maxquant_patterns[[quantity]]) %>% which()
    expr_rows  <- 2:nrow(dt)
    expr_cols  <- sid_cols
    fvar_rows  <- 1
    fvar_cols  <- match(fvars, names(dt))
    fdata_rows <- 2:nrow(dt)
    fdata_cols <- fvar_cols

    # Read sumexp
    object <- file %>% read_omics(fid_rows   = fid_rows,     fid_cols   = fid_cols,
                                  sid_rows   = sid_rows,     sid_cols   = sid_cols,
                                  expr_rows  = expr_rows,    expr_cols  = expr_cols,
                                  fvar_rows  = fvar_rows,    fvar_cols  = fvar_cols,
                                  fdata_rows = fdata_rows,   fdata_cols = fdata_cols,
                                  transpose  = FALSE,
                                  verbose    = verbose)
    # Clean sdata
    if (standardize_snames) object %<>% standardize_maxquant_snames(verbose = verbose)
    if (demultiplex_snames) object %<>% demultiplex_snames(verbose = verbose)
    object$subgroup  <- object$sample_id %>% guess_subgroup_values(verbose = verbose)
    object$replicate <- object$sample_id %>% guess_subgroup_values(invert = TRUE, verbose = FALSE)
    #object$block    <- object$sample_id %>% guess_subject_values( verbose = TRUE)

    # Clean fdata
    contaminant_var <- c('Contaminant', 'Potential contaminant') %>% intersect(fvars(object))
    fdata(object)[[contaminant_var]] %<>% (function(x){x[is.na(x)] <- ''; x})
    fdata(object)[['Reverse'      ]] %<>% (function(x){x[is.na(x)] <- ''; x})
    fdata(object)$feature_name    <- fdata(object)$`Gene names`
    fdata(object)$feature_uniprot <- fdata(object)$`Majority protein IDs`
    fdata(object) %<>% autonomics.support::pull_columns(c('feature_id', 'feature_name', 'feature_uniprot'))

    # Return
    object

}

#' phosphosites fvars
#' @export
phosphosite_fvars <- c('id', 'Protein group IDs', 'Positions within proteins', 'Localization prob')


#' Read phosphosites
#' @param file                string: phosphosites filepath
#' @param proteingroups_file  string: proteingroups filepath
#' @param quantity            NULL or value in names(maxquant_patterns)
#' @param fvars               string vector
#' @param standardize_snames  logical
#' @param demultiplex_snames  logical
#' @param verbose             logical
#' @examples
#' \dontrun{
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    file <- 'extdata/stemdiff/maxquant/phospho (STY)Sites.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% read_phosphosites()
#' }
#' }
#' @importFrom magrittr %>%
#' @export
read_phosphosites <- function(
    file,
    proteingroups_file = dirname(file) %>% paste0('/phospho (STY)Sites.txt'),
    quantity           = guess_maxquant_quantity(file),
    fvars              = phosphosite_fvars,
    standardize_snames = TRUE,
    demultiplex_snames = TRUE,
    verbose            = FALSE
){

    # Check
    `Protein group IDs` <- `Localization prob` <- NULL

    # Assert
    assertive.files::assert_all_are_existing_files(file)
    assertive.sets::assert_is_subset(quantity, names(maxquant_patterns))
    assertive.types::assert_is_character(fvars)

    # Initial Read
    dt <- data.table::fread(file, integer64 = 'numeric', header = TRUE)
    pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
    value_cols <- names(dt) %>% (function(x) stringi::stri_detect_regex(x,pattern) & !stringi::stri_detect_regex(x, '___[1-3]')) %>% which()
    fvar_cols  <- which(names(dt) %in% fvars)
    #dt %<>% magrittr::extract(, c(fvars, value_cols), with = FALSE)

    # Read phosphosites
    fid_rows   <- 2:nrow(dt)
    fid_cols   <- which(names(dt) == 'id')
    sid_rows   <- 1
    sid_cols   <- value_cols
    expr_rows  <- 2:nrow(dt)
    expr_cols  <- value_cols
    fvar_rows  <- 1
    fvar_cols  <- fvar_cols
    fdata_rows <- 2:nrow(dt)
    fdata_cols <- fvar_cols
    phosphosites  <- file  %>% read_omics(fid_rows   = fid_rows,  fid_cols   = fid_cols,
                                          sid_rows   = sid_rows,     sid_cols   = sid_cols,
                                          expr_rows  = expr_rows,    expr_cols  = expr_cols,
                                          fvar_rows  = fvar_rows,    fvar_cols  = fvar_cols,
                                          fdata_rows = fdata_rows,   fdata_cols = fdata_cols,
                                          transpose  = FALSE,
                                          verbose    = verbose)

    # Calculate occupancies (allows to disentangle phosphorylation and protein expression)
    autonomics.support::cmessage('\t\toccupancies(phosphosites) = exprs(phosphosites) - exprs(proteingroups)')
    proteingroups <- read_proteingroups(proteingroups_file, quantity = quantity, verbose = FALSE) %>%
        magrittr::extract(phosphosites %>% fvalues("Protein group IDs"), ) %>%
        magrittr::extract(, phosphosites %>% snames())
    occupancies(phosphosites) <- exprs(phosphosites) - exprs(proteingroups)

    # Simplify snames
    if (standardize_snames) phosphosites %<>% standardize_maxquant_snames(verbose = verbose)
    if (demultiplex_snames) phosphosites %<>% demultiplex_snames(verbose = verbose)

    # Return
    phosphosites

}

