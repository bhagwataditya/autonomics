

#==============================================================================
#
#                       is_excel_file
#                       is_fixed_col_file
#                       nrows
#                       ncols
#
#==============================================================================

is_excel_file <- function(file){
    stri_detect_fixed(
            file_ext(file), 'xls', case_insensitive = TRUE)
}

nrows <- function(x, sheet=1){
    if (is_excel_file(x)){
        nrow(read_excel(
                x, sheet = sheet, .name_repair = 'minimal', col_names = FALSE))
    } else {
        length(readLines(x, warn=FALSE))
    }
}


ncols <- function(x, sheet=1){
    if (is_excel_file(x)){
        ncol(read_excel(
            x, sheet=sheet, .name_repair = 'minimal', col_names = FALSE))
    } else {
        max(count.fields(x, quote = "", sep = '\t'))
    }
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
#' file <- download_data('halama18.metabolon.xlsx')
#' is_fixed_col_file(file)
#' @noRd
is_fixed_col_file <- function(file){
    if (is_excel_file(file)){ TRUE
    } else {
        fields <- count.fields(file, quote = "", sep = '\t')
        all(fields==fields[1])
    }
}

#==============================================================================
#
#                           extract_rectangle
#
#==============================================================================

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
#'    x <- download_data('halama18.metabolon.xlsx')
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
#'    x <-download_data('halama18.metabolon.xlsx') %>%
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
#' @export
extract_rectangle <- function(x, ...){
    UseMethod('extract_rectangle')
}

#' @rdname extract_rectangle
#' @export
extract_rectangle.character <- function(x, rows=seq_len(nrows(x, sheet=sheet)),
    cols = seq_len(ncols(x, sheet=sheet)), verbose = FALSE, transpose = FALSE,
    drop = FALSE, sheet = 1, ...
){
    # Assert
    assert_all_are_existing_files(x)
    assert_is_numeric(rows)
    assert_is_numeric(cols)
    row1 <- min(rows); rown <- max(rows)
    col1 <- min(cols); coln <- max(cols)

    # Read file
    if (verbose) message('\t\tRead ',
                if (sheet==1) '' else paste0("sheet '", sheet, "' of"), file)
    dt <- if (is_excel_file(x)){
            data.table(read_excel(
                x,
                sheet       = sheet,
                col_names   = FALSE,
                range       = sprintf('R%dC%d:R%dC%d', row1, col1, rown, coln),
                .name_repair = 'minimal'))
        } else {
            fread(
                x,
                na.strings = "",
                header     = FALSE,
                integer64  = 'numeric',
                skip       = row1-1,
                nrows      = 1+rown-row1,
                select     = col1:coln) }

    # Extract rectangle
    extract_rectangle.data.table(dt, transpose = transpose, drop = drop)

}

# Extract row
extract_dt_row <- function(dt, i) unname(as.matrix(dt[i,])[1,])
extract_dt_col <- function(dt, i) dt[[i]]


#' @rdname extract_rectangle
#' @export
extract_rectangle.data.table <- function(
    x,
    rows = seq_len(nrow(x)),
    cols = seq_len(ncol(x)),
    transpose = FALSE,
    drop = FALSE,
    ...
){
    extract(x, rows, cols, with = FALSE) %>%
    as.matrix() %>%
    extract_rectangle.matrix(transpose = transpose, drop = drop)
}


#' @rdname extract_rectangle
#' @export
extract_rectangle.matrix <- function(
    x,
    rows      = seq_len(nrow(x)),
    cols      = seq_len(ncol(x)),
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

#==============================================================================
#
#                         extract_fdata
#                         extract_sdata
#
#==============================================================================

# Leave rownames(fdata1) empty: fids1 may contain non-valid values
# This happens in MaxQuant files, which sometimes contain missing rows
# (probably after opening in excell)
extract_fdata <- function(
    x, sheet, fids, fvar_rows, fvar_cols, fdata_rows, fdata_cols, transpose
){
    fdata1 <- data.frame(feature_id = fids, stringsAsFactors = FALSE)#,
                        #row.names  = fids) doesn't work with NA fids
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
        sdata1 %<>% numerify()
    }
    sdata1
}


is_numeric <- function(x)   all(grepl('^[0-9.]+$', x) | is.na(x) | x=='NA')
numerify   <- function(df){
    for (i in names(df)){
        if (is_numeric(df[[i]]))  df[[i]] %<>% as.numeric()
    }
    df
}


#==============================================================================
#
#                         read_omics
#
#==============================================================================


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
#' @return SummarizedExperiment
#' @examples
#' # RNASEQ
#'    file <- download_data('billing16.rnacounts.txt')
#'    read_omics(file,fid_rows   = 2:58736,   fid_cols   = 1,
#'                    sid_rows   = 1,         sid_cols   = 4:14,
#'                    expr_rows  = 2:58736,   expr_cols  = 4:14,
#'                    fvar_rows  = 1,         fvar_cols  = 1:3,
#'                    fdata_rows = 2:58736,   fdata_cols = 1:3,
#'                    transpose  = FALSE)
#' # LCMSMS PROTEINGROUPS
#'    file <- download_data('billing19.proteingroups.txt')
#'    read_omics(file,fid_rows   = 2:9783,  fid_cols   = 383,
#'                    sid_rows   = 1,       sid_cols   = seq(124, 316, by = 6),
#'                    expr_rows  = 2:9783,  expr_cols  = seq(124, 316, by = 6),
#'                    fvar_rows  = 1,       fvar_cols  = c(2, 6, 7, 383),
#'                    fdata_rows = 2:9783,  fdata_cols = c(2, 6, 7, 383),
#'                    transpose  = FALSE)
#' # SOMASCAN
#'    file <- download_data('billing16.somascan.adat')
#'    read_omics(file,fid_rows   = 21,       fid_cols   = 19:1146,
#'                    sid_rows   = 30:40,    sid_cols   = 4,
#'                    expr_rows  = 30:40,    expr_cols  = 19:1146,
#'                    fvar_rows  = 21:28,    fvar_cols  = 18,
#'                    svar_rows  = 29,       svar_cols  = 1:17,
#'                    fdata_rows = 21:28,    fdata_cols = 19:1146,
#'                    sdata_rows = 30:40,    sdata_cols = 1:17,
#'                    transpose  = TRUE)
#' # METABOLON
#'    file <- download_data('halama18.metabolon.xlsx')
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
    sdata_cols = NULL, transpose  = FALSE, verbose    = TRUE){
# Assert
    assert_all_are_existing_files(file);     assert_is_a_bool(transpose)
# Read (in one go if fixed col file)
    if (verbose) message("\tRead: ", file)
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
    suppressWarnings(class(exprs1) <- 'numeric') # prevent NA warning
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
    metadata(object)$file <- file
    object
}


#=============================================================================
#
#                        add_designvars
#
#=============================================================================


split_values <- function(x){
    sep <- guess_sep(x)
    dt <- data.table::data.table(x = x)
    dt[, data.table::tstrsplit(x, sep) ]
}

#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' split_subgroup_values(object)
#' @noRd
split_subgroup_values <- function(object){
    subgroupvalues <- subgroup_values(object)
    cbind(subgroup = subgroupvalues, split_values(subgroupvalues))
}



#==============================================================================
#
#                           read_affymetrix
#
#==============================================================================

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
#' if (!requireNamespace("BiocManager", quietly = TRUE)){
#'     install.packages('BiocManager')
#' }
#' BiocManager::install('hgu95av2.db')
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


#==============================================================================
#
#                           read_genex
#
#==============================================================================


#' Read genex file
#' @param file string: path to exiqon genex file
#' @return SummarizedExperiment
#' @noRd
read_genex <- function(file){
    assert_all_are_existing_files(file)
    dt <- extract_rectangle(file, sheet=1)
    read_omics(
        file,
        sheet = 1,
        fid_rows   = 1,                     fid_cols   = 2:(ncol(dt)-2),
        sid_rows   = 2:(nrow(dt)-3),        sid_cols   = 1,
        expr_rows  = 2:(nrow(dt)-3),        expr_cols  = 2:(ncol(dt)-2),
        fvar_rows  = (nrow(dt)-2):nrow(dt), fvar_cols  = 1,
        svar_rows  = 1,                     svar_cols  = (ncol(dt)-1):ncol(dt),
        fdata_rows = (nrow(dt)-2):nrow(dt), fdata_cols = 2:(ncol(dt)-2),
        sdata_rows = 2:(nrow(dt)-3),        sdata_cols = (ncol(dt)-1):ncol(dt),
        transpose  = TRUE,
        verbose    = TRUE)
}


