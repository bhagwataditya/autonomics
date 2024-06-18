#---------------------------------------------------------------------------
# 
#                       guess_maxquant_quantity
#
#---------------------------------------------------------------------------


#' maxquant quantity patterns
#' @examples
#' MAXQUANT_PATTERNS
#' @export
MAXQUANT_PATTERNS <- c(
    `normalizedratio`            = '^Ratio ([HM]/[ML]) normalized (.+)$',
    `ratio`                      = '^Ratio ([HM]/[ML]) (?!count|type|variability|iso-count|normalized)(.+)',
    `correctedreporterintensity` = '^Reporter intensity corrected ([0-9]+) (.+)$',
    `reporterintensity`          = '^Reporter intensity ([0-9]+) (.+)$',
    `maxlfq`                     = '^LFQ intensity ([HML] )? ?(.+)$',
    `labeledintensity`           = '^Intensity ([HML]) (.+)$',
    `intensity`                  = '^Intensity (.+)$'
)


#' Guess maxquant quantity from snames
#'
#' @param x character vector
#' @return  string: value from names(MAXQUANT_PATTERNS)
#' @examples
#' # file
#'     file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'     guess_maxquant_quantity(file)
#'
#' # character vector
#'     x <- "Ratio M/L normalized STD(L)_E00(M)_E01(H)_R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "Ratio M/L STD(L)_E00(M)_E01(H)_R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "LFQ intensity E00.R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "Reporter intensity corrected 0 STD(0)E00(1)E01(2)_R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "Reporter intensity 0 STD(0)E00(1)E01(2)_R1"
#'     guess_maxquant_quantity(x)
#'
#'     x <- "Intensity H STD(L)_E00(M)_E01(H)_R1"
#'     guess_maxquant_quantity(x)
#'
#' @export
guess_maxquant_quantity <- function(x){
# Assert
    assert_is_character(x)
# read if x is filename
    if (is_scalar(x)){
        if (is_existing_file(x)){  
            x <- names(fread(x, header = TRUE, nrows = 1))
        }
    }
# guess from character vector
    for (quantity in names(MAXQUANT_PATTERNS)){
        pattern <- MAXQUANT_PATTERNS[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#---------------------------------------------------------------------------
# 
#                       .read_maxquant_proteingroups
#                       .read_maxquant_phosphosites
#
#---------------------------------------------------------------------------

utils::globalVariables('where')
un_int64 <- function(x) {
    dplyr::mutate(x, dplyr::across(where(bit64::is.integer64), as.numeric))
}

#' Read proteingroups/phosphosites as-is
#' @param file         proteingroups / phosphosites file
#' @param profile  proteingroups file
#' @param quantity     string
#' @param verbose      TRUE / FALSE
#' @return data.table
#' @examples 
#' profile <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
#' fosfile <- system.file('extdata/billing19.phosphosites.txt',  package = 'autonomics')
#' prodt <- .read_maxquant_proteingroups(file = profile)
#' fosdt <- .read_maxquant_phosphosites( file = fosfile, profile = profile)
#' @export
.read_maxquant_proteingroups <- function(file, quantity = guess_maxquant_quantity(file), verbose = TRUE){
# Assert
    assert_maxquant_proteingroups(file)
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS))
# Read
    if (verbose)  cmessage('%sRead%sproteingroups   %s', spaces(4), spaces(6), file)
    prodt <- fread(file, colClasses = c(id = 'character'), integer64 = 'numeric')
    prodt %<>% un_int64()
  # prodt[Reverse == '+', `Majority protein IDs` := split_extract_fixed(`Majority protein IDs`, ';', 1)]
            # In FOS `Proteins` column is empty for Reverse, so `Protein` (singular!) is copied over
            # In PRO therefore take leading protein only for naming identity with FOS.
            # Upon further insight naming parity is probably not so important
            # So dont unnecessarily modify
# Extract relevant columns
    pattern <- MAXQUANT_PATTERNS[[quantity]]   # the proteins column is empty in (only) reverse
    anncols <- c('id', 'Majority protein IDs', 'Reverse', 
                 'Potential contaminant', 'Contaminant', 'Fasta headers')#, 'Phospho (STY) site IDs')
    anncols %<>% intersect(names(prodt))
    valcols <- grep(pattern, names(prodt), value = TRUE)
    pepcols <- grep('Razor + unique peptides ', names(prodt), fixed = TRUE, value = TRUE)
    prodt %<>% extract(, c(anncols, pepcols, valcols), with = FALSE)
    digits <- ceiling(log10(nrow(prodt)))
  # if (verbose)  message('\t\t', nrow(prodt), ' proteins, contaminants, reverse')
# Return
    names(prodt) %<>% stri_replace_first_fixed(          'Reverse',     'reverse')
    names(prodt) %<>% stri_replace_first_fixed(          'Contaminant', 'contaminant') # older MaxQuant
    names(prodt) %<>% stri_replace_first_fixed('Potential contaminant', 'contaminant') # newer MaxQuant
    setnames(prodt, 'id',                   'proId')
    setnames(prodt, 'Majority protein IDs', 'uniprot')
    prodt[]
}


#' @rdname dot-read_maxquant_proteingroups
#' @export
.read_maxquant_phosphosites <- function(
    file, profile, quantity = guess_maxquant_quantity(file), verbose = TRUE
){
# Assert
    assert_maxquant_proteingroups(profile)
    assert_maxquant_phosphosites(file)
    `Protein group IDs` <- Reverse <- Proteins <- Protein <- NULL
# Read    
    if (verbose)  cmessage('%sphosphosites    %s', spaces(14), file)
    colclasses <- c(id = 'character', `Protein group IDs` = 'character')
    fosdt <- fread(file, colClasses = colclasses, integer64 = 'numeric')
    fosdt %<>% un_int64()
# Extract relevant columns
    fosdt[Reverse == '+', Proteins := Protein]  # the proteins column is empty in (only) reverse
    pattern <- MAXQUANT_PATTERNS[[quantity]]    # this fails the feature_id naming
    anncols <- c('id', 'Protein group IDs', 'Proteins', 
                 'Positions within proteins', 'Amino acid', 
                 'Reverse', 'Contaminant', 'Potential contaminant', 'Fasta headers')
    anncols %<>% intersect(names(fosdt))
    valcols <- grep(pattern, names(fosdt), value = TRUE)
    valcols %<>% extract(stri_detect_fixed(., '___1'))
    fosdt %<>% extract(, c(anncols, valcols), with = FALSE)
    digits <- ceiling(log10(nrow(fosdt)))
    if (verbose)  message(spaces(38), 'Read   ', formatC(nrow(fosdt), digits = digits), 
                          ' phosphosites in proteins, contaminants, reverse')
    names(fosdt) %<>% stri_replace_first_fixed(          'Reverse',     'reverse')
    names(fosdt) %<>% stri_replace_first_fixed(          'Contaminant', 'contaminant') # older MaxQuant
    names(fosdt) %<>% stri_replace_first_fixed('Potential contaminant', 'contaminant') # newer MaxQuant
# Filter
    fosdt <- fosdt[stri_count_fixed(`Protein group IDs`, ';') == 0]
    if (verbose)  message(spaces(38), 'Retain ', formatC(nrow(fosdt), digits = digits), 
                          ' phosphosites in single proteingroup')
    idx <- rowSums(fosdt[, valcols, with = FALSE], na.rm = TRUE) > 0
    fosdt <- fosdt[which(idx)]
    if (verbose)  message(spaces(38), 'Retain ', formatC(nrow(fosdt), digits = digits), 
                          ' phosphosites with signal in some sample')
# Rename 
    names(fosdt) %<>% stri_replace_first_fixed('___1', '')
    setnames(fosdt, 'id',                'fosId')
    setnames(fosdt, 'Protein group IDs', 'proId')
    setnames(fosdt, 'Proteins',          'uniprot')
    fosdt[]
}


#---------------------------------------------------------------------------
# 
#                      drop_differing_uniprots
#
#---------------------------------------------------------------------------


drop_differing_uniprots <- function(fosdt, prodt, verbose){
# Assert
    contaminant <- fosId <- `Positions within proteins` <- proId <- NULL
    reverse <- uniprot <- NULL
# Avoid changing the global env
    prodt %<>% copy()
    fosdt %<>% copy()
# Extract annotation cols
    if (verbose)  message(spaces(38), 'Keep proteingroup uniprots')
    proanncols <- c('proId', 'uniprot', 'reverse', 'contaminant')
    fosanncols <- c('fosId', 'proId', 'uniprot', 'Positions within proteins', 'reverse', 'contaminant')
    proanndt <- prodt[, proanncols, with = FALSE]
    fosanndt <- fosdt[, fosanncols, with = FALSE]
# Separate contaminants-reverse from other proteins.
    conrevdt <- fosanndt[reverse == '+' | contaminant == '+']
    fosanndt <- fosanndt[reverse == ''  & contaminant == '']
    proanndt <- proanndt[reverse == ''  & contaminant == '']
    fosanndt[, c('reverse', 'contaminant') := NULL]
    proanndt[, c('reverse', 'contaminant') := NULL]
    conrevdt[, c('reverse', 'contaminant') := NULL]
# Merge
    fosanndt %<>% uncollapse(uniprot, `Positions within proteins`, sep = ';')
    proanndt %<>% uncollapse(uniprot,                              sep = ';')
    fosanndt %<>% merge(proanndt, by = c('proId', 'uniprot'))
    fosanndt %<>% extract(order(as.integer(fosId)))
    fosanndt %<>% extract(, lapply(.SD, paste0, collapse = ';'), by = 'fosId')
# Add back contaminants/reverse
    fosanndt %<>% rbind(conrevdt)
    fosanndt[, proId := NULL]
    fosdt[, c('uniprot', 'Positions within proteins') := NULL]
    fosdt %<>% merge(fosanndt, by = 'fosId', all.x = TRUE, sort = FALSE)
    fosdt %<>% pull_columns(names(fosanndt))
    fosdt    
}


#---------------------------------------------------------------------------
# 
#                 mqdt_to_mat
#                 dequantify
#                     label2index
#
#---------------------------------------------------------------------------


#' Convert maxquant data.table to matrix
#' @param dt       data.table
#' @param pattern  string
#' @param logbase  base for value logaritmization
#' @param verbose  TRUE / FALSE
#' @return matrix
#' @examples 
#' profile <- system.file('extdata/billing19.proteingroups.txt',  package = 'autonomics')
#' fastafile <- system.file('extdata/uniprot_hsa_20140515.fasta', package = 'autonomics')
#' prodt <- .read_maxquant_proteingroups(file = profile)
#'     uniprothdrs <- read_uniprotdt(fastafile)
#' contaminanthdrs <- read_contaminantdt()
#'    maxquanthdrs <- parse_maxquant_hdrs(prodt$`Fasta headers`)
#' prodt %<>% annotate_maxquant(uniprothdrs, contaminanthdrs, maxquanthdrs)
#' quantity <- guess_maxquant_quantity(profile)
#' pattern <- MAXQUANT_PATTERNS[[quantity]]
#' mqdt_to_mat(prodt, pattern = pattern)[1:2, 1:2]
#' @noRd
mqdt_to_mat <- function(dt, pattern, logbase = 2, verbose = TRUE){
    mat <- dt[, .SD, .SDcols = patterns(pattern)]
    mat %<>% data.matrix()
    rownames(mat) <- dt$feature_id
    mat %<>% zero_to_na(verbose = verbose) 
    mat %<>% nan_to_na( verbose = verbose)
    if (!is.null(logbase)) mat <- log(mat, base = logbase)
    mat
}


#' Dequantify maxquant snames
#' 
#' Drop quantity ('Reporter intensity'). \cr
#' Encode `{channel}` as suffix.
#' 
#' `               Ratio H/L WT(L).KD(H).R1  ->  WT(L).KD(H).R1{H/L}`
#' `                    LFQ intensity WT.R1  ->  WT.R1`
#' `Reporter intensity 0 WT(126).KD(127).R1  ->  WT(1).KD(2).R1{1}`
#' @param x        `character`
#' @param quantity `'ratio',              'normalizedratio'`,  \cr
#'                 `'LFQ intensity'`, \cr
#'                 `'intensity',          'labeledintensity'`
#'                 `'reporterintensity', 'correctedreporterintensity'`
#' @param verbose  `TRUE` or `FALSE`
#' @return `character`
#' @examples
#' dequantify(c('Ratio H/L WT(L).KD(M).OE(H).R1',             # Ratios
#'              'Ratio M/L WT(L).KD(M).OE(H).R1'))
#' dequantify(c('Ratio H/L normalized WT(L).KD(M).OE(H).R1',  # Norm. Ratios
#'              'Ratio M/L normalized WT(L).KD(M).OE(H).R1'))
#' dequantify(c('LFQ intensity WT.R1',                        # LFQ intensity
#'              'LFQ intensity KD.R1'))
#' dequantify(c('Reporter intensity 1 WT(126).KD(127).R1',    # Rep.intensities
#'              'Reporter intensity 2 WT(126).KD(127).R1'))
#' @md
#' @export
dequantify <- function(
    x, quantity = guess_maxquant_quantity(x), verbose  = FALSE
){
# x = multiplex + channel. Return multiplex if single channel.
# Decompose multiplex and channel
    pattern <- MAXQUANT_PATTERNS[[quantity]]
    if (quantity == 'intensity'){
        multiplex <- stri_replace_first_regex(x, pattern, '$1')
        channel   <- rep('', length(multiplex))
    } else {
        multiplex <- stri_replace_first_regex(x, pattern, '$2')
        channel   <- stri_replace_first_regex(x, pattern, '$1')
        channel %<>% stri_replace_first_fixed(' ', '')
    }
# Reporter intensity
    if (quantity %in% c('reporterintensity', 'correctedreporterintensity')){
        multiplex %<>% lapply(label2index)
        channel %<>% as.numeric()
        if (0 %in% channel){                           # pre-2018 mq is 0-based
            channel %<>% as.numeric() %>% add(1) %>% as.character()
        }
    }
# Standardize
    if (all(channel=='')){
        cleanx <- multiplex
    } else {
        cleanx <- sprintf('%s{%s}', multiplex, channel)
    }
    if (verbose)  cmessage('%sStandardize snames: %s  ->  %s', spaces(14), x[1], cleanx[1])
    return(cleanx)
}


#' Convert labels into indices
#' @param x  `character`
#' @examples
#' label2index(x = 'Reporter intensity 0 WT(0).KD(1).OE(2).R1')
#' label2index(x = 'Reporter intensity 1 WT(1).KD(2).OE(3).R1')
#' label2index(x = 'Reporter intensity 0 WT(126).KD(127).OE(128).R1')
#' label2index(x = 'Reporter intensity 1 WT(126).KD(127).OE(128).R1')
#' label2index(x = 'Reporter intensity 1 Mix1')
#' @export
label2index <- function(x){
    labels <- unlist(stri_extract_all_regex(x, '\\(.+?\\)'))
    if (any(is.na(labels)))  return(x)
    for (i in rev(seq_along(labels))){ 
        # important to do this in rev order!!! otherwise ...
        #                 Reporter intensity 0 WT(0).KD(1).OE(2).R1
        # i=1: (0) -> 1 : Reporter intensity 0 WT(1).KD(1).OE(2).R1
        # i=2: (1) -> 2 : Reporter intensity 0 WT(2).KD(1).OE(2).R1
        # not what we want!!!
        x %<>% stri_replace_first_fixed(labels[[i]], paste0('(',i,')'))
    }
    x
}


#---------------------------------------------------------------------------
# 
#                   read_maxquant_proteingroups
#                   read_maxquant_phosphosites
#
#---------------------------------------------------------------------------



#' Is a file?
#'
#' Is a file (and not a dir)
#'
#' This function distinguishies between dir and file.
#' Others dont: is.file, fs::file_exists, assertive::is_existing_file
#' @param file filepath
#' @examples
#' dir  <- tempdir();  dir.create(dir, showWarnings = FALSE)
#' file <- tempfile(); invisible(file.create(file))
#' is_file(dir)
#' is_file(file)
#' @export
is_file <- function(file){
    file.exists(file) & !dir.exists(file)
}


#' Read maxquant proteingroups
#' @param dir           proteingroups directory
#' @param file          proteingroups file
#' @param fastafile     uniprot fastafile
#' @param restapi       TRUE or FALSE : use uniprot restapi to annotate uniprots not in fastadt ?
#' @param quantity     'normalizedratio', 'ratio', 'correctedreporterintensity', 
#'                     'reporterintensity', 'maxlfq', 'labeledintensity', 
#'                     'intensity' or NULL
#' @param subgroups     NULL or string vector : subgroups to retain
#' @param contaminants  TRUE or FALSE : retain contaminants ?
#' @param reverse       TRUE or FALSE : include reverse hits ?
#' @param rm_missing_in_all_samples TRUE or FALSE
#' @param invert        string vector : subgroups which require inversion
#' @param impute        TRUE or FALSE: impute group-specific NA values?
#' @param plot          TRUE or FALSE: plot ?
#' @param label         fvar
#' @param pca           TRUE or FALSE: run pca ?
#' @param pls           TRUE or FALSE: run pls ?
#' @param fit           model engine: 'limma', 'lm', 'lme(r)', 'wilcoxon' or NULL
#' @param formula       model formula
#' @param block         model blockvar: string or NULL
#' @param coefs         model coefficients    of interest: character vector or NULL
#' @param contrasts     coefficient contrasts of interest: character vector or NULL
#' @param palette       color palette : named character vector
#' @param verbose       TRUE or FALSE : message ?
#' @param ...           maintain deprecated functions
#' @return SummarizedExperiment
#' @examples
#' # fukuda20 - LFQ
#'     file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'     pro <- read_maxquant_proteingroups(file = file)
#'     
#' # billing19 - Normalized Ratios
#'     file <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
#'     fastafile <- system.file('extdata/uniprot_hsa_20140515.fasta', package = 'autonomics')
#'     subgroups <- sprintf('%s_STD', c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'))
#'     pro <- read_maxquant_proteingroups(file = file, subgroups = subgroups)
#'     pro <- read_maxquant_proteingroups(file = file, fastafile = fastafile, subgroups = subgroups)
#' @export
read_maxquant_proteingroups <- function(
                          dir = getwd(), 
                         file = if (is_file(dir)) dir else file.path(dir, 'proteinGroups.txt'), 
                    fastafile = NULL,
                      restapi = FALSE,
                     quantity = NULL,
                    subgroups = NULL,
                       invert = character(0),
                 contaminants = FALSE,
                      reverse = FALSE,
    rm_missing_in_all_samples = TRUE, 
                       impute = FALSE,
                         plot = FALSE,
                        label = 'feature_id',
                          pca = plot,
                          pls = plot, 
                          fit = if (plot) 'limma' else NULL,
                      formula = ~ subgroup, 
                        block = NULL, 
                        coefs = NULL, 
                    contrasts = NULL,
                      palette = NULL, 
                      verbose = TRUE
){
# Assert
    assert_maxquant_proteingroups(file)
    if (is.null(quantity))  quantity <- guess_maxquant_quantity(file)
    assert_is_a_string(quantity)
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS))
    assert_is_a_bool(verbose)
    `Fasta headers` <- NULL
# Read/Annotate
    prodt <- .read_maxquant_proteingroups(file = file, quantity = quantity, verbose = verbose)
        uniprothdrs <- read_uniprotdt(fastafile)
    contaminanthdrs <- read_contaminantdt()
       maxquanthdrs <- parse_maxquant_hdrs(prodt$`Fasta headers`); prodt[, `Fasta headers` := NULL ]
    prodt %<>% annotate_maxquant( uniprothdrs = uniprothdrs, 
                              contaminanthdrs = contaminanthdrs, 
                                 maxquanthdrs = maxquanthdrs, 
                                      restapi = restapi,
                                      verbose = verbose )
# SumExp
    if (verbose)  cmessage('%sSumExp', spaces(4))
    pattern <- MAXQUANT_PATTERNS[[quantity]]
    promat <- mqdt_to_mat(prodt, pattern, verbose = verbose)
    pepcols <- names(prodt) %>% extract(stri_detect_fixed(., 'eptides'))
    pepdt <- prodt[, pepcols, with = FALSE]
    prodt %<>% extract(, names(prodt) %>% setdiff(colnames(promat)) %>% setdiff(names(pepdt)), with = FALSE)
    object <- list(promat)
    names(object) <- paste0('log2', quantity)
    object %<>% SummarizedExperiment(rowData = prodt)
# Dequantify. Add pepcounts
    object$mqcol <- colnames(object)
    colnames(object) %<>% dequantify(quantity = quantity, verbose = verbose)
    pepcols <- paste0('Razor + unique peptides ', gsub('\\{.+\\}', '', colnames(object)))
    pepmat <- pepdt[, pepcols, with = FALSE ]
    pepmat %<>% data.matrix()
    dimnames(pepmat) <- dimnames(object)
    assays(object)$pepcounts <- pepmat
    object %<>% process_maxquant( subgroups = subgroups,      
                                     invert = invert,
                                    reverse = reverse,  
                               contaminants = contaminants,   
                  rm_missing_in_all_samples = rm_missing_in_all_samples, 
                                     impute = impute,
                                    verbose = verbose )
    assays <- c(assayNames(object)[1], 'pepcounts')
    for (assay in assays)  object %<>% add_assay_means(assay)
    object %<>% analyze(
        pca          = pca,           pls       = pls,
        fit          = fit,           formula   = formula,
        block        = block,         coefs     = coefs,
        contrasts    = contrasts,     plot      = plot,
        label        = label,
        palette      = palette,       verbose   = verbose )
    object
}


#' @rdname read_maxquant_proteingroups
#' @export
read_proteingroups <- function(...){
    .Deprecated('read_maxquant_proteingroups')
    read_maxquant_proteingroups(...)
}


#' Read maxquant phosphosites
#' @param dir           proteingroups directory
#' @param fosfile       phosphosites  file
#' @param profile       proteingroups file
#' @param fastafile     uniprot fastafile
#' @param restapi       TRUE or FALSE : annotate non-fastadt uniprots using uniprot restapi
#' @param quantity     'normalizedratio', 'ratio', 'correctedreporterintensity', 
#'                     'reporterintensity', 'maxlfq', 'labeledintensity', 
#'                     'intensity' or NULL
#' @param subgroups     NULL or string vector : subgroups to retain
#' @param contaminants  TRUE or FALSE: retain contaminants ?
#' @param reverse       TRUE or FALSE: include reverse hits 
#' @param rm_missing_in_all_samples TRUE or FALSE
#' @param localization  number: min localization probability (for phosphosites)
#' @param invert        string vector: subgroups which require inversion
#' @param impute        TRUE or FALSE: impute group-specific NA values?
#' @param plot          TRUE or FALSE
#' @param label         fvar
#' @param pca           TRUE or FALSE: run pca ?
#' @param pls           TRUE or FALSE: run pls ?
#' @param fit           model engine: 'limma', 'lm', 'lme(r)', 'wilcoxon' or NULL
#' @param formula       model formula
#' @param block         model blockvar: string or NULL
#' @param coefs         model coefficients          of interest: string vector or NULL
#' @param contrasts     model coefficient contrasts of interest: string vector or NULL
#' @param palette       color palette: named string vector
#' @param verbose       TRUE or FALSE: message ?
#' @param ...           maintain deprecated functions
#' @return SummarizedExperiment
#' @examples
#'   profile <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
#'   fosfile <- system.file('extdata/billing19.phosphosites.txt',  package = 'autonomics')
#' fastafile <- system.file('extdata/uniprot_hsa_20140515.fasta', package = 'autonomics')
#' subgroups <- sprintf('%s_STD', c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'))
#' pro <- read_maxquant_proteingroups(file = profile, subgroups = subgroups)
#' fos <- read_maxquant_phosphosites(
#'               fosfile = fosfile, profile = profile, subgroups = subgroups)
#' fos <- read_maxquant_phosphosites(
#'               fosfile = fosfile, profile = profile, fastafile = fastafile, subgroups = subgroups)
#' @export
read_maxquant_phosphosites <- function(
                      dir = getwd(), 
                  fosfile = if (is_file(dir)) dir else file.path(dir, 'phospho (STY)Sites.txt'), 
                  profile = file.path(dirname(fosfile), 'proteinGroups.txt'), 
                fastafile = NULL,
                  restapi = FALSE,
                 quantity = NULL,   
                subgroups = NULL,
                   invert = character(0),
             contaminants = FALSE,
                  reverse = FALSE,
rm_missing_in_all_samples = TRUE,
             localization = 0.75,
                   impute = FALSE,
                     plot = FALSE,
                    label = 'feature_id',
                      pca = plot, 
                      pls = plot,
                      fit = if (plot) 'limma' else NULL,
                  formula = ~ subgroup,
                    block = NULL,
                    coefs = NULL,
                contrasts = NULL, 
                  palette = NULL,
                  verbose = TRUE
){
# Assert
    assert_all_are_existing_files(c(fosfile, profile))
    if (is.null(quantity))  quantity <- guess_maxquant_quantity(fosfile)
    assert_is_a_string(quantity)
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS))
    assert_is_a_bool(verbose)
    `Fasta headers` <- NULL
# Read
    prodt <- .read_maxquant_proteingroups(file = profile, quantity = quantity, verbose = verbose)
    fosdt <- .read_maxquant_phosphosites( file = fosfile, quantity = quantity, profile = profile, verbose = verbose)
    fosdt %<>% drop_differing_uniprots(prodt, verbose = verbose)
        uniprothdrs <- read_uniprotdt(fastafile)
    contaminanthdrs <- read_contaminantdt()
       maxquanthdrs <- parse_maxquant_hdrs(prodt$`Fasta headers`); fosdt[, `Fasta headers` := NULL ]
    fosdt %<>% annotate_maxquant( uniprothdrs = uniprothdrs, 
                              contaminanthdrs = contaminanthdrs, 
                                 maxquanthdrs = maxquanthdrs, 
                                      restapi = restapi,
                                      verbose = verbose )
    prodt %<>% extract(fosdt$proId, on = 'proId')
# SumExp
    if (verbose)  cmessage('%sSumExp', spaces(4))
    pattern <- MAXQUANT_PATTERNS[[quantity]]
    promat <- mqdt_to_mat(prodt, pattern, verbose = verbose)
    pepcols <- names(prodt) %>% extract(stri_detect_fixed(., 'eptides'))
    pepdt <- prodt[, pepcols, with = FALSE]
    prodt %<>% extract(, names(prodt) %>% setdiff(colnames(promat)) %>% setdiff(names(pepdt)), with = FALSE)
    fosmat <- mqdt_to_mat(fosdt, pattern, verbose = verbose)
    fosdt <- fosdt %>% extract(, setdiff(names(.), colnames(fosmat)), with = FALSE)
    object <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(log2sites    = fosmat,
                       log2proteins = promat,
                       log2diffs = fosmat - promat),
        rowData = fosdt)
# Dequantify. Add pepcounts
    object$mqcol <- colnames(object)
    colnames(object) %<>% dequantify()
    pepcols <- paste0('Razor + unique peptides ', gsub('\\{.+\\}', '', colnames(object)))
    pepmat <- pepdt[, pepcols, with = FALSE ]
    pepmat %<>% data.matrix()
    dimnames(pepmat) <- dimnames(object)
    assays(object)$pepcounts <- pepmat
# Process / Analyze
    object %<>% process_maxquant( subgroups = subgroups,     
                                     invert = invert,      
                                    reverse = reverse, 
                               contaminants = contaminants,  
                  rm_missing_in_all_samples = rm_missing_in_all_samples,
                               localization = localization, 
                                     impute = impute,
                                    verbose = verbose )
    assays <- c('log2sites', 'pepcounts')
    for (assay in assays)  object %<>% add_assay_means(assay)
    object %<>% analyze(
        pca       = pca,           pls       = pls,
        fit       = fit,           formula   = formula,
        block     = block,         coefs     = coefs,
        contrasts = contrasts,     plot      = plot,
        label     = label,
        palette   = palette,       verbose   = verbose )
    object
}


#' @rdname read_maxquant_phosphosites
#' @export
read_phosphosites <- function(...){
    .Deprecated('read_maxquant_phosphosites')
   read_maxquant_phosphosites(...) 
}


#-----------------------------------------------------------------------------
#
#                            invert_subgroups
#
#-----------------------------------------------------------------------------

# x <- c('Ctrl_A', 'Ctrl_B')
# .invert_subgroups(x)
.invert_subgroups <- function(x, sep = guess_sep(x)){
    stri_split_fixed(x, sep)                      %>%
    lapply(rev)                                   %>%
    vapply(paste, character(1), collapse = sep)
}


#' Invert subgroups
#'
#' Invert expressions , subgroups, and sample ids
#'
#' @param  object     SummarizedExperiment
#' @param  subgroups  character vector: subgroup levels to be inversed
#' @param  sep        string: collapsed string separator
#' @return character vector or SummarizedExperiment
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' invert_subgroups(object)
#' @export
invert_subgroups <- function(
    object, 
    subgroups = slevels(object, 'subgroup'), 
    sep = guess_sep(object, 'subgroup')
){
# Assert / Msg
    if (length(subgroups)==0) return(object)
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset('subgroup', svars(object))
    assert_is_subset(subgroups, subgroup_levels(object))
    idx <- which(object$subgroup %in% subgroups)
    first <- values(object)[, idx[1]] %>% (function(y) which(!is.na(y))[[1]])
    oldvalue <- values(object)[first, idx[1]] %>% round(2) %>% as.character()
    cmessage('%sInvert subgroups %s', spaces(14), paste0(subgroups, collapse = ', '))
# Invert (log) ratios
    if (all(values(object)>0, na.rm=TRUE)){ values(object)[, idx] %<>% (function(object){1/object})
    } else {                                values(object)[, idx] %<>% (function(object){ -object})}
    newvalue <- as.character(round(values(object)[first, idx[1]], 2))
    cmessage('%sexprs    : %s -> %s', spaces(21), as.character(oldvalue), 
            as.character(newvalue))
# Invert subgroup and sampleid values
    oldsubgroups <- sdata(object)$subgroup[idx]
    newsubgroups <- vapply(oldsubgroups, .invert_subgroups, character(1),sep=sep)
    oldsampleids <- sdata(object)$sample_id[idx]
    for (i in seq_along(idx)){
        sdata(object)$subgroup[ idx[i]] <- newsubgroups[i]
        sdata(object)$sample_id[idx[i]] %<>% stri_replace_first_fixed(oldsubgroups[i], newsubgroups[i])
        snames(object)[         idx[i]] %<>% stri_replace_first_fixed(oldsubgroups[i], newsubgroups[i])
    }
    newsampleids <- sdata(object)$sample_id[idx]
    cmessage('%ssubgroups: %s -> %s', spaces(21), oldsubgroups[1], newsubgroups[1])
    cmessage('%ssampleids: %s -> %s', spaces(21), oldsampleids[1], newsampleids[1])
# Order on subgroup
    object %<>% arrange_samples_('subgroup')
    cmessage('%sOrder on subgroup', spaces(14))
    object$subgroup %<>% droplevels()

# Return
    return(object)
}

arrange_samples_ <- function(object, svars){
    idx <- do.call(order, sdata(object)[, svars, drop = FALSE])
    object %<>% extract(, idx)
    return(object)
}


#---------------------------------------------------------------------------
#
#                     demultiplex
#                         .is_multiplexed
#                         .demultiplex
#
#---------------------------------------------------------------------------


#' Demultiplex snames
#' 
#' Demultiplex maxquant samplenames
#' 
#' `WT(L).KD(H).R1{H/L}  -> KD_WT.R1`
#' `WT(1).KD(2).R1{1}    -> WT.R1`
#' `         WT.R1       -> WT.R1`
#' @param x       character vector
#' @param verbose TRUE or FALSE
#' @return character
#' @examples
#' # uniplexed / intensity / ratio
#'    demultiplex(c('KD.R1','OE.R1'))
#'    demultiplex(c('WT(L).KD(M).OE(H).R1{M}',  'WT(L).KD(M).OE(H).R1{H}'))
#'    demultiplex(c('WT(L).KD(M).OE(H).R1{M/L}','WT(L).KD(M).OE(H).R1{H/L}'))
#' # run / replicate
#'    demultiplex(c('WT(L).OE(H).R1{L}',    'WT(L).OE(H).R1{H}'))     # run
#'    demultiplex(c('WT.R1(L).OE.R1(H){L}', 'WT.R1(L).OE.R1(H){H}'))  # repl
#' # label / index
#'    demultiplex(c('WT(L).OE(H).R1{L}',    'WT(L).OE(H).R1{H}'))     # label
#'    demultiplex(c('WT(1).OE(2).R1{1}',    'WT(1).OE(2).R1{2}'))     # index
#' # with unused channels
#'    demultiplex('WT(1).KD(2).OE(3).R1{6}')
#' @md
#' @export
demultiplex <- function(x, verbose = FALSE){
    assert_is_character(x)
    assert_is_a_bool(verbose)
    . <- NULL
    if (!.is_multiplexed(x)) return(x)
    y <- unname(vapply(x, .demultiplex, character(1)))
    make.unique(y, sep = '-')
}


.is_multiplexed <- function(x){
# Components
    pattern <- '(.+)\\{(.+)\\}'
    n_open   <- stri_count_fixed(x, '(')
    n_closed <- stri_count_fixed(x, ')')
    channel <- gsub(pattern, '\\2', x)
# Multiplexed consistently ?
    y <- all(stri_detect_regex(x, pattern) &
                 n_open > 0                    &  
                 n_open == n_closed)
# All channels defined (TMT) ? 
    if (all(is_numeric_string(channel))){
        channel %<>% as.numeric()
        y %<>% and(all(as.numeric(channel) <= n_open))
    }
    y
}

.demultiplex <- function(y){
    y0 <- y
    channel <- split_extract_fixed(y, '{', 2) %>% substr(1, nchar(.)-1)
    y %<>% stri_replace_last_fixed(paste0('{', channel, '}'), '')
    
    labels     <- y %>% stri_extract_all_regex('\\([^()]+\\)') %>% unlist() %>% substr(2, nchar(.)-1)
    run        <- y %>% stri_split_regex('\\([^[()]]+\\)')     %>% unlist() %>% extract(length(.))
    biosamples <- y %>% stri_split_regex('\\([^[()]]+\\)')     %>% unlist() %>% extract(-length(.))
    biosamples %<>% stri_replace_first_regex('^[._ ]', '')
    assert_are_same_length(biosamples, labels)
    names(biosamples) <- labels
    channel %<>% stri_split_fixed('/') %>% unlist()
    
    if (!all(channel %in% labels)) return(y0)
    
    result <- paste0(biosamples[channel], collapse = '_')
    result %<>% paste0(run)
    result
}


#----------------------------------------------------------------------------
#
#                 process_maxquant
#
#----------------------------------------------------------------------------

process_maxquant <- function(
    object, subgroups, invert, contaminants, reverse, rm_missing_in_all_samples = TRUE, 
    localization = 0.75, impute, verbose
){
# Demultiplex. Infer Subgroup
    contaminant <- `Localization prob` <- NULL
    colnames(object) %<>% demultiplex(verbose = verbose)
    object$sample_id <- colnames(object)
    object %<>% add_subgroup(verbose = verbose)
    cols <- c('sample_id', 'subgroup', 'replicate', 'mqcol')
    cols %<>% intersect(svars(object))
    sdt(object) %<>% pull_columns(cols)
# Samples
    object %<>% filter_samples_available_for_some_feature(verbose = verbose)
    if (!is.null(subgroups)){
        assert_is_subset(subgroups, as.character(object$subgroup))
        object %<>% filter_samples(subgroup %in% subgroups, verbose = verbose)
    }
    object %<>% invert_subgroups(invert)
# Features
    analysis(object)$nfeatures <- c(nrow(object))
    if (!reverse)             object %<>% filter_features(   reverse == '', verbose = verbose)
    if (!contaminants)        object %<>% filter_features(contaminant== '', verbose = verbose)
    if ({rm_missing_in_all_samples})  object %<>% rm_missing_in_all_samples(verbose = verbose)
    #object %<>% filter_exprs_replicated_in_some_subgroup(verbose = verbose) # doesnt work for single-instance subgroups
    if ('Localization prob' %in% fvars(object)){                             # subgroup could be increasing concentrations or so
        object %<>% filter_features(
            `Localization prob` >= localization, verbose = verbose)  }
# Impute
    if ({{impute}})  object %<>% impute(plot = FALSE)
    object
}


#---------------------------------------------------------------------------
#
#                       log2proteins
#
#---------------------------------------------------------------------------


#' Get/Set log2proteins
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' log2proteins(object)[1:3, 1:3]
#' @rdname log2proteins
#' @export
setGeneric('log2proteins', function(object)   standardGeneric("log2proteins"))

#' @rdname log2proteins
setMethod("log2proteins", signature("SummarizedExperiment"),
function(object)   assays(object)$log2proteins)


#' @rdname log2proteins
#' @export
setGeneric('log2proteins<-',
function(object, value) standardGeneric("log2proteins<-"))

#' @rdname log2proteins
setReplaceMethod("log2proteins", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2proteins <- value
    object })

#' @rdname log2proteins
setReplaceMethod("log2proteins", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2proteins[] <- value
    object })


#---------------------------------------------------------------------------
#
#                       log2sites
#
#---------------------------------------------------------------------------


#' Get/Set log2sites
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' log2sites(object)[1:3, 1:3]
#' @rdname log2sites
#' @export
setGeneric('log2sites', function(object)   standardGeneric("log2sites"))

#' @rdname log2sites
setMethod("log2sites", signature("SummarizedExperiment"),
function(object)   assays(object)$log2sites)


#' @rdname log2sites
#' @export
setGeneric('log2sites<-',
function(object, value) standardGeneric("log2sites<-"))

#' @rdname log2sites
setReplaceMethod("log2sites", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2sites <- value
    object })

#' @rdname log2sites
setReplaceMethod("log2sites", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2sites[] <- value
    object })


#---------------------------------------------------------------------------
#
#                       log2diffs
#
#---------------------------------------------------------------------------


#' Get/Set log2diffs
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' log2diffs(object)[1:3, 1:3]
#' @rdname log2diffs
#' @export
setGeneric('log2diffs', function(object)   standardGeneric("log2diffs"))

#' @rdname log2diffs
setMethod("log2diffs", signature("SummarizedExperiment"),
function(object)   assays(object)$log2diffs)


#' @rdname log2diffs
#' @export
setGeneric('log2diffs<-',
function(object, value) standardGeneric("log2diffs<-"))

#' @rdname log2diffs
setReplaceMethod("log2diffs", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2diffs <- value
    object })

#' @rdname log2diffs
setReplaceMethod("log2diffs", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2diffs[] <- value
    object })
