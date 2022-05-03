#---------------------------------------------------------------------------
# 
#                       guess_maxquant_quantity
#
#---------------------------------------------------------------------------


#' maxquant quantity patterns
#' @examples
#' MAXQUANT_PATTERNS_QUANTITY
#' @export
MAXQUANT_PATTERNS_QUANTITY <- c(
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
#' @param x character vector
#' @return  string: value from names(MAXQUANT_PATTERNS_QUANTITY)
#' @examples
#' # file
#'     file <- download_data('fukuda20.proteingroups.txt')
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
    for (quantity in names(MAXQUANT_PATTERNS_QUANTITY)){
        pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#---------------------------------------------------------------------------
# 
#                       .read_proteingroups
#                       .read_phosphosites
#
#---------------------------------------------------------------------------

.read_proteingroups <- function(proteinfile, verbose){
# Assert
    assert_all_are_existing_files(proteinfile)
    assert_all_are_matching_fixed(proteinfile, 'roups.txt')
# Read
    if (verbose)  message('\tRead ', proteinfile)
    prodt <- fread(proteinfile, colClasses = c(id = 'character'), integer64 = 'numeric')
    n0 <- nrow(prodt)
    pattern <- MAXQUANT_PATTERNS_QUANTITY[[guess_maxquant_quantity(proteinfile)]]
    anncols <- c('id', 'Majority protein IDs', 'Reverse', 
                 'Potential contaminant', 'Contaminant', 'Fasta headers')#, 'Phospho (STY) site IDs')
    anncols %<>% intersect(names(prodt))
    valcols <- grep(pattern, names(prodt), value = TRUE)
    pepcols <- grep('Razor + unique peptides ', names(prodt), fixed = TRUE, value = TRUE)
    prodt %<>% extract(, c(anncols, pepcols, valcols), with = FALSE)
    digits <- ceiling(log10(nrow(prodt)))
    if (verbose)  message('\t\t', nrow(prodt), 
                          ' proteins, contaminants, reverse')
# Return
    names(prodt) %<>% stri_replace_first_fixed('Contaminant', 'Potential contaminant') # older MaxQuant
    setnames(prodt, c('id', 'Majority protein IDs'), c('proId', 'Uniprot'))
    prodt
}


.read_phosphosites <- function(phosphofile, proteinfile, verbose){
# Assert
    assert_all_are_existing_files(c(proteinfile, phosphofile))
    assert_all_are_matching_fixed(proteinfile, 'roups.txt')
    assert_all_are_matching_fixed(phosphofile, 'ites.txt')
    if (!is.null(fastafile))  assert_all_are_matching_fixed(fastafile, 'fasta')
# Read    
    if (verbose)  message('\tRead ', phosphofile)
    fosdt <- fread(phosphofile, colClasses = c(id = 'character'), integer64 = 'numeric')
    pattern <- MAXQUANT_PATTERNS_QUANTITY[[guess_maxquant_quantity(phosphofile)]]
    anncols <- c('id', 'Protein group IDs', 'Proteins', 
                 'Positions within proteins', 'Amino acid', 
                 'Reverse', 'Contaminant', 'Potential contaminant', 'Fasta headers')
    anncols %<>% intersect(names(fosdt))
    valcols <- grep(pattern, names(fosdt), value = TRUE)
    valcols %<>% extract(stri_detect_fixed(., '___1'))
    fosdt %<>% extract(, c(anncols, valcols), with = FALSE)
    digits <- ceiling(log10(nrow(fosdt)))
    if (verbose)  message('\t\tRead   ', formatC(nrow(fosdt), digits = digits), 
                          ' phosphosites in proteins, contaminants, reverse')
    names(fosdt) %<>% stri_replace_first_fixed('Contaminant', 'Potential contaminant')
        # Older MaxQuant files use 'Contaminant'
# Filter
    fosdt <- fosdt[stri_count_fixed(`Protein group IDs`, ';') == 0]
    if (verbose)  message('\t\tRetain ', formatC(nrow(fosdt), digits = digits), 
                          ' phosphosites in single proteingroup')
    idx <- rowSums(fosdt[, valcols, with = FALSE], na.rm = TRUE) > 0
    fosdt <- fosdt[which(idx)]
    if (verbose)  message('\t\tRetain ', formatC(nrow(fosdt), digits = digits), 
                          ' phosphosites with signal in some sample')
# Rename 
    names(fosdt) %<>% stri_replace_first_fixed('___1', '')
    setnames(fosdt, c('id',    'Protein group IDs', 'Proteins'), 
                 c('fosId', 'proId',             'Uniprot'))
}


#---------------------------------------------------------------------------
# 
#                      drop_differing_uniprots
#                      add_feature_id
#                      process_maxquant
#
#---------------------------------------------------------------------------


drop_differing_uniprots <- function(fosdt, prodt, verbose){
# Extract annotation cols
    if (verbose)  message('\t\t\tKeep proteingroup Uniprots')
    proanncols <- c('proId', 'Uniprot', 'Reverse', 'Potential contaminant')
    fosanncols <- c('fosId', 'proId', 'Uniprot', 'Positions within proteins', 
                    'Reverse', 'Potential contaminant')
    proanndt <- prodt[, proanncols, with = FALSE]
    fosanndt <- fosdt[, fosanncols, with = FALSE]
# Separate contaminants-reverse from other proteins.
    conrevdt <- fosanndt[Reverse == '+' | `Potential contaminant` == '+']
    fosanndt <- fosanndt[Reverse == ''  & `Potential contaminant` == '']
    proanndt <- proanndt[Reverse == ''  & `Potential contaminant` == '']
    fosanndt[, c('Reverse', 'Potential contaminant') := NULL]
    proanndt[, c('Reverse', 'Potential contaminant') := NULL]
    conrevdt[, c('Reverse', 'Potential contaminant') := NULL]
# Merge
    fosanndt %<>% uncollapse(Uniprot, `Positions within proteins`, sep = ';')
    proanndt %<>% uncollapse(Uniprot,                              sep = ';')
    fosanndt %<>% merge(proanndt, by = c('proId', 'Uniprot'))
    fosanndt %<>% extract(order(as.integer(fosId)))
    fosanndt %<>% extract(, lapply(.SD, paste0, collapse = ';'), by = 'fosId')
# Add back contaminants/reverse
    fosanndt %<>% rbind(conrevdt)
    fosanndt[, proId := NULL]
    fosdt[, c('Uniprot', 'Positions within proteins') := NULL]
    fosdt %<>% merge(fosanndt, by = 'fosId', all.x = TRUE, sort = FALSE)
    fosdt %<>% pull_columns(names(fosanndt))
    fosdt    
}


add_feature_id <- function(dt){
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    dt[, feature_id := `Entry names`]
    dt[Isoform!=1, feature_id := paste0(feature_id, '(', stri_replace_all_fixed(Isoform, ';', ''), ')')]
    if (idcol=='fosId'){
        dt$feature_id %<>% paste0('-',dt$`Amino acid`)
        dt$feature_id %<>% paste0('-', split_extract_fixed(dt$`Positions within proteins`, ';', 1))
    }
    dt[                Reverse == '+', feature_id := paste0('REV_', idcol)]
    dt[`Potential contaminant` == '+', feature_id := paste0('CON_', idcol)]
        # Contaminants and actual proteins get mixed up. For naming: make sure to use contaminant.
    dt %<>% pull_columns(c(idcol, 'feature_id'))
    dt[]
}


#-----------------------------------------------------------------------------
#
#                       read_fastahdrs
#                           parse_fastahdrs
#
#-----------------------------------------------------------------------------

uncollapse <- function(dt, ..., sep = ';'){
    dt %>% 
    separate_rows(..., sep = sep) %>%
    data.table()
}

nastring_to_nachar <- function(x){ x[x=='NA'] <- NA_character_;  x }

nastring_to_1 <- function(x){ x[x=='NA'] <- 1; x}

extract_reviewed <- function(`Fasta headers`){
    `Fasta headers` %>% 
    substr(2, 3)    %>%
    equals('sp')    %>%
    as.integer()
}

extract_entryname <- function(`Fasta headers`){
    `Fasta headers` %>% 
    split_extract_fixed(' ', 1) %>% 
    split_extract_fixed('|', 3) %>% 
    split_extract_fixed('_', 1)
}

extract_genename <- function(`Fasta headers`){
    `Fasta headers` %>% 
    split_extract_fixed('GN=', 2)      %>% 
    split_extract_fixed(' ', 1)}

extract_uniprot <- function(`Fasta headers`){
    `Fasta headers` %>%
    split_extract_fixed(' ', 1) %>% 
    split_extract_fixed('|', 2) 
}

extract_canonical <- function(`Fasta headers`){
    `Fasta headers`   %>%
    extract_uniprot() %>%
    split_extract_fixed('-', 1) 
}

extract_isoform <- function(`Fasta headers`){
    `Fasta headers`   %>%
    extract_uniprot() %>%
    split_extract_fixed('-', 2) %>%
    nastring_to_1() %>%
    as.integer()
}

extract_proteinname <- function(`Fasta headers`){
    `Fasta headers` %>%
        split_extract_regex('_[A-Z]+ ', 2) %>% 
        split_extract_regex(' OS=', 1)
}


extract_fragment <- function(`Fasta headers`){
    `Fasta headers`                %>%
        extract_proteinname()          %>%
        stri_detect_fixed(  'ragment') %>%
        as.integer()
}

extract_existence <- function(`Fasta headers`){
    `Fasta headers`      %>%
    split_extract_fixed('PE=', 2) %>%
    split_extract_fixed(' ', 1) %>%
    nastring_to_nachar() %>%
    as.integer()
}

parse_fastahdrs <- function(`Fasta headers`){
    dt <- data.table(                                           # Reviewed
        Reviewed        = extract_reviewed(   `Fasta headers`), #   0 tr
        `Entry names`   = extract_entryname(  `Fasta headers`), #   1 sp
        `Gene names`    = extract_genename(   `Fasta headers`),
        Uniprot         = extract_uniprot(    `Fasta headers`), # Existence 
        Canonical       = extract_canonical(  `Fasta headers`), #   1 protein 
        Isoform         = extract_isoform(    `Fasta headers`), #   1 transcript
        `Protein names` = extract_proteinname(`Fasta headers`), #   2 homolog
        Fragment        = extract_fragment(   `Fasta headers`), #   3 prediction
        Existence       = extract_existence(  `Fasta headers`)) #   4 uncertain
    dt[, Existence := unique(.SD)[, Existence[!is.na(Existence)]], by = 'Canonical']
    # `unique`: for phosphosites the Fasta headers are sometimes
    #  replicated (when protein has multiple phosphosites)
    #  This duplication needs to be eliminated before proceeding.
    dt[]
}


#' Read headers from uniprot fastafile
#'
#' @param fastafile    path to fasta file
#' @param fastafields  character vector
#' @param verbose      bool
#' @examples
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' read_fastafile_headers(fastafile)
#' @return data.table(uniprot, genename, proteinname, reviewed, existence)
#' @note Existence values are always those of the canonical isoform
#'       (no isoform-level resolution for this field)
#' @export
read_fastahdrs <- function(fastafile, verbose = TRUE){
# Assert
    if (is.null(fastafile)) return(NULL)
    if (!requireNamespace('seqinr', quietly = TRUE)){
        stop("BiocManager::install('seqinr'). Then re-run.") }
    assert_all_are_existing_files(fastafile)
# Read
    if (verbose) message('\tRead ', fastafile)
    fastahdrs <- seqinr::read.fasta(fastafile)
    fastahdrs %<>% vapply(attr, character(1), 'Annot') %>% unname()
# Parse
    parse_fastahdrs(fasta)
}


#-----------------------------------------------------------------------------
#
#                       curate_uniprots
#
#-----------------------------------------------------------------------------


# tests/testthat/test_3_curate_uniprots.R
#' Curate annotations
#' @param dt data.table
#' @return data.table
#' @examples
#' proteinfile <- download_data('billing19.proteingroups.txt')
#' phosphofile <- download_data('billing19.phosphosites.txt')
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' prodt <- .read_pro(proteinfile)
#' fosdt <- .read_fos(phosphofile)
#' curate_uniprots(prodt)[, 1:8]
#' curate_uniprots(fosdt)[, 1:8]
#' curate_uniprots(prodt, fastafile = fastafile)[, 1:8]
#' curate_uniprots(fosdt, fastafile = fastafile)[, 1:8]
#' @export
curate_uniprots <- function(dt, fastafile = NULL, verbose = TRUE){
    # Split
    if (verbose)  message('\tCurate uniprots')
    if (verbose)  message('\t\tRm Contaminants/Reverse')
    idx <- dt$Reverse == '+' | dt$`Potential contaminant` == '+' 
    conrevdt <- dt[idx]
    dt %<>% extract(!idx)
    # Curate using MaxQuant
    if (verbose)  message('\t\tMaxQuant Fastahdrs')
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    anndt <- dt[, c(idcol, 'Fasta headers'), with = FALSE]
    if (verbose)  message('\t\t\tUncollapse, Drop truncated, Parse')
    anndt %<>% uncollapse(`Fasta headers`, sep = ';')
    anndt %<>% extract( stri_count_fixed(`Fasta headers`, '|') == 2)
    anndt %<>% cbind(parse_fastahdrs(anndt$`Fasta headers`))
    anndt$`Fasta headers` <- NULL
    anndt %<>% drop_inferior()
    
    if (verbose)  message('\t\t\tCollapse')
    anndt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = idcol)
    # Curate using fastafile
    if (!is.null(fastafile)){
        if (verbose)  message('\t\tFastafile Fastahdrs')
        if (verbose)  message('\t\t\tUncollapse Uniprot accessions')
        anndt2 <- dt[, c(idcol, 'Uniprot'), with = FALSE]
        anndt2 %<>% uncollapse(`Uniprot`, sep = ';')
        if (verbose)  message('\t\t\tRead/Parse Fastafile ', fastafile)
        fastahdrs <- read_fastahdrs(fastafile)
        anndt2 %<>% merge(fastadt, by = 'Uniprot', sort = FALSE)
        anndt2 %<>% drop_inferior()
        if (verbose)  message('\t\t\tCollapse')
        anndt2 %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = idcol)
        missingInFasta  <- setdiff(  dt[[idcol]], anndt2[[idcol]])
        if (verbose)  message('\t\t\t', length(missingInFasta), 
                              ' proteins missing in fastafile - revert to MaxQuant headers for these')
        anndt %<>% extract(missingInFasta,  on = idcol)
        anndt %<>% rbind(anndt2)
        anndt %<>% extract(order(as.integer(get(idcol))))
    }
    # Merge back annotations
    if (verbose)  message('\t\tBring back Contaminants/Reverse')
    setnames(anndt, 'Uniprot', 'Curated')
    dt %<>% merge(anndt, by = idcol, sort = FALSE)
    newcols <- setdiff(names(dt), names(conrevdt))   
    conrevdt[, (newcols) := '']
    dt %<>% rbind(conrevdt)
    dt %<>% pull_columns(names(anndt))
    dt$`Fasta headers` <- NULL
    dt[order(as.integer(get(idcol)))]
}

drop_inferior <- function(anndt, verbose = TRUE){
    idcol <- if ('fosId' %in% names(anndt)) 'fosId' else 'proId'
    
    if (verbose)  message('\t\t\tDrop proteins with NA `Entry names` / `Protein names`')
    anndt <- anndt[`Entry names` != 'NA']
    anndt <- anndt[`Protein names` != 'NA']
    
    if (verbose)  message('\t\t\tWithin ', idcol, ': drop trembl    in favour of swissprot')
    anndt <- anndt[, .SD[ Reviewed == max(Reviewed) ], by = idcol]
    
    if (verbose)  message('\t\t\t       ','     ','  drop fragments in favour of full proteins')
    anndt <- anndt[, .SD[ Fragment == min(Fragment) ], by = idcol]
    
    if (verbose)  message('\t\t\t       ','     ','  drop worse existences', 
                          ' (1=protein, 2=transcript, 3=homolog, 4=prediction, 5=uncertain)')
    if (!any(is.na(anndt$Existence))){
        anndt <- anndt[, .SD[Existence == min(Existence)], by = idcol] }
    anndt[, c('Reviewed', 'Fragment', 'Existence') := NULL]
    
    if (verbose)  message("\t\t\tDrop 'Isoform x of ' in `Protein names`")
    anndt$`Protein names` %<>% stri_replace_first_regex('Isoform [0-9A-Z]+ of ', '')
    
    if (verbose)  message("\t\t\tDrop '(Fragment)'    in `Protein names`")
    anndt$`Protein names` %<>% stri_replace_first_fixed(' (Fragment)', '')
    anndt[order(as.integer(get(idcol)), Canonical, Isoform)]
}


#---------------------------------------------------------------------------
# 
#                 mqdt_to_mat
#                 dequantify
#                     label2index
#
#---------------------------------------------------------------------------


mqdt_to_mat <- function(dt, pattern, verbose){
    mat <- dt[, .SD, .SDcols = patterns(pattern)]
    mat %<>% data.matrix()
    rownames(mat) <- dt$feature_id
    mat %<>% zero_to_na(verbose = verbose) 
    mat %<>% nan_to_na( verbose = verbose)
    mat <- log2(1 + mat)
    mat
}


#' Dequantify maxquant snames
#' 
#' Drop quantity ('Reporter intensity'). \cr
#' Encode {channel} as suffix.
#' 
#' `               Ratio H/L WT(L).KD(H).R1  ->  WT(L).KD(H).R1{H/L}`
#' `                    LFQ intensity WT.R1  ->  WT.R1`
#' `Reporter intensity 0 WT(126).KD(127).R1  ->  WT(1).KD(2).R1{1}`
#' @param x        `character`
#' @param quantity `'Ratio',              'Ratio normalized'`,  \cr
#'                 `'LFQ intensity'`, \cr
#'                 `'Intensity',          'Intensity labeled'`
#'                 `'Reporter intensity', 'Reporter intensity corrected'`
#' @param verbose  `TRUE` or `FALSE`
#' @return `character`
#' @examples
#' dequantify(c('Ratio H/L WT(L).KD(M).OE(H).R1',             # Ratios
#'              'Ratio M/L WT(L).KD(M).OE(H).R1'))
#' 
#' dequantify(c('Ratio H/L normalized WT(L).KD(M).OE(H).R1',  # Norm. Ratios
#'              'Ratio M/L normalized WT(L).KD(M).OE(H).R1'))
#' 
#' dequantify(c('LFQ intensity WT.R1',                        # LFQ intensity
#'              'LFQ intensity KD.R1'))
#' 
#' dequantify(c('Reporter intensity 1 WT(126).KD(127).R1',    # Rep.intensities
#'              'Reporter intensity 2 WT(126).KD(127).R1'))
#' @md
#' @export
dequantify <- function(
    x, quantity = guess_maxquant_quantity(x), verbose  = FALSE
){
# x = multiplex + channel. Return multiplex if single channel.
# Decompose multiplex and channel
    pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
    if (quantity == 'Intensity'){
        multiplex <- stri_replace_first_regex(x, pattern, '$1')
        channel   <- rep('', length(multiplex))
    } else {
        multiplex <- stri_replace_first_regex(x, pattern, '$2')
        channel   <- stri_replace_first_regex(x, pattern, '$1')
    }
# Reporter intensity
    if (quantity %in% c('Reporter intensity', 'Reporter intensity corrected')){
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
    if (verbose) message('\t\tStandardize snames: ', x[1], '  ->  ', cleanx[1])
    return(cleanx)
}


#' Convert labels into indices
#' @param x  `character`
#' @examples
#' label2index('Reporter intensity 0 WT(0).KD(1).OE(2).R1')
#' label2index('Reporter intensity 1 WT(1).KD(2).OE(3).R1')
#' label2index('Reporter intensity 0 WT(126).KD(127).OE(128).R1')
#' label2index('Reporter intensity 1 WT(126).KD(127).OE(128).R1')
#' label2index(x)
#' @export
label2index <- function(x){
    labels <- unlist(stri_extract_all_regex(x, '\\(.+?\\)'))
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
#                   read_proteingroups
#                   read_phosphosites
#
#---------------------------------------------------------------------------


#' Read proteingroups / phosphosites
#'
#' @param proteinfile   proteingroups file
#' @param phosphofile   phosphosites  file
#' @param fastafile     fastafile or NULL
#' @param curate        whether to curate uniprots
#' @param subgroups     character / NULL : subset of subgroups to retain
#' @param contaminants  whether to include contaminants
#' @param reverse       whether to include reverse hits
#' @param localization  min site localization probability
#' @param invert        character : subset of subgroups to invert
#' @param impute        whether to impute
#' @param plot          whether to plot
#' @param pca           whether to pca
#' @param fit           model fit engine: 'limma', 'lm', 'lmer', 'lme'
#' @param formula       model formula
#' @param block         block var (sdt)
#' @param coefs         character: coefficients to test
#' @param contrastdefs  character: coefficient contrasts to test
#' @param feature_id    string: summary plot feature
#' @param sample_id     string: summary plot sample
#' @param palette       character: color palette
#' @param verbose       whether to msg
#' @return SummarizedExperiment
#' @examples 
#' proteinfile <- download_data('billing19.proteingroups.txt')
#' phosphofile <- download_data('billing19.phosphosites.txt')
#'   fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' subgroups <- sprintf('%s_STD', c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'))
#' pro <- read_proteingroups(proteinfile,              subgroups = subgroups)
#' fos <- read_phosphosites( phosphofile, proteinfile, subgroups = subgroups)
#' fos <- read_phosphosites( phosphofile, proteinfile, fastafile = fastafile, subgroups = subgroups)
#' @export
read_proteingroups <- function(
    proteinfile, fastafile = NULL, curate = TRUE, 
    subgroups = NULL, invert = character(0),
    contaminants = FALSE, reverse = FALSE, impute = TRUE,
    plot = FALSE, pca = plot, fit = if (plot) 'limma' else NULL,
    formula = NULL, block = NULL, coefs = NULL, contrastdefs = NULL, 
    feature_id = NULL, sample_id = NULL, palette = NULL, verbose = TRUE
){
# Assert
    assert_is_a_string(quantity)
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS_QUANTITY))
    assert_is_a_bool(verbose)
# Read
    prodt <- .read_proteingroups(proteinfile = proteinfile, verbose = verbose)
    if (curate)  prodt %<>% curate_uniprots(fastafile = fastafile)
    prodt %<>% add_feature_id()
# SumExp
    if (verbose)  message('\tCreate SummarizedExperiment')
    pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
    promat <- mqdt_to_mat(prodt, pattern, verbose = verbose)
    pepcols <- names(prodt) %>% extract(stri_detect_fixed(., 'eptides'))
    pepdt <- prodt[, pepcols, with = FALSE]
    prodt %<>% extract(, names(prodt) %>% setdiff(colnames(promat)) %>% setdiff(names(pepdt)), with = FALSE)
    object <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(log2pro = promat), 
        rowData = prodt)
# Dequantify. Add pepcounts
    object$mqcol <- colnames(object)
    colnames(object) %<>% dequantify()
    pepcols <- paste0('Razor + unique peptides ', gsub('\\{.+\\}', '', colnames(object)))
    pepmat <- pepdt[, pepcols, with = FALSE ]
    pepmat %<>% data.matrix()
    dimnames(pepmat) <- dimnames(object)
    assays(object)$pepcounts <- pepmat
    object %<>% process_maxquant(
        subgroups = subgroups,      invert = invert,       reverse = reverse,
     contaminants = contaminants,   impute = impute,       verbose = verbose)
    object %<>% analyze(
              pca = pca,               fit = fit,          formula = formula,
            block = block,           coefs = coefs,   contrastdefs = contrastdefs, 
          verbose = verbose,          plot = plot,      feature_id = feature_id, 
        sample_id = sample_id,     palette = palette )
    object
}


#' @rdname read_proteingroups
#' @export
read_phosphosites <- function(
    phosphofile, proteinfile, fastafile = NULL, curate = TRUE, 
    subgroups = NULL, invert = character(0), 
    contaminants = FALSE, reverse = FALSE, localization = 0.75, 
    plot = FALSE, pca = plot, fit = if (plot) 'limma' else NULL,  
    formula = NULL, block = NULL, coefs = NULL, contrastdefs = NULL, 
    feature_id = NULL, sample_id = NULL, palette = NULL, verbose = TRUE
){
# Assert
    assert_is_a_string(quantity)
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS_QUANTITY))
    assert_is_a_bool(verbose)
# Read
    prodt <- .read_proteingroups(proteinfile = proteinfile, verbose = verbose)
    fosdt <- .read_phosphosites(phosphofile = phosphofile, proteinfile = proteinfile, verbose = verbose)
    fosdt %<>% drop_differing_uniprots(prodt, verbose = verbose)
    if (curate)  fosdt %<>% curate_uniprots(fastafile = fastafile)
    fosdt %<>% add_feature_id()
    prodt %<>% extract(fosdt$proId, on = 'proId')
# SumExp
    if (verbose)  message('\tCreate SummarizedExperiment')
    pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
    promat <- mqdt_to_mat(prodt, pattern, verbose = verbose)
    pepcols <- names(prodt) %>% extract(stri_detect_fixed(., 'eptides'))
    pepdt <- prodt[, pepcols, with = FALSE]
    prodt %<>% extract(, names(prodt) %>% setdiff(colnames(promat)) %>% setdiff(names(pepdt)), with = FALSE)
    fosmat <- mqdt_to_mat(fosdt, pattern, verbose = verbose)
    fosdt <- fosdt %>% extract(, setdiff(names(.), colnames(fosmat)), with = FALSE)
    object <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(log2fos = fosmat,
                       log2pro = promat,
                       log2dif = fosmat - promat),
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
    object %<>% process_maxquant(
        subgroups = subgroups,     invert = invert,      
          reverse = reverse, contaminants = contaminants,  localization = localization, 
           impute = impute,      verbose = verbose)
    object %<>% analyze(
              pca = pca,              fit = fit,         formula = formula, 
            block = block,          coefs = coefs,  contrastdefs = contrastdefs, 
          verbose = verbose,         plot = plot,     feature_id = feature_id, 
        sample_id = sample_id,    palette = palette)
    object
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
#' @param  x          character vector or SummarizedExperiment
#' @param  sep        string: collapsed string separator
#' @param  subgroups  character vector: subgroup levels to be inversed
#' @param  ... to enable S3 method dispatch
#' @return character vector or SummarizedExperiment
#' @examples
#'
#' # SummarizedExperiment
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(file, plot=FALSE)
#'     invert_subgroups(object)
#' @export
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
    message('\t\tInvert subgroups ', paste0(subgroups, collapse = ', '))
# Invert (log) ratios
    if (all(values(object)>0, na.rm=TRUE)){ values(object)[, idx] %<>% (function(object){1/object})
    } else {                                values(object)[, idx] %<>% (function(object){ -object})}
    newvalue <- as.character(round(values(object)[first, idx[1]], 2))
    message('\t\t\texprs    : ', as.character(oldvalue), ' -> ', 
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
    message('\t\t\tsubgroups: ', oldsubgroups[1], ' -> ', newsubgroups[1])
    message('\t\t\tsampleids: ', oldsampleids[1], ' -> ', newsampleids[1])
# Order on subgroup
    object %<>% arrange_samples_('subgroup')
    message('\t\tOrder on subgroup')
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
#' @param channels    character: c('M', 'H')
#' @param biosamples  list: list(c(L='t0', M='t1', H='t2'), 
#'                               c(L='t0', M='t1', H='t2'))
#' @param replicate   character: c('_R1', '_R2')
#' @return character
#' @examples
#' # uniplexed / intensity / ratio
#'    demultiplex(c('KD.R1','OE.R1'))
#'    demultiplex(c('WT(L).KD(M).OE(H).R1{M}',  'WT(L).KD(M).OE(H).R1{H}'))
#'    demultiplex(c('WT(L).KD(M).OE(H).R1{M/L}','WT(L).KD(M).OE(H).R1{H/L}'))
#'
#' # run / replicate
#'    demultiplex(c('WT(L).OE(H).R1{L}',    'WT(L).OE(H).R1{H}'))     # run
#'    demultiplex(c('WT.R1(L).OE.R1(H){L}', 'WT.R1(L).OE.R1(H){H}'))  # repl
#'
#' # label / index
#'    demultiplex(c('WT(L).OE(H).R1{L}',    'WT(L).OE(H).R1{H}'))     # label
#'    demultiplex(c('WT(1).OE(2).R1{1}',    'WT(1).OE(2).R1{2}'))     # index
#' 
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
    if (all(assertive::is_numeric_string(channel))){
        channel %<>% as.numeric()
        y %<>% `&`(all(as.numeric(channel) <= n_open))
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
    assertive::assert_are_same_length(biosamples, labels)
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
    object, subgroups, invert, contaminants, reverse, localization = 0.75, 
    impute, verbose
){
# Demultiplex. Infer Subgroup
    colnames(object) %<>% demultiplex(verbose = verbose)
    object$sample_id <- colnames(object)
    object %<>% add_subgroup(verbose = verbose)
    sdt(object) %<>% pull_columns(c('sample_id', 'subgroup', 'replicate', 'mqcol'))
# Samples
    object %<>% filter_samples_available_for_some_feature(verbose = verbose)
    object %<>% filter_samples(subgroup %in% subgroups, verbose = verbose)
    object %<>% invert_subgroups(invert)
# Features
    if (verbose) message('\tFilter features')
    if (!reverse)      object %<>% filter_features(               Reverse == '', verbose = verbose)
    if (!contaminants) object %<>% filter_features(`Potential contaminant`== '', verbose = verbose)
    object %<>% rm_missing_in_all_samples(verbose = verbose)
    object %<>% filter_exprs_replicated_in_some_subgroup(verbose = verbose)
    if (!'Localization prob' %in% fvars(object)){
        object %<>% filter_features(
            `Localization prob` >= localization, verbose = verbose)  }
# Impute
    if (impute)  object %<>% impute_systematic_nondetects(plot = FALSE)
    object
}



#---------------------------------------------------------------------------
#
#                       annotate_uniprot_ws
#
#---------------------------------------------------------------------------

paste_unique <- function(x, collapse) paste0(unique(x), collapse=collapse)

.annotate_uniprot_ws <- function(fdt, upws, columns='ENSEMBL', collapse=';'){
    # Map uniprot -> columns
        resdt <- data.table(UniProt.ws::select(upws, fdt$uniprot, columns))
    # Trim whitespace
        resdt <- resdt[, lapply(.SD, trimws),                          by = 'UNIPROTKB']
    # Collapse per uniprot
        resdt <- resdt[, lapply(.SD, paste_unique, collapse=collapse), by = 'UNIPROTKB']
    # Ensure original order
        merge.data.table(
            fdt, resdt, by.x='uniprot', by.y='UNIPROTKB', all.x = TRUE, sort = FALSE)
}

#' Annotate uniprotids using UniProt.ws
#'
#' Annotate uniprotids in data.table or SummarizedExperiment
#'
#' data.table: column "uniprot" with uniprotid values (e.g. "Q5JTY5"). \cr
#'
#' SummarizedExperiment:  svar "feature_id" with collapsed uniprot values for 
#' a single protein ("A0A068F9P7") or a proteingroup ("A0A068F9P7;F1Q7I4").
#' On these the mapping is performed. 1:many mappings are collapsed and only then returned.
#' @param x          data.table/SummarizedExperiment with (s)var "uniprot" \cr
#'                   with single (data.table) or collapsed (SummarizedExperiment) uniprot values
#' @param upws       return value of Uniprot.ws::Uniprot.ws()
#' @param columns    subset of UniProt.ws::columns(up)
#' @param collapse   string: used in paste(collapse = .)
#' @param ...        used for S3 generic definition
#' @return SummarizedExperiment/data.table
#' @examples
#' # data.table
#'     x <- data.table::data.table(uniprot = c('A0A068F9P7', 'Q7ZVA2'))
#'     # upws <- UniProt.ws::UniProt.ws(taxId=7955)
#'     # annotate_uniprot_ws(x, upws)
#'
#' # SummarizedExperiment
#'     require(magrittr)
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     x <- read_proteingroups(file, plot=FALSE)
#'     x %<>% extract(1:10, )
#'     fdata(x)[1:3, ]
#'     # x %<>% annotate_uniprot_ws(upws)
#'     # fdata(x)[1:3, ]
#' @export
annotate_uniprot_ws <- function(x, ...)  UseMethod('annotate_uniprot_ws')


#' @rdname annotate_uniprot_ws
#' @export
annotate_uniprot_ws.data.table <- function(
    x, upws, columns=c('ENSEMBL'), collapse=';', ...
){
    # Assert valid inputs
        if (!requireNamespace('UniProt.ws', quietly = TRUE)){
            message("`BiocManager::install('UniProt.ws')`. Then re-run.")
            return(x)
        }
        assert_is_data.table(x)
        assert_is_subset('uniprot', names(x))
        assert_is_identical_to_true(all(!is.na(x$uniprot)))
        assert_is_all_of(upws, 'UniProt.ws')
        assert_is_character(columns)
        assert_is_subset(columns, UniProt.ws::columns(upws))
        assert_is_a_string(collapse)
        columns %<>% setdiff(names(x))
    # Chunk data (to avoid UniProt.ws warning)
        n <- ceiling(nrow(x)/99)
        chunks <- rep(seq_len(n), each=99)
        chunks %<>% extract(seq_len(nrow(x)))
        chunk <- NULL
        x[, chunk := chunks]
    # Call backend for each chunk
        x <- x[, .annotate_uniprot_ws(.SD, upws, columns=columns, collapse=collapse), by='chunk']
    # Return
        x$chunk <- NULL
        x[]
}

#' @rdname annotate_uniprot_ws
#' @export
annotate_uniprot_ws.SummarizedExperiment <- function(
    x, upws, columns = c('ENSEMBL'), collapse=';', ...
){
    # Split proteingroups into proteins and extract uniprot
        fdt <- data.table(fdata(x))
        fdt %<>% uncollapse(uniprot, sep=collapse)
        fdt[, uniprot := split_extract_fixed(uniprot, '-', 1)]
    # Annotate
        uniprot <- NULL
        fdt %<>% extract(uniprot %in% UniProt.ws::keys(upws, keytype = 'UNIPROTKB'))
        fdt %<>% annotate_uniprot_ws.data.table(upws, columns=columns, collapse=collapse)
    # Collapse
        fdt <- fdt[, lapply(.SD, trimws),                          by='feature_id'] # trim whitespace
        fdt <- fdt[, lapply(.SD, paste_unique, collapse=collapse), by='feature_id'] # collapse
        x %<>% merge_fdata(fdt)
        x
}



#---------------------------------------------------------------------------
#
#                       log2pro
#
#---------------------------------------------------------------------------


#' Get/Set log2pro
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' log2pro(object)[1:3, 1:3]
#' @rdname log2pro
#' @export
setGeneric('log2pro', function(object)   standardGeneric("log2pro"))

#' @rdname log2pro
setMethod("log2pro", signature("SummarizedExperiment"),
function(object)   assays(object)$log2pro)


#' @rdname log2pro
#' @export
setGeneric('log2pro<-',
function(object, value) standardGeneric("log2pro<-"))

#' @rdname log2pro
setReplaceMethod("log2pro", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2pro <- value
    object })

#' @rdname log2pro
setReplaceMethod("log2pro", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2pro[] <- value
    object })


#---------------------------------------------------------------------------
#
#                       log2fos
#
#---------------------------------------------------------------------------


#' Get/Set log2fos
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' log2fos(object)[1:3, 1:3]
#' @rdname log2fos
#' @export
setGeneric('log2fos', function(object)   standardGeneric("log2fos"))

#' @rdname log2fos
setMethod("log2fos", signature("SummarizedExperiment"),
function(object)   assays(object)$log2fos)


#' @rdname log2fos
#' @export
setGeneric('log2fos<-',
function(object, value) standardGeneric("log2fos<-"))

#' @rdname log2fos
setReplaceMethod("log2fos", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2fos <- value
    object })

#' @rdname log2fos
setReplaceMethod("log2fos", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2fos[] <- value
    object })


#---------------------------------------------------------------------------
#
#                       log2dif
#
#---------------------------------------------------------------------------


#' Get/Set log2dif
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' log2dif(object)[1:3, 1:3]
#' @rdname log2dif
#' @export
setGeneric('log2dif', function(object)   standardGeneric("log2dif"))

#' @rdname log2dif
setMethod("log2dif", signature("SummarizedExperiment"),
function(object)   assays(object)$log2dif)


#' @rdname log2dif
#' @export
setGeneric('log2dif<-',
function(object, value) standardGeneric("log2dif<-"))

#' @rdname log2dif
setReplaceMethod("log2dif", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$log2dif <- value
    object })

#' @rdname log2dif
setReplaceMethod("log2dif", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$log2dif[] <- value
    object })
