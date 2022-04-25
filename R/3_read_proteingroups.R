#===========================================================================
#
#                            Getters/Setters
#
#===========================================================================

#' @title Get/Set occupancies
#' @description Get / Set phosphosite occupancies matrix
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' occupancies(object)
#' occupancies(object) <- values(object)
#' occupancies(object)[1:3, 1:3]
#' @rdname occupancies
#' @export
setGeneric('occupancies', function(object)   standardGeneric("occupancies"))

#' @rdname occupancies
setMethod("occupancies", signature("SummarizedExperiment"),
function(object)   assays(object)$occupancies)


#' @rdname occupancies
#' @export
setGeneric('occupancies<-',
function(object, value) standardGeneric("occupancies<-"))

#' @rdname occupancies
setReplaceMethod("occupancies", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$occupancies <- value
    object })

#' @rdname occupancies
setReplaceMethod("occupancies", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$occupancies[] <- value
    object })



#' @title Get/Set proteingroups
#' @description Get / Set proteingroups matrix
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' proteingroups(object)[1:3, 1:3]
#' @rdname proteingroups
#' @export
setGeneric('proteingroups', function(object)   standardGeneric("proteingroups"))

#' @rdname proteingroups
setMethod("proteingroups", signature("SummarizedExperiment"),
function(object)   assays(object)$proteingroups)


#' @rdname proteingroups
#' @export
setGeneric('proteingroups<-',
function(object, value) standardGeneric("proteingroups<-"))

#' @rdname proteingroups
setReplaceMethod("proteingroups", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)$proteingroups <- value
    object })

#' @rdname proteingroups
setReplaceMethod("proteingroups", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)$proteingroups[] <- value
    object })



#============================================================================
#
#                            Constants
#
#============================================================================


#' proteingroup fvars
#' @noRd
PROTEINGROUP_FVARS <- c(
    'id', 'Majority protein IDs', 'Protein names', 'Gene names',
    'Contaminant', 'Potential contaminant', 'Reverse', 'Phospho (STY) site IDs')


#' phosphosites fvars
#' @noRd
PHOSPHOSITE_FVARS <- c('id', 'Protein group IDs', 'Proteins', 'Protein names',
    'Gene names', 'Positions within proteins', 'Localization prob', 'Reverse',
    'Potential contaminant', 'Contaminant')


#==============================================================================
#
#                         guess_maxquant_quantity
#
#==============================================================================

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


#==============================================================================
#
#                       dequantify
#                       demultiplex
#
#==============================================================================

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


is_multiplexed <- function(x){
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
        y %<>% `&`(all(as.numeric(channel) < n_open))
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
    if (!is_multiplexed(x)) return(x)
    y <- unname(vapply(x, .demultiplex, character(1)))
    make.unique(y, sep = '-')
}


#==============================================================================
#
#                       simplify_proteingroups
#                            read_uniprot_fasta
#                                extract_from_name
#                                extract_from_annot
#                                rm_from_annot
#
#=============================================================================

#' Load uniprot annotations from protein sequence fastafile
#'
#' @param fastafile    path to fasta file
#' @param fastafields  character vector
#' @param verbose      bool
#' @examples
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' read_uniprot_fasta(fastafile)
#' @return data.table(uniprot, genename, proteinname, reviewed, existence)
#' @note Existence values are always those of the canonical isoform
#'       (no isoform-level resolution for this field)
#' @export
read_uniprot_fasta <- function(fastafile, verbose = TRUE, use_cache = TRUE){
# Assert
    if (!requireNamespace('seqinr', quietly = TRUE)){
        stop("BiocManager::install('seqinr'). Then re-run.") }
    assert_all_are_existing_files(fastafile)
    Canonical <- Reviewed <- `Entry names` <- Version <- Existence <- `Gene names` <- NULL
    Orgid <- Orgname <- `Protein names` <- Uniprot <- Annotation <- NULL
# Return fasta tsv if available
    dir <- R_user_dir('autonomics', 'cache')
    tsvfile <- basename(fastafile)
    tsvfile %<>% stri_replace_first_fixed('.fasta', '.tsv')
    tsvfile %<>% file.path(dir, .)
    if (use_cache & file.exists(tsvfile)){
        if (verbose) message('\t\tRead ', tsvfile)
        return(fread(tsvfile))
    }
# Read fasta
    if (verbose) message('\tRead ', fastafile)
    fasta <- seqinr::read.fasta(fastafile)
    all_accessions <- extract_from_name(names(fasta), 2)
    fastadt <- data.table(Uniprot  = extract_from_name(names(fasta), 2),
                          Annotation = unname(vapply(fasta, attr, character(1), 'Annot')))
    fastadt[, Canonical := split_extract_fixed(Uniprot, '-', 1)]
    idx <- which(stri_detect_fixed(fastadt$Uniprot, '-'))
    fastadt[ idx, Isoform  := split_extract_fixed(Uniprot, '-', 2)]
    fastadt[!idx, Isoform  := 1]
# Extract components
    # Reviewed: trembl/swissprot
        message('\t\t\tExtract Reviewed: 0=trembl, 1=swissprot')
        fastadt[, Reviewed := as.numeric(extract_from_name(names(fasta),1)=='sp')]
    # `Entry names`
        message('\t\t\tExtract `Entry names`')
        fastadt [, `Entry names` := extract_from_name(names(fasta), 3)]
        fastadt [, `Entry names` := split_extract_fixed(`Entry names`, '_', 1)]
    # Sequence version
        # message('\t\t\tExtract (sequence) Version')
        pattern <- ' SV=[0-9]'
        #fastadt[, Version    := as.numeric(extract_from_annot(Annotation, pattern,5))]
        fastadt [, Annotation := rm_from_annot(Annotation, pattern)]
    # Existence
        message('\t\t\tExtract Existence: 1=protein, 2=transcript, ',
                '3=homology, 4=prediction, 5=uncertain, NA=isoform')
        pattern <- ' PE=[0-9]'
        fastadt [, Existence  := as.numeric(extract_from_annot(Annotation, pattern, 5))]
        fastadt [, Annotation := rm_from_annot(Annotation, pattern) ]
        fastadt [, Existence  := Existence[Uniprot == Canonical], by = 'Canonical']
    # Genes names (canonical only)
        message('\t\t\tExtract `Gene names`')
        pattern <- ' GN=.+$'
        fastadt [,`Gene names` := extract_from_annot(Annotation, pattern, 5)]
        fastadt [, Annotation   := rm_from_annot(Annotation, pattern)]
    # Orgid (canonical only)
        # message('\t\t\tExtract Orgid')
        pattern <- ' OX=[0-9]+'
        #fastadt[, Orgid        := extract_from_annot(Annotation, pattern, 5)]
        fastadt [, Annotation   := rm_from_annot(Annotation, pattern)]
    # Orgname (canonical only)
        # message('\t\t\tExtract Orgname')
        pattern <- ' OS=.+$'
        #fastadt[, Orgname      := extract_from_annot(Annotation, pattern, 5)]
        fastadt [, Annotation   := rm_from_annot(Annotation, pattern)]
    # Protein names (canonical only)
        message('\t\t\tExtract `Protein names`')
        pattern <- ' .+$'
        fastadt [,`Protein names` := extract_from_annot(Annotation, pattern, 2)]
        fastadt [, Annotation     := rm_from_annot(Annotation, pattern)]
        fastadt [,`Protein names` := `Protein names`[Uniprot == Canonical], by = 'Canonical']
        fastadt[, c('Annotation', 'Canonical') := NULL]
# Write / Return
    message('\t\t', tsvfile)
    fwrite(fastadt, tsvfile, sep = '\t')
    return(fastadt)
    
}


extract_from_name <- function(fastanames, i){
    fastanames                          %>%
    stri_split_fixed('|')               %>%
    vapply(extract, character(1), i)
}


extract_from_annot <- function(annotation, pattern, i){
    . <- NULL
    stri_extract_last_regex(annotation, pattern) %>%
    substr(i, nchar(.))
}


rm_from_annot <- function(annotation, pattern){
    stri_replace_last_regex(annotation, pattern, '')
}


#============================================================================
#
#                        simplify_proteingroups()
#
#============================================================================


#' Simplify proteingroups
#'
#' @param object SummarizedExperiment
#' @param fastafile string
#' @param verbose TRUE (default) or FALSE
#' @return data.table
#' @examples
#' require(magrittr)
#' file <- download_data('billing19.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' fdata(object)[1:5, ]
#' object %<>% simplify_proteingroups(fastafile)
#' fdata(object)[1:5, ]
#' @noRd
simplify_proteingroups <- function(
    object, fastafile, verbose = TRUE
){
# Return if NULL
    if (is.null(fastafile)) return(object)
# Uncollapse fdata annotations
    fasta_dt <- read_uniprot_fasta(fastafile)
    feature_dt  <-  fdata(object)[, c('feature_id', 'uniprot')] %>%
                    separate_rows('uniprot', sep=';') %>%
                    data.table()
# Merge in fasta annotations
    feature_dt %<>% merge(fasta_dt,  by.x = 'uniprot', by.y='Uniprot',
                        sort=FALSE, all.x=TRUE)
# Simplify
    if (verbose) message('\t\tSimplify proteingroups')
    feature_dt %<>% prefer_best_existence()
    feature_dt %<>% prefer_swissprot_over_trembl()
    feature_dt %<>% drop_fragments()
    feature_dt %<>% collapse_isoforms_paralogs()
# Merge into sumexp
    fdata(object)$uniprot <- fdata(object)$feature_name <- NULL
    fdata(object)$`Protein names` <- NULL
    fdata(object) %<>% merge(
        feature_dt, by = 'feature_id', sort = FALSE, all.x = TRUE)
    fdata(object) %<>% pull_columns(c('feature_id', 'feature_name',
                        'uniprot', 'canonical', 'Protein names'))
    object
}


prefer_best_existence <- function(feature_dt, verbose=TRUE){
    Existence <- NULL
    feature_dt[is.na(Existence), Existence:=5]
    feature_dt %<>% extract(, .SD[Existence == min(Existence)], 
                            by = 'feature_id')
    feature_dt[, Existence := NULL]
    if (verbose) message('\t\t\tDrop inferior existences')
    feature_dt
}

prefer_swissprot_over_trembl <- function(feature_dt, verbose=TRUE){
    Reviewed <- NULL
    feature_dt[is.na(Reviewed), Reviewed:=0]
    feature_dt %<>% extract(, .SD[Reviewed == max(Reviewed)], by = 'feature_id')
    if (verbose) message(
                '\t\t\tDrop trembl when swissprot available')
    feature_dt[, Reviewed := NULL]
    feature_dt
}

drop_fragments <- function(feature_dt, verbose=TRUE){
    `Protein names` <- IS.FRAGMENT <- NULL
    feature_dt[is.na(`Protein names`), `Protein names`:='']
    feature_dt[, IS.FRAGMENT := 
                as.numeric(stri_detect_fixed( `Protein names`, '(Fragment)'))]
    feature_dt %<>% extract(, .SD[IS.FRAGMENT == min(IS.FRAGMENT)], 
                            by = c('feature_id', 'Gene names'))
    feature_dt[, IS.FRAGMENT     := NULL]
    if (verbose) message(
            '\t\t\tDrop fragments when full seqs available')
    feature_dt
}

collapse_isoforms_paralogs <- function(feature_dt, verbose=TRUE){
    if (nrow(feature_dt)==0) return(feature_dt)
    `Gene names` <- uniprot <- Canonical <- `Protein names`   <- NULL
    
    groupby <- 'feature_id'
    feature_dt[is.na(`Gene names`), `Gene names` := '']
    feature_dt[, uniprot  := paste0(unique(uniprot),  collapse=';'), by=groupby]
    feature_dt[, Canonical:= paste0(unique(Canonical),collapse=';'), by=groupby]
    feature_dt[, `Gene names`   := paste0(unique(`Gene names`),    collapse=';'), by=groupby]
    feature_dt[,`Protein names` :=
                        commonify_strings(unique(`Protein names`)), by=groupby]          
    feature_dt %<>% unique()
    setnames(feature_dt, 'Gene names',     'feature_name')
    if (verbose) message('\t\t\tCollapse isoforms and paralogs')
    feature_dt
}


#' Commonify strings
#' @param x character vector
#' @examples
#' # NO DIFFERENCES
#'    x <- c( 'Retrotransposon Gag-like protein 8B',
#'            'Retrotransposon Gag-like protein 8B')
#'    commonify_strings(x)
#' # TAILS DIFFER
#'    x <- c( 'Histone H2B type 1-K',
#'            'Histone H2B type 1-C/E/F/G/I')
#'    commonify_strings(x)
#'    x <- c("Small nuclear ribonucleoprotein-associated proteins B and B'",
#'           "Small nuclear ribonucleoprotein-associated protein N")
#'    commonify_strings(x)
#' # MORE COMPLEX DIFFERENCES
#'    x <- c( 'Fatty acid binding protein, isoform 3',
#'            'Fatty acid binding protein',
#'            'heart-specific Fatty acid binding protein',
#'            'heart-specific Fatty acid binding protein, isoform 3')
#'    commonify_strings(x)
#' # NOTHING IN COMMON
#'    x <- c('ABC1', 'DEF2')
#'    commonify_strings(x)
#' @noRd
commonify_strings <- function(x){
    . <- NULL
    common <- Reduce(extract_common_substr, x)
    alternate  <- if (common==''){  x
                } else {            stri_replace_first_fixed(x, common, '') %>%
                                    stri_replace_first_fixed(', ', '')      %>%
                                    trimws()
                }
    if (all(alternate == '')) return(common)

    alternate                          %>%
    unique()                           %>%
    (function(s){s[s==''] <- '.'; s})  %>%
    sort()                             %>%
    #magrittr::extract(.!='')          %>%
    paste0(collapse=' | ')             %>%
    paste0('( ', ., ' )')              %>%
    paste0(common, ' ', .)
}


#' Extract common substring
#' @param a first string
#' @param b second string
#' @return  string
#' @examples
#' # Sequences
#'   a <- "heart-specific Fatty acid binding protein"
#'   b <- "Fatty acid binding protein, isoform 3"
#'   extract_common_substr(a, b)
#'
#'   a <- "Small nuclear ribonucleoprotein-associated proteins B and B'"
#'   b <- "Small nuclear ribonucleoprotein-associated protein N"
#'   extract_common_substr(a, b)
#' @references https://stackoverflow.com/questions/28261825
#' @noRd
extract_common_substr <- function(a, b){

    tt <- drop(attr(adist(a, b, counts=TRUE), "trafos"))

    # Nothing in common
    if (!stri_detect_regex(tt, 'M+')) return('')

    # Something in common
    aa  <-  stri_sub(tt, stri_locate_all_regex(tt, '[DM]+')[[1]]) %>%
            paste0(collapse = '') %>% trimws()
            # paste is required because multiple substrings can be found
    #bb <- tt %>%
    # stri_sub(stri_locate_all_regex(tt, '[IM]+')[[1]]) %>%
    # paste0(collapse = '')

    stri_sub(a, stri_locate_all_regex(aa, 'M+')[[1]]) %>%
    paste0(collapse = '') %>%
    trimws()

    # different  = c(a %>%
    #    stri_sub(stri_locate_all_regex(aa, 'D+')[[1]]) %>%
    #    trimws(),
    # b %>%
    #    stri_sub(stri_locate_all_regex(bb, 'I+')[[1]])) %>%
    #    trimws())

}

#==============================================================================
#
#                            invert
#
#==============================================================================

#' Invert
#'
#' For character vectors: invert collapsed strings.
#' For SummarizedExperiments: invert expressions , subgroups, and sample ids
#'
#' @param  x          character vector or SummarizedExperiment
#' @param  sep        string: collapsed string separator
#' @param  subgroups  character vector: subgroup levels to be inversed
#' @param  ... to enable S3 method dispatch
#' @return character vector or SummarizedExperiment
#' @examples
#' # character
#'     x <- c('Ctrl_A', 'Ctrl_B')
#'     invert(x)
#'
#' # SummarizedExperiment
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(file, plot=FALSE)
#'     invert(object)
#' @export
invert <- function(x, ...)    UseMethod('invert', x)


#' @rdname invert
#' @export
invert.character <- function(x, sep = guess_sep(x), ... ){
    stri_split_fixed(x, sep)                      %>%
    lapply(rev)                                   %>%
    vapply(paste, character(1), collapse = sep)
}


#' @rdname invert
#' @export
invert.SummarizedExperiment <- function(
    x, subgroups = slevels(x, 'subgroup'), sep = guess_sep(x, 'subgroup'), ...
){
# Assert
    if (length(subgroups)==0) return(x)
    assert_is_all_of(x, 'SummarizedExperiment')
    assert_is_subset('subgroup', svars(x))
    assert_is_subset(subgroups, subgroup_levels(x))
# Initialize message
    idx <- which(x$subgroup %in% subgroups)
    first <- values(x)[, idx[1]] %>% (function(y) which(!is.na(y))[[1]])
    oldvalue <- values(x)[first, idx[1]] %>% round(2) %>% as.character()
    message('\t\tInvert subgroups ', paste0(subgroups, collapse = ', '))

# Invert (log) ratios
    if (all(values(x)>0, na.rm=TRUE)){ values(x)[, idx] %<>% (function(x){1/x})
    } else {                           values(x)[, idx] %<>% (function(x){ -x})}
    newvalue <- as.character(round(values(x)[first, idx[1]], 2))
    message('\t\t\texprs    : ', as.character(oldvalue), ' -> ', 
            as.character(newvalue))
# Invert subgroup and sampleid values
    oldsubgroups <- sdata(x)$subgroup[idx]
    newsubgroups <- vapply(oldsubgroups, invert.character, character(1),sep=sep)
    oldsampleids <- sdata(x)$sample_id[idx]
    for (i in seq_along(idx)){
        sdata(x)$subgroup[ idx[i]] <- newsubgroups[i]
        sdata(x)$sample_id[idx[i]] %<>% stri_replace_first_fixed(
                                            oldsubgroups[i], newsubgroups[i])
        snames(x)[         idx[i]] %<>% stri_replace_first_fixed(
                                            oldsubgroups[i], newsubgroups[i])
    }
    newsampleids <- sdata(x)$sample_id[idx]
    message('\t\t\tsubgroups: ', oldsubgroups[1], ' -> ', newsubgroups[1])
    message('\t\t\tsampleids: ', oldsampleids[1], ' -> ', newsampleids[1])

# Order on subgroup
    x %<>% arrange_samples_('subgroup')
    message('\t\tOrder on subgroup')
    x$subgroup %<>% droplevels()

# Return
    return(x)
}

arrange_samples_ <- function(x, svars){
    idx <- do.call(order, sdata(x)[, svars, drop = FALSE])
    x %<>% extract(, idx)
    return(x)
}

#==============================================================================
#
#                       create_sfile
#
#==============================================================================


#' Create sfile
#' @param object      SummarizedExperiment
#' @param sfile  sample file
#' @param verbose     TRUE/FALSE
#' @return sample file path
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' create_sfile(object, paste0(tempfile(), '.tsv'))
#' @export
create_sfile <- function(object, sfile, verbose = TRUE){
    if (verbose) message('\tCreate sfile: ', sfile)
    assert_all_are_dirs(dirname(sfile))
    fwrite(sdata(object), sfile, sep = '\t', row.names = FALSE)
    return(sfile)
}


#==============================================================================
#
#               filter_maxquant_features
#
#==============================================================================


rm_reverse <- function(object, verbose){
    Reverse <- NULL
    fdata(object)$Reverse %<>% na_to_string()
    object %<>% filter_features(Reverse != '+', verbose=verbose)
    fdata(object)$Reverse <- NULL
    object
}


#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <-  read_proteingroups(file, contaminants=TRUE, pca=FALSE,
#'                             fit='limma', plot=FALSE) # read + analyze
#' rm_contaminants(object, verbose=TRUE)
#' @noRd
rm_contaminants <- function(object, verbose){
    contaminant_var <- c('Potential contaminant', 'Contaminant')
    contaminant_var %<>% intersect(fvars(object))
    fdata(object)[[contaminant_var]] %<>% na_to_string()

    idx <- fdata(object)[[contaminant_var]]==''
    if (verbose) message('\t\tRetain ',
        sum(idx), '/', length(idx), " features: contaminant != '+'")
    object %<>% extract_features(idx)
    object

}


rm_unlocalized <- function(object, min_localization_prob, verbose){
    `Localization prob` <- NULL
    if (!'Localization prob' %in% fvars(object)) return(object)
    assert_all_are_in_range(min_localization_prob, 0, 1)
    object %<>% filter_features(`Localization prob` >= min_localization_prob,
                                verbose = verbose)
    object
}

filter_maxquant_features <- function(
    object, reverse, contaminants, min_localization_prob = 0.75,
    verbose
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_a_bool(reverse)
    assert_is_a_bool(contaminants)
    assert_all_are_in_range(min_localization_prob, 0, 1)
    assert_is_a_bool(verbose)
# Filter
    if (verbose) message('\tFilter features')
    if (!reverse)      object %<>% rm_reverse(     verbose = verbose)
    if (!contaminants) object %<>% rm_contaminants(verbose = verbose)
    object %<>% rm_missing_in_all_samples(verbose = verbose)
    object %<>% filter_exprs_replicated_in_some_subgroup(verbose = verbose)
    object %<>% rm_unlocalized(min_localization_prob, verbose = verbose)
# Return
    object
}

rename_proteingroup_fvars <- function(object){
    stri_rep <- stri_replace_first_fixed
    names(fdata(object)) %<>% stri_rep('Gene names', 'feature_name')
    names(fdata(object)) %<>% stri_rep('Majority protein IDs', 'uniprot')
    fdata(object) %<>% pull_columns(c('feature_id', 'feature_name', 'uniprot'))
    object
}

rename_phospho_fvars <- function(object){
    stri_rep <- stri_replace_first_fixed
    names(fdata(object)) %<>% stri_rep('Gene names', 'feature_name')
    names(fdata(object)) %<>% stri_rep('Proteins', 'uniprot')
    names(fdata(object)) %<>% stri_rep('Positions within proteins', 'position')
    fdata(object) %<>% pull_columns(c('feature_id', 'feature_name', 'uniprot'))
    object
}




#==============================================================================
#
#                     filter_maxquant samples
#
#==============================================================================

filter_maxquant_samples <- function(object, select_subgroups, verbose){
    object %<>% filter_samples_available_for_some_feature(verbose = verbose)
    if (!is.null(select_subgroups))  object %<>%
        filter_samples(subgroup %in% select_subgroups, verbose = verbose)
    object
}


#==============================================================================
#
#                     transform_maxquant
#
#==============================================================================

transform_maxquant <- function(object, impute, verbose, plot){
    if (verbose && impute) message('\tTransform exprs')
# Impute
    if (impute) object %<>% impute_systematic_nondetects(plot = FALSE)
    object
}

#==============================================================================
#
#                  subtract_proteingroups
#
#==============================================================================

#' Which are the phospho expr columns?
#' @param x phosphosites colnames
#' @examples
#' phosphosites <- download_data('billing19.phosphosites.txt')
#' x <- names(data.table::fread(phosphosites))
#' quantity <- guess_maxquant_quantity(phosphosites)
#' phospho_expr_columns(x, quantity)
#' @noRd
phospho_expr_columns <- function(x, quantity){
    pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
    which(stri_detect_regex(x,pattern) & stri_detect_regex(x, '___1'))
}


subtract_proteingroups <- function(phosphosites, proteingroups, verbose){
# Initialize
    . <- phospho <- protein <- occupancy <- NULL
    `Protein group IDs` <- NULL
# Report
    if (verbose) message(
        '\tAdd occupancies(phospho) = values(phospho) - values(proteins)')
# phospho datatable
    cols <- c('feature_id', 'Protein group IDs')
    fosdt <- sumexp_to_wide_dt(phosphosites, fvars = 'Protein group IDs')
    fosdt %<>% separate_rows(`Protein group IDs`, sep = ';' )
    fosdt %<>% data.table()
    fosdt %<>% melt.data.table(id.vars = c('feature_id', 'Protein group IDs'))
    setnames(fosdt,
        c('feature_id', 'Protein group IDs', 'variable',  'value'),
        c('phospho_id', 'protein_id',        'sample_id', 'phospho'))
# proteingroups datatable
    protein_dt <- sumexp_to_long_dt(proteingroups)
    setnames(
        protein_dt,c('feature_id', 'value'), c('protein_id', 'protein'))
# merge
    fosdt %<>% merge(  protein_dt,
                            by = c('protein_id', 'sample_id'), all.x = TRUE)
    fosdt %<>% extract(,
                        .(  phospho  = unique(phospho),
                            protein = median(protein, na.rm = TRUE),
                            subgroup      = unique(subgroup)),
                        by = c('phospho_id', 'sample_id') )
# compute occupancy
    fosdt[, occupancy := phospho - protein]
    occ <- dt2mat(dcast(fosdt, phospho_id ~ sample_id, value.var = 'occupancy'))
    pgs <- dt2mat(dcast(fosdt, phospho_id ~ sample_id, value.var = 'protein'))
    assert_are_identical(rownames(phosphosites), rownames(occ))
    assert_are_identical(rownames(phosphosites), rownames(pgs))
    assert_are_identical(colnames(phosphosites), colnames(occ))
    assert_are_identical(colnames(phosphosites), colnames(pgs))
    proteingroups(phosphosites) <- pgs
    occupancies(phosphosites)   <- occ
    phosphosites
}

#==============================================================================
#
#               extract peptide count information
#
#==============================================================================
#' maxquant peptide count patterns
#' @examples
#' MAXQUANT_PATTERNS_PEPCOUNTS
#' @export
MAXQUANT_PATTERNS_PEPCOUNTS <- c(
    uniquepeps         = '^Unique peptides (.+)$',
    razoranduniquepeps = '^Razor \\+ unique peptides (.+)$',
    peps               = '^Peptides (.+)$')

#==============================================================================
#
#                  .read_maxquant
#                       add_pepcounts
#                   read_proteingroups
#                   read_phosphosites
#                       phospho_expr_columns
#                       subtract_proteingroups
#
#==============================================================================

add_pepcounts <- function(object, file, pepcountpattern, quantity){
    . <- NULL
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_all_are_existing_files(file)
    assert_is_a_string(pepcountpattern)
    assert_has_names(pepcountpattern)
    assert_is_subset(pepcountpattern, MAXQUANT_PATTERNS_PEPCOUNTS)
    assert_is_subset(quantity,  names(MAXQUANT_PATTERNS_QUANTITY))
    
    names1 <- names(fread(file, integer64 = 'numeric', nrows=1))
    pepcselect <- names1 %>% extract(stri_detect_regex(., pepcountpattern))
    pepcs <- as.matrix(fread(file, select=pepcselect, integer64='numeric'))
    rownames(pepcs) <- rownames(object)
    colnames(pepcs) %<>% stri_replace_first_regex(pepcountpattern, "$1")
    if (quantity == 'Intensity'){
        exprsselect <- stri_replace_first_regex(
            colnames(object), MAXQUANT_PATTERNS_QUANTITY[quantity], "$1")
    } else {
        exprsselect <- stri_replace_first_regex(
            colnames(object), MAXQUANT_PATTERNS_QUANTITY[quantity], "$2")
    }
    pepcs %<>% extract(, exprsselect)
    colnames(pepcs) <- colnames(object)
    assays(object)$pepcounts <- pepcs
    assayNames(object) %<>% stri_replace_first_fixed(
                            'pepcounts', names(pepcountpattern))
    object
}


#' @rdname read_proteingroups
#' @export
.read_maxquant <- function(file, quantity = guess_maxquant_quantity(file),
    sfile = NULL, by.x = 'sample_id', by.y = 'sample_id', 
    subgroupvar = 'subgroup', select_subgroups = NULL, 
    invert_subgroups = character(0),
    pepcountpattern = MAXQUANT_PATTERNS_PEPCOUNTS[1], verbose = TRUE){
# Read
    . <- NULL
    assert_all_are_existing_files(file)
    assert_is_subset(quantity,  names(MAXQUANT_PATTERNS_QUANTITY))
    assert_is_a_bool(verbose)
    phospho <- stri_detect_fixed(basename(file),'phospho',case_insensitive=TRUE)
    if (verbose)  message('\tRead ', file)
    names1 <- names(fread(file, integer64 = 'numeric', nrows=1))
    fids1  <- fread(file, select = 'id', colClasses = 'character')[[1]]
    FVARS <- if (phospho) PHOSPHOSITE_FVARS else PROTEINGROUP_FVARS
    fdata1 <- fread(file, select = intersect(FVARS, names1))
    fdata1$id %<>% as.character()
    setnames(fdata1, 'id', 'feature_id')
    if ('Gene names' %in% names(fdata1)) setnames(fdata1, 'Gene names',
                                                        'feature_name')
# Exprs/pepcounts
    select <- names1 %>%
            extract(stri_detect_regex(.,MAXQUANT_PATTERNS_QUANTITY[[quantity]]))
    if (phospho)  select %<>% extract(stri_endswith_fixed(., '___1'))
    exprs1 <- fread(file, select = select, integer64 = 'numeric') %>%
        un_int64() %>%
        as.matrix()
    if (is.character(exprs1)){
        message("File contains strings rather than numbers - convert")
        storage.mode(exprs1) <- 'numeric'}
    assert_all_are_non_missing_nor_empty_character(fids1)
    rownames(exprs1) <- fids1
    object <- SummarizedExperiment(list(exprs = exprs1), rowData=fdata1)
    if (!phospho)  object %<>% add_pepcounts(file, pepcountpattern, quantity)
# Process
    if (phospho)  colnames(object) %<>% stri_replace_last_fixed('___1', '')
    colnames(object) %<>% dequantify(quantity)
    colnames(object) %<>% demultiplex()
    object$sample_id <- colnames(object)
    object %<>% merge_sfile(sfile = sfile, by.x = by.x, by.y = by.y)
    object %<>% add_subgroup(subgroupvar=subgroupvar, verbose=verbose)
    object %<>% filter_maxquant_samples(
                    select_subgroups = select_subgroups, verbose)
    values(object) %<>% zero_to_na(verbose = verbose)
    values(object) %<>% nan_to_na( verbose = verbose)
    object %<>% log2transform(verbose=verbose)
    object %<>% invert(invert_subgroups)
# Return
    assayNames(object)[1] <- 'maxquant'
    metadata(object)$quantity <- quantity
    metadata(object)$file <- file
    object
}

utils::globalVariables('where')
un_int64 <- function(x) {
    dplyr::mutate(x, dplyr::across(where(bit64::is.integer64), as.numeric))
}

#' Read/Analyze proteingroups/phosphosites
#'
#' @param file         proteingroups/phosphosites file
#' @param proteinfile  proteingroups file
#' @param fastafile    NULL or fastafile (to deconvolute proteingroups)
#' @param quantity     string: "Ratio normalized",
#'                             "Ratio",
#'                             "LFQ intensity",
#'                             "Reporter intensity corrected",
#'                             "Reporter intensity",
#'                             "Intensity labeled",
#'                             "Intensity"
#' @param sfile         sample file
#' @param by.x          `file`  column to merge sdata: string
#' @param by.y          `sfile` column to merge sdata: string
#' @param subgroupvar   subgroup svar
#' @param select_subgroups  subgroups to be selected (character vector)
#' @param invert_subgroups  subgroups to be inverted (character vector)
#' @param contaminants  whether to return contaminants
#' @param reverse       whether to return reverse peptides
#' @param min_localization_prob min site localization probability (number)
#' @param impute        whether to impute consistent nondetects
#' @param pepcountpattern value in MAXQUANT_PATTERNS_PEPCOUNTS
#' @param formula       desgnmat formula
#' @param block         block svar
#' @param coefs        NULL or character vector: model coefficients to test
#' @param contrastdefs NULL or character vector: coefficient contrasts to test
#' @param pca           whether to pca
#' @param fit           fit model: NULL, 'limma', 'lm', 'lme', 'lmer','wilcoxon'
#' @param verbose       whether to message
#' @param plot          whether to plot
#' @return SummarizedExperiment
#' @examples
#' # Fukuda 2020
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(file, pca=TRUE, fit='limma')
#' # Billing 2019
#'     profile <- download_data('billing19.proteingroups.txt')
#'     subgroups <- c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00')
#'     subgroups %<>% paste0('_STD')
#'     palette <- make_colors(c(subgroups[1], 'E00.8_STD', subgroups[-1]))
#'     pro <- read_proteingroups(
#'         profile, select_subgroups = subgroups, palette = palette, 
#'         coefs = 'subgroupE05_STD', sample_id = 'E05_STD.R2')
#' @export
read_proteingroups <- function(
    file, quantity = guess_maxquant_quantity(file), sfile = NULL,
    by.x = 'sample_id', by.y = NULL, select_subgroups = NULL, 
    contaminants = FALSE,
    reverse = FALSE, fastafile = NULL, invert_subgroups = character(0),
    impute = stri_detect_regex(quantity, "[Ii]ntensity"),
    pepcountpattern = MAXQUANT_PATTERNS_PEPCOUNTS[1], subgroupvar = NULL,
    formula = NULL, block = NULL, coefs = NULL, contrastdefs = NULL,
    pca = TRUE, fit = 'limma', verbose = TRUE, plot = pca & !is.null(fit), 
    feature_id = NULL, sample_id = NULL, palette = NULL
){
# Assert
    . <- NULL
    assert_all_are_existing_files(file)
    if (!is.null(fastafile)) assert_all_are_existing_files(fastafile)
# Read
    object <- .read_maxquant(file, quantity,
        sfile = sfile, by.x = by.x, by.y = by.y, 
        select_subgroups = select_subgroups,
        invert_subgroups = invert_subgroups,
        pepcountpattern = pepcountpattern, verbose = verbose)
    assayNames(object)[1] %<>% gsub('maxquant', 'proteingroups', .)
# Preprocess
    object %<>% filter_maxquant_features(reverse = reverse,
                    contaminants = contaminants, verbose = verbose)
    object %<>% rename_proteingroup_fvars()
    object %<>% simplify_proteingroups(fastafile)
    object %<>% transform_maxquant(impute=impute, verbose=verbose, plot=plot)
# Analyze
    object %<>% analyze(
        pca = pca, fit = fit, subgroupvar = subgroupvar, 
        formula = formula, block = block, 
        coefs = coefs, contrastdefs = contrastdefs, 
        verbose = verbose, plot = plot, 
        feature_id = feature_id, sample_id = sample_id, palette = palette)
# Return
    object
}

#' @rdname read_proteingroups
#' @export
read_phosphosites <- function(
    file, proteinfile = paste0(dirname(file), '/proteinGroups.txt'),
    quantity = guess_maxquant_quantity(file),
    sfile = NULL, by.x = 'sample_id', by.y = 'sample_id', 
    select_subgroups = NULL, contaminants = FALSE,
    reverse = FALSE, min_localization_prob = 0.75, fastafile = NULL,
    invert_subgroups = character(0), pca = TRUE,
    fit = 'limma', subgroupvar = NULL, formula = NULL, block = NULL, 
    coefs = NULL, contrastdefs = NULL, 
    verbose = TRUE, plot = pca & !is.null(fit)
){
# Assert
    . <- NULL
    `Protein group IDs` <- `Localization prob` <- NULL
    assert_all_are_existing_files(c(file, proteinfile))
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS_QUANTITY))
# Read
    prot <- .read_maxquant(file=proteinfile, quantity = quantity,
        sfile = sfile, select_subgroups = select_subgroups,
        invert_subgroups = invert_subgroups, verbose=verbose)
    object  <- .read_maxquant(file = file, quantity = quantity,
        sfile = sfile, by.x = by.x, by.y = by.y,
        invert_subgroups = invert_subgroups, verbose = verbose)
    assayNames(prot)[1]   %<>% gsub('maxquant', 'proteingroups', ., fixed=TRUE)
    assayNames(object)[1] %<>% gsub('maxquant', 'phosphosites',  ., fixed=TRUE)
    object %<>% filter_maxquant_features(
                    reverse = reverse, contaminants = contaminants,
                    min_localization_prob = min_localization_prob,
                    verbose = verbose)
    object %<>% subtract_proteingroups(prot, verbose)
# Prepare
    object %<>% rename_phospho_fvars()
    object %<>% simplify_proteingroups(fastafile)
    object %<>% transform_maxquant(impute=FALSE, verbose=verbose, plot=plot)
# Analyze
    object %<>% analyze(pca=pca, fit=fit, subgroupvar=subgroupvar, 
                        formula=formula, block=block, 
                        coefs=coefs, contrastdefs=contrastdefs,
                        verbose=verbose, plot=plot)
# Return
    object
}


#==========================================================================
#
#   ANNOTATE_UNIPROT
#
#==========================================================================

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
        fdt %<>% tidyr::separate_rows(uniprot, sep=collapse)
        fdt %<>% data.table::as.data.table()
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



