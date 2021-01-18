#==============================================================================
#
#        guess_maxquant_quantity
#
#==============================================================================

#' maxquant patterns
#' @examples
#' MAXQUANT_PATTERNS
#' @export
MAXQUANT_PATTERNS <- c(
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
#' @return  string: value from names(MAXQUANT_PATTERNS)
#' @examples
#' # file
#'     x <- download_data('billing16.proteingroups.txt')
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
#'     file <- download_data( 'billing16.proteingroups.txt')
#'     x <- data.table::fread(file)
#'     guess_maxquant_quantity(x)
#'
#' # SummarizedExperiment
#'      file <-download_data( 'billing16.proteingroups.txt')
#'      # x <- .read_proteingroups(file)
#'      # guess_maxquant_quantity(x)
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
    for (quantity in names(MAXQUANT_PATTERNS)){
        pattern <- MAXQUANT_PATTERNS[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.data.frame <- function(x, ...){
    x <- names(x)
    for (quantity in names(MAXQUANT_PATTERNS)){
        pattern <- MAXQUANT_PATTERNS[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.SummarizedExperiment <- function(x, ...){
    x <- snames(x)
    for (quantity in names(MAXQUANT_PATTERNS)){
        pattern <- MAXQUANT_PATTERNS[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}

#==============================================================================
#
#                       standardize_maxquant_snames
#
#==============================================================================


#' Standardize maxquant snames
#'
#' Standardize maxquant sample names
#'
#' Drop "Ratio normalized", "LFQ intensity" etc from maxquant sample names
#'
#' @param x        character vector or SummarizedExperiment
#' @param quantity maxquant quantity
#' @param verbose  TRUE (default) or FALSE
#' @param ...      allow for proper S3 method dispatch
#' @return character vector or SummarizedExperiment
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
#'    file <- download_data('billing16.proteingroups.txt')
#'    # x <- .read_proteingroups(file)
#'    # standardize_maxquant_snames(x)
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
    # x = multiplexes + channels. Return multiplexes if single channel.
    pattern <- MAXQUANT_PATTERNS[[quantity]]

    # Decompose multiplexes and channels
    if (quantity == 'Intensity'){
        multiplexes     <- stri_replace_first_regex(x, pattern, '$1')
        channel <- rep('', length(multiplexes))
    } else {
        multiplexes     <- stri_replace_first_regex(x, pattern, '$2')
        channel <- stri_replace_first_regex(x, pattern, '$1')
    }

    # Standardize
    if (all(channel=='')){
        cleanx <- multiplexes
    } else {
        cleanx <- sprintf('%s{%s}', multiplexes, channel)
    }
    if (verbose) message('\t\tStandardize snames: ', x[1], '  ->  ', cleanx[1])
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

#==============================================================================
#
#                         demultiplex
#
#==============================================================================



#' Demultiplex snames
#'
#' Demultiplex sample names
#'
#' @param x        character vector or SummarizedExperiment
#' @param verbose  logical
#' @param ...      allow for S3 dispatch
#' @return character vector or SummarizedExperiment
#' @examples
#' # character vector
#'    # Alternate multiplexing forms supported
#'    demultiplex("STD(L)_EM00(M)_EM01(H)_R1{M/L}") # Label Ratio
#'    demultiplex('A(0)_B(1)_C(2)_D(3)_R1{0}'     ) # Reporter intensity
#'    demultiplex('STD(L)_EM00(M)_EM01(H)_R1{L}')   # Label Intensity
#'
#'    # Alternate separators supported
#'    demultiplex('STD(L)_EM00(M)_EM01(H)_R1{L}')   # underscore
#'    demultiplex('STD(L).EM00(M).EM01(H).R1{L}')   # dot
#'    demultiplex('STD(L)EM00(M)EM01(H).R1{L}')     # no separator
#'
#'    # Composite snames supported
#'    demultiplex("WT.t0(L)_WT.t1(M)_WT.t2(H)_R1{H/L}")
#'
#'    # Uniqueness ensured by appending labels when necessary
#'    demultiplex(c("STD(L).BM00(M).BM00(H).R10{M/L}",
#'                        "STD(L).BM00(M).BM00(H).R10{H/L}"))
#'
#'    # Uniplexed snames are returned unchanged
#'    demultiplex(c('STD_R1', 'EM0_R1'))
#'
#' # SummarizedExperiment
#'    require(magrittr)
#'    file <- download_data('billing16.proteingroups.txt')
#'    # x <- .read_proteingroups(file)
#'    # x %<>% standardize_maxquant_snames()
#'    # demultiplex(x, verbose = TRUE)
#' @export
demultiplex <- function (x, ...) {
    UseMethod("demultiplex", x)
}

is_multiplexed <- function(x){
    pattern <- '(.+)\\{(.+)\\}'
    n_open   <- stri_count_fixed(x, '(')
    n_closed <- stri_count_fixed(x, ')')
    all(stri_detect_regex(x, pattern) & (n_open==n_closed) & (n_open>0))
}

are_all_identical <- function(y){
    if (length(y)==1) TRUE else all(y[-1] == y[1])
}

#' Separate multiplexes from channels
#'
#'     STD(L)_EM00(M)_EM01(H)_R1{M/L}
#'     ------------------------- ...
#'           multiplex         channel
#'
#'
#' @param x character vector: e.g. "STD(L)_EM00(M)_EM01(H)_R1{M/L}"
#' @examples
#' x <- "STD(L)_EM00(M)_EM01(H)_R1{M/L}"
#' extract_multiplexes(x)     # STD(L)_EM00(M)_EM01(H)_R1
#' extract_channels(x) #                           M/L
#' @noRd
extract_multiplexes <- function(x){
    pattern <- '(.+)\\{(.+)\\}'
    stri_replace_first_regex(x, pattern, '$1')
}
extract_channels <- function(x){
    pattern <- '(.+)\\{(.+)\\}'
    stri_replace_first_regex(x, pattern, '$2')
}

#' Separate samples from labels
#' @param multiplexes character vector: e.g. "STD(L)_EM00(M)_EM01(H)_R1"
#' @examples
#' x <- "STD(L)_EM00(M)_EM01(H)_R1"
#' extract_labels(x)
#' extract_samples(x)
#' @noRd
extract_labels <- function(multiplexes){
    pattern <- '\\(.+?\\)'
    stri_extract_all_regex(multiplexes, pattern)       %>%
    lapply(stri_replace_first_fixed, '(', '')  %>%
    lapply(stri_replace_first_fixed, ')', '')
}
extract_biosamples <- function(multiplexes, labels){
    pattern <- '\\(.+?\\)'
    stri_split_regex(multiplexes, pattern) %>%
    lapply(function(y){
            y[seq_len(length(labels[[1]]))] %<>%
            stri_replace_first_regex('^[_. ]', ''); y})
    # rm sep from samples (but not from replicate - needed to glue back later!)
}


#' Are number of samples/labels equally multiplexed across mixes?
#'
#' @param samples list of char vectors (see examples)
#' @param labels  list of char vectors (see examples)
#' @return TRUE or FALSE
#' @examples
#' samples <- list(c("STD", "EM00", "EM01", "_R1" ),
#'                 c("STD", "EM00", "EM01", "_R2"))
#' labels  <- list(c("L", "M", "H"),
#'                 c("L", "M", "H"))
#' are_equally_multiplexed(samples, labels)
#' @noRd
are_equally_multiplexed <- function(samples, labels){
    n_samples <- vapply(samples,  length, integer(1))
    n_labels  <- vapply(labels, length, integer(1))
    are_all_identical(n_samples) & are_all_identical(n_labels)
}


#' Extract/Drop replicate
#' @param samples list of char vectors: list(c("STD", "EM00", "EM01", "_R1" ))
#' @param labels  list of char vectors: list(c("L", "M", "H"))
#' @return character vector: "_R1"
#' @examples
#' samples <- list(c("STD", "EM00", "EM01", "_R1" ))
#' labels  <- list(c("L", "M", "H"))
#' extract_replicates(samples, labels) # "_R1"
#' drop_replicates(samples, labels)
#' @noRd
extract_replicates <- function(multiplexes, biosamples, labels){

    pattern <- '\\(.+?\\)'
    n_biosamples <- unique(vapply(biosamples, length, integer(1)))
    n_labels  <- unique(vapply(labels,  length, integer(1)))
    if (n_biosamples > n_labels){
        stri_split_regex(multiplexes, pattern) %>%
        vapply((function(y) extract(y, length(y))), character(1))
    } else {
        rep('', length(biosamples))
    }
}
drop_replicates <- function(biosamples, labels){
    n_samples <- unique(vapply(biosamples, length, integer(1)))
    n_labels  <- unique(vapply(labels,  length, integer(1)))
    if (n_samples > n_labels){
        biosamples %<>% lapply(extract, seq_len(n_samples-1))
    }
    biosamples
}

#' First stri_split_fixed. Then extract component i
#' @param x character vector
#' @param split   string character on which to split
#' @param i number: which component to extract
#' @return character vector
#' @noRd
stri_split_fixed_extract <- function(x, split, i){
    stri_split_fixed(x, split) %>%
    vapply(extract, character(1), i)
}

#' Forge channel samples
#' @param channels    character vector
#' @param biosamples    character vector
#' @param replicate  character vector
#' @return character vector
#' @examples
#' channels <- c('M/L', 'H/L')
#' biosamples <- list(c(  L = 'STD', M = 'EM00', H = 'EM01'),
#'                 c(  L = 'STD', M = 'EM00', H = 'EM01'))
#' replicate <- c('_R1', '_R2')
#' extract_channel_samples(channels, biosamples, replicate)
#' @noRd
extract_channelsamples <- function(channels, biosamples, replicate){
    . <- NULL
    is_ratio <- all(stri_detect_fixed(channels, '/'))
    if (is_ratio){
        num_label <- stri_split_fixed_extract(channels, '/', 1)
        den_label <- stri_split_fixed_extract(channels, '/', 2)
        den_samples <- mapply(extract, biosamples, den_label)
        num_samples <- mapply(extract, biosamples, num_label)
        sprintf('%s_%s%s', num_samples, den_samples, replicate)
    } else {
        biosamples %<>% mapply(extract, ., channels)
        sprintf('%s%s', biosamples, replicate)
    }
}

uniquify_channelsamples <- function(channelsamples, channels, verbose){
    idx <- which(cduplicated(channelsamples))
    if (length(idx)>0){
        label_tags <- channels[idx] %>% stri_replace_first_fixed('/', '')
        if (verbose)   message('\t\tUniquify snames: ',
                                channelsamples[idx][1], ' -> ',
                                channelsamples[idx][1],
                                label_tags[1],
                                ' (for ',
                                length(idx),
                                '/',
                                length(channelsamples),
                                ' snames)')
        channelsamples[idx] %<>% paste0(label_tags)
    }
    channelsamples
}

#' @rdname demultiplex
#' @export
demultiplex.character <- function(x, verbose = FALSE, ...){

    # Return unchanged if not multiplexed
    . <- NULL
    if (!is_multiplexed(x)) return(x)

    # Separate channels from multiplexes.
    channels    <- extract_channels(x)
    multiplexes <- extract_multiplexes(x)

    # Separate labels from biosamples. Return if inconsistent no of components.
    labels      <- extract_labels(multiplexes)
    biosamples  <- extract_biosamples(multiplexes, labels)
    if (!are_equally_multiplexed(biosamples, labels)){
        message('\t\tDemultiplexing failed: no of biosamples or labels differs')
        return(x)
    }
    biosamples %<>% mapply(set_names, ., labels, SIMPLIFY = FALSE)
    replicates  <- extract_replicates(multiplexes, biosamples, labels)
    biosamples %<>% drop_replicates(labels)

    # Extract channelsamples from multiplexes
    channelsamples <- extract_channelsamples(channels, biosamples, replicates)
    if (verbose) message(
        '\t\tDemultiplex snames: ', x[1], '  ->  ', channelsamples[1])

    # Ensure uniqueness. Add labels if required.
    uniquify_channelsamples(channelsamples, channels, verbose)

}


#' @rdname demultiplex
#' @export
demultiplex.SummarizedExperiment <- function(x,verbose  = FALSE, ...){
    newsnames <- demultiplex(snames(x), verbose = verbose)
    snames(x) <- sdata(x)$sample_id <- newsnames
    x
}

#==============================================================================
#
#              load_uniprot_fasta
#
#=============================================================================

# DEFAULT_FASTAFIELDS <- c('GENES', 'PROTEIN-NAMES', 'REVIEWED',
#                     'EXISTENCE', 'ENTRYNAME', 'ORGID', 'ORGNAME', 'VERSION')
DEFAULT_FASTAFIELDS <- c('GENES', 'EXISTENCE', 'REVIEWED', 'PROTEIN-NAMES')


#' Load uniprot annotations from protein sequence fastafile
#'
#' @param fastafile    path to fasta file
#' @param fastafields  character vector
#' @examples
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' load_uniprot_fasta(fastafile)
#' @return data.table(uniprot, genename, proteinname, reviewed, existence)
#' @note EXISTENCE values are always those of the canonical isoform
#'       (no isoform-level resolution for this field)
#' @noRd
load_uniprot_fasta <- function(
    fastafile, fastafields = DEFAULT_FASTAFIELDS, verbose = TRUE){
# Assert, Read, Extract
    assert_all_are_existing_files(fastafile)
    CANONICAL <- REVIEWED <- ENTRYNAME <- VERSION <- EXISTENCE <- GENES <- NULL
    ORGID <- ORGNAME <- `PROTEIN-NAMES` <- UNIPROTKB <- annotation <- NULL
    if (verbose) message('\t\tLoad fasta file')
    fasta <- read.fasta(fastafile)
    all_accessions <- extract_from_name(names(fasta), 2)
    dt <- data.table(UNIPROTKB = extract_from_name(names(fasta), 2),
            annotation = unname(vapply(fasta, attr, character(1), 'Annot')))
    dt[, CANONICAL := stri_replace_last_regex(UNIPROTKB, '[-][0-9]+','')]
# Extract name components
    message('\t\t\tExtract REVIEWED: 0=trembl, 1=swissprot')
    dt[, REVIEWED := as.numeric(extract_from_name(names(fasta),1)=='sp')]
    message('\t\t\tExtract ENTRYNAME')
    dt [, ENTRYNAME := extract_from_name(names(fasta), 3)]
# Extract annotation components
message('\t\t\tExtract (sequence) VERSION')
    pattern <- ' SV=[0-9]'  # VERSION
    dt [, VERSION := as.numeric(extract_from_annot(annotation, pattern,5))]
    dt [, annotation   := rm_from_annot(annotation, pattern)]
message('\t\t\tExtract EXISTENCE: 1=protein, 2=transcript, ',
        '3=homology, 4=prediction, 5=uncertain, NA=isoform')
    pattern <- ' PE=[0-9]'  # EXISTENCE
    dt [, EXISTENCE  := as.numeric(extract_from_annot(annotation, pattern, 5))]
    dt [, annotation := rm_from_annot(annotation, pattern) ]
    dt [, EXISTENCE  := EXISTENCE[UNIPROTKB==CANONICAL], by = 'CANONICAL']
    # only canonical have this
message('\t\t\tExtract GENES')
    pattern <- ' GN=.+$'    # GENES
    dt [, GENES       := extract_from_annot(annotation, pattern, 5)]
    dt [, annotation  := rm_from_annot(annotation, pattern)]
message('\t\t\tExtract ORGID')
    pattern <- ' OX=[0-9]+' # ORGID
    dt [, ORGID        := extract_from_annot(annotation, pattern, 5)]
    dt [, annotation   := rm_from_annot(annotation, pattern)]
message('\t\t\tExtract ORGNAME')
    pattern <- ' OS=.+$'    # ORGNAME
    dt [, ORGNAME      := extract_from_annot(annotation, pattern, 5)]
    dt [, annotation   := rm_from_annot(annotation, pattern)]
message('\t\t\tExtract PROTEIN-NAMES')
    pattern <- ' .+$'       # PROTEIN-NAMES
    dt [, `PROTEIN-NAMES` := extract_from_annot(annotation, pattern, 2)]
    dt [,  annotation     := rm_from_annot(annotation, pattern)]
    dt [, `PROTEIN-NAMES` := `PROTEIN-NAMES`[UNIPROTKB==CANONICAL],
        by = 'CANONICAL'] # only canonical have this
# Order & Return
    dt[, c('UNIPROTKB', 'CANONICAL', fastafields), with = FALSE]
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
#             simplify_proteingroups()
#
#============================================================================


#' Simplify proteingroups
#'
#' @param object SummarizedExperiment
#' @param fastafile string
#' @param fastafields character vector
#' @param verbose TRUE (default) or FALSE
#' @return data.table
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' fdata(object)[1:5, ]
#' object %<>% simplify_proteingroups(fastafile)
#' fdata(object)[1:5, ]
#' @noRd
simplify_proteingroups <- function(
    object, fastafile, fastafields = DEFAULT_FASTAFIELDS, verbose = TRUE
){
# Return if NULL
    if (is.null(fastafile)) return(object)
# Uncollapse fdata annotations
    fasta_dt <- load_uniprot_fasta(fastafile, fastafields)
    feature_dt  <-  fdata(object)[, c('feature_id', 'uniprot')] %>%
                    separate_rows('uniprot', sep=';') %>%
                    data.table()
# Merge in fasta annotations
    feature_dt %<>% merge(fasta_dt, by.x = 'uniprot', by.y='UNIPROTKB',
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
    EXISTENCE <- NULL
    feature_dt[is.na(EXISTENCE), EXISTENCE:=5]
    feature_dt %<>% extract(, .SD[EXISTENCE == min(EXISTENCE)],
                            by = 'feature_id')
    feature_dt[, EXISTENCE := NULL]
    if (verbose) message('\t\t\tDrop inferior existences')
    feature_dt
}


prefer_swissprot_over_trembl <- function(feature_dt, verbose=TRUE){
    REVIEWED <- NULL
    feature_dt[is.na(REVIEWED), REVIEWED:=0]
    feature_dt %<>% extract(, .SD[REVIEWED == max(REVIEWED)], by = 'feature_id')
    if (verbose) message(
                '\t\t\tDrop trembl when swissprot available')
    feature_dt[, REVIEWED := NULL]
    feature_dt
}


drop_fragments <- function(feature_dt, verbose=TRUE){
    `PROTEIN-NAMES` <- IS.FRAGMENT <- NULL
    feature_dt[is.na(`PROTEIN-NAMES`), `PROTEIN-NAMES`:='']
    feature_dt[, IS.FRAGMENT :=
                as.numeric(stri_detect_fixed( `PROTEIN-NAMES`, '(Fragment)'))]
    feature_dt %<>% extract(, .SD[IS.FRAGMENT == min(IS.FRAGMENT)],
                            by = c('feature_id', 'GENES'))
    feature_dt[, IS.FRAGMENT     := NULL]
    if (verbose) message(
            '\t\t\tDrop fragments when full seqs available')
    feature_dt
}



collapse_isoforms_paralogs <- function(feature_dt, verbose=TRUE){
    if (nrow(feature_dt)==0) return(feature_dt)
    GENES <- uniprot <- CANONICAL <- GENES <- `PROTEIN-NAMES` <- NULL

    groupby <- c('feature_id')
    feature_dt[is.na(GENES), GENES := '']
    feature_dt[, uniprot  := paste0(unique(uniprot),  collapse=';'), by=groupby]
    feature_dt[, CANONICAL:= paste0(unique(CANONICAL),collapse=';'), by=groupby]
    feature_dt[, GENES    := paste0(unique(GENES),    collapse=';'), by=groupby]
    feature_dt[,`PROTEIN-NAMES` :=
                        commonify_strings(unique(`PROTEIN-NAMES`)), by=groupby]
    feature_dt %<>% unique()
    setnames(feature_dt, 'GENES',     'feature_name')
    setnames(feature_dt, 'CANONICAL', 'canonical')
    setnames(feature_dt, 'PROTEIN-NAMES', 'Protein names')
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
#' #----------
#'    x <- c('Ctrl_A', 'Ctrl_B')
#'    invert(x)
#'
#' # SummarizedExperiment
#' #---------------------
#'    file <- download_data('billing16.proteingroups.txt')
#'    x <- read_proteingroups(file)
#'    invert(x, subgroups = c('E_EM', 'E_BM', 'EM_BM'))
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
    first <- exprs(x)[, idx[1]] %>% (function(y) which(!is.na(y))[[1]])
    oldvalue <- exprs(x)[first, idx[1]] %>% round(2) %>% as.character()
    cmessage('\t\tInvert subgroups %s', paste0(subgroups, collapse = ', '))

# Invert (log) ratios
    if (all(exprs(x)>0, na.rm = TRUE)){ exprs(x)[, idx] %<>% (function(x){1/x})
    } else {                            exprs(x)[, idx] %<>% (function(x){ -x})}
    newvalue <- as.character(round(exprs(x)[first, idx[1]], 2))
    cmessage('\t\t\texprs    : %s -> %s',
            as.character(oldvalue), as.character(newvalue))
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
    cmessage('\t\t\tsubgroups: %s -> %s', oldsubgroups[1], newsubgroups[1])
    cmessage('\t\t\tsampleids: %s -> %s', oldsampleids[1], newsampleids[1])

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
#                       create_maxquant_samplefile
#
#==============================================================================


#' Create maxquant samplefile
#' @param file        proteingroups/phosphosite file
#' @param samplefile  sample file
#' @param verbose     TRUE/FALSE
#' @param quantity    maxquant quantity
#' @return samplefile
#' @export
create_samplefile <- function(object, samplefile, verbose = TRUE){
    if (verbose) message('\tCreate samplefile: ', samplefile)
    assert_all_are_dirs(dirname(samplefile))
    fwrite(sdata(object), samplefile, sep = '\t', row.names = FALSE)
    return(samplefile)
}

#==============================================================================
#
#                            .read_proteingroups
#                            .read_phosphosites
#
#==============================================================================

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

#' Read maxquant file
#' @param file              proteingroups/phosphosites file
#' @param quantity          value in MAXQUANT_PATTERNS
#' @param samplefile        sample file
#' @param sampleidvar       sampleid var
#' @param subgroupvar       subgroup var
#' @param select_subgroups  character vector
#' @param invert_subgroups  character vector
#' @param verbose           TRUE/FALSE
#' @examples
#' # LFQ INTENSITIES
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     object <- .read_maxquant(file)
#'     samplefile <- create_samplefile(object, default_samplefile(file))
#'     object <- .read_maxquant(file, samplefile=samplefile)
#'
#' # LABEL RATIOS
#'     file <- download_data('billing16.proteingroups.txt')
#'     inv <- c('EM_E','BM_E','BM_EM')
#'     object <- .read_maxquant(file, invert_subgroups = inv)
#'
#' # LABEL INTENSITIES
#'     object <- .read_maxquant(file, quantity = 'Intensity labeled')
#'
#' # SELECTED LABEL RATIOS
#'     file <-  download_data('billing19.proteingroups.txt')
#'     select_subgroups <-  c(sprintf(
#'         '%s_STD', c('EM00','EM01', 'EM02','EM05','EM15','EM30', 'BM00')))
#'     object <- .read_maxquant(file, select_subgroups = select_subgroups)
#'
#' # PHOSPHO RATIOS
#'     file <- download_data('billing19.phosphosites.txt')
#'     .read_maxquant(file, select_subgroups = select_subgroups)
#'
#' @export
.read_maxquant <- function(file, quantity = guess_maxquant_quantity(file),
    samplefile = NULL, sampleidvar = 'sample_id', subgroupvar = 'subgroup',
    select_subgroups = NULL, invert_subgroups = character(0),
    verbose = TRUE){
# Read
    assert_all_are_existing_files(file)
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS))
    assert_is_a_bool(verbose)
    phospho <- stri_detect_fixed(basename(file), 'phospho')
    if (verbose) message('\tRead ', file)
    names1 <- names(fread(file, integer64 = 'numeric', nrows=1))
    fids1  <- fread(file, select = 'id')[[1]]
    FVARS <- if (phospho) PHOSPHOSITE_FVARS else PROTEINGROUP_FVARS
    fdata1 <- fread(file, select = intersect(FVARS, names1))
    setnames(fdata1, 'id', 'feature_id')
    if ('Gene names' %in% names(fdata1)) setnames(fdata1, 'Gene names',
                                                        'feature_name')
    select <- names1 %>%
            extract(stri_detect_regex(., MAXQUANT_PATTERNS[[quantity]]))
    if (phospho)  select %<>% extract(stri_endswith_fixed(., '___1'))
    exprs1 <- data.matrix(fread(file, select = select, integer64 = 'numeric'))
# Rm NA fids
    idx <- !is.na(fids1)
    fids1 <- fids1[idx]; fdata1 <- fdata1[idx, ]; exprs1 <- exprs1[idx, ]
    rownames(exprs1) <- fids1  # first na need to be removed!
# Simplify maxquant snames
    sids1 <- colnames(exprs1)
    if (phospho) sids1 %<>% stri_replace_last_fixed('___1', '')
    sids1 %<>% standardize_maxquant_snames(quantity)
    sids1 %<>% demultiplex()
    colnames(exprs1) <- sids1
# Create sumexp
    object <- matrix2sumexp(exprs1, featuredata = fdata1)
    object$subgroup <- NULL
    object %<>% merge_samplefile(samplefile = samplefile,
                        sampleidvar = sampleidvar, subgroupvar  = subgroupvar)
    object %<>% add_subgroup(verbose=verbose)
# Filter/Transform
    object %<>% filter_maxquant_samples(
                    select_subgroups = select_subgroups, verbose)
    exprs(object) %<>% zero_to_na(verbose = verbose)
    exprs(object) %<>% nan_to_na( verbose = verbose)
    object %<>% log2transform(verbose=verbose)
    object %<>% invert(invert_subgroups)
# Return
    metadata(object)$quantity <- quantity
    metadata(object)$platform <- 'maxquant'
    metadata(object)$file     <- file
    object
}


#==============================================================================
#
#               prepare maxquant features
#
#==============================================================================


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


rm_reverse <- function(object, verbose){
    Reverse <- NULL
    fdata(object)$Reverse %<>% na_to_string()
    object %<>% filter_features(Reverse != '+', verbose=verbose)
    fdata(object)$Reverse <- NULL
    object
}


#' @examples
#' file <- download_data('billing19.proteingroups.txt')
#' object <- read_proteingroups(file, contaminants=TRUE)
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
#                     prepare_maxquant samples
#
#==============================================================================

add_maxquant_sdata <- function(
    object, samplefile = default_samplefile(object), verbose
){
    snames(object) %<>% stri_replace_last_fixed('___1', '') # PHOSPHOSITES
    object %<>% standardize_maxquant_snames(verbose = verbose)
    object %<>% demultiplex(verbose = verbose)
    object %<>% merge_samplefile(samplefile = samplefile, verbose = verbose)
    object
}


filter_maxquant_samples <- function(object, select_subgroups, verbose){
    object %<>% filter_samples_available_for_some_feature(verbose = verbose)
    if (!is.null(select_subgroups)){
        object %<>% filter_samples(subgroup %in% select_subgroups,
                                    verbose = verbose)
        object$subgroup %<>% factor(select_subgroups)
    } else {
        object$subgroup %<>% factor()
    }
    object
}


#==============================================================================
#
#                     transform_maxquant
#
#==============================================================================

transform_maxquant <- function(object, impute, verbose, plot){
# Remove batch effect
    if (verbose) message('\tTransform exprs')
    if (grepl('Reporter intensity', metadata(object)$quantity)){
        message('\t\tTMT: rm run effect')
        suppressWarnings(exprs(object) %<>% limma::removeBatchEffect(
                                                batch = object$replicate))
    }
# Impute
    if (impute) object %<>% impute_systematic_nondetects(plot = FALSE)
    object
}


#==============================================================================
#
#                   read_proteingroups
#                   read_phosphosites
#                       phospho_expr_columns
#                       add_occupancies
#==============================================================================


#' Read proteingroups/phosphosites
#'
#' @param proteinfile       proteingroups file path
#' @param phosphofile       phosphosites  file path
#' @param fastafile         NULL or fastafile (to deconvolute proteingroups)
#' @param quantity          string: "Ratio normalized",
#'                                   "Ratio",
#'                                   "LFQ intensity",
#'                                   "Reporter intensity corrected",
#'                                   "Reporter intensity",
#'                                   "Intensity labeled",
#'                                   "Intensity"
#' @param contaminants      whether to return contaminants (TRUE/FALSE)
#' @param reverse           whether to return reverse peptides (TRUE/FALSE)
#' @param min_localization_prob min site localization probability (number)
#' @param select_subgroups  subgroups to be selected (character vector)
#' @param invert_subgroups  subgroups to be inverted (character vector)
#' @param samplefile       samplefile path
#' @param impute            whether to impute consistent nondetects (TRUE/FALSE)
#' @param formula           formula to create design matrix (using svars)
#' @param contrastdefs      contrast definition vector/matrix/list
#' @param verbose           TRUE/FALSE
#' @param plot              TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' # LFQ intensities
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(file)
#'
#' # Ratios: invert if required
#'     file <- download_data('billing16.proteingroups.txt')
#'     inv <- c('EM_E','BM_E','BM_EM')
#'     object <- read_proteingroups(file, invert_subgroups = inv,
#'                                  contrastdefs = c('E_EM', 'E_BM', 'EM_BM'))
#'
#' # Unnormalized Intensities: specify quantity
#'     object <- read_proteingroups(file, quantity = 'Intensity labeled')
#'
#' # Internal Standard Ratios: rm meaningless ratios
#'     file <-  download_data('billing19.proteingroups.txt')
#'     select_subgroups <-  c(sprintf(
#'         '%s_STD', c('EM00','EM01', 'EM02','EM05','EM15','EM30', 'BM00')))
#'     object <- read_proteingroups(file, select_subgroups = select_subgroups)
#'
#' # Phosphosites
#'     phosphofile <- download_data('billing19.phosphosites.txt')
#'     proteinfile <- download_data('billing19.proteingroups.txt')
#'     read_phosphosites(phosphofile, proteinfile,
#'                         select_subgroups = select_subgroups)
#'
#' # Deconvolute proteingroups
#'     fdata(object)[1:3, 1:4]
#'     object <- read_proteingroups(file,quantity = 'Intensity labeled',
#'         fastafile=download_data('uniprot_hsa_20140515.fasta'), plot=FALSE)
#'     fdata(object)[1:3, 1:4]
#' @export
read_proteingroups <- function(
    proteinfile, quantity = guess_maxquant_quantity(proteinfile),
    samplefile = default_samplefile(proteinfile),
    select_subgroups = NULL, contaminants = FALSE,
    reverse = FALSE, fastafile = NULL, invert_subgroups = character(0),
    impute = stri_detect_regex(quantity, "[Ii]ntensity"),
    formula      = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    contrastdefs = contrast_subgroups(object),
    verbose      = TRUE, plot = TRUE
){
# Assert
    assert_all_are_existing_files(proteinfile)
    if (!is.null(fastafile)) assert_all_are_existing_files(fastafile)
    formula      <- enexpr(formula)
    contrastdefs <- enexpr(contrastdefs)
# Read
    object <- .read_maxquant(proteinfile, quantity,
        samplefile = samplefile, select_subgroups = select_subgroups,
        invert_subgroups = invert_subgroups, verbose = verbose)
# Prepare
    object %<>% filter_maxquant_features(reverse = reverse,
                    contaminants = contaminants, verbose = verbose)
    object %<>% rename_proteingroup_fvars()
    object %<>% simplify_proteingroups(fastafile)
    object %<>% transform_maxquant(impute=impute, verbose=verbose, plot=plot)
# Analyze
    object %<>% pca()
    object %<>% add_limma(formula = eval_tidy(formula),
                    contrastdefs  = eval_tidy(contrastdefs), plot = FALSE)
# Return
    if (plot)  plot_samples(object)
    object
}

#' @rdname read_proteingroups
#' @export
read_phosphosites <- function(
    phosphofile,
    proteinfile = paste0(dirname(phosphofile), '/proteinGroups.txt'),
    quantity = guess_maxquant_quantity(phosphofile),
    samplefile = default_samplefile(proteinfile),
    select_subgroups = NULL,
    contaminants = FALSE, reverse = FALSE, min_localization_prob = 0.75,
    fastafile = NULL, invert_subgroups = character(0),
    formula      = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup,
    contrastdefs = contrast_subgroups(object),
    verbose      = TRUE, plot = TRUE
){
# Assert
    `Protein group IDs` <- `Localization prob` <- NULL
    assert_all_are_existing_files(c(phosphofile, proteinfile))
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS))
    formula      <- enexpr(formula)
    contrastdefs <- enexpr(contrastdefs)
# Read
    proteingroups <- .read_maxquant(file=proteinfile, quantity = quantity,
        samplefile = samplefile, select_subgroups = select_subgroups,
        invert_subgroups = invert_subgroups, verbose=verbose)
    object  <- .read_maxquant(file = phosphofile, quantity = quantity,
        samplefile = samplefile, select_subgroups = select_subgroups,
        invert_subgroups = invert_subgroups, verbose = verbose)
    object %<>% filter_maxquant_features(
                    reverse = reverse, contaminants = contaminants,
                    min_localization_prob = min_localization_prob,
                    verbose = verbose)
    object %<>% add_occupancies(proteingroups, verbose)
# Prepare
    object %<>% rename_phospho_fvars()
    object %<>% simplify_proteingroups(fastafile)
    object %<>% transform_maxquant(impute=FALSE,verbose=verbose,plot=plot)
# Contrast
    object %<>% pca()
    object %<>% add_limma(formula = eval_tidy(formula),
                    contrastdefs  = eval_tidy(contrastdefs), plot = FALSE)
# Return
    if (plot)  plot_samples(object)
    object
}




#' Which are the phospho expr columns?
#' @param x phosphosites colnames
#' @examples
#' phosphosites <- download_data('billing19.phosphosites.txt')
#' x <- names(data.table::fread(phosphosites))
#' quantity <- guess_maxquant_quantity(phosphosites)
#' phospho_expr_columns(x, quantity)
#' @noRd
phospho_expr_columns <- function(x, quantity){
    pattern <- MAXQUANT_PATTERNS[[quantity]]
    which(stri_detect_regex(x,pattern) & stri_detect_regex(x, '___1'))
}


add_occupancies <- function(phosphosites, proteingroups, verbose){
# Initialize
    . <- phospho_value <- protein_value <- occupancy <- NULL
    `Protein group IDs` <- NULL
# Report
    if (verbose) message(
        '\t\tAdd occupancies(phospho) = exprs(phospho) - exprs(proteins)')
# phospho datatable
    cols <- c('feature_id', 'Protein group IDs')
    phospho_dt <- sumexp_to_wide_dt(phosphosites, fvars = 'Protein group IDs')
    phospho_dt %<>% separate_rows(`Protein group IDs`, sep = ';' )
    phospho_dt %<>% data.table()
    phospho_dt %<>% data.table::melt.data.table(
        id.vars = c('feature_id', 'Protein group IDs'))
    setnames(phospho_dt,
        c('feature_id', 'Protein group IDs', 'variable',  'value'),
        c('phospho_id', 'protein_id',        'sample_id', 'phospho_value'))
# proteingroups datatable
    protein_dt <- sumexp_to_long_dt(proteingroups)
    setnames(
        protein_dt,c('feature_id', 'value'), c('protein_id', 'protein_value'))
# merge
    phospho_dt %<>% merge(  protein_dt,
                            by = c('protein_id', 'sample_id'), all.x = TRUE)
    phospho_dt %<>% extract(,
                        .(  phospho_value = unique(phospho_value),
                            protein_value = median(protein_value, na.rm = TRUE),
                            subgroup      = unique(subgroup)),
                        by = c('phospho_id', 'sample_id') )
    phospho_dt[ is.na(protein_value), protein_value := 0 ]
# compute occupancy
    phospho_dt[, occupancy := phospho_value - protein_value]
    phospho_dt %<>% data.table::dcast(
                        phospho_id ~ sample_id, value.var = 'occupancy')
    occupancy_mat <- data.matrix(phospho_dt[, -1])
    rownames(occupancy_mat) <- phospho_dt[[1]]
    occupancy_mat %<>% extract(, snames(phosphosites))
    occupancy_mat %<>% extract(fnames(phosphosites), )
    occupancies(phosphosites) <- occupancy_mat

    phosphosites
}




#=========================================================================

#' @title Get/Set occupancies
#' @description Get / Set occupancies matrix
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' proteinfile <- download_data('billing19.proteingroups.txt')
#' phosphofile <- download_data('billing19.phosphosites.txt')
#' object <- read_phosphosites(phosphofile, proteinfile)
#' exprs(object)[1:3, 1:3]
#' occupancies(object)[1:3, 1:3]
#' @rdname occupancies
#' @export
setGeneric('occupancies', function(object)   standardGeneric("occupancies"))

#' @rdname occupancies
setMethod(
    "occupancies",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$occupancies)


#' @rdname occupancies
#' @export
setGeneric(
    'occupancies<-',
    function(object, value) standardGeneric("occupancies<-"))

#' @rdname occupancies
setReplaceMethod(
    "occupancies",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$occupancies <- value
        object })

#' @rdname occupancies
setReplaceMethod(
    "occupancies",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        assays(object)$occupancies[] <- value
        object })
