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
#' character vector, dataframe, or SummarizedExperiment.
#'
#' @param x character vector, dataframe, or SummarizedExperiment
#' @param ... used for proper S3 method dispatch
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
#' # dataframe
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     x <- data.table::fread(file)
#'     guess_maxquant_quantity(x)
#'
#' # SummarizedExperiment
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(file, plot=FALSE)
#'     guess_maxquant_quantity(file)
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
    for (quantity in names(MAXQUANT_PATTERNS_QUANTITY)){
        pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.data.frame <- function(x, ...){
    x <- names(x)
    for (quantity in names(MAXQUANT_PATTERNS_QUANTITY)){
        pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.SummarizedExperiment <- function(x, ...){
    x <- snames(x)
    for (quantity in names(MAXQUANT_PATTERNS_QUANTITY)){
        pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
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
#'    x <- "Ratio M/L normalized STD(L)_E00(M)_E01(H)_R1"
#'    standardize_maxquant_snames(x)
#'
#'    x <- "Ratio M/L STD(L)_E00(M)_E01(H)_R1"
#'    standardize_maxquant_snames(x)
#'
#'    x <-'LFQ intensity STD_R1'
#'    standardize_maxquant_snames(x)
#'
#'    x <- 'LFQ intensity L STD(L)_E00(M)_E01(H)_R1'
#'    standardize_maxquant_snames(x)
#'
#'    x <-'Reporter intensity 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1'
#'    standardize_maxquant_snames(x)
#'
#'    x <- 'Reporter intensity corrected 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1'
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
    # x = multiplexes + channels. Return multiplexes if single channel.
    pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]

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
extract_reviewed <- function(fastahdrs){
    fastahdrs                    %>%   # seqinr::read.fasta          gives '>sp|tr' fastahdrs
    split_extract_fixed('|',1)         %>%   # Biostrings::readAAStringSet gives  'sp|tr' fastahdrs
    substr(nchar(.)-1, nchar(.))       %>%   # this function now works with both scenarios
    equals('sp')                       %>%
    as.integer()
}

#==============================================================================
#
#                               demultiplex
#
#==============================================================================

extract_protein <- function(fastahdrs){
    y <- fastahdrs %>% split_extract_fixed(' ', 1)
    idx <- stri_detect_fixed(y, '|')         # >IL2-FAT1  (645 aa)
    y[idx] %<>% split_extract_fixed('|', 3)  # >sp|P22234|PUR6_HUMAN
    #split_extract_fixed('_', 1)             # drop _HUMAN: not robust in multi-organism db!
    y
}

extract_gene <- function(fastahdrs){
    fastahdrs %>% 
    split_extract_fixed('GN=', 2)      %>% 
    split_extract_fixed(' ', 1)}

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
#'    demultiplex("STD(L)_E00(M)_E01(H)_R1{M/L}") # Label Ratio
#'    demultiplex('A(0)_B(1)_C(2)_D(3)_R1{0}'     ) # Reporter intensity
#'    demultiplex('STD(L)_E00(M)_E01(H)_R1{L}')   # Label Intensity
#'
#'    # Alternate separators supported
#'    demultiplex('STD(L)_E00(M)_E01(H)_R1{L}')   # underscore
#'    demultiplex('STD(L).E00(M).E01(H).R1{L}')   # dot
#'    demultiplex('STD(L)E00(M)E01(H).R1{L}')     # no separator
#'
#'    # Composite snames supported
#'    demultiplex("WT.t0(L)_WT.t1(M)_WT.t2(H)_R1{H/L}")
#'
#'    # Uniqueness ensured by appending labels when necessary
#'    demultiplex(c("STD(L).M00(M).M00(H).R10{M/L}",
#'                        "STD(L).M00(M).M00(H).R10{H/L}"))
#'
#'    # Uniplexed snames are returned unchanged
#'    demultiplex(c('STD_R1', 'EM0_R1'))
#'
#' # SummarizedExperiment
#'    require(magrittr)
#'    file <- download_data('billing19.proteingroups.txt')
#'    # x <- .read_proteingroups(file)
#'    # x %<>% standardize_maxquant_snames()
#'    # demultiplex(x, verbose = TRUE)
#' @noRd
demultiplex <- function (x, ...) {
    UseMethod("demultiplex", x)
extract_uniprot <- function(fastahdrs){
    fastahdrs                    %>%
    split_extract_fixed(' ', 1)        %>% 
    split_extract_fixed('|', 2) 
}

is_multiplexed <- function(x){
    pattern <- '(.+)\\{(.+)\\}'
    n_open   <- stri_count_fixed(x, '(')
    n_closed <- stri_count_fixed(x, ')')
    all(stri_detect_regex(x, pattern) & (n_open==n_closed) & (n_open>0))
extract_canonical <- function(fastahdrs){
    fastahdrs                    %>%
    extract_uniprot()                  %>%
    split_extract_fixed('-', 1) 
}

are_all_identical <- function(y){
    if (length(y)==1) TRUE else all(y[-1] == y[1])
extract_isoform <- function(fastahdrs){
    fastahdrs                    %>%
    extract_uniprot()                  %>%
    split_extract_fixed('-', 2)        %>%
    nastring_to_0()                    %>%
    as.integer()
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
#' x <- "STD(L)_E00(M)_E01(H)_R1{M/L}"
#' extract_multiplexes(x)     # STD(L)_E00(M)_E01(H)_R1
#' extract_channels(x) #                           M/L
#' @noRd
extract_multiplexes <- function(x){
    pattern <- '(.+)\\{(.+)\\}'
    stri_replace_first_regex(x, pattern, '$1')
extract_description <- function(fastahdrs){
    fastahdrs                    %>%
    split_extract_regex('_[A-Z]+ ', 2) %>% 
    split_extract_regex(' OS=', 1)
}
extract_channels <- function(x){
    pattern <- '(.+)\\{(.+)\\}'
    stri_replace_first_regex(x, pattern, '$2')

extract_fragment <- function(fastahdrs){
    fastahdrs                    %>%
    extract_description()              %>%
    stri_detect_fixed(  'ragment')     %>%
    as.integer()
}

#' Separate samples from labels
#' @param multiplexes character vector: e.g. "STD(L)_E00(M)_E01(H)_R1"
#' @examples
#' x <- "STD(L)_E00(M)_E01(H)_R1"
#' extract_labels(x)
#' extract_samples(x)
#' @noRd
extract_labels <- function(multiplexes){
    pattern <- '\\(.+?\\)'
    stri_extract_all_regex(multiplexes, pattern)       %>%
    lapply(stri_replace_first_fixed, '(', '')  %>%
    lapply(stri_replace_first_fixed, ')', '')
extract_existence <- function(fastahdrs){
    fastahdrs                   %>%
    split_extract_fixed('PE=', 2)     %>%
    split_extract_fixed(' ', 1)       %>%
    nastring_to_nachar()              %>%
    as.integer()
}
extract_biosamples <- function(multiplexes, labels){
    pattern <- '\\(.+?\\)'
    stri_split_regex(multiplexes, pattern) %>%
    lapply(function(y){
            y[seq_len(length(labels[[1]]))] %<>%
            stri_replace_first_regex('^[_. ]', ''); y})
    # rm sep from samples (but not from replicate - needed to glue back later!)

FASTAFIELDS <- c('reviewed', 'protein', 'gene', 'canonical', 
                 'isoform', 'fragment', 'existence', 'organism', 'description')

parse_fastahdrs <- function(fastahdrs, fastafields = setdiff(FASTAFIELDS, 'description')){
    reviewed <- protein <- gene <- canonical <- isoform <- fragment <- NULL
    existence <- organism <- description <- NULL
    dt <- data.table(uniprot = extract_uniprot(    fastahdrs))                              #   0 tr
    if ('reviewed'    %in% fastafields)  dt[, reviewed    := extract_reviewed( fastahdrs)]  #   1 sp
    if ('protein'     %in% fastafields)  dt[, protein     := extract_protein(  fastahdrs)]  # existence
    if ('gene'        %in% fastafields)  dt[, gene        := extract_gene(     fastahdrs)]  #   1 protein
    if ('canonical'   %in% fastafields)  dt[, canonical   := extract_canonical(fastahdrs)]  #   2 transcript
    if ('isoform'     %in% fastafields)  dt[, isoform     := extract_isoform(  fastahdrs)]  #   3 homolog
    if ('fragment'    %in% fastafields)  dt[, fragment    := extract_fragment( fastahdrs)]  #   4 prediction
    if ('existence'   %in% fastafields)  dt[, existence   := extract_existence(fastahdrs)]  #   5 uncertain
    if ('organism'    %in% fastafields)  dt[, organism    := split_extract_fixed(protein, '_', 2)]
    if ('description' %in% fastafields)  dt[, description := extract_description(fastahdrs)]
    if ('existence'   %in% fastafields)  dt[, existence   := unique(.SD)[, existence[!is.na(existence)]], by = 'canonical'] 
        # `unique`: for phosphosites the fastahdrs are
        #  replicated when protein has multiple phosphosites
        #  This duplication needs to be eliminated before proceeding.
    dt[]
}


#' Are number of samples/labels equally multiplexed across mixes?
#' Read headers from uniprot fastafile
#'
#' @param samples list of char vectors (see examples)
#' @param labels  list of char vectors (see examples)
#' @return TRUE or FALSE
#' @param fastafile    string (or charactervector)
#' @param fastafields  charactervector : which fastahdr fields to extract ?
#' @param verbose      bool
#' @examples
#' samples <- list(c("STD", "E00", "E01", "_R1" ),
#'                 c("STD", "E00", "E01", "_R2"))
#' labels  <- list(c("L", "M", "H"),
#'                 c("L", "M", "H"))
#' are_equally_multiplexed(samples, labels)
#' @noRd
are_equally_multiplexed <- function(samples, labels){
    n_samples <- vapply(samples,  length, integer(1))
    n_labels  <- vapply(labels, length, integer(1))
    are_all_identical(n_samples) & are_all_identical(n_labels)
#' # Single fastafile
#'    fastafile <- download_data('uniprot_hsa_20140515.fasta')
#'    read_fastahdrs(fastafile)
#'    
#' # Multiple fastafiles
#'    # dir <- R_user_dir('autonomics', 'cache')
#'    # dir %<>% file.path('uniprot', 'uniprot_sprot-only2014_01')
#'    # fastafile <- file.path(dir, c("uniprot_sprot_human.fasta", "uniprot_sprot_mouse.fasta"))
#'    # read_fastahdrs(fastafile)
#' @return data.table(uniprot, protein, gene, uniprot, reviewed, existence)
#' @note existence values are always those of the canonical isoform
#'       (no isoform-level resolution for this field)
#' @export
read_fastahdrs <- function(
    fastafile, fastafields = setdiff(FASTAFIELDS, 'description'), verbose = TRUE
){
# Assert
    if (is.null(fastafile)) return(NULL)
    assert_all_are_existing_files(fastafile)
# Read
    if (verbose) message('\tRead ', paste0(fastafile, collapse = '\n\t     '))
    fastahdrs <- mapply(.read_fastahdrs, fastafile, MoreArgs = list(fastafields = fastafields), SIMPLIFY = FALSE)
    fastahdrs %<>% data.table::rbindlist()
    fastahdrs
}


#' Extract/Drop replicate
#' @param samples list of char vectors: list(c("STD", "E00", "E01", "_R1" ))
#' @param labels  list of char vectors: list(c("L", "M", "H"))
#' @return character vector: "_R1"
#' @examples
#' samples <- list(c("STD", "E00", "E01", "_R1" ))
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
.read_fastahdrs <- function(
    fastafile, fastafields = setdiff(FASTAFIELDS, 'description'), verbose = TRUE
){
    if (!requireNamespace('Biostrings', quietly = TRUE)){
        stop("BiocManager::install('Biostrings'). Then re-run.") }
    fastahdrs <- Biostrings::readAAStringSet(fastafile)
    fastahdrs %<>% names()
    fastahdrs %<>% parse_fastahdrs(fastafields = fastafields)
    fastahdrs
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
#' biosamples <- list(c(  L = 'STD', M = 'E00', H = 'E01'),
#'                 c(  L = 'STD', M = 'E00', H = 'E01'))
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

#' @noRd
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


#' @noRd
demultiplex.SummarizedExperiment <- function(x,verbose  = FALSE, ...){
    newsnames <- demultiplex(snames(x), verbose = verbose)
    snames(x) <- sdata(x)$sample_id <- newsnames
    x
}

#==============================================================================
#
#                            load_uniprot_fasta
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
#' @param fastafields character vector
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
#---------------------------------------------------------------------------
#
#                      add_psp
#
#---------------------------------------------------------------------------

#' Add psp
#' 
#' Add PhosphoSitePlus literature counts
#' 
#' Go to www.phosphosite.org                   \cr
#' Register and Login.                         \cr
#' Download Phosphorylation_site_dataset.gz'.  \cr
#' Save into: file.path(R_user_dir('autonomics','cache'),'phosphositeplus')
#' @param object     SummarizedExperiment
#' @param pspfile    phosphositeplus file
#' @return  SummarizedExperiment
#' @examples 
#' phosphofile <- download_data('billing19.phosphosites.txt')
#' proteinfile <- download_data('billing19.proteingroups.txt')
#' object <- read_maxquant_phosphosites(phosphofile = phosphofile, proteinfile = proteinfile)
#' fdt(object)
#' object %<>% add_psp()
#' fdt(object)
#' @export
add_psp <- function(
    object, 
    pspfile = file.path(R_user_dir('autonomics', 'cache'), 
            'phosphositeplus', 'Phosphorylation_site_dataset.gz')
){
# Initialize
    AA <- ACC_ID <- `Amino acid` <- LT_LIT <- MOD_RSD <- MS_CST <- MS_LIT <- psp <- NULL
# Read
    assert_is_all_of(object, 'SummarizedExperiment')
    if (!file.exists(pspfile))  return(object)
    dt <- data.table::fread(pspfile)
    dt[is.na(LT_LIT), LT_LIT := 0]
    dt[is.na(MS_LIT), MS_LIT := 0]
    dt[is.na(MS_CST), MS_CST := 0]
    dt[, psp := LT_LIT + MS_LIT + MS_CST]
    dt$AA <- dt$MOD_RSD %>% substr(1, 1)
    dt$MOD_RSD %<>% substr(2, nchar(.)-2)
# Rm duplicates
    idx <- cduplicated(dt[, .(ACC_ID, MOD_RSD)])
    # dt[idx]  # Seems like a (single) psiteplus bug: two different rows for the same site
    #          # On the website it says 'under review' for this particular site
    #          # Rm these
    dt %<>% extract(!idx)
    dt %<>% extract(, .(ACC_ID, MOD_RSD, psp, AA))
    # Look at example
    # dt %<>% extract('ZZZ3',   on = 'PROTEIN')
    # dt %<>% extract(c('S606-p', 'S613-p'), on = 'MOD_RSD')
# Merge
    fdt(object)$uniprot1  <- split_extract_fixed(fdt(object)$Curated, ';', 1)
    fdt(object)$Position1 <- split_extract_fixed(fdt(object)$`Positions within proteins`, ';', 1)
    object %<>% merge_fdt(dt, by.x = c('uniprot1', 'Position1'), 
                                by.y = c('ACC_ID', 'MOD_RSD'))
    idx <- fdt(object)[, !is.na(psp) &  AA != `Amino acid`]
    fdt(object)$psp[idx] <- NA
    fdt(object)$AA <- NULL
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
    if (verbose) message('\tTransform exprs')
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
    fosdt[ is.na(protein), protein := 0 ]
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
    sfile = NULL, sfileby = NULL, subgroupvar = 'subgroup',
    select_subgroups = NULL, invert_subgroups = character(0),
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
    object <- SummarizedExperiment(exprs1, rowData=fdata1)
    if (!phospho)  object %<>% add_pepcounts(file, pepcountpattern, quantity)
# Process
    if (phospho)  colnames(object) %<>% stri_replace_last_fixed('___1', '')
    colnames(object) %<>% standardize_maxquant_snames(quantity)
    colnames(object) %<>% demultiplex()
    object$sample_id <- colnames(object)
    object %<>% merge_sfile(sfile = sfile, by.x = 'sample_id', by.y = sfileby)
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

un_int64 <- function(x) {
    dplyr::mutate(x, dplyr::across(dplyr::where(bit64::is.integer64), as.numeric))
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
#' @param sfileby       sample file mergeby column
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
#' @param contrastdefs  contrastdef vector/matrix/list
#' @param pca           whether to pca
#' @param fit           fit model: NULL, 'limma', 'lm', 'lme', 'lmer','wilcoxon'
#' @param verbose       whether to message
#' @param plot          whether to plot
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, pca=TRUE, fit='limma')
#' @export
read_proteingroups <- function(
    file, quantity = guess_maxquant_quantity(file), sfile = NULL,
    sfileby = NULL, select_subgroups = NULL, contaminants = FALSE,
    reverse = FALSE, fastafile = NULL, invert_subgroups = character(0),
    impute = stri_detect_regex(quantity, "[Ii]ntensity"),
    pepcountpattern = MAXQUANT_PATTERNS_PEPCOUNTS[1], subgroupvar = NULL,
    formula = NULL, block = NULL, contrastdefs = NULL,
    pca = FALSE, fit = NULL, verbose = TRUE, plot = TRUE
){
# Assert
    . <- NULL
    assert_all_are_existing_files(file)
    if (!is.null(fastafile)) assert_all_are_existing_files(fastafile)
# Read
    object <- .read_maxquant(file, quantity,
        sfile = sfile, sfileby = sfileby, select_subgroups = select_subgroups,
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
    object %<>% analyze(pca=pca, fit=fit, subgroupvar = subgroupvar, 
                        formula = formula, block = block, 
                        contrastdefs = contrastdefs, 
                        verbose = verbose, plot = plot)
# Return
    object
}

#' @rdname read_proteingroups
#' @export
read_phosphosites <- function(
    file, proteinfile = paste0(dirname(file), '/proteinGroups.txt'),
    quantity = guess_maxquant_quantity(file),
    sfile = NULL, sfileby = NULL, select_subgroups = NULL, contaminants = FALSE,
    reverse = FALSE, min_localization_prob = 0.75, fastafile = NULL,
    invert_subgroups = character(0), pca = FALSE,
    fit = NULL, subgroupvar = NULL, formula = NULL, block = NULL, 
    contrastdefs = NULL, verbose = TRUE, plot = TRUE
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
        sfile = sfile, sfileby = sfileby, select_subgroups = select_subgroups,
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
                        formula=formula, block=block, contrastdefs=contrastdefs,
                        verbose=verbose, plot=plot)
# Return
    object
}

