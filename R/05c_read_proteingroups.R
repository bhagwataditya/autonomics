#====================================
# LCMS PROTEINGROUPS & PHOSPHOSITES
#====================================

#' maxquant patterns
#' @examples
#' maxquant_patterns
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
#'     x <- download_data('stemcells.proteinGroups.txt')
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
#'     file <- download_data( 'stemcells.proteinGroups.txt')
#'     x <- data.table::fread(file)
#'     guess_maxquant_quantity(x)
#'
#' # SummarizedExperiment
#'      file <-download_data( 'stemcells.proteinGroups.txt')
#'      x <- read_proteingroups(file, demultiplex_snames = FALSE)
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
#' @noRd
proteingroups_fvars <- c(
    'id', 'Majority protein IDs', 'Protein names', 'Gene names',
    'Contaminant', 'Potential contaminant', 'Reverse', 'Phospho (STY) site IDs')


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
#'    file <- download_data('stemcells.proteinGroups.txt')
#'    x <- read_proteingroups(file, demultiplex_snames = FALSE)
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
    pattern <- maxquant_patterns[[quantity]]

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
#'    file <- download_data('stemcells.proteinGroups.txt')
#'    x <- read_proteingroups(file, demultiplex = FALSE)
#'    x %<>% standardize_maxquant_snames()
#'    demultiplex(x, verbose = TRUE)
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
    is_ratio <- channels %>% stri_detect_fixed('/') %>% all()
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


#' Read proteingroups
#' @param file        string: path to 'proteinGroups.txt'
#' @param quantity    string: any value in \code{maxquant_patterns}
#' @param fvars       character vector: annotation columns to extract
#' @param demultiplex_snames  TRUE (default) or FALSE: whether to demultiplex
#'    maxquant snames: Ratio normalized H/L WT(L).KD(H).R1 -> KD(H)_WT(L).R1 ?
#' @param verbose     TRUE (default) or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
#' read_proteingroups(file)
#' @export
read_proteingroups <- function(file, quantity = guess_maxquant_quantity(file),
    fvars = proteingroups_fvars, demultiplex_snames = TRUE, verbose = TRUE){
# Assert
    assert_all_are_existing_files(file)
    assert_is_subset(quantity, names(maxquant_patterns))
    assert_is_a_bool(verbose)
# Initial Read
    assert_all_are_existing_files(file)
    dt <- fread(file, integer64 = 'numeric', header = TRUE)
    fvars %<>% intersect(names(dt))
# Define components
    fid_rows   <- 2:nrow(dt)
    fid_cols   <- which(names(dt) == 'id')
    sid_rows   <- 1
    sid_cols   <- which(
                    stri_detect_regex(names(dt), maxquant_patterns[[quantity]]))
    expr_rows  <- 2:nrow(dt)
    expr_cols  <- sid_cols
    fvar_rows  <- 1
    fvar_cols  <- match(fvars, names(dt))
    fdata_rows <- 2:nrow(dt)
    fdata_cols <- fvar_cols
# Read sumexp
    object <- read_omics(file, fid_rows = fid_rows,     fid_cols   = fid_cols,
                               sid_rows   = sid_rows,   sid_cols   = sid_cols,
                               expr_rows  = expr_rows,  expr_cols  = expr_cols,
                               fvar_rows  = fvar_rows,  fvar_cols  = fvar_cols,
                               fdata_rows = fdata_rows, fdata_cols = fdata_cols,
                               transpose  = FALSE,      verbose    = verbose)
# Clean sdata
    if (demultiplex_snames){
        object %<>% standardize_maxquant_snames(verbose = verbose)
        object %<>% demultiplex(verbose = verbose) }
    object$subgroup <- guess_subgroup_values(object$sample_id, verbose=verbose)
    object$replicate<- guess_subgroup_values(object$sample_id, invert = TRUE,
                                            verbose = FALSE)
    #object$block  <- object$sample_id %>% guess_subject_values( verbose = TRUE)
# Clean fdata
    contaminant_var <- intersect(
        fvars(object), c('Contaminant', 'Potential contaminant'))
    fdata(object)[[contaminant_var]] %<>% (function(x){x[is.na(x)] <- ''; x})
    fdata(object)[['Reverse'      ]] %<>% (function(x){x[is.na(x)] <- ''; x})
    fdata(object)$feature_name    <- fdata(object)$`Gene names`
    fdata(object)$feature_uniprot <- fdata(object)$`Majority protein IDs`
    fdata(object) %<>% pull_columns(
                        c('feature_id', 'feature_name', 'feature_uniprot'))
# Return
    object
}


#' phosphosites fvars
#' @noRd
phosphosite_fvars <- c(
  'id', 'Protein group IDs', 'Positions within proteins', 'Localization prob')

#' Which are the phospho expr columns?
#' @param x phosphosites colnames
#' @examples
#' phosphofile <- download_data('differentiation.phosphoSites.txt')
#' x <- names(data.table::fread(phosphofile))
#' quantity <- guess_maxquant_quantity(phosphofile)
#' phospho_expr_columns(x, quantity)
#' @noRd
phospho_expr_columns <- function(x, quantity){
    pattern <- maxquant_patterns[[quantity]]
    which(stri_detect_regex(x,pattern) & stri_detect_regex(x, '___1'))
}


filter_single_proteingroup <- function(phosphosites, verbose){
    idx <- stri_detect_fixed(
                rowData(phosphosites)$`Protein group IDs`, ';', negate = TRUE)
    phosphosites %<>% subset(idx, )
    analysis(phosphosites)$nfeatures %<>% c(singlePG = nrow(phosphosites))
    if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx),
                        ' phosphosites mapping to a single proteingroup')
    phosphosites
}

add_occupancies <- function(phosphosites, proteingroups, verbose){
    #if (verbose)  message("\t\tRead: ", proteinfile)
  proteingroups %<>% extract(rowData(phosphosites)$`Protein group IDs`, )
    if (verbose) message(
        '\t\tCalculate occupancies(phospho) = exprs(phospho) - exprs(proteins)')
    occupancies(phosphosites) <- exprs(phosphosites) - exprs(proteingroups)
    phosphosites
}

#' Read phosphosites
#' @param phosphofile         string: phosphosites filepath
#' @param proteinfile         string: proteingroups filepath
#' @param quantity            NULL or value in names(maxquant_patterns)
#' @param fvars               string vector
#' @param demultiplex_snames  TRUE (default) or FALSE
#' @param verbose             TRUE (default) or FALSE
#' @return SummarizedExperiment
#' @examples
#' phosphofile <- download_data('differentiation.phosphoSites.txt')
#' proteinfile <- download_data('differentiation.proteinGroups.txt')
#' read_phosphosites(phosphofile, proteinfile)
#' @export
read_phosphosites <- function(
    phosphofile,
    proteinfile        = paste0(dirname(phosphofile), '/proteinGroups.txt'),
    quantity           = guess_maxquant_quantity(phosphofile),
    fvars              = phosphosite_fvars,
    demultiplex_snames = TRUE,
    verbose            = TRUE
){
# Check. Assert
    `Protein group IDs` <- `Localization prob` <- NULL
    assert_all_are_existing_files(c(phosphofile, proteinfile))
    assert_is_subset(quantity, names(maxquant_patterns))
    assert_is_character(fvars)
# Initial Read
    dt <- fread(phosphofile, integer64 = 'numeric', header = TRUE)
    value_cols <- phospho_expr_columns(names(dt), quantity)
    fvar_cols  <- which(names(dt) %in% fvars)
# Read phosphosites
    fid_rows   <- 2:nrow(dt);     fid_cols   <- which(names(dt) == 'id')
    sid_rows   <- 1;              sid_cols   <- value_cols
    expr_rows  <- 2:nrow(dt);     expr_cols  <- value_cols
    fvar_rows  <- 1;              fvar_cols  <- fvar_cols
    fdata_rows <- 2:nrow(dt);     fdata_cols <- fvar_cols
    if (verbose) message("\t\tRead: ", phosphofile)
    phosphosites  <- read_omics(phosphofile,
                                fid_rows   = fid_rows,    fid_cols = fid_cols,
                                sid_rows   = sid_rows,    sid_cols = sid_cols,
                                expr_rows  = expr_rows,  expr_cols = expr_cols,
                                fvar_rows  = fvar_rows,  fvar_cols = fvar_cols,
                                fdata_rows = fdata_rows,fdata_cols = fdata_cols,
                                transpose  = FALSE, verbose = verbose)
    phosphosites %<>% filter_single_proteingroup(verbose = verbose)
# Add occupancies
    proteingroups <- read_proteingroups(proteinfile, quantity = quantity,
                                        verbose = FALSE)
    phosphosites %<>% add_occupancies(proteingroups, verbose = verbose)
# Demultiplex snames and return
    if (demultiplex_snames){
        phosphosites %<>% standardize_maxquant_snames(verbose = verbose)
        phosphosites %<>% demultiplex(verbose = verbose)
    }
    phosphosites
}

