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

utils::globalVariables('where')
un_int64 <- function(x) {
    dplyr::mutate(x, dplyr::across(where(bit64::is.integer64), as.numeric))
}

#' Read proteingroups/phosphosites as-is
#' @param proteinfile  file
#' @param phosphofile  file
#' @param verbose  TRUE / FALSE
#' @return data.table
#' @examples 
#' proteinfile <- download_data('billing19.proteingroups.txt')
#' phosphofile <- download_data('billing19.phosphosites.txt')
#' prodt <- .read_proteingroups(proteinfile = proteinfile)
#' fosdt <- .read_phosphosites(phosphofile = phosphofile, proteinfile = proteinfile)
#' @export
.read_proteingroups <- function(proteinfile, verbose = TRUE){
# Assert
    assert_all_are_existing_files(proteinfile)
    assert_all_are_matching_fixed(proteinfile, 'roups.txt')
# Read
    if (verbose)  message('\tRead ', proteinfile)
    prodt <- fread(proteinfile, colClasses = c(id = 'character'), integer64 = 'numeric')
    prodt %<>% un_int64()
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
    setnames(prodt, c('id', 'Majority protein IDs'), c('proId', 'uniprot'))
    prodt
}


.read_phosphosites <- function(phosphofile, proteinfile, verbose = TRUE){
# Assert
    assert_all_are_existing_files(c(proteinfile, phosphofile))
    assert_all_are_matching_fixed(proteinfile, 'roups.txt')
    assert_all_are_matching_fixed(phosphofile, 'ites.txt')
# Read    
    if (verbose)  message('\tRead ', phosphofile)
    fosdt <- fread(phosphofile, colClasses = c(id = 'character'), integer64 = 'numeric')
    fosdt %<>% un_int64()
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
                    c('fosId', 'proId',             'uniprot'))
}


#---------------------------------------------------------------------------
# 
#                      drop_differing_uniprots
#
#---------------------------------------------------------------------------


drop_differing_uniprots <- function(fosdt, prodt, verbose){
# Avoid changing the global env
    prodt %<>% copy()
    fosdt %<>% copy()
# Extract annotation cols
    if (verbose)  message('\t\t\tKeep proteingroup uniprots')
    proanncols <- c('proId', 'uniprot', 'Reverse', 'Potential contaminant')
    fosanncols <- c('fosId', 'proId', 'uniprot', 'Positions within proteins', 
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

nastring_to_0 <- function(x){ 
    x[x=='NA'] <- 0; x
        # Sometimes the canonical isoform is NOT isoform-1 !
        # https://www.uniprot.org/uniprotkb/Q9H0P0/entry#sequences
}

extract_reviewed <- function(`Fasta headers`){
    `Fasta headers`                    %>% 
    substr(2, 3)                       %>%
    equals('sp')                       %>%
    as.integer()
}

extract_entry <- function(`Fasta headers`){
    `Fasta headers` %>% 
    split_extract_fixed(' ', 1)        %>% 
    split_extract_fixed('|', 3)        %>% 
    split_extract_fixed('_', 1)
}

extract_genesymbol <- function(`Fasta headers`){
    `Fasta headers` %>% 
    split_extract_fixed('GN=', 2)      %>% 
    split_extract_fixed(' ', 1)}

extract_uniprot <- function(`Fasta headers`){
    `Fasta headers`                    %>%
    split_extract_fixed(' ', 1)        %>% 
    split_extract_fixed('|', 2) 
}

extract_canonical <- function(`Fasta headers`){
    `Fasta headers`                    %>%
    extract_uniprot()                  %>%
    split_extract_fixed('-', 1) 
}

extract_isoform <- function(`Fasta headers`){
    `Fasta headers`                    %>%
    extract_uniprot()                  %>%
    split_extract_fixed('-', 2)        %>%
    nastring_to_0()                    %>%
    as.integer()
}

extract_proteinname <- function(`Fasta headers`){
    `Fasta headers`                    %>%
    split_extract_regex('_[A-Z]+ ', 2) %>% 
    split_extract_regex(' OS=', 1)
}

extract_fragment <- function(`Fasta headers`){
    `Fasta headers`                    %>%
    extract_proteinname()              %>%
    stri_detect_fixed(  'ragment')     %>%
    as.integer()
}

extract_existence <- function(`Fasta headers`){
    `Fasta headers`                   %>%
    split_extract_fixed('PE=', 2)     %>%
    split_extract_fixed(' ', 1)       %>%
    nastring_to_nachar()              %>%
    as.integer()
}

parse_fastahdrs <- function(`Fasta headers`){
    dt <- data.table(                                       # reviewed
        reviewed    = extract_reviewed(   `Fasta headers`), #   0 tr
        entry       = extract_entry(      `Fasta headers`), #   1 sp
        genesymbol  = extract_genesymbol( `Fasta headers`),
        uniprot     = extract_uniprot(    `Fasta headers`), # existence 
        canonical   = extract_canonical(  `Fasta headers`), #   1 protein 
        isoform     = extract_isoform(    `Fasta headers`), #   1 transcript
        proteinname = extract_proteinname(`Fasta headers`), #   2 homolog
        fragment    = extract_fragment(   `Fasta headers`), #   3 prediction
        existence   = extract_existence(  `Fasta headers`)) #   4 uncertain
    dt[, existence := unique(.SD)[, existence[!is.na(existence)]], by = 'canonical']
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
#' @return data.table(uniprot, entry, gene, protein, reviewed, existence)
#' @note existence values are always those of the canonical isoform
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
    parse_fastahdrs(fastahdrs)
}


#-----------------------------------------------------------------------------
#
#                   maxquant_curate
#                   fasta_curate
#                       .drop_inferior
#
#-----------------------------------------------------------------------------

drop_inferior <- function(anndt, verbose = TRUE){
    anndt %<>% copy()
    idcol <- if ('fosId' %in% names(anndt)) 'fosId' else 'proId'
    
    if (verbose)  message('\t\t\tDrop proteins with NA entry / proteinname')
    anndt <- anndt[    entry   != 'NA']
    anndt <- anndt[proteinname != 'NA']
    
    if (verbose)  message('\t\t\tWithin ', idcol, ': drop trembl    in favour of swissprot')
    anndt <- anndt[, .SD[ reviewed == max(reviewed) ], by = idcol]
    
    if (verbose)  message('\t\t\t       ','     ','  drop fragments in favour of full proteins')
    anndt <- anndt[, .SD[ fragment == min(fragment) ], by = idcol]
    
    if (verbose)  message('\t\t\t       ','     ','  drop worse existences', 
                          ' (1=protein, 2=transcript, 3=homolog, 4=prediction, 5=uncertain)')
    if (!any(is.na(anndt$existence))){
        anndt <- anndt[, .SD[existence == min(existence)], by = idcol] }
    anndt[, c('reviewed', 'fragment', 'existence') := NULL]
    
    if (verbose)  message("\t\t\tDrop 'Isoform x of ' in Protein names")
    anndt$proteinname %<>% stri_replace_first_regex('Isoform [0-9A-Z]+ of ', '')
    
    if (verbose)  message("\t\t\tDrop '(Fragment)'    in Protein names")
    anndt$proteinname %<>% stri_replace_first_fixed(' (Fragment)', '')
    if (assertive::is_numeric_string(anndt[[idcol]][1])){
        anndt[order(as.integer(get(idcol)), canonical, isoform)]
    } else {
        anndt[order(get(idcol), canonical, isoform)]
    }
}


.merge <- function(dt1, dt2, by){
    dt1 %<>% copy()
    dt2 %<>% copy()
    commoncols <- intersect(names(dt1), names(dt2))
    commoncols %<>% setdiff(by)
    if (length(commoncols) > 0)  dt1[, (commoncols) := NULL]
    dt1 %<>% merge(dt2, by = by, sort = FALSE, all.x = TRUE)
    dt1 %<>% pull_columns(names(dt2))
    dt1[]
}

.rbind <- function(dt1, dt2, sortby, as.integer = FALSE){
    dt1only <- setdiff(names(dt1), names(dt2))
    dt2only <- setdiff(names(dt2), names(dt1))
    if (length(dt1only)>0)  dt2[, (dt1only) := '']
    if (length(dt2only)>0)  dt1[, (dt2only) := '']
    dt <- rbind(dt1, dt2)
    sorter <- dt[[sortby]]
    if ({{as.integer}})  sorter %<>% as.integer()
    sorter %<>% order()
    dt[sorter]
}

CURATEDCOLS <- c('entry', 'isoform', 'uniprot', 'canonical', 'proteinname', 'genesymbol')

#' Curate and Annotate.
#'
#' Using Fastafile/MaxQuant fastahdrs
#' 
#' `curate_annotate_maxquant`: MaxQuant  fastahdrs \cr
#' `curate_annotate_fastafile` Fastafile fastahdrs \cr
#' `curate`: Fastafile + (for missing entries) MaxQuant fastahdrs \cr
#' Steps: \cr
#' Within proteingroup
#'   1. Uncollapse
#'   2. Drop lower-quality uniprots
#'         * `sp > tr`
#'         * `fullseqs > fragments`
#'         * `pro > rna > hom > pred`
#'   3. Annotate: entry, isoform, genesymbol, proteinname
#'         * `maxquant_curate :  MaxQuant `
#'         * `   fasta_curate :  Fastafile`
#'   4. Recollapse into `Curated`
#' @param dt      `data.table`
#' @param fastadt `data.table`
#' @param verbose `TRUE / FALSE`
#' @return data.table
#' @examples
#' # Fukuda 2020: MaxQuant 
#'     proteinfile <- download_data('fukuda20.proteingroups.txt')
#'     dt <- .read_proteingroups(proteinfile)
#'     curate_annotate_maxquant(dt)[, 1:7]
#'     curate_annotate_fastafile(dt, fastadt = NULL)[, 1:2]
#'     curate_annotate(dt, fastadt = NULL)[, 1:7]
#'     
#' # Billing 2019: Fastafile + MaxQuant
#'     proteinfile <- download_data('billing19.proteingroups.txt')
#'     phosphofile <- download_data('billing19.phosphosites.txt')
#'     fastafile   <- download_data('uniprot_hsa_20140515.fasta')
#'     fastadt <-  read_fastahdrs(fastafile)
#'     dt <- .read_proteingroups(proteinfile, verbose = TRUE)
#'     dt[, 1:4]
#'     curate_annotate_maxquant( dt)[, 1:7]
#'     curate_annotate_fastafile(dt, fastadt = NULL)[, 1:7]
#'     curate_annotate_fastafile(dt, fastadt)[, 1:7]
#'     curate_annotate(dt, fastadt)[, 1:7]
#' @md
#' @export
curate_annotate <- function(dt, fastadt = NULL, verbose = TRUE){
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    dt %<>% copy()
    dt[, entry := NA_character_]
    dt %<>% curate_annotate_fastafile(fastadt, verbose = verbose)
    idx <- dt$Reverse == '' & dt$`Potential contaminant` == '' & is.na(dt$entry)
    dt1 <- dt[ idx] %>% curate_annotate_maxquant(verbose = verbose)
    dt  <- dt[!idx]
    dt %<>% rbind(dt1)
    dt %<>% extract(order(as.integer(get(idcol))))
    dt[, `Fasta headers` := NULL]
    dt[]
}


#' @rdname curate_annotate
#' @export
curate_annotate_fastafile <- function(dt, fastadt, verbose = TRUE){
# Return if NULL fastadt
    dt %<>% copy()
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    if (is.null(fastadt)){
        dt[, (CURATEDCOLS) := NA_character_]
        dt %<>% pull_columns(c(idcol, CURATEDCOLS))
        return(dt[])
    }
# Drop Curated / Contaminants / Reverse
    idxdrop <- dt$`Potential contaminant`=='+'  |  dt$Reverse=='+'
    if (sum(idxdrop)==length(idxdrop))  return(dt)
    dropdt <- dt[ idxdrop]
    dt     <- dt[!idxdrop]
# Curate
    if (verbose)  message('\t\tFasta data.table')
    if (verbose)  message('\t\t\tUncollapse uniprot accessions')
    anndt <- dt[, c(idcol, 'uniprot'), with = FALSE]
    anndt %<>% extract(, c(idcol, 'uniprot'), with = FALSE)
    anndt %<>% uncollapse(uniprot, sep = ';')
    anndt %<>% merge(fastadt, by = 'uniprot', sort = FALSE)
    anndt %<>% drop_inferior()
    #setnames(anndt, 'uniprot', 'Curated')
    setnames(dt, 'uniprot', 'Original')
    if (verbose)  message('\t\t\tCollapse')
    anndt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = idcol)
    dt %<>% .merge(anndt, by = idcol)
# Pickup
    dt %<>% .rbind(dropdt, sortby = idcol, as.integer = TRUE)
    dt[is.na(uniprot), uniprot := Original]
    dt[, Original := NULL]
    dt %<>% pull_columns(c(idcol, intersect(CURATEDCOLS, names(.))))
    dt[]
}


#' @rdname curate_annotate
#' @export
curate_annotate_maxquant <- function(dt, verbose = TRUE){
# Drop Curated / Contaminants / Reverse
    dt %<>% copy()
    idxdrop <- dt$`Potential contaminant`=='+'  |  dt$Reverse=='+'
    if (sum(idxdrop)==length(idxdrop))  return(dt)
    dropdt <- dt[ idxdrop]
    dt     <- dt[!idxdrop]
# Curate
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    if (verbose)  message('\t\tMaxQuant Fastahdrs')
    if (verbose)  message('\t\t\tUncollapse, Drop truncated, Parse')
    anndt <- dt[, c(idcol, 'Fasta headers'), with = FALSE]
    anndt %<>% uncollapse(`Fasta headers`, sep = ';')
    anndt %<>% extract( stri_count_fixed(`Fasta headers`, '|') == 2)
    anndt %<>% cbind(parse_fastahdrs(anndt$`Fasta headers`))
    anndt %<>% drop_inferior()
    setnames(dt, 'uniprot', 'Original')
    if (verbose)  message('\t\t\tCollapse')
    anndt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = idcol)
    dt %<>% .merge(anndt, by = idcol)
    dt %<>% .rbind(dropdt, sortby = idcol, as.integer = TRUE)
    dt[is.na(uniprot), uniprot := Original]
    dt[, Original := NULL]
    dt %<>% pull_columns(c(idcol, intersect(CURATEDCOLS, names(dt))))
    dt[]
}



#--------------------------------------------------------------------------
#
#                      add_feature_id
#                      process_maxquant
#
#---------------------------------------------------------------------------

add_feature_id <- function(dt){
    dt %<>% copy()
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    dt[, feature_id := entry]
    dt[isoform!=0, feature_id := paste0(feature_id, '(', stri_replace_all_fixed(isoform, ';', ''), ')')]
    if (idcol=='fosId'){
        dt$feature_id %<>% paste0('-',dt$`Amino acid`)
        dt$feature_id %<>% paste0(split_extract_fixed(dt$`Positions within proteins`, ';', 1))
    }
    dt[                Reverse == '+', feature_id := paste0('REV_', get(idcol))]
    dt[`Potential contaminant` == '+', feature_id := paste0('CON_', get(idcol))]
        # Contaminants and actual proteins get mixed up. For naming: make sure to use contaminant.
    dt %<>% pull_columns(c(idcol, 'feature_id'))
    dt[]
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
#' @param verbose  TRUE / FALSE
#' @return matrix
#' @examples 
#' proteinfile <- download_data('billing19.proteingroups.txt')
#' phosphofile <- download_data('billing19.phosphosites.txt')
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' prodt <- .read_proteingroups(proteinfile = proteinfile)
#' fosdt <- .read_phosphosites(phosphofile = phosphofile, proteinfile = proteinfile)
#' prodt %<>% maxquant_curate()
#' prodt %<>% add_feature_id()
#' quantity <- guess_maxquant_quantity(proteinfile)
#' pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
#' mqdt_to_mat(prodt, pattern = pattern)[1:2, 1:2]
#' @export
mqdt_to_mat <- function(dt, pattern, verbose = TRUE){
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
#' @param fastadt       fasta data.table or NULL
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
#' # fukuda20
#'     proteinfile <- download_data('fukuda20.proteingroups.txt')
#'     pro <- read_proteingroups(proteinfile, plot = TRUE)
#' # billing19
#'     proteinfile <- download_data('billing19.proteingroups.txt')
#'     phosphofile <- download_data('billing19.phosphosites.txt')
#'     fastafile <- download_data('uniprot_hsa_20140515.fasta')
#'     fastadt <- read_fastahdrs(fastafile)
#'     subgroups <- sprintf('%s_STD', c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'))
#'     pro <- read_proteingroups(proteinfile, subgroups = subgroups, plot = TRUE)
#'     pro <- read_proteingroups(proteinfile, subgroups = subgroups)
#'     fos <- read_phosphosites( phosphofile, proteinfile, subgroups = subgroups)
#'     fos <- read_phosphosites( phosphofile, proteinfile, fastadt = fastadt, subgroups = subgroups)
#' @export
read_proteingroups <- function(
    proteinfile, fastadt = NULL, 
    quantity = guess_maxquant_quantity(proteinfile), curate = TRUE,
    subgroups = NULL, invert = character(0),
    contaminants = FALSE, reverse = FALSE, impute = FALSE,
    plot = FALSE, pca = plot, fit = if (plot) 'limma' else NULL,
    formula = NULL, block = NULL, coefs = NULL, contrastdefs = NULL,
    feature_id = NULL, sample_id = NULL, palette = NULL, verbose = TRUE
){
# Assert
    assert_is_a_string(quantity)
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS_QUANTITY))
    assert_is_a_bool(verbose)
# Read/Curate
    prodt <- .read_proteingroups(proteinfile = proteinfile, verbose = verbose)
    if (curate)  prodt %<>% curate_annotate(fastadt = fastadt, verbose = verbose)
    prodt %<>% add_feature_id()
# SumExp
    if (verbose)  message('\tCreate SummarizedExperiment')
    pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
    promat <- mqdt_to_mat(prodt, pattern, verbose = verbose)
    pepcols <- names(prodt) %>% extract(stri_detect_fixed(., 'eptides'))
    pepdt <- prodt[, pepcols, with = FALSE]
    prodt %<>% extract(, names(prodt) %>% setdiff(colnames(promat)) %>% setdiff(names(pepdt)), with = FALSE)
    object <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(log2proteins = promat), 
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
    phosphofile, proteinfile, fastadt = NULL, 
    quantity = guess_maxquant_quantity(proteinfile), curate = TRUE, 
    subgroups = NULL, invert = character(0), 
    contaminants = FALSE, reverse = FALSE, localization = 0.75, 
    impute = FALSE, plot = FALSE, pca = plot, fit = if (plot) 'limma' else NULL,  
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
    if (curate)  fosdt %<>% curate_annotate(fastadt = fastadt, verbose = verbose)
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
#' @param x character vector
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
    if (!is.null(subgroups)){
        assert_is_subset(subgroups, as.character(object$subgroup))
        object %<>% filter_samples(subgroup %in% subgroups, verbose = verbose)
    }
    object %<>% invert_subgroups(invert)
# Features
    if (verbose) message('\tFilter features')
    if (!reverse)      object %<>% filter_features(               Reverse == '', verbose = verbose)
    if (!contaminants) object %<>% filter_features(`Potential contaminant`== '', verbose = verbose)
    object %<>% rm_missing_in_all_samples(verbose = verbose)
    object %<>% filter_exprs_replicated_in_some_subgroup(verbose = verbose)
    if ('Localization prob' %in% fvars(object)){
        object %<>% filter_features(
            `Localization prob` >= localization, verbose = verbose)  }
# Impute
    if (impute)  object %<>% impute_consistent_nas(plot = FALSE)
    object
}


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
#' object <- read_phosphosites(phosphofile, proteinfile)
#' fdt(object)
#' object %<>% add_psp()
#' fdt(object)
#' @export
add_psp <- function(
    object, 
    pspfile = file.path(R_user_dir('autonomics', 'cache'), 
            'phosphositeplus', 'Phosphorylation_site_dataset.gz')
){
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
    object %<>% merge_fdata(dt, by.x = c('uniprot1', 'Position1'), 
                                by.y = c('ACC_ID', 'MOD_RSD'))
    idx <- fdt(object)[, !is.na(psp) &  AA != `Amino acid`]
    fdt(object)$psp[idx] <- NA
    fdt(object)$AA <- NULL
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
#                       log2proteins
#
#---------------------------------------------------------------------------


#' Get/Set log2proteins
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
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
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
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
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
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
