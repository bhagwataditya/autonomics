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
    for (quantity in names(MAXQUANT_PATTERNS)){
        pattern <- MAXQUANT_PATTERNS[[quantity]]
        if (any(stri_detect_regex(x, pattern)))   return(quantity)
    }
    stop('quantity could not be infered')
}


#' @export
}


#' @export
}


#---------------------------------------------------------------------------
# 
#                      drop_differing_uniprots
#


#' @examples
#' @export
}

#' @export


}

extract_reviewed <- function(fastahdrs){
    fastahdrs                    %>%   # seqinr::read.fasta          gives '>sp|tr' fastahdrs
    split_extract_fixed('|',1)         %>%   # Biostrings::readAAStringSet gives  'sp|tr' fastahdrs
    substr(nchar(.)-1, nchar(.))       %>%   # this function now works with both scenarios
    equals('sp')                       %>%
    as.integer()
}

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

extract_uniprot <- function(fastahdrs){
    fastahdrs                    %>%
    split_extract_fixed(' ', 1)        %>% 
    split_extract_fixed('|', 2) 
}

extract_canonical <- function(fastahdrs){
    fastahdrs                    %>%
    extract_uniprot()                  %>%
    split_extract_fixed('-', 1) 
}

extract_isoform <- function(fastahdrs){
    fastahdrs                    %>%
    extract_uniprot()                  %>%
    split_extract_fixed('-', 2)        %>%
    nastring_to_0()                    %>%
    as.integer()
}

extract_description <- function(fastahdrs){
    fastahdrs                    %>%
    split_extract_regex('_[A-Z]+ ', 2) %>% 
    split_extract_regex(' OS=', 1)
}

extract_fragment <- function(fastahdrs){
    fastahdrs                    %>%
    extract_description()              %>%
    stri_detect_fixed(  'ragment')     %>%
    as.integer()
}

extract_existence <- function(fastahdrs){
    fastahdrs                   %>%
    split_extract_fixed('PE=', 2)     %>%
    split_extract_fixed(' ', 1)       %>%
    nastring_to_nachar()              %>%
    as.integer()
}

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


#' Read headers from uniprot fastafile
#'
#' @param fastafile    string (or charactervector)
#' @param fastafields  charactervector : which fastahdr fields to extract ?
#' @param verbose      bool
#' @examples
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


#-----------------------------------------------------------------------------
#
#                   maxquant_curate
#                   fasta_curate
#                       .drop_inferior
#
#-----------------------------------------------------------------------------

drop_inferior <- function(anndt, verbose = TRUE){
    # Assert
    canonical <- existence <- fragment <- isoform <- protein <- reviewed <- NULL
    
    anndt %<>% copy()
    idcol <- if ('fosId' %in% names(anndt)) 'fosId' else 'proId'
    
    if (verbose)  message('\t\t\tDrop instances with NA protein entries') #/ proteinname')
    anndt <- anndt[    protein   != 'NA']
    #anndt <- anndt[proteinname != 'NA']
    
    if (verbose)  message('\t\t\tWithin ', idcol, ': drop trembl    in favour of swissprot')
    anndt <- anndt[, .SD[ reviewed == max(reviewed) ], by = idcol]
    
    if (verbose)  message('\t\t\t       ','     ','  drop fragments in favour of full proteins')
    anndt <- anndt[, .SD[ fragment == min(fragment) ], by = idcol]
    
    if (verbose)  message('\t\t\t       ','     ','  drop worse existences', 
                          ' (1=protein, 2=transcript, 3=homolog, 4=prediction, 5=uncertain)')
    if (!any(is.na(anndt$existence))){
        anndt <- anndt[, .SD[existence == min(existence)], by = idcol] }
    anndt[, c('reviewed', 'fragment', 'existence') := NULL]
    
    #if (verbose)  message("\t\t\tDrop 'Isoform x of ' in Protein names")
    #anndt$proteinname %<>% stri_replace_first_regex('Isoform [0-9A-Z]+ of ', '')
    
    #if (verbose)  message("\t\t\tDrop '(Fragment)'    in Protein names")
    #anndt$proteinname %<>% stri_replace_first_fixed(' (Fragment)', '')
    if (is_numeric_string(anndt[[idcol]][1])){
        anndt[order(as.integer(get(idcol)), canonical, isoform)]
    } else {
        anndt[order(get(idcol), canonical, isoform)]
    }
}


}

}

CURATEDCOLS <- c('protein', 'isoform', 'uniprot', 'canonical', 'gene', 'organism')

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
#'   3. Annotate: protein, isoform, gene
#'         * `maxquant_curate :  MaxQuant `
#'         * `   fasta_curate :  Fastafile`
#'   4. Recollapse into `Curated`
#' @param dt      `data.table`
#' @param fastadt `data.table`
#' @param verbose `TRUE / FALSE`
#' @return data.table
#' @examples
#' # Fukuda 2020: MaxQuant 
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     dt <- .read_maxquant_proteingroups(file)
#'     curate_annotate_maxquant(dt)[, 1:7]
#'     curate_annotate_fastafile(dt, fastadt = NULL)[, 1:2]
#'     curate_annotate(dt, fastadt = NULL)[, 1:7]
#'     
#' # Billing 2019: Fastafile + MaxQuant
#'     file <- download_data('billing19.proteingroups.txt')
#'     phosphofile <- download_data('billing19.phosphosites.txt')
#'     fastafile   <- download_data('uniprot_hsa_20140515.fasta')
#'     fastadt <-  read_fastahdrs(fastafile)
#'     dt <- .read_maxquant_proteingroups(file, verbose = TRUE)
#'     dt[, 1:4]
#'     curate_annotate_maxquant( dt)[, 1:7]
#'     curate_annotate_fastafile(dt, fastadt = NULL)[, 1:7]
#'     curate_annotate_fastafile(dt, fastadt)[, 1:7]
#'     curate_annotate(dt, fastadt)[, 1:7]
#' @md
#' @export
curate_annotate <- function(dt, fastadt = NULL, verbose = TRUE){
    contaminant <- fastahdrs <- NULL
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    dt %<>% copy()
    dt %<>% curate_annotate_fastafile(fastadt, verbose = verbose)
    dt[, contaminant := as.character(contaminant)]  # deal with 'all contaminants NA' case
    dt[is.na(contaminant), contaminant := '']
    idx <- dt$reverse == '' & dt$contaminant == '' & is.na(dt$protein)
    dt1 <- dt[ idx] %>% curate_annotate_maxquant(verbose = verbose)
    dt  <- dt[!idx]
    dt %<>% rbind(dt1)
    dt %<>% extract(order(as.integer(get(idcol))))
    dt[, fastahdrs := NULL]
    dt[]
}


#' @rdname curate_annotate
#' @export
curate_annotate_fastafile <- function(dt, fastadt, verbose = TRUE){
# Assert
    assert_fastadt_or_null(fastadt)
    dt %<>% copy()
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    if (is.null(fastadt)){
        dt[, (CURATEDCOLS) := NA_character_]
        dt %<>% pull_columns(c(idcol, CURATEDCOLS))
        return(dt[])
    }
    uniprot <- Original <- NULL
# Drop Curated / contaminants / reverse
    idxdrop <- dt$contaminant=='+'  |  dt$reverse=='+'
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
    anndt %<>% drop_inferior(verbose = verbose)
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
# Drop Curated / contaminants / reverse
    dt %<>% copy()
    idxdrop <- dt$contaminant=='+'  |  dt$reverse=='+'
    if (sum(idxdrop)==length(idxdrop))  return(dt)
    dropdt <- dt[ idxdrop]
    dt     <- dt[!idxdrop]
    fastahdrs <- ok <- Original <- uniprot <- NULL
# Curate
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'           # Maxquant doesnt provide protein entries
    if (verbose)  message('\t\tMaxQuant Fastahdrs')                     # They can be inferred from the fastahdrs
    if (verbose)  message('\t\t\tUncollapse, Drop truncated, Parse')    # These fastahdrs are sometimes truncated
    anndt <- dt[, c(idcol, 'fastahdrs'), with = FALSE]
    anndt %<>% uncollapse(fastahdrs, sep = ';')
    
    anndt[, ok := stri_count_fixed(fastahdrs, '|')]                            # prefer non-truncated fastahdrs: starting with 'tr|F1Q6Z9|F1Q6Z9_DANRE'
    anndt <- anndt[, .SD[ok==max(ok)], by = idcol]
    anndt[, ok := NULL]
    
    anndt[, ok := fastahdrs %>% stri_count_regex('.+SV=[0-9]+'), by = idcol]   # prefer non-trnuncated fastahdrs: ending with SV=1
    anndt <- anndt[, .SD[ok==max(ok)], by = idcol]
    anndt[, ok := NULL]
    
    anndt %<>% cbind(parse_fastahdrs(anndt$fastahdrs))
    anndt %<>% drop_inferior(verbose = verbose)
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

#' Add feature_id
#' @param dt       data.table
#' @return data.table
#' @examples 
#' file <- download_data('billing19.proteingroups.txt')
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' prodt <- .read_maxquant_proteingroups(file = file)
#' prodt %<>% curate_annotate_maxquant()
#' prodt %<>% add_feature_id()
#' @export
add_feature_id <- function(dt){
# Initialize
    `Amino acid`                <- contaminant <- isoform <- protein <- NULL
    `Positions within proteins` <- reverse     <- NULL
# Add
    dt %<>% copy()
    dt[, contaminant := as.character(contaminant)]   # one dataset had NA_logical values
    dt[is.na(contaminant), contaminant := '']        # which made all features being filtered out in the next lines
    dt1 <- dt[reverse ==''  & contaminant == '']
    dt2 <- dt[reverse =='+' & contaminant == '' ]
    dt3 <- dt[reverse ==''  & contaminant == '+']
    dt4 <- dt[reverse =='+' & contaminant == '+']
    
    idcol <- if ('fosId' %in% names(dt1)) 'fosId' else 'proId'
    dt1[, feature_id := protein]
    if (length(unique(dt1$organism)) == 1)  dt1$feature_id %<>% stri_replace_all_regex('_[^;]+', '')
    dt1[isoform!=0, feature_id := paste0(feature_id, '(', stri_replace_all_fixed(isoform, ';', ''), ')')]
    if (idcol=='fosId'){
        dt1[, feature_id := paste0(feature_id, '-', `Amino acid`) ]
        dt1[, feature_id := paste0(feature_id, split_extract_fixed(`Positions within proteins`, ';', 1)) ]
    }
    dt2[, feature_id := paste0('REV_',     get(idcol))]
    dt3[, feature_id := paste0('CON_',     get(idcol))]
    dt4[, feature_id := paste0('REV_CON_', get(idcol))]
        # Contaminants and actual proteins get mixed up. For naming: make sure to use contaminant.
    dt <- dt1 %>% rbind(dt2) %>% rbind(dt3) %>% rbind(dt4)
    dt %<>% pull_columns(c(idcol, 'feature_id'))
    assert_all_are_non_missing_nor_empty_character(dt$feature_id)
    assert_has_no_duplicates(dt$feature_id)
    dt[]
}


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
    if (verbose) message('\t\tStandardize snames: ', x[1], '  ->  ', cleanx[1])
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
#' @param fastadt       NULL or data.table
#' @param quantity     'normalizedratio', 'ratio', 'correctedreporterintensity', 
#'                     'reporterintensity', 'maxlfq', 'labeledintensity', 
#'                     'intensity' or NULL
#' @param subgroups     NULL or string vector : subgroups to retain
#' @param contaminants  TRUE or FALSE : retain contaminants ?
#' @param reverse       TRUE or FALSE : include reverse hits ?
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
    object, subgroups, invert, contaminants, reverse, localization = 0.75, 
    impute, verbose
){
# Demultiplex. Infer Subgroup
    contaminant <- `Localization prob` <- NULL
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
    analysis(object)$nfeatures <- c(nrow(object))
    if (!reverse)       object %<>% filter_features(reverse == '', verbose = verbose)
    if (!contaminants)  object %<>% filter_features(contaminant== '', verbose = verbose)
    object %<>% rm_missing_in_all_samples(verbose = verbose)
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
        resdt <- resdt[, lapply(.SD, trimws),                          by = 'From']
    # Collapse per uniprot
        resdt <- resdt[, lapply(.SD, paste_unique, collapse=collapse), by = 'From']
    # Ensure original order
        merge.data.table(
            fdt, resdt, by.x = 'uniprot', by.y = 'From', all.x = TRUE, sort = FALSE)
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
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     x <- read_maxquant_proteingroups(file, plot = FALSE)
#'     x %<>% extract(1:10, )
#'     fdt(x)[1:3, ]
#'     # x %<>% annotate_uniprot_ws(upws)
#'     # fdt(x)[1:3, ]
#' @export
annotate_uniprot_ws <- function(x, ...)  UseMethod('annotate_uniprot_ws')

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
#' @rdname annotate_uniprot_ws
#' @export
annotate_uniprot_ws.data.table <- function(
    x, upws, columns=c('xref_ensembl'), collapse=';', ...
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
        x <- x[, .annotate_uniprot_ws(.SD, upws, columns = columns, collapse = collapse), by = 'chunk']
    # Return
        x$chunk <- NULL
        x[]
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
#' @rdname annotate_uniprot_ws
#' @export
annotate_uniprot_ws.SummarizedExperiment <- function(
    x, upws, columns = c('ENSEMBL'), collapse=';', ...
){
    # Assert
        assert_is_all_of(upws, 'UniProt.ws')
    # Split proteingroups into proteins and extract uniprot
        fdt <- data.table(fdata(x))
        fdt %<>% uncollapse(uniprot, sep=collapse)
        fdt[, uniprot := split_extract_fixed(uniprot, '-', 1)]
    # Annotate
        uniprot <- NULL
        # fdt %<>% extract(uniprot %in% UniProt.ws::keys(upws, keytype = 'UniProtKB'))  # latest version returns only small number of keys!
        fdt %<>% annotate_uniprot_ws.data.table(upws, columns = columns, collapse = collapse)
    # Collapse
        fdt <- fdt[, lapply(.SD, trimws),                          by = 'feature_id'] # trim whitespace
        fdt <- fdt[, lapply(.SD, paste_unique, collapse=collapse), by = 'feature_id'] # collapse
        x %<>% merge_fdt(fdt)
        x
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

