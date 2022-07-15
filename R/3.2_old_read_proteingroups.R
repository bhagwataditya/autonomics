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



#============================================================================
#
#                        old_simplify_proteingroups()
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
#' object %<>% old_simplify_proteingroups(fastafile)
#' fdata(object)[1:5, ]
#' @noRd
old_simplify_proteingroups <- function(
    object, fastafile, verbose = TRUE
){
# Return if NULL
    if (is.null(fastafile)) return(object)
# Uncollapse fdata annotations
    fasta_dt <- read_fastafile_headers(fastafile)
    feature_dt  <-  fdata(object)[, c('feature_id', 'uniprot')] %>%
                    uncollapse('uniprot', sep=';')
# Merge in fasta annotations
    feature_dt %<>% merge(fasta_dt,  by.x = 'uniprot', by.y='Uniprot',
                        sort=FALSE, all.x=TRUE)
# Simplify
    if (verbose) message('\t\tSimplify proteingroups')
    feature_dt %<>% old_prefer_best_existence()
    feature_dt %<>% old_prefer_swissprot_over_trembl()
    feature_dt %<>% old_drop_fragments()
    feature_dt %<>% old_collapse_isoforms_paralogs()
# Merge into sumexp
    fdata(object)$uniprot <- fdata(object)$feature_name <- NULL
    fdata(object)$`Protein names` <- NULL
    fdata(object) %<>% merge(
        feature_dt, by = 'feature_id', sort = FALSE, all.x = TRUE)
    fdata(object) %<>% pull_columns(c('feature_id', 'feature_name',
                        'uniprot', 'canonical', 'Protein names'))
    object
}


old_prefer_best_existence <- function(feature_dt, verbose=TRUE){
    Existence <- NULL
    feature_dt[is.na(Existence), Existence:=5]
    feature_dt %<>% extract(, .SD[Existence == min(Existence)], 
                            by = 'feature_id')
    feature_dt[, Existence := NULL]
    if (verbose) message('\t\t\tDrop inferior existences')
    feature_dt
}

old_prefer_swissprot_over_trembl <- function(feature_dt, verbose=TRUE){
    Reviewed <- NULL
    feature_dt[is.na(Reviewed), Reviewed:=0]
    feature_dt %<>% extract(, .SD[Reviewed == max(Reviewed)], by = 'feature_id')
    if (verbose) message(
                '\t\t\tDrop trembl when swissprot available')
    feature_dt[, Reviewed := NULL]
    feature_dt
}

old_drop_fragments <- function(feature_dt, verbose=TRUE){
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

old_collapse_isoforms_paralogs <- function(feature_dt, verbose=TRUE){
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
#               old_filter_maxquant_features
#
#==============================================================================


old_rm_reverse <- function(object, verbose){
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
#' old_rm_contaminants(object, verbose=TRUE)
#' @noRd
old_rm_contaminants <- function(object, verbose){
    contaminant_var <- c('Potential contaminant', 'Contaminant')
    contaminant_var %<>% intersect(fvars(object))
    fdata(object)[[contaminant_var]] %<>% na_to_string()

    idx <- fdata(object)[[contaminant_var]]==''
    if (verbose) message('\t\tRetain ',
        sum(idx), '/', length(idx), " features: contaminant != '+'")
    object %<>% extract_features(idx)
    object

}


old_rm_unlocalized <- function(object, localization, verbose){
    `Localization prob` <- NULL
    if (!'Localization prob' %in% fvars(object)) return(object)
    assert_all_are_in_range(localization, 0, 1)
    object %<>% filter_features(`Localization prob` >= localization,
                                verbose = verbose)
    object
}

old_filter_maxquant_features <- function(
    object, reverse, contaminants, localization = 0.75,
    verbose
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_a_bool(reverse)
    assert_is_a_bool(contaminants)
    assert_all_are_in_range(localization, 0, 1)
    assert_is_a_bool(verbose)
# Filter
    if (verbose) message('\tFilter features')
    if (!reverse)      object %<>% old_rm_reverse(     verbose = verbose)
    if (!contaminants) object %<>% old_rm_contaminants(verbose = verbose)
    object %<>% rm_missing_in_all_samples(verbose = verbose)
    object %<>% filter_exprs_replicated_in_some_subgroup(verbose = verbose)
    object %<>% old_rm_unlocalized(localization, verbose = verbose)
# Return
    object
}

old_rename_proteingroup_fvars <- function(object){
    stri_rep <- stri_replace_first_fixed
    names(fdata(object)) %<>% stri_rep('Gene names', 'feature_name')
    names(fdata(object)) %<>% stri_rep('Majority protein IDs', 'uniprot')
    fdata(object) %<>% pull_columns(c('feature_id', 'feature_name', 'uniprot'))
    object
}

old_rename_phospho_fvars <- function(object){
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

filter_maxquant_samples <- function(object, subgroups, verbose){
    object %<>% filter_samples_available_for_some_feature(verbose = verbose)
    if (!is.null(subgroups))  object %<>%
        filter_samples(subgroup %in% subgroups, verbose = verbose)
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
    if (impute) object %<>% impute_consistent_nas(plot = FALSE)
    object
}

#==============================================================================
#
#                  old_subtract_proteingroups
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


old_subtract_proteingroups <- function(phosphosites, proteingroups, verbose){
# Initialize
    . <- phospho <- protein <- occupancy <- NULL
    `Protein group IDs` <- NULL
# Report
    if (verbose) message(
        '\tAdd occupancies(phospho) = values(phospho) - values(proteins)')
# phospho datatable
    cols <- c('feature_id', 'Protein group IDs')
    fosdt <- sumexp_to_wide_dt(phosphosites, fvars = 'Protein group IDs')
    fosdt %<>% uncollapse(`Protein group IDs`, sep = ';' )
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
#                       old_add_pepcounts
#                   read_proteingroups
#                   read_phosphosites
#                       phospho_expr_columns
#                       old_subtract_proteingroups
#
#==============================================================================

old_add_pepcounts <- function(object, file, pepcountpattern, quantity){
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


old_read_maxquant <- function(file, quantity = guess_maxquant_quantity(file),
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
    if (!phospho)  object %<>% old_add_pepcounts(file, pepcountpattern, quantity)
# Process
    if (phospho)  colnames(object) %<>% stri_replace_last_fixed('___1', '')
    colnames(object) %<>% dequantify(quantity)
    colnames(object) %<>% demultiplex()
    object$sample_id <- colnames(object)
    object %<>% merge_sfile(sfile = sfile, by.x = by.x, by.y = by.y)
    object %<>% add_subgroup(subgroupvar=subgroupvar, verbose=verbose)
    object %<>% filter_maxquant_samples(
                    subgroups = select_subgroups, verbose)
    values(object) %<>% zero_to_na(verbose = verbose)
    values(object) %<>% nan_to_na( verbose = verbose)
    object %<>% log2transform(verbose=verbose)
    object %<>% invert_subgroups({{invert_subgroups}})
# Return
    assayNames(object)[1] <- 'maxquant'
    metadata(object)$quantity <- quantity
    metadata(object)$file <- file
    object
}


#' Read/Analyze proteingroups/phosphosites
#'
#' @param proteinfile  proteingroups file
#' @param phosphofile  proteingroups file
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
#' @param localization min site localization probability (number)
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
#'     proteinfile <- download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(proteinfile, pca=TRUE, fit='limma')
#' # Billing 2019
#'     proteinfile <- download_data('billing19.proteingroups.txt')
#'     phosphofile <- download_data('billing19.phosphosites.txt')
#'     subgroups <- c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00')
#'     subgroups %<>% paste0('_STD')
#'     palette <- make_colors(c(subgroups[1], 'E00.8_STD', subgroups[-1]))
#'     pro <- read_proteingroups(
#'                proteinfile, select_subgroups = subgroups, palette = palette, 
#'                coefs = 'subgroupE05_STD', sample_id = 'E05_STD.R2')
#'     fos <- read_phosphosites(phosphofile, proteinfile, 
#'                select_subgroups = subgroups)
#' @export
old_read_proteingroups <- function(
    proteinfile, fastafile = NULL, 
    quantity = guess_maxquant_quantity(proteinfile), 
    select_subgroups = NULL, invert_subgroups = character(0),
    contaminants = FALSE, reverse = FALSE, 
    impute = stri_detect_regex(quantity, "[Ii]ntensity"),
    subgroupvar = NULL,
    formula = NULL, block = NULL, coefs = NULL, contrastdefs = NULL,
    pca = TRUE, fit = 'limma', verbose = TRUE, plot = pca & !is.null(fit), 
    feature_id = NULL, sample_id = NULL, palette = NULL
){
# Assert
    . <- NULL
    assert_all_are_existing_files(proteinfile)
    if (!is.null(fastafile)) assert_all_are_existing_files(fastafile)
# Read
    object <- read_maxquant(proteinfile, fastafile = fastafile, 
        select_subgroups = select_subgroups, invert_subgroups = invert_subgroups,
        quantity = quantity, verbose = verbose)
# Preprocess
    object %<>% old_filter_maxquant_features(reverse = reverse,
                    contaminants = contaminants, verbose = verbose)
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
old_read_phosphosites <- function(
    phosphofile, proteinfile = paste0(dirname(phosphofile), '/proteinGroups.txt'),
    quantity = guess_maxquant_quantity(phosphofile),
    select_subgroups = NULL, contaminants = FALSE,
    reverse = FALSE, localization = 0.75, fastafile = NULL,
    invert_subgroups = character(0), pca = TRUE,
    fit = 'limma', subgroupvar = NULL, formula = NULL, block = NULL, 
    coefs = NULL, contrastdefs = NULL, 
    verbose = TRUE, plot = pca & !is.null(fit)
){
# Assert
    . <- NULL
    `Protein group IDs` <- `Localization prob` <- NULL
    assert_all_are_existing_files(c(phosphofile, proteinfile))
    assert_is_subset(quantity, names(MAXQUANT_PATTERNS_QUANTITY))
# Read / Prepare
    object <- read_maxquant(
        proteinfile = proteinfile, phosphofile = phosphofile, 
        fastafile = fastafile, quantity = quantity, 
        select_subgroups = select_subgroups, invert_subgroups = invert_subgroups, 
        verbose = verbose)
    object %<>% old_filter_maxquant_features(
        reverse = reverse, contaminants = contaminants,
        localization = localization,
        verbose = verbose)
    object %<>% transform_maxquant(impute=FALSE, verbose=verbose, plot=plot)
# Analyze / Return
    object %<>% analyze(pca = pca, fit = fit, subgroupvar = subgroupvar, 
                        formula = formula, block = block, 
                        coefs = coefs, contrastdefs = contrastdefs,
                        verbose = verbose, plot = plot)
    object
}


