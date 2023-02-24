
#=================
# FILTER FEATURES
#=================

#' Filter features on condition
#' @param object     SummarizedExperiment
#' @param condition  expression
#' @param verbose    logical
#' @param ...        required for s4 dispatch
#' @return filtered eSet
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' filter_features(object, SUPER_PATHWAY == 'Lipid')
#' @export
setGeneric('filter_features', function(object, ...)  standardGeneric('filter_features'))

#' @rdname filter_features
#' @export
setMethod('filter_features', signature(object = 'SummarizedExperiment'), 
function(object, condition, verbose = TRUE){
    . <- NULL
    condition <- enquo(condition)
    idx <- eval_tidy(condition, fdata(object))
    idx <- idx & !is.na(idx)
    if (verbose & sum(idx) < length(idx)){
        message('\t\tRetain ', sum(idx), '/', length(idx), 
                ' features: ', expr_text(condition) %>% substr(1, min(120, nchar(.))))}
    object %<>% extract(idx,)
    fdata(object) %<>% droplevels()
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>%
            c(structure(sum(idx), names = quo_name(condition)))
    }
    object
})

#' @rdname filter_features
#' @export
setMethod('filter_features', signature(object = 'MultiAssayExperiment'), 
function(object, condition, verbose = TRUE){
    condition <- enquo(condition)
    for (ass in names(object)){
        object[[ass]] %<>% filter_features(!!condition, verbose = verbose)
    }
    object
})

#' Rm features missing in all samples
#' @param object  SummarizedExperiment
#' @param verbose TRUE (default) or FALSE
#' @return  filtered SummarizedExperiment
#' @export
rm_missing_in_all_samples <- function(object, verbose = TRUE){
    # . != 0 needed due to stupid behaviour of rowAnys
    # https://github.com/HenrikBengtsson/matrixStats/issues/89
    selector <- rowAnys(values(object) != 0, na.rm = TRUE)
    if (verbose && sum(selector)<length(selector)){
        message('\t\tRetain ', sum(selector), '/', length(selector),
                ' features: non-zero, non-NA, and non-NaN for some sample')
        object %<>% extract(selector, )
        if (!is.null(analysis(object))) {
            analysis(object)$nfeatures %<>% c(structure(sum(selector),
                names = "non-zero, non-NA, and non-NaN for some sample"))
        }
    }
    object
}

is_available_in_all_samples <- function(object)  rowAlls(!is.na(values(object)))


#' Keep features that are available in all samples
#' @param object SummarizedExperiment
#' @param verbose TRUE (default) or FALSE
#' @return updated object
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' rm_missing_in_some_samples(object)
#' @noRd
rm_missing_in_some_samples <- function(object, verbose = TRUE){

    # Restrict to available values
    selector <- is_available_in_all_samples(object)
    if (verbose)  message('\t\t\tUse ', sum(selector), '/', length(selector),
                            ' features with available value for each sample')
    object %<>% extract(selector, )
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>% c(structure(sum(selector),
            names = "available value for each sample"))
    }
    object
}


#==================
# FILTER EXPRS
#==================

#' Filter features with replicated expression in some subgroup
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup svar
#' @param comparator   '>' or '!='
#' @param lod          number: limit of detection
#' @param verbose      TRUE or FALSE
#' @return Filtered SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% filter_exprs_replicated_in_some_subgroup()
#' filter_exprs_replicated_in_some_subgroup(object, character(0))
#' filter_exprs_replicated_in_some_subgroup(object, NULL)
#' @export
filter_exprs_replicated_in_some_subgroup <- function(
    object, subgroupvar = 'subgroup', assay = assayNames(object)[1],
    comparator = if (contains_ratios(object)) '!=' else '>',
    lod = 0, nsample = 2, nsubgroup = 1, verbose = TRUE
){
# Assert
    assert_is_subset(subgroupvar, svars(object))
# Datatablify
    replicated_in_its_subgroup <- replicated_in_any_subgroup <- value <- NULL
    dt <- sumexp_to_longdt(object, svars = subgroupvar, assay = assay)
# Find replicated features
    exceeds_lod <- if (comparator == '>'){ function(value, lod) value >  lod
            } else if (comparator == '!=') function(value, lod) value != lod
    
    condition <- sprintf('value %s %s', comparator, lod)
    repfeatures <- dt[,  .(replicated = .(sum(eval(parse(text = condition)), na.rm = TRUE) >= nsample)) , by = c('feature_id', subgroupvar)]
    repfeatures %<>% extract( , .(replicated = sum(as.numeric(replicated)) >= 1), by = 'feature_id')
    repfeatures %<>% extract(replicated == TRUE)
    repfeatures %<>% extract2('feature_id')
    repfeatures %<>% as.character()
# Keep only replicated features
    idx <- fid_values(object) %in% repfeatures
    if (verbose)  if (any(!idx))  cmessage('\t\tFilter %d/%d features: %s %s %s for at least %d samples in %d %s', sum(idx), 
            length(idx), assay, comparator, as.character(lod), nsample, nsubgroup, subgroupvar)
    object %<>% extract(idx, )
# Update analysis log
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>% c(structure(
                sum(idx),
                names = sprintf(
                    "expr %s %s, for at least two samples in some %s",
                    comparator, as.character(lod), subgroupvar)))
    }
    object
}

#' Filter for replicated features
#' @param object     SummarizedExperiment
#' @param comparator string
#' @param lod        number: limit of detection
#' @param n          number: number of replicates above lod
#' @param verbose    TRUE/FALSE
#' @return  SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file, plot = FALSE)
#' object %<>% filter_replicated()
#' @export
filter_replicated  <- function(
    object, comparator = `>`, lod=0, n=2, verbose=TRUE
){
    assert_is_all_of(object, "SummarizedExperiment")
    assert_is_function(comparator)
    assert_is_a_number(lod)
    assert_is_a_number(n)

    nreplicates <- rowSums(comparator(values(object), lod), na.rm=TRUE)
    idx <- nreplicates >= n
    if (verbose)  if (!any(idx))  message(
        '\t\t\tRetain ', sum(idx), '/', length(idx), 
        'features replicated in at least ', n, ' samples')
    object[idx, ]
}




#=======================
# FILTER SAMPLES
#=======================


#' Filter samples on condition
#' @param object    SummarizedExperiment
#' @param condition filter condition
#' @param verbose   TRUE/FALSE 
#' @param record    TRUE/FALSE 
#' @return filtered SummarizedExperiment
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' filter_samples(object, subgroup != 't0', verbose = TRUE)
#' @export
filter_samples <- function(object, condition, verbose = TRUE, record = TRUE){
    . <- NULL
    condition <- enquo(condition)
    idx <- eval_tidy(condition, sdata(object))
    idx <- idx & !is.na(idx)
    if (verbose & sum(idx)<length(idx)){
        message('\t\t\tRetain ', sum(idx), '/', length(idx), ' samples: ', 
                expr_text(condition) %>% substr(1, min(120, nchar(.))))}
    object %<>% extract(, idx)
    sdata(object) %<>% droplevels()
    if (record && !is.null(analysis(object))) {
        analysis(object)$nsamples %<>%  
            c(structure(sum(idx), names = expr_text(condition)))
    }
    object
}


#' Filter samples available for some feature
#' @param object SummarizedExperiment
#' @param verbose TRUE or FALSE
#' @return SummarizedExperiment
#' @noRd
filter_samples_available_for_some_feature <- function(object, verbose = FALSE){
    subsetter <- is_available_for_some_feature(object)
    if (any(!subsetter)){
        if (verbose)  message('\t\t\tRetain ', sum(subsetter), '/', 
                            length(subsetter),
                            ' samples with a value available for some feature')
    object %<>% extract(, subsetter)
    }
    object
}

is_available_for_some_feature <- function(object){
    subsetter <- (!is.na(values(object))) & (values(object) != 0)
    set_names(colAnys(subsetter), snames(object))
}



