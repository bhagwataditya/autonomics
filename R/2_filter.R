
#=================
# FILTER FEATURES
#=================

#' Filter features on condition
#' @param object SummarizedExperiment
#' @param condition filter condition
#' @param verbose logical
#' @return filtered eSet
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' filter_features(object, SUPER_PATHWAY == 'Lipid')
#' @export
filter_features <- function(object, condition, verbose = TRUE){
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
}


#' @rdname rm_missing_in_some_samples
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


#' Rm features missing in some samples
#' @param object SummarizedExperiment
#' @param verbose TRUE (default) or FALSE
#' @return updated object
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' rm_missing_in_all_samples( object)
#' rm_missing_in_some_samples(object)
#' @export
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
#' @param assay        string
#' @param comparator   '>' or '!='
#' @param lod          number: limit of detection
#' @param nsample      number
#' @param nsubgroup    number
#' @param verbose      TRUE or FALSE
#' @return Filtered SummarizedExperiment
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
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
    replicated <- NULL
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


#' Keep replicated features
#'
#' Keep features replicated for each slevel
#'
#' @param object   SummarizedExperiment
#' @param formula  formula
#' @param n        min replications required
#' @param verbose  TRUE or FALSE
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% keep_replicated_features()
#' object %<>% keep_replicated_features(~ subgroup)
#' @export
keep_replicated_features <- function(
    object, formula = ~ 1, n = 3, verbose = TRUE
){
# Drop NA values
    dt <- sumexp_to_longdt(object, svars = c(all.vars(formula)), assay = assayNames(object)[1])
    n0 <- length(unique(dt$feature_id))
    value <- NULL
    dt %<>% extract(!is.na(value))
    dt %<>% extract(, .SD[.N>=n], by = 'feature_id')
    n1 <- length(unique(dt$feature_id))
    if (n1<n0 & verbose)  cmessage('\t\tKeep %d/%d features with %d+ values', n1, n0, n)
# Feature covers each slevel
    for (var in all.vars(formula)){
        # must span all slevels
        nlevels <- length(unique(dt[[var]]))
        n0 <- length(unique(dt$feature_id))
        dt %<>% extract(, .SD[length(unique(get(var))) == nlevels], by = 'feature_id')
        n1 <- length(unique(dt$feature_id))
        if (n1<n0 & verbose)  cmessage('\t\tKeep %d/%d features spanning all %s levels', n1, n0, var)

        # must have n+ obs per slevel
        n0 <- length(unique(dt$feature_id))
        dt %<>% extract(, .SD[.N>=n], by = c('feature_id', var))
        n1 <- length(unique(dt$feature_id))
        if (n1<n0 & verbose)  cmessage('\t\tKeep %d/%d features with %d+ values per %s', n1, n0, n, var)
    }
# Return
    idx <- fnames(object) %in% as.character(unique(dt$feature_id))
    object[idx, ]
}




#=======================
# FILTER SAMPLES
#=======================


#' Filter samples on condition
#' @param object    SummarizedExperiment
#' @param condition filter condition
#' @param verbose   TRUE or FALSE (default)
#' @param record    TRUE (default) or FALSE
#' @return filtered SummarizedExperiment
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' filter_samples(object, Group != 't0', verbose = TRUE)
#' @export
filter_samples <- function(object, condition, verbose = FALSE, record = TRUE){
    condition <- enquo(condition)
    idx <- eval_tidy(condition, sdata(object))
    idx <- idx & !is.na(idx)
    if (verbose & sum(idx)<length(idx))  message('\t\t\tRetain ', sum(idx), '/',
                                length(idx), ' samples: ', expr_text(condition))
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


#' Rm singleton samples
#'
#' @param object  SummarizedExperiment
#' @param svar    sample var
#' @param verbose TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' object %<>% filter_samples(SampleGroup %in% c('t1', 't2'), verbose = TRUE)
#' rm_singleton_samples(object, 'Subject_ID')
#' @export
rm_singleton_samples <- function(object, svar = 'subgroup', verbose = TRUE){
    selectedsamples <-
        data.table(sdata(object))[, .SD[.N>1], by = svar][['sample_id']]
    if (verbose)   message('\t\t\tRetain ', length(selectedsamples), '/', 
                            ncol(object), ' samples with replicated ', svar)
    object[, selectedsamples]
}


