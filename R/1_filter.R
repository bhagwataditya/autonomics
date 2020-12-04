
#===================================================
# EXTRACT FEATURES
#===================================================

#' Extract features
#' @param   object SummarizedExperiment
#' @param   extractor logical/numeric vector
#' @return  SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' (object %<>% extract_features(c(5,4)))
#' @export
extract_features <- function(object, extractor){
    object %<>% extract(extractor, )
    if (!is.null(limma(object))){
        limma(object) %<>% extract(fnames(object), , , drop = FALSE)
    }
    object
}

#' Extract first from collapsed values
#' @param x charactervector or factorvector
#' @param sep collapsed string separator, e.g. ';'
#' @param ... to allow for S3 method dispatch
#' @return Updated x
#' @examples
#' x <- c('a;b;c', '1;2;3', 'alpha;beta;gamma')
#' extract_first_from_collapsed(x, sep = ';')
#' @noRd
extract_first_from_collapsed <- function (x, ...) {
    UseMethod("extract_first_from_collapsed", x)
}

extract_first_from_collapsed.character <- function(x, sep = guess_sep(x), ...){
    if (is.null(sep)) return(x)

    stringi::stri_split_fixed(x, sep) %>%
    vapply(extract, character(1), 1)
}

extract_first_from_collapsed.factor <- function(x, sep = guess_sep(x), ...){
    levels(x) %<>% extract_first_from_collapsed.character(sep=sep)
    x
}



#=================
# FILTER FEATURES
#=================

#' Filter features on condition
#' @param object SummarizedExperiment
#' @param condition filter condition
#' @param verbose logical
#' @return filtered eSet
#' @examples
#' # GLUTAMINASE
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     filter_features(object,   SUPER_PATHWAY=='Lipid',  verbose = TRUE)
#' @export
filter_features <- function(object, condition, verbose = FALSE){
    condition <- enquo(condition)
    idx <- eval_tidy(condition, fdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) message('\t\tRetain ',
        sum(idx), '/', length(idx), ' features: ', expr_text(condition))
    object %<>% extract(idx,)
    fdata(object) %<>% droplevels()
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>%
            c(structure(sum(idx), names = quo_name(condition)))
    }
    object
}


#' Rm features missing in all samples
#' @param object  SummarizedExperiment
#' @param verbose TRUE (default) or FALSE
#' @return  filtered SummarizedExperiment
#' @noRd
rm_missing_in_all_samples <- function(object, verbose = TRUE){
    # . != 0 needed due to stupid behaviour of rowAnys
    # https://github.com/HenrikBengtsson/matrixStats/issues/89
    selector <- rowAnys(exprs(object) != 0, na.rm = TRUE)
    if (verbose) message(
                    '\t\tRetain ', sum(selector), '/', length(selector),
                    ' features: non-zero, non-NA, and non-NaN for some sample')
    object %<>% extract(selector, )
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>% c(structure(sum(selector),
            names = "non-zero, non-NA, and non-NaN for some sample"))
    }
    object
}

is_available_in_all_samples <- function(object)  rowAlls(!is.na(exprs(object)))


#' Keep features that are available in all samples
#' @param object SummarizedExperiment
#' @param verbose TRUE (default) or FALSE
#' @return updated object
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' rm_missing_in_some_samples(object)
#' @noRd
rm_missing_in_some_samples <- function(object, verbose = TRUE){

    # Restrict to available values
    selector <- is_available_in_all_samples(object)
    if (verbose)  message('\t\tUse ', sum(selector), '/', length(selector),
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
#' @param object      SummarizedExperiment
#' @param group       group svar
#' @param comparator  '>' or '!='
#' @param lod         number: limit of detection
#' @param verbose     TRUE or FALSE
#' @return Filtered SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' invert_subgroups <- c('E_EM', 'BM_EM', 'E_BM')
#' object <- read_proteingroups(file, invert_subgroups=invert_subgroups)
#' object %<>% filter_exprs_replicated_in_some_subgroup()
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% filter_exprs_replicated_in_some_subgroup()
#'
#' filter_exprs_replicated_in_some_subgroup(object, character(0))
#' @export
filter_exprs_replicated_in_some_subgroup <- function(
    object, group = 'subgroup',
    comparator = if (contains_ratios(object)) '!=' else '>',
    lod = 0, verbose = TRUE
){
# Assert
    assert_is_subset(group, svars(object))
# Datatablify
    replicated_in_its_subgroup <- replicated_in_any_subgroup <- value <- NULL
    dt <- sumexp_to_long_dt(object, svars = group)
# Find replicated features
    exceeds_lod <- if (comparator == '>'){ function(value, lod) value >  lod
            } else if (comparator == '!=') function(value, lod) value != lod
    V1 <- dt[,.I[sum(exceeds_lod(value, lod), na.rm=TRUE)>1],
            by = c('feature_id', group)]$V1
            #https://stackoverflow.com/questions/16573995
# Keep only replicated features
    replicated_features <- dt[V1]$feature_id
    idx <- fid_values(object) %in% replicated_features
    if (verbose)   message('\t\tFilter ', sum(idx), '/', length(idx),
            ' features: expr ', comparator, ' ', as.character(lod),
            ' for at least two samples in some ', group)
    object %<>% extract_features(idx) # also handles limma in metadata
# Update analysis log
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>% c(structure(
                sum(idx),
                names = sprintf(
                    "expr %s %s, for at least two samples in some %s",
                    comparator, as.character(lod), group)))
    }
    object
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
#' # GLUTAMINASE
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     filter_samples(object,   TIME_POINT=='h10',  verbose = TRUE)
#' @export
filter_samples <- function(object, condition, verbose = FALSE, record = TRUE){
    condition <- enquo(condition)
    idx <- eval_tidy(condition, sdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx),
                                        ' samples: ', expr_text(condition))
    object %<>% extract(, idx)
    sdata(object) %<>% droplevels()
    if (record && !is.null(analysis(object))) {
        analysis(object)$nsamples   %<>%  c(structure(sum(idx),
                                            names = expr_text(condition)))
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
        if (verbose) cmessage(
            '\t\tRetain %d/%d samples with a value available for some feature',
            sum(subsetter), length(subsetter))
    object %<>% extract(, subsetter)
    }
    object
}

is_available_for_some_feature <- function(object){
    subsetter <- (!is.na(exprs(object))) & (exprs(object) != 0)
    set_names(colAnys(subsetter), snames(object))
}







