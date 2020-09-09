
#===================================================
# EXTRACT FEATURES
#===================================================

#' Extract features
#' @param object SummarizedExperiment
#' @param extractor logical/numeric vector
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
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
   vapply(magrittr::extract, character(1), 1)
}

extract_first_from_collapsed.factor <- function(x, sep = guess_sep(x), ...){
   levels(x) %<>% extract_first_from_collapsed.character(sep=sep)
   x
}



#=================
# FILTER FEATURES
#=================

#' @rdname filter_features
#' @export
filter_features_ <- function(object, condition, verbose = FALSE){
    if (is.null(condition)) return(object)
    idx <- lazy_eval(condition, fdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) message(
          '\t\tRetain ', sum(idx), '/', length(idx), ' features: ',
          if (class(condition)=='lazy') deparse(condition$expr) else condition)
    object %<>% extract_features(idx)
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>%
            c(structure(sum(idx), names = condition))
    }
    object
}

#' Filter features on condition
#' @param object SummarizedExperiment
#' @param condition filter condition
#' @param verbose logical
#' @return filtered eSet
#' @examples
#' # GLUTAMINASE
#'     file <- download_data('glutaminase.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     filter_features(object,   SUPER_PATHWAY=='Lipid',  verbose = TRUE)
#'     filter_features_(object, "SUPER_PATHWAY=='Lipid'", verbose = TRUE)
#' @export
filter_features <- function(object, condition, verbose = FALSE){
    # filter_features_(object, lazyeval::lazy(condition), verbose = verbose)
        # older version using -now deprecated- lazyeval package
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


#' Keep features that are non-zero, non-NA, and non-NaN for some sample
#' @param object  SummarizedExperiment
#' @param verbose logical
#' @return  filtered eSet
#' @export
filter_features_nonzero_in_some_sample <- function(object, verbose = TRUE){
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
#' @return updated object
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' filter_features_available_in_all_samples(object)
#' @export
filter_features_available_in_all_samples <- function(object){

    # Restrict to available values
    selector <- is_available_in_all_samples(object)
    message('\t\tUse ', sum(selector), '/', length(selector),
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
#' @param comparator  '>' or '!='
#' @param lod         number
#' @return Filtered SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' object %<>% filter_exprs_replicated_in_some_subgroup()
#'
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% filter_exprs_replicated_in_some_subgroup()
#' @export
filter_exprs_replicated_in_some_subgroup <- function(
    object,
    comparator = if (contains_ratios(object)) '!=' else '>',
    lod = 0
){
    # Return if no subgroups or replicates
    if (!'subgroup' %in% svars(object))             return(object)
    if (all(!duplicated(sdata(object)$subgroup)))   return(object)

    # Datatablify
    replicated_in_its_subgroup <- replicated_in_any_subgroup <- value <- NULL
    dt <- sumexp_to_long_dt(object, svars = 'subgroup')

    # Find replicated features
    exceeds_lod <- if (comparator == '>'){ function(value, lod) value >  lod
            } else if (comparator == '!=') function(value, lod) value != lod
    V1 <- dt[,.I[sum(exceeds_lod(value, lod), na.rm=TRUE)>1],
              by = c('feature_id', 'subgroup')]$V1
            #https://stackoverflow.com/questions/16573995

    # Keep only replicated features
    replicated_features <- dt[V1]$feature_id
    idx <- fid_values(object) %in% replicated_features
    message('\t\tFilter ', sum(idx), '/', length(idx), ' features: expr ',
            comparator, ' ', as.character(lod),
            ' for at least two samples in some subgroup')
    object %<>% extract_features(idx) # also handles limma in metadata

    # Update analysis log
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>% c(structure(
                sum(idx),
                names = sprintf(
                    "expr %s %s, for at least two samples in some subgroup",
                    comparator, as.character(lod))))
    }
    object
}

#=======================
# FILTER SAMPLES
#=======================


#' @rdname filter_samples
#' @export
filter_samples_ <- function(object, condition, verbose = FALSE, record = TRUE){
    if (is.null(condition)) return(object)
    idx <- lazy_eval(condition, sdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) if (verbose) message(
        '\t\tRetain ', sum(idx), '/', length(idx), ' samples: ',
        if (class(condition)=='lazy') deparse(condition$expr) else condition)
    object %<>% extract(, idx)
    sdata(object) %<>% droplevels()
    if (record && !is.null(analysis(object))) {
        analysis(object)$nsamples %<>% c(structure(sum(idx), names =
        if (class(condition)=='lazy') deparse(condition$expr) else condition))
    }
    object
}


#' Filter samples on condition
#' @param object    SummarizedExperiment
#' @param condition filter condition
#' @param verbose   TRUE or FALSE (default)
#' @param record    TRUE (default) or FALSE
#' @return filtered SummarizedExperiment
#' @examples
#' # GLUTAMINASE
#'     file <- download_data('glutaminase.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     filter_samples(object,   TIME_POINT=='h10',  verbose = TRUE)
#'     filter_samples_(object, "TIME_POINT=='h10'", verbose = TRUE)
#' @export
filter_samples <- function(object, condition, verbose = FALSE, record = TRUE){
    # filter_samples_(object, lazyeval::lazy(condition), verbose = verbose)
          # old version using deprecated lazyeval package
    condition <- enquo(condition)
    idx <- eval_tidy(condition, sdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx),
                                        ' samples: ', expr_text(condition))
    object %<>% extract(, idx)
    sdata(object) %<>% droplevels()
    if (record && !is.null(analysis(object))) {
        analysis(object)$nsamples %<>%  c(structure(sum(idx),
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







