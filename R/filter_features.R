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


