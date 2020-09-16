
#=====================================================================================
# has/contains

is_summarized_experiment <- function(x, .xname = get_name_in_parent(x)){
   if (!methods::is(x, 'SummarizedExperiment')){
      return(false('%s is not a SummarizedExperiment', .xname))
   }
   TRUE
}


has_valid_featureNames <- function(x, .xname = get_name_in_parent(x)){
   # SummarizedExperiments do not allow row naming of fdata
   # if (!all(fnames(x) == rownames(fdata(x)))){
   #   return(false('fnames(%s) differ from rownames(fdata(%s))',.xname,.xname))
   # }
   if (!all(fnames(x) == rownames(exprs(x)))){
      return(false('fnames(%s) differ from rownames(exprs(%s))',.xname,.xname))
   }
   TRUE
}


has_valid_sampleNames <- function(x, .xname = get_name_in_parent(x)){
   if (!all(snames(x) == rownames(sdata(x)))){
      return(false('snames(%s) differ from rownames(sdata(%s))',.xname,.xname))
   }
   if (!all(snames(x) == colnames(exprs(x)))){
      return(false('snames(%s) differ from colnames(exprs(%s))',.xname,.xname))
   }
   TRUE
}


#' Does object have complete svalues
#' @param object SummarizedExperiment
#' @param svar   sample var
#' @return logical
#' @examples
#' file <- load_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' has_complete_subgroup_values(object)
#' has_complete_block_values(object)
#' @noRd
has_complete_svalues <- function(object, svar){

    # svar missing
    var_present <- svar %in% svars(object)
    if (!var_present) return(FALSE)

    # svalues missing
    values_present <-
    if (values_present) return(FALSE)

    # svar and svalues both present
    return(TRUE)
}


#' @rdname has_complete_svalues
#' @noRd
has_complete_subgroup_values <- function(object){
    has_complete_svalues(object, 'subgroup')
}


#' @rdname has_complete_svalues
#' @noRd
has_complete_block_values <- function(object){
    has_complete_svalues(object, 'block')
}


#' @title Does object contain prepro?
#' @description Does object contain preprocessing info?
#' @param object SummarizedExperiment
#' @return logical
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' contains_prepro(object)
#' @noRd
contains_prepro <- function(object){
    assert_is_valid_object(object)
    length(prepro(object))!=0
}

#' Is max quant object
#' @param object SummarizedExperiment
#' @return TRUE or FALSE
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' is_maxquant_eset(object)
#' @noRd
is_maxquant_eset <- function(object){
   assert_is_valid_object(object)
   all(c('Contaminant', 'Reverse') %in% fvars(object))
}


#' Does object contain ratio values?
#' @param object SummarizedExperiment
#' @return logical
#' @examples
#' # STEM CELL COMPARISON
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' contains_ratios(object)
#'
#' # GLUTAMINASE
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' contains_ratios(object)
#' @noRd
contains_ratios <- function(object){
    quantity <- metadata(object)$quantity
    if (is.null(quantity)) return(FALSE)
    return(stri_detect_fixed(quantity, 'Ratio'))
}


#=============================================================================
# Assert

#' Is valid SummarizedExperiment
#' @param x SummarizedExperiment
#' @param .xname see assertive.base::get_name_in_parent
#' @return TRUE or FALSE
#' @export
is_valid_object <- function(x, .xname = get_name_in_parent(x)){
    if (!(ok <- is_summarized_experiment(x, .xname = .xname)))  {   return(ok)}
    if (!(ok <- has_valid_featureNames(x, .xname = .xname)))    {   return(ok)}
    if (!(ok <- has_valid_sampleNames(x,  .xname = .xname)))    {   return(ok)}
    TRUE
}

#' Assert that x is a valid SummarizedExperiment
#' @param x SummarizedExperiment
#' @return error if not true
#' @export
assert_is_valid_object <- function(x){
    assert_engine(is_valid_object, x, .xname = get_name_in_parent(x))
}

#' Assert that features are valid
#' @param features \code{\link{numeric}} index or \code{\link{character}} names
#' of \code{\link{fdata}} entries
#' @param object SummarizedExperiment
#' @noRd
assert_all_are_valid_features <- function(features, object){
    if (is.character(features)) {
        assert_is_subset(features, fdata(object)[["feature_id"]])
    } else if (is.numeric(features)) {
        assert_all_are_whole_numbers(features)
        assert_all_are_in_closed_range(features, 1, nrow(object))
    } else {
        stop("Features must be numeric indexes or character objects.")
    }
    invisible(features)
}


