
#===============================================================================
# has/contains


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

