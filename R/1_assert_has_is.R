
#==============================================================================
# has/contains


#' Does object have complete svalues
#' @param object SummarizedExperiment
#' @param svar   sample var
#' @return logical
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' has_complete_subgroup_values(object)
#' has_complete_block_values(object)
#' @noRd
has_complete_svalues <- function(object, svar){

    # svar missing
    var_present <- svar %in% svars(object)
    if (!var_present) return(FALSE)

    # svalues missing
    values_present <- any(is_empty_character(svalues(object, svar)))
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


#' Does object contain ratio values?
#' @param object SummarizedExperiment
#' @return logical
#' @examples
#' # STEM CELL COMPARISON
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' contains_ratios(object)
#'
#' # GLUTAMINASE
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' contains_ratios(object)
#' @noRd
contains_ratios <- function(object){
    quantity <- metadata(object)$quantity
    if (is.null(quantity)) return(FALSE)
    return(stri_detect_fixed(quantity, 'Ratio'))
}


#==============================================================================
#
#                        assert_is_valid_sumexp
#
#==============================================================================


has_valid_fnames <- function(x, .xname = get_name_in_parent(x)){

    if (is.null(fnames(x))){
        return(false('fnames(%s) are NULL', .xname))}

    if (!all(fnames(x) == fdata(x)$feature_id)){
        return(false('fnames(%s) != fdata(%s)$feature_id', .xname, .xname))}

    #if (!all(fnames(x) == rownames(exprs(x)))){
    #    return(false('fnames(%s) != rownames(exprs(%s))', .xname, .xname))}

    #if (!all(fnames(x) == rownames(fdata(x)))){
    #    return(false('fnames(%s) != rownames(fdata(%s))', .xname, .xname))}

    TRUE
}


has_valid_snames <- function(x, .xname = get_name_in_parent(x)){

    if (is.null(snames(x))){
        return(false('snames(%s) are NULL', .xname))}

    if (!all(snames(x) == sdata(x)$sample_id)){
        return(false('snames(%s) != sdata(%s)$sample_id', .xname, .xname))}

    #if (!all(snames(x) == colnames(exprs(x)))){
    #    return(false('snames(%s) != colnames(exprs(%s))', .xname, .xname))}

    #if (!all(snames(x) == rownames(sdata(x)))){
    #    return(false('snames(%s) != colnames(sdata(%s))', .xname, .xname))}

    TRUE
}




#' Is valid SummarizedExperiment
#' @param x SummarizedExperiment
#' @param .xname see assertive.base::get_name_in_parent
#' @return TRUE or FALSE
#' @noRd
is_valid_sumexp <- function(x, .xname = get_name_in_parent(x)){
    if (!(ok <- assertive::is2(x, "SummarizedExperiment")))  return(ok)
    if (!(ok <- has_valid_fnames(x, .xname = .xname)))       return(ok)
    if (!(ok <- has_valid_snames(x, .xname = .xname)))       return(ok)
    TRUE
}


#' Assert that x is a valid SummarizedExperiment
#'
#' @param x SummarizedExperiment
#' @param .xname see assertive.base::get_name_in_parent
#' @return TRUE or FALSE
#' @examples
#' # VALID
#'     file <- download_data('halama18.metabolon.xlsx')
#'     x <- read_metabolon(file, plot = FALSE)
#'     assert_is_valid_sumexp(x)
#' # NOT VALID
#'     rownames(SummarizedExperiment::colData(x)) <- NULL
#'     # assert_is_valid_sumexp(x)
#' @export
assert_is_valid_sumexp <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_valid_sumexp, x, .xname = get_name_in_parent(x))
}







