
#==============================================================================
# has/contains


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



