#' Create design matrix for statistical analysis
#' @param object      SummarizedExperiment
#' @param intercept   TRUE or FALSE: include an intercept in the design?
#' @param confounders confounder svars (character)
#' @return design matrix
#' @examples
#' # STEM CELL COMPARISON
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' create_design_matrix(object)
#'
#' # GLUTAMINASE
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' create_design_matrix(object)[1:10, 1:10]
#' @export
create_design_matrix <- function(
    object,
    intercept = length(unique(sdata(object)$subgroup)) == 1,
    confounders = character(0)
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset('subgroup', svars(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    sdata1 <- sdata(object)
    if(is.character(sdata1$subgroup)) sdata1$subgroup %<>% factorify()
# Create formula
    formula <- if (intercept) '~ 1' else '~ 0'
    if (length(unique(sdata(object)$subgroup)) > 1) formula %<>%
        paste0(' + subgroup')
    if (length(confounders)>0){
        formula %<>% sprintf('%s + %s', ., paste0(confounders, collapse=' + '))
    }
    formula %<>% as.formula()
# Create design matrix
    myDesign <- model.matrix(formula,  data = sdata1)
# Rename coefficients
    # ~ 1 + subgroup
    if (intercept){
        subgroup1 <- as.character(unique(sdata(object)$subgroup)[1])
        colnames(myDesign) %<>% gsub('(Intercept)', subgroup1, ., fixed = TRUE)
        colnames(myDesign) %<>% stri_replace_first_regex('subgroup(.+)',
                                                    paste0('$1_', subgroup1))
    # ~ 0 + subgroup
    } else { colnames(myDesign) %<>% gsub('subgroup', '', ., fixed = TRUE) }
    # Rename confounders
    if (length(confounders) > 0){
        numeric_confounders <-  sdata(object)[, confounders, drop = FALSE] %>%
                                vapply(is.numeric, logical(1)) %>%
                                which() %>% names()
        factor_confounders  <- confounders %>% setdiff(numeric_confounders)
    }
    # Validify names
    colnames(myDesign) %<>% gsub(':', '..', ., fixed = TRUE)
    colnames(myDesign) %<>% make.names()
# Return
    return(myDesign)
}
