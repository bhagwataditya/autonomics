

single_subgroup <- function(object){
    assert_is_subset('subgroup', svars(object))
    length(unique(object$subgroup))==1
}


are_factor <- function(df) vapply(df, is.factor, logical(1))


#' Create design
#'
#'  Create design matrix  for statistical analysis
#' @param object      SummarizedExperiment
#' @param formula     formula with svars
#' @param ...         backward compatibility
#' @return design matrix
#' @examples
#' file <- download_data('hypoglycemia.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' create_design(object)
#' create_design(object, ~ 0 + subgroup + Sex + T2D + age + bmi)
#' @export
create_design <- function(
    object,
    formula = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(all.vars(formula), svars(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    for (var in all.vars(formula)){
        if (is.character(sdata(object)[[var]])){
            sdata(object)[[var]] %<>% factor()
        }
    }
# Create design matrix
    myDesign <- model.matrix(formula,  data = sdata(object))
# Rename: intercept -> factor1level1
    factors <- svars(object)[are_factor(sdata(object))]
    factor1 <- factors[1]
    level1  <- levels(sdata(object)[[factor1]])[1]
    colnames(myDesign) %<>% gsub('(Intercept)', level1, ., fixed = TRUE)
# Rename regressors
    for (var in factors) colnames(myDesign) %<>% gsub(var, '', ., fixed = TRUE)
        # Fails for e.g. T2D = YES/NO: a meaningless column "YES" is created
        # For other cases it works wonderfully, so I keep it for now.
        # If it gives too many issues, roll back to doing the dropping only
        # for "subgroup" levels:
        #colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
# Validify names
    colnames(myDesign) %<>% gsub(':', '..', ., fixed = TRUE)
    colnames(myDesign) %<>% make.names()
# Return
    return(myDesign)
}


#' @rdname create_design
#' @export
create_design_matrix <- function(...){
    .Deprecated('create_design')
    create_design(...)
}






