#' @rdname sumexp_to_long_dt
#' @export
sumexp_to_wide_dt <- function(object, fid = 'feature_id', fvars = character(0)){

    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fid,   fvars(object))
    assert_is_subset(fvars, fvars(object))

    # Extract
    fdata1 <- data.table(fdata(object)[, unique(c(fid, fvars)), drop = FALSE])
    exprs1 <- data.table(exprs(object))
    wide1  <- cbind(fdata1, exprs1)

    # Return
    wide1
}


#' Convert SummarizedExperiment into data.table
#' @details
#' \itemize{
#'    \item \code{sumexp_to_wide_dt}:   feature          x sample
#'    \item \code{sumexp_to_subrep_dt}: feature.subgroup x replicate
#'    \item \code{sumexp_to_long_dt}:   feature.sample
#' }
#' @param object sumexp
#' @param fid    fvar carrying feature id
#' @param fvars  additional fvars to include in table
#' @param sid    svar carrying sample id
#' @param svars  additional svars to include in table
#' @return data.table
#' @examples
#' # Stem cell differentiation
#'     file <- download_data('stemcells.proteinGroups.txt')
#'     object <- read_proteingroups(file)
#'     sumexp_to_wide_dt(object)
#'     sumexp_to_long_dt(object)
#'     sumexp_to_subrep_dt(object)
#'
#' # Glutaminase
#'    file <- download_data('glutaminase.metabolon.xlsx')
#'    object <- read_metabolon(file)
#'    object %<>% add_pca()
#'    sumexp_to_wide_dt(object)
#'    sumexp_to_long_dt(object)
#'    sumexp_to_subrep_dt(object)
#' @export
sumexp_to_long_dt <- function(
    object,
    fid = 'feature_id',
    fvars = character(0),
    sid = 'sample_id',
    svars = if ('subgroup' %in% importomics::svars(object)){ 'subgroup'
            } else {                                         character(0) }
){
    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fid,   importomics::fvars(object))
    assert_is_subset(sid,   importomics::svars(object))
    assert_is_subset(fvars, importomics::fvars(object))
    assert_is_subset(svars, importomics::svars(object))
    common <- intersect(svars, fvars)
    if (length(common) > 0){
        message('\t\tRemove clashing svars/fvars: ', paste0(common,collapse = ', '))
        fvars %<>% setdiff(common) # Avoid name clashes
        svars %<>% setdiff(common)
    }

    # Extract & Return
    sdata1 <- sdata(object)[, c('sample_id', svars), drop = FALSE]
    sumexp_to_wide_dt(object, fid, fvars) %>%
    data.table::melt.data.table(
        id.vars = unique(c(fid, fvars)),
        variable.name = sid, value.name = 'value') %>%
        # Note: unique is to avoid duplication of same fields in fid and fvars
    merge(sdata1, by = sid) %>%
    extract(, unique(c(fid, fvars, sid, svars, 'value')), with = FALSE)
}

#' @export
#' @rdname sumexp_to_long_dt
sumexp_to_subrep_dt <- function(object){

    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset('subgroup', svars(object))
    sample_id <- subgroup <- . <- NULL

    # Melt
    dt <- sumexp_to_long_dt(object, svars = 'subgroup')
    sep <- guess_sep(object)
    dt[, replicate := stri_replace_first_fixed(sample_id, subgroup, '') %>%
                        substr(2, nchar(.))]

    # Partially recast
    subrepdt <- data.table::dcast(dt, feature_id + subgroup ~ replicate)
    subrepdt$feature_id %<>% factor(fid_values(object))
    setorderv(subrepdt, 'feature_id')

    # Return
    subrepdt
}
