#' @rdname sumexp_to_long_dt
#' @export
sumexp_to_wide_dt <- function(
    object,
    fid = 'feature_id',
    fvars = character(0)
){
    fdata1 <- data.table(fdata(object)[, unique(c(fid, fvars)), drop = FALSE])
    exprs1 <- data.table(exprs(object))
    wide1  <- cbind(fdata1, exprs1)
    wide1
}


#' Convert SummarizedExperiment into data.table
#' @param object sumexp
#' @param fid    fvar carrying feature id
#' @param fvars  additional fvars to include in table
#' @param sid    svar carrying sample id
#' @param svars  additional svars to include in table
#' @param ...    (backward compatibility)
#' @return data.table
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' sumexp_to_wide_dt(object)
#' sumexp_to_long_dt(object)
#' sumexp_to_long_dt(object, svars = 'subgroup')
#' @export
sumexp_to_long_dt <- function(object, fid = 'feature_id', fvars = character(0),
    sid = 'sample_id', svars = character(0)
){
    sdata1 <- sdata(object)[, c('sample_id', svars), drop = FALSE]

    # Note: unique is to avoid duplication of same fields in fid and fvars
    sumexp_to_wide_dt(object, fid, fvars) %>%
    data.table::melt.data.table(id.vars = unique(c(fid, fvars)),
                                variable.name = sid, value.name = 'value') %>%
    merge(sdata1, by = sid) %>%
    extract(, unique(c(fid, fvars, sid, svars, 'value')), with = FALSE)
}
