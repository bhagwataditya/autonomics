#' @rdname sumexp_to_long_dt
#' @export
sumexp_to_wide_dt <- function(
    object,
    fid   = 'feature_id',
    fvars = 'feature_name',
    assay = 'exprs'
){

    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fid,   fvars(object))
    assert_is_subset(fvars, fvars(object))

    # Extract
    fdata1 <- data.table(fdata(object)[, unique(c(fid, fvars)), drop = FALSE])
    exprs1 <- data.table(assays(object)[[assay]])
    wide1  <- cbind(fdata1, exprs1)
    wide1[, (fid) := factor(get(fid), unique(fdata(object)[[fid]]))]

    # Return
    wide1[]
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
#' @param assay  matrix in assays(object) to be used
#' @return data.table
#' @examples
#' # Stem cell comparison
#'     file <- download_data('billing16.proteingroups.txt')
#'     object <- read_proteingroups(file)
#'     sumexp_to_wide_dt(object)
#'     sumexp_to_long_dt(object)
#'     sumexp_to_subrep_dt(object)
#'
#' # Glutaminase
#'    require(magrittr)
#'    file <- download_data('halama18.metabolon.xlsx')
#'    object <- read_metabolon(file)
#'    object %<>% pca()
#'    sumexp_to_wide_dt(object)
#'    sumexp_to_long_dt(object)
#'    sumexp_to_subrep_dt(object)
#' @export
sumexp_to_long_dt <- function(
    object,
    fid   = 'feature_id',
    fvars = 'feature_name',
    sid   = 'sample_id',
    svars = if ('subgroup' %in% importomics::svars(object)){ 'subgroup'
            } else {                                         character(0) },
    assay = 'exprs'
){
    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fid,   importomics::fvars(object))
    assert_is_subset(sid,   importomics::svars(object))
    assert_is_subset(fvars, importomics::fvars(object))
    assert_is_subset(svars, importomics::svars(object))
    common <- intersect(svars, fvars)
    if (length(common) > 0){
        message('\t\tRemove clashing svars/fvars: ',
                paste0(common,collapse = ', '))
        fvars %<>% setdiff(common) # Avoid name clashes
        svars %<>% setdiff(common)
    }

    # Extract & Return
    sdata1 <- sdata(object)[, c('sample_id', svars), drop = FALSE]
    dt  <-  sumexp_to_wide_dt(object, fid, fvars, assay = assay) %>%
            data.table::melt.data.table(
                id.vars       = unique(c(fid, fvars)),
                variable.name = sid,
                value.name    = 'value') %>%
            # Note: unique is to avoid duplication of same fields in fid and fvars
            merge(sdata1, by = sid) %>%
            extract(, unique(c(fid, fvars, sid, svars, 'value')), with = FALSE)

    # Encode fid/sid order
    dt[, (fid) := factor(get(fid), unique(fdata(object)[[fid]]))]
    dt[, (sid) := factor(get(sid), unique(sdata(object)[[sid]]))]
    dt[]
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

dt2mat    <- function(dt) dt[,-1] %>% as.matrix() %>% set_rownames(dt[[1]])
dt2DF     <- function(dt) DataFrame(dt, row.names = dt[[1]])
dt2exprs  <- function(dt) dt2mat(data.table::dcast(
                            dt, feature_id ~ sample_id, value.var = 'value'))
dt2sumexp  <- function(
    dt,
    fvars = character(0),
    svars = setdiff(names(dt), c('feature_id', fvars, 'sample_id', 'value'))
){
    exprs1 <- dt2exprs(dt)
    sdata1 <- dt2DF(dt[, .SD[1], by = 'sample_id' ])[, c('sample_id',  svars),
                                                     drop = FALSE]
    fdata1 <- dt2DF(dt[, .SD[1], by = 'feature_id'])[, c('feature_id', fvars),
                                                     drop = FALSE]
    SummarizedExperiment(
        assays = list(exprs = exprs1),
        rowData = fdata1,
        colData = sdata1)
}


#' Convert matrix into SummarizedExperiment
#' @param x matrix
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' x <- exprs(read_metabolon(file, plot=FALSE))
#' object <- matrix2sumexp(x)
#' biplot(object, nloadings=0)
#' @export
matrix2sumexp <- function(x){
    object <- SummarizedExperiment(list(exprs = x))
    fdata(object)$feature_id <- rownames(object)
    fdata(object)$feature_name <- rownames(object)
    sdata(object)$sample_id    <- colnames(object)
    object$subgroup <- 'group1'
    #object %<>% add_designvars(designfile = NULL) # too slow for large matrices
    object
}
