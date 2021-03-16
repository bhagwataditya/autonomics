#' @rdname sumexp_to_long_dt
#' @export
sumexp_to_wide_dt <- function(
    object,
    fid   = 'feature_id',
    fvars = intersect('feature_name', autonomics::fvars(object)),
    assay = assayNames(object)[1]
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
#' @param subgroup subgroup (sym)
#' @param fid    fvar carrying feature id
#' @param fvars  additional fvars to include in table
#' @param sid    svar carrying sample id
#' @param svars  additional svars to include in table
#' @param assay  matrix in assays(object) to be used
#' @return data.table
#' @examples
#' # Stem cell comparison
#'     file <- download_data('billing16.proteingroups.txt')
#'     invert_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
#'     object <- read_proteingroups(file, invert_subgroups = invert_subgroups,
#'                   plot=FALSE)
#'     sumexp_to_wide_dt(object)
#'     sumexp_to_long_dt(object)
#'     sumexp_to_subrep_dt(object)
#'
#' # Glutaminase
#'    require(magrittr)
#'    file <- download_data('atkin18.metabolon.xlsx')
#'    object <- read_metabolon(file, plot=FALSE)
#'    sumexp_to_wide_dt(object)
#'    sumexp_to_long_dt(object)
#'    sumexp_to_subrep_dt(object, Group)
#'
#' # Fukuda
#'    require(magrittr)
#'    file <- download_data('fukuda20.proteingroups.txt')
#'    object <- read_proteingroups(file, impute=FALSE, plot=FALSE)
#'    sumexp_to_long_dt(object)
#'    object %<>% impute_systematic_nondetects(plot=FALSE)
#'    sumexp_to_long_dt(object)
#' @export
sumexp_to_long_dt <- function(
    object,
    fid   = 'feature_id',
    fvars = intersect('feature_name', autonomics::fvars(object)),
    sid   = 'sample_id',
    svars = intersect('subgroup', autonomics::svars(object)),
    assay = assayNames(object) %>% intersect(c(.[1], 'is_imputed'))
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fid,   autonomics::fvars(object))
    assert_is_subset(sid,   autonomics::svars(object))
    assert_is_subset(fvars, autonomics::fvars(object))
    assert_is_subset(svars, autonomics::svars(object))
    common <- intersect(svars, fvars)
    if (length(common) > 0){
        message('\t\tRemove clashing svars/fvars: ',
                paste0(common,collapse = ', '))
        fvars %<>% setdiff(common) # Avoid name clashes
        svars %<>% setdiff(common)
    }
# Melt
    melt <- data.table::melt.data.table
    dt <- sumexp_to_wide_dt(object, fid, fvars, assay = assay[1])
    dt %<>% melt(id.vars = unique(c(fid, fvars)), variable.name = sid,
                value.name = 'value')
# Merge
    if (length(assay)>1){
        for (ass in assay[-1]){
            assdt <- sumexp_to_wide_dt(
                        object, fid=fid, fvars = character(0), assay = ass)
            assdt %<>% melt(id.vars = fid, variable.name = sid, value.name=ass)
            dt %<>% merge(assdt, by = c(fid, sid))
        }
    }
    sdata1 <- sdata(object)[, c(sid, svars), drop = FALSE]
    dt %<>% merge(sdata1, by = sid)
    cols <- intersect(unique(
                c(fid, fvars, sid, svars, 'value', assay[-1])), names(dt))
    dt %<>% extract(, cols, with = FALSE)
        # Note: unique is to avoid duplication of same fields in fid and fvars
# Order
    dt[, (fid) := factor(get(fid), unique(fdata(object)[[fid]]))]
    dt[, (sid) := factor(get(sid), unique(sdata(object)[[sid]]))]
# Return
    dt[]
}


#' @export
#' @rdname sumexp_to_long_dt
sumexp_to_subrep_dt <- function(object, subgroup=subgroup){
    subgroup <- enquo(subgroup)
    subgroupvar <- as_name(subgroup)

    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(subgroupvar, svars(object))
    sample_id <- subgroup <- . <- NULL

    # Melt
    dt <- sumexp_to_long_dt(object, svars = subgroupvar)
    sep <- guess_sep(object)
    dt[, replicate := stri_replace_first_fixed(
                        sample_id, dt[[subgroupvar]], '') %>%
                    substr(2, nchar(.))]

    # Partially recast
    formula <- as.formula(sprintf('feature_id + %s ~ replicate', subgroupvar))
    subrepdt <- data.table::dcast(dt, formula=formula)
    subrepdt$feature_id %<>% factor(fid_values(object))
    setorderv(subrepdt, 'feature_id')

    # Return
    subrepdt
}


#' `data.table` to `matrix`
#'
#' Convert between `data.table` and `matrix`
#'
#' @param x  data.table / matrix
#' @param idvar idvar string
#' @return matrix / data.table
#' @examples
#' x <- data.table::data.table(
#'         gene    = c('ENSG001', 'ENSG002', 'ENSG003'),
#'         sampleA = c(1787, 10, 432),
#'         sampleB = c(1143,  3, 268))
#' dt2mat(x)
#' mat2dt(dt2mat(x), 'gene')
#' @export
dt2mat    <- function(x) x[,-1] %>% as.matrix() %>% set_rownames(x[[1]])


#' @rdname dt2mat
#' @export
mat2dt    <- function(x, idvar){
                x %<>% data.table(keep.rownames = TRUE)
                names(x)  %<>% gsub('rn', idvar, .,  fixed = TRUE)
                x }

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
                                                        drop = FALSE ]
    fdata1 <- dt2DF(dt[, .SD[1], by = 'feature_id'])[, c('feature_id', fvars),
                                                        drop = FALSE ]
    SummarizedExperiment(
        assays = list(exprs = exprs1),
        rowData = fdata1,
        colData = sdata1)
}


#' Convert matrix into SummarizedExperiment
#' @param x             matrix
#' @param sdt           sample data.table / data.frame / DataFrame
#' @param sdtby         sample data mergeby column
#' @param subgroupvar   string / NULL
#' @param fdt           feature data.table / data.frame / DataFrame
#' @param fdtby         feature data mergeby column
#' @param fnamevar      string / NULL
#' @param verbose       TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' require(magrittr) 
#' file <- download_data('atkin18.metabolon.xlsx')
#' x <- values(read_metabolon(file, plot=FALSE))
#' object <- matrix2sumexp(x)
#' object %<>% pca()
#' biplot(object, nloadings=0, color=subgroup)
#' @export
matrix2sumexp <- function(
    x,
    sdt         = NULL,
    sdtby       = if (is.null(sdt))   NULL else names(sdt)[1],
    subgroupvar = NULL,
    fdt         = NULL,
    fdtby       = if (is.null(fdt)) NULL else names(fdt)[1],
    fnamevar    = NULL,
    verbose     = TRUE
){
# exprs
    object <- SummarizedExperiment(list(exprs = x))
    fdata(object)$feature_id   <- rownames(object)
    fdata(object)$feature_name <- rownames(object)
    sdata(object)$sample_id    <- colnames(object)
# sdata
    object %<>% merge_sdata(sdt, by.x = 'sample_id', by.y = sdtby)
    object %<>% merge_fdata(fdt, by.x ='feature_id', by.y = fdtby)
    object %<>% add_subgroup(subgroupvar, verbose = verbose)
# return
    object
}


#' Create MultiAssayExperiment from SummarizedExperiment list
#' @param experiments named list of SummarizedExperiments
#' @return MultiAssayExperiment
#' @examples
#' require(magrittr)
#' somascanfile  <- download_data('atkin18.somascan.adat')
#' metabolonfile <- download_data('atkin18.metabolon.xlsx')
#' somascan <- read_somascan(somascanfile,   plot=FALSE)
#' metabolon<- read_metabolon(metabolonfile, plot=FALSE)
#' svars(somascan)  %<>% stringi::stri_replace_first_fixed(
#'                         'SampleGroup', 'subgroup')
#' svars(metabolon) %<>% stringi::stri_replace_first_fixed(
#'                         'Group',       'subgroup')
#' metabolon$replicate <- NULL
#' object   <- sumexp2mae(list(somascan=somascan, metabolon=metabolon))
#' @export
sumexp2mae <- function(experiments){
    assert_is_list(experiments)
    assert_has_names(experiments)
    for (experiment in experiments){
        assert_is_all_of(experiment, 'SummarizedExperiment')
        assert_is_subset(c('sample_id', 'subgroup'), svars(experiment))
    }
    for (i in seq_along(experiments))  experiments[[i]] %<>%
                                            extract(, order(colnames(.)))
    extract_sdata <- function(sumexp){
        extractvars <- c('sample_id', 'subgroup', 'replicate')
        extractvars %<>% intersect(svars(sumexp))
        sdata(sumexp)[, extractvars, drop=FALSE]
    }
    sdata1 <- unique(Reduce(rbind, lapply(experiments, extract_sdata)))
    sdata1 %<>% extract(order(.$sample_id), )
    assert_all_are_true(table(sdata1$sample_id)==1)
    MultiAssayExperiment(experiments = experiments, colData = sdata1)
}


