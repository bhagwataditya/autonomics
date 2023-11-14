#' @rdname sumexp_to_longdt
#' @export
sumexp_to_widedt <- function(
    object,
    fvars = autonomics::fvars(object),
    assay = assayNames(object)[1]
){

    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fvars, fvars(object))

    # Extract
    wide1 <- fdt(object)[, unique(c('feature_id', fvars)), with = FALSE]
    if (is.numeric(assay))  assay <- assayNames(object)[assay]
    for (ass in assay){
        exprs1 <- data.table(assays(object)[[ass]])
        if (length(assay) > 1)  colnames(exprs1) %<>% paste0(ass, '.', .)
        wide1  %<>% cbind(exprs1)
    }
    wide1[, feature_id := factor(feature_id, unique(fdt(object)$feature_id))]

    # Return
    wide1[]
}



#' SummarizedExperiment to data.table
#'
#' @details
#' \itemize{
#'    \item \code{sumexp_to_widedt}:    feature          x sample
#'    \item \code{sumexp_to_subrep_dt}: feature.subgroup x replicate
#'    \item \code{sumexp_to_longdt}:    feature.sample
#' }
#' @param object sumexp
#' @param subgroup subgroup (sym)
#' @param fvars  additional fvars to include in table
#' @param svars  additional svars to include in table
#' @param assay  matrix in assays(object) to be used
#' @return data.table
#' @examples
#' # Atkin Hypoglycemia
#'    file <- download_data('atkin.metabolon.xlsx')
#'    object <- read_metabolon(file)
#'    sumexp_to_widedt(object)
#'    sumexp_to_longdt(object)
#'    sumexp_to_subrep_dt(object)
#'
#' # Stem cell comparison
#'     file <- download_data('billing16.proteingroups.txt')
#'     invert_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
#'     object <- read_maxquant_proteingroups(file, invert = invert_subgroups, plot = FALSE)
#'     sumexp_to_widedt(object)
#'     sumexp_to_longdt(object)
#'     sumexp_to_subrep_dt(object)
#'
#' # Fukuda
#'    file <- download_data('fukuda20.proteingroups.txt')
#'    object <- read_maxquant_proteingroups(file, impute = FALSE)
#'    sumexp_to_longdt(object)
#'    object %<>% impute(plot = FALSE)
#'    sumexp_to_widedt(object)
#'    sumexp_to_longdt(object)
#' @export
sumexp_to_longdt <- function(
    object,
    fvars = intersect('feature_name', autonomics::fvars(object)),
    svars = intersect('subgroup',     autonomics::svars(object)),
    assay = assayNames(object) %>% intersect(c(.[1], 'is_imputed'))
){
# Assert
    . <- sample_id <- NULL
    assert_is_valid_sumexp(object)
    assert_is_subset(fvars, autonomics::fvars(object))
    assert_is_subset(svars, autonomics::svars(object))
    svars %<>% setdiff('sample_id')
    fvars %<>% setdiff('feature_id')
    common <- intersect(svars, fvars)
    if (length(common) > 0){
        message('\t\tRemove clashing svars/fvars: ', paste0(common,collapse = ', '))
        fvars %<>% setdiff(common) # Avoid name clashes
        svars %<>% setdiff(common)
    }
# Melt
    melt <- data.table::melt.data.table
    dt <- sumexp_to_widedt(object, fvars = fvars, assay = assay[1])
    dt %<>% melt(id.vars = unique(c('feature_id', fvars)), variable.name = 'sample_id', value.name = 'value')
# Merge
    if (length(assay)>1){
        for (ass in assay[-1]){
            assdt <- sumexp_to_widedt(object, fvars = character(0), assay = ass)
            assdt %<>% melt(id.vars = 'feature_id', variable.name = 'sample_id', value.name = ass)
            dt %<>% merge(assdt, by = c('feature_id', 'sample_id'))
        }
    }
    sdt1 <- sdt(object)[, c('sample_id', svars), with = FALSE]
    dt %<>% merge(sdt1, by = 'sample_id')
    cols <- c('feature_id', fvars, 'sample_id', svars, 'value', assay[-1])
    cols %<>% unique() # avoid duplication of same fields in fid and fvars
    cols %<>% intersect(names(dt))
    dt %<>% extract(, cols, with = FALSE)
# Order and Return
    dt[, feature_id := factor(feature_id, unique(fdt(object)$feature_id))]
    dt[, sample_id  := factor(sample_id,  unique(sdt(object)$sample_id ))]
    dt[]
}


#' @export
#' @rdname sumexp_to_longdt
sumexp_to_subrep_dt <- function(object, subgroup=subgroup){
    subgroup <- enquo(subgroup)
    subgroupvar <- as_name(subgroup)

    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(subgroupvar, svars(object))
    sample_id <- subgroup <- . <- NULL

    # Melt
    dt <- sumexp_to_longdt(object, svars = subgroupvar)
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
                . <- NULL
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
    . <- NULL
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


