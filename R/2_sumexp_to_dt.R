#' @rdname sumexp_to_longdt
#' @export
sumexp_to_widedt <- function(
    object,
    fid   = 'feature_id',
    fvars = autonomics::fvars(object), # intersect('feature_name', autonomics::fvars(object)),
    assay = assayNames(object)[1]
){

    # Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fid,   fvars(object))
    assert_is_subset(fvars, fvars(object))

    # Extract
    wide1 <- fdt(object)[, unique(c(fid, fvars)), with = FALSE]
    if (is.numeric(assay))  assay <- assayNames(object)[assay]
    for (ass in assay){
        exprs1 <- data.table(assays(object)[[ass]])
        if (length(assay) > 1)  colnames(exprs1) %<>% paste0(ass, '.', .)
        wide1  %<>% cbind(exprs1)
    }
    wide1[, (fid) := factor(get(fid), unique(fdata(object)[[fid]]))]

    # Return
    wide1[]
}



#' Convert SummarizedExperiment into data.table
#'
#' @details
#' \itemize{
#'    \item \code{sumexp_to_widedt}:   feature          x sample
#'    \item \code{sumexp_to_subrep_dt}: feature.subgroup x replicate
#'    \item \code{sumexp_to_longdt}:   feature.sample
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
#'     object <- read_maxquant_proteingroups(file, invert = invert_subgroups, plot = FALSE)
#'     sumexp_to_widedt(object)
#'     sumexp_to_longdt(object)
#'     sumexp_to_subrep_dt(object)
#'
#' # Glutaminase
#'    file <- download_data('atkin18.metabolon.xlsx')
#'    object <- read_metabolon(file)
#'    sumexp_to_widedt(object)
#'    sumexp_to_longdt(object)
#'    sumexp_to_subrep_dt(object)
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
    fid   = 'feature_id',
    fvars = intersect('feature_name', autonomics::fvars(object)),
    sid   = 'sample_id',
    svars = intersect('subgroup', autonomics::svars(object)),
    assay = assayNames(object) %>% intersect(c(.[1], 'is_imputed'))
){
    # Assert
    . <- NULL
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
    dt <- sumexp_to_widedt(object, fid, fvars, assay = assay[1])
    dt %<>% melt(id.vars = unique(c(fid, fvars)), variable.name = sid,
                value.name = 'value')
    # Merge
    if (length(assay)>1){
        for (ass in assay[-1]){
            assdt <- sumexp_to_widedt(
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

#' Write sumexp to tsv
#' @param object SummarizedExperiment
#' @param assay string
#' @param file filename
#' @return NULL
#' @examples 
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file, fit = 'limma')
#' tsv <- file.path(tempdir(), 'fukuda20.proteingroups.tsv')
#' sumexp_to_tsv(object, tsv)
#' @export
sumexp_to_tsv <- function(object, assay = assayNames(object)[1], file){
    widedt <- sumexp_to_widedt(object, assay = assay)
    fwrite(widedt, file, sep = '\t')
}

#' Get fit vars/dt
#' @param object SummarizedExperimenmt
#' @return string vector
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma')
#' fitvars(object)
#' fitdt(object)
#' @export
fitvars <- function(object){
    fvars(object) %>% extract(stri_detect_fixed(., FITSEP))
}

#' @rdname fitvars
#' @export
fitdt <- function(object){
    cols <- c('feature_id', fitvars(object))
    fdt(object)[, cols, with = FALSE]
}

#' fitcoefs
#' @param object SummarizedExperiment
#' @return string vector
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma')
#' fitcoefs(object)
#' default_fitc
#' @export
fitcoefs <- function(object){
    fitvars(object)  %>%  split_extract_fixed(FITSEP, 2:3) %>%  unique()
}


extract_contrast_fdt <- function(object, fitcoef){
# Order on p
       pvar   <- paste0(     'p', FITSEP, fitcoef)
       fdrvar <- paste0(   'fdr', FITSEP, fitcoef)
    effectvar <- paste0('effect', FITSEP, fitcoef)
    idx <- order(fdt(object)[[pvar]])
    object %<>% extract(idx, )
# Extract results
    cols <- fvars(object) %>% setdiff('feature_id') %>% setdiff(fitvars(object))
    cols %<>% c('feature_id', pvar, fdrvar, effectvar, .)
    
    fdt0 <- fdt(object)[, cols, with = FALSE]
    names(fdt0) %<>% stri_replace_first_fixed(paste0(FITSEP, fitcoef), '')
    fdt0
}


#' Write xl/ods
#' @param object  SummarizedExperiment
#' @param xlfile  file
#' @return filepath
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma')
#' xlfile  <- file.path(tempdir(), 'fukuda20.proteingroups.fdt.xlsx')
#' odsfile <- file.path(tempdir(), 'fukuda20.proteingroups.fdt.ods')
#' # write_xl(object,  xlfile)
#' # write_ods(object, odsfile)
#' @export
write_xl <- function(object, xlfile){
# Assert
    if (!requireNamespace('writexl', quietly = TRUE)){
        message("`BiocManager::install('AnnotationDbi')`. Then re-run.")
        return(NULL)
    }
    assert_is_valid_sumexp(object)
    assert_all_are_dirs(dirname(xlfile))
# Write    
    list0 <- mapply(extract_contrast_fdt, fitcoef = fitcoefs(object), 
                    MoreArgs = list(object = object), SIMPLIFY = FALSE)
    writexl::write_xlsx(list0, path = xlfile)
# Return
    return(xlfile)
}


#' @rdname write_xl
#' @export
write_ods <- function(object, odsfile){
# Assert
    if (!requireNamespace('readODS', quietly = TRUE)){
        message("`BiocManager::install('readODS')`. Then re-run.")
        return(NULL)
    }
    assert_is_valid_sumexp(object)
    assert_all_are_dirs(dirname(odsfile))
# Write    
    list0 <- mapply(extract_contrast_fdt, fitcoef = fitcoefs(object), 
                    MoreArgs = list(object = object), SIMPLIFY = FALSE)
    if (file.exists(odsfile))  unlink(odsfile)
    mapply(readODS::write_ods, x = list0, sheet = names(list0), 
           MoreArgs = list(path = odsfile, append = TRUE), 
           SIMPLIFY = FALSE)
    
# Return
    return(xlfile)
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
#' @param verbose       TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' require(magrittr) 
#' file <- download_data('atkin18.metabolon.xlsx')
#' x <- values(read_metabolon(file))
#' object <- matrix2sumexp(x)
#' object %<>% pca()
#' biplot(object, color = subgroup)
#' @export
matrix2sumexp <- function(x, verbose = TRUE){
    object <- SummarizedExperiment(list(exprs = x))
    fdata(object)$feature_id   <- rownames(object)
    sdata(object)$sample_id    <- colnames(object)
    object %<>% add_subgroup(verbose = verbose)
    object
}


