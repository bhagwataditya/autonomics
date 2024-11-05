#==============================================================================
#
#                           rm_unmatched_samples
#                           rm_singleton_samples
#
#==============================================================================

#' rm unmatched/singleton samples
#'
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup variable (string)
#' @param subgroupctr  control subgroup (string)
#' @param block        block variable (string)
#' @param verbose      TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
#' object <- read_somascan(file)
#' object %<>% filter_samples(subgroup %in% c('t1', 't2'), verbose = TRUE)
#' rm_singleton_samples(object, subgroupvar = 'Subject')
#' rm_unmatched_samples(object, subgroupvar = 'SampleGroup', block = 'Subject')
#' @export
rm_unmatched_samples <- function(
    object,
    subgroupvar = 'subgroup',
    subgroupctr = slevels(object, subgroupvar)[1],
    block,
    verbose = TRUE
){
    snames1 <- sdt(object)[,
                .SD[(sum(get(subgroupvar)==subgroupctr)==1) &
                    (sum(get(subgroupvar)!=subgroupctr) >0)],  by = block]$sample_id
    n <- length(snames1)
    if (verbose & n < ncol(object)){
        message('\t\tRetain ', n, '/', ncol(object), ' samples with matching ', subgroupctr) }
    object %<>% extract(, snames1)
    object
}

#' @rdname rm_unmatched_samples
#' @export
rm_singleton_samples <- function(object, subgroupvar = 'subgroup', verbose = TRUE){
    selectedsamples <- sdt(object)[, .SD[.N>1], by = subgroupvar][['sample_id']]
    n <- length(selectedsamples)
    if (verbose & n < ncol(object)){
        message('\t\tRetain ', length(selectedsamples), '/',
                ncol(object), ' samples with replicated ', subgroupvar) }
    object[, selectedsamples]
}


#==============================================================================
#
#                           log2transform()
#                           zscore()
#                           quantnorm()
#                           invnorm()
#
#==============================================================================

# mat <- cbind(s1=c(-1,-1), s2=c(-1,1), s3=c(1,-1), s4=c(0.1,0.1))
# which.medoid(mat)
which.medoid <- function(mat){
    idx <- matrixStats::rowAnyNAs(mat)
    assert_all_are_not_na(idx)
    if (any(idx))  cmessage('\t\t\t\tUse %d/%d non-NA rows to compute spatial median', 
                           sum(idx), length(idx))
    mat %<>% extract(!idx, )
    spatmed <- ICSNP::spatial.median(t(mat))
    which.min(sqrt(colSums((sweep(mat, 1, spatmed))^2)))
}

.filter_medoid <- function(object, verbose = FALSE){
    if (ncol(object)==1)  return(object)
    medoid <- which.medoid(values(object))
    object %<>% extract(, medoid)
    if (verbose)  message('\t\t\t\t', object$sample_id)
    object
}

#' Filter medoid sample
#' @param object SummarizedExperiment
#' @param by svar
#' @param verbose  whether to message
#' @return SummarizedExperiment
#' @examples 
#' file <- system.file('extdata/billing19.rnacounts.txt', package = 'autonomics')
#' object <- read_rnaseq_counts(file, plot = FALSE)
#' object %<>% filter_medoid(by = 'subgroup', verbose=TRUE)
#' @export
filter_medoid <- function(object, by = NULL, verbose = FALSE){
    if (!requireNamespace('ICSNP', quietly = TRUE)){
        message("`BiocManager::install('ICSNP')`. Then re-run.")
        return(object) }
    if (is.null(by))  return(.filter_medoid(object, verbose=verbose))
    object %<>% split_samples(by)
    if (verbose)  message('\t\t\tRetain medoid sample')
    object %<>% lapply(.filter_medoid, verbose=verbose)
    do.call(BiocGenerics::cbind, object)
}

.subtract_baseline <- function(
    object, subgroupvar, subgroupctr, assaynames
){
    idx <- object[[subgroupvar]] == subgroupctr
    controls <- object[, idx]
    perturbs <- object[, !idx]
    controls %<>% filter_medoid()
    for (ass in assaynames){
        assays(perturbs)[[ass]] %<>% sweep(1, assays(controls)[[ass]]) 
        # Confirm that sweep works as thought (it does)
        # controlmat <- assays(controls)[[ass]]
        # assert_is_identical_to_true(ncol(controlmat)==1)
        # controlmat %<>% extract(, 1)
        # controlmat %<>% matrix(byrow = FALSE, nrow = length(controlvec), ncol = ncol(perturbs))
        # assays(perturbs)[[ass]] %<>% subtract(controlmat)
    }
    perturbs
}

#' Subtract baseline
#' 
#' Subtract baseline level within block
#' 
#' \code{subtract_baseline} subtracts baseline levels within block, using the 
#' medoid baseline sample if multiple exist. \cr
#' 
#' \code{subtract_pairs} also subtracts baseline level within block. 
#' It cannot handle multiple baseline samples, but has instead been optimized
#' for many blocks \cr
#' 
#' \code{subtract_differences} subtracts differences between subsequent levels, 
#' again within block
#' @param  object       SummarizedExperiment
#' @param  subgroupvar  subgroup svar
#' @param  subgroupctr  control subgroup
#' @param  block        block svar (within which subtraction is performed)
#' @param  assaynames   which assays to subtract for
#' @param  verbose      TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' # read 
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object0 <- read_metabolon(file)
#'     pca(object0, plot = TRUE, color = 'Time')
#' 
#' # subtract_baseline: takes medoid of baseline samples if multiple
#'     object <- subtract_baseline(object0, block = 'Subject', subgroupvar = 'Time')
#'     pca(object, plot = TRUE, color = 'Time')
#' 
#' # subtract_pairs: optimized for many blocks
#'     object <- subtract_pairs(object0, block = 'Subject', subgroupvar = 'Time')
#'     pca(object, plot = TRUE, color = 'Time')
#' 
#' # subtract_differences
#'     object <- subtract_differences(object0, block = 'Subject', subgroupvar = 'Time')
#'     values(object) %<>% na_to_zero()
#'     pca(object, plot = TRUE, color = 'Time')
#' @export 
subtract_baseline <- function(
    object, subgroupvar, subgroupctr = slevels(object, subgroupvar)[1], 
    block = NULL, assaynames = setdiff(assayNames(object), c('weights', 'pepcounts')), 
    verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(subgroupvar, svars(object))
    assert_scalar_subset(subgroupctr, unique(object[[subgroupvar]]))
    if (!is.null(block))  assert_scalar_subset(block, svars(object))
    assert_is_subset(assaynames, assayNames(object))
    assert_is_a_bool(verbose)
# Subtract
    if (verbose){ 
        message("\t\tSubtract controls")
        message("\t\t\tcontrols  : ", subgroupvar, "=", subgroupctr," (medoid)")
        if (!is.null(block))  message("\t\t\tin block  : ", block)
        message("\t\t\tfor assays: ", paste0(assaynames, collapse = ', '))
    }
    objlist <- if (is.null(block)) list(object) else split_samples(object, block)
    objlist %<>% lapply(.subtract_baseline, 
                        subgroupvar = subgroupvar, subgroupctr = subgroupctr, 
                        assaynames = assaynames)
# Return
    objlist %<>% do.call(BiocGenerics::cbind, .)
    objlist
}


#' @rdname subtract_baseline
#' @export
subtract_pairs <- function(
    object, 
    subgroupvar = 'subgroup', 
    subgroupctr = slevels(object, subgroupvar)[1], 
    block,
    assaynames = assayNames(object)[1], verbose = TRUE
){
# Report
    # PRO: optimized for many block levels
    # CON: works only with exactly one ref per block 
    . <- NULL
    if (verbose){ 
        txt <- paste0("\t\tSubtract ", subgroupvar, "==", subgroupctr, ' ')
        txt %<>% paste0(paste0(assaynames, collapse = '/'))
        if (!is.null(block))  txt %<>% paste0(" per ", block)
        message(txt)
    }
# Ensure single ref per block
    sdt1 <- sdt(object)[, c('sample_id', subgroupvar, block), with = FALSE]
    singlerefperblock <- sdt1[, sum(get(subgroupvar)==subgroupctr)==1, by=block]$V1
    assert_is_identical_to_true(all(singlerefperblock))
# Subtract ref
    splitobjects <- split_samples(object, subgroupvar)
    refobj <- splitobjects[[subgroupctr]]
    splitobjects %<>% extract(setdiff(names(splitobjects), subgroupctr))
    splitobjects %<>% lapply(function(obj){
        idx <- match(obj[[block]], refobj[[block]])
        refobj %<>% extract(, idx)
        assert_is_identical_to_true(all(obj[[block]] == refobj[[block]]))
        for (assayname in assaynames){
            assays(obj)[[assayname]] %<>% subtract(assays(refobj)[[assayname]])
            assayNames(obj) %<>% stri_replace_first_fixed(
                                    assayname, paste0(assayname, 'ratios'))
        }
        obj        
    })
    splitobjects %<>% do.call(S4Vectors::cbind, .)
    idx <- na.exclude(match(object$sample_id, splitobjects$sample_id))
    splitobjects[, idx]
}


#' @rdname subtract_baseline
#' @export
subtract_differences <- function(object, block, subgroupvar, verbose=TRUE){
    # PRO: robust (works with missing levels, or multiple replicates per level)
    # CON: performance not optimized for many block levels
    #      currently only performed on first assay (can off course be updated)
    sample_id <- NULL
    if (verbose){ 
        message("\t\tSubtract differences")
        if (!is.null(block))  message("\t\t\tin block  : ", block)
        message("\t\t\tfor assays: ", assayNames(object)[1])
    }
    fvars0 <- intersect(c('feature_id', 'feature_name'), fvars(object))
    dt <- sumexp_to_longdt(object, fvars=fvars0, svars=c(subgroupvar, block))
    subgroups <- slevels(object, subgroupvar)
    n <- length(subgroups)
    formula <- paste0(c(fvars0, block), collapse = ' + ')
    formula %<>% paste0(' ~ ', subgroupvar)
    formula %<>% as.formula()
    dt %<>% dcast(formula, value.var ='value')
    
    newdt  <- dt[, c(fvars0, block), with = FALSE]
    diffdt <- dt[, setdiff(names(dt), c(fvars0, block)), with=FALSE]
    diffdt  <-  diffdt[, subgroups[-1], with=FALSE] - 
                diffdt[, subgroups[-n], with=FALSE]
    names(diffdt) %<>% paste0('_', subgroups[-n])
    newdt %<>% cbind(diffdt)
    
    newdt %<>% melt.data.table(
                id.vars =  c(fvars0, block), variable.name = subgroupvar)
    data.table::setorderv(newdt, c('feature_id', block, subgroupvar))
    newdt[, sample_id := paste0(get(block), '.', get(subgroupvar))]
    newobject <- dt2sumexp(newdt)
    assayNames(newobject)[1] <- assayNames(object)[1]
    newobject
}


#' Transform values
#' 
#' @param  object   SummarizedExperiment
#' @param  mat      matrix
#' @param  assay    character vector : assays for which to perform transformation
#' @param  pseudo   number           : pseudo value to be added prior to transformation
#' @param  verbose  TRUE or FALSE    : whether to msg
#' @param  delog    TRUE or FALSE (vsn)
#' @return Transformed sumexp
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#'
#' object                       %>% plot_sample_densities()
#' invnorm(object)              %>% plot_sample_densities()
#'
#' object                       %>% plot_sample_densities()
#' quantnorm(object)            %>% plot_sample_densities()
#'
#' object                       %>% plot_sample_densities()
#'#vsn(object)                  %>% plot_sample_densities()  # dataset too small
#'
#' object                       %>% plot_sample_densities()
#' zscore(object)               %>% plot_sample_densities()
#'
#' object                       %>% plot_sample_densities()
#' exp2(object)                 %>% plot_sample_densities()
#' log2transform(exp2(object))  %>% plot_sample_densities()
#' @export
log2transform <- function(
    object, 
    assay   = assayNames(object)[1], 
    pseudo  = 0,
    verbose = FALSE
){
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_not_null(assayNames(object))
    assert_is_subset(assay, assayNames(object))
    . <- NULL
    for (ass in assay){
        i <- match(ass, assayNames(object))
        if (verbose)  if (pseudo != 0)  message('\t\tAdd pseudo value ', pseudo)
        assays(object)[[i]] %<>% add(pseudo) 
         
        if (verbose)  cmessage('%slog2 %s', spaces(14), ass)
        assays(object)[[i]] %<>% log2()
        assayNames(object)[i] %<>% paste0('log2', .)
    }
    object
}


#' @rdname log2transform
#' @export
exp2 <- function(object, verbose = FALSE){
    . <- NULL
    if (verbose)  message('\t\tExp2 transform')
    values(object) %<>% magrittr::raise_to_power(2, .)
    object
}


#' @rdname log2transform
#' @export
zscore <- function(object, verbose = FALSE){
    values(object) %<>% sscale(verbose = verbose)
    object
}

#' @rdname log2transform
#' @export
sscale <- function(mat, verbose = FALSE){
    assert_is_matrix(mat)
    if (verbose)  message('\t\tZscore samples')
    scale(mat)
}

#' @rdname log2transform
#' @export
fscale <- function(mat, verbose = FALSE){
    assert_is_matrix(mat)
    if (verbose)  message('\t\tZscore features')
    t(scale(t(mat)))
}

#' Center samples
#' @param object   SummarizedExperiment
#' @param selector logical vector (length = nrow(object))
#' @param fun      aggregation function (string)
#' @param verbose  TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' require(matrixStats)
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' fdt(object)$housekeeping <- FALSE
#' fdt(object)$housekeeping[order(rowVars(values(object)))[1:5]] <- TRUE
#' values(object)[, object$subgroup=='Adult'] %<>% magrittr::add(5)
#' plot_sample_densities(object)
#' plot_sample_densities(center(object))
#' plot_sample_densities(center(object, housekeeping))
#' @export
center <- function(
    object, selector = rep(TRUE, nrow(object))==TRUE, fun = 'median', verbose = TRUE
){
    selector <- enexpr(selector)
    selector <- rlang::eval_tidy(selector, data = fdata(object))
    if (verbose)  message('\t\t', fun, ' center samples on ', 
                            nrow(object[selector, ]), ' features')
    correction_factors <- apply(values(object[selector, ]), 2, fun, na.rm=TRUE)
    correction_factors[is.na(correction_factors)] <- 0
    values(object) %<>% sweep(2, correction_factors)
    object
}


#' @rdname log2transform
#' @export
quantnorm <- function(object, verbose = FALSE){
    if (verbose)  message('\t\tQuantnorm')
    values(object) %<>% limma::normalizeBetweenArrays()
    object
}


#' @rdname log2transform
#' @export
invnorm <- function(object, verbose = FALSE){
    if (verbose)  message('Invnorm')
    values(object) %<>% apply(2, transform_to_fitting_normal)
    object
}

#' @rdname log2transform
#' @export
vsn <- function(object, verbose = FALSE, delog = TRUE){
    if (verbose) message('\t\tVSN')
    if (delog) object %<>% exp2()
    values(object) %<>% vsn::justvsn(verbose = FALSE)
    object
}

#' Transform vector to fitting normal distribution
#' @param x numeric vector
#' @return transformed vector
#' @noRd
transform_to_fitting_normal <- function(x){
    pars <- estimate_mean_sd(x)
    transform_to_normal(x, mean = pars[['mean']], sd = pars[['sd']])
}

estimate_mean_sd <- function(x){
    if (!requireNamespace('MASS', quietly = TRUE)){
        stop("BiocManager::install('MASS'). Then re-run.")
    }
    . <- NULL
    x %<>% extract(!is.na(.) & !is.infinite(.))
    MASS::fitdistr(x, 'normal')[['estimate']]
}



#' Transform vector to normal distribution
#' @param x numeric vector
#' @param mean  mean
#' @param sd    standard deviation
#' @return transformed vector
#' @noRd
transform_to_normal <- function(x, mean, sd){
    selector <- !is.na(x) & !is.nan(x) & !is.infinite(x)
    pvals <- rank(x[selector]) / (length(x[selector]) + 1)
    y <- x
    y[selector] <- qnorm(pvals, mean = mean, sd = sd)
    y
}


#' Transform vector to standard normal distribution
#' @param x numeric vector
#' @return transformed vector
#' @noRd
transform_to_standard_normal <- function(x){
    transform_to_normal(x, mean = 0, sd = 1)
}


#==============================================================================
#
#                        plot_transformation_biplots
#                        plot_transformation_densities
#                        plot_transformation_violins
#
#==============================================================================

gglegend<-function(p){
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(vapply(
                    tmp$grobs, function(x) x$name, character(1))=="guide-box")
    if (length(leg)==0)  grid::nullGrob()  else  tmp$grobs[[leg]]
}


plot_transformation_densities <- function(
    object,
    subgroupvar = 'subgroup',
    transformations = c('quantnorm', 'vsn' , 'zscore', 'invnorm'),
    ...,
    fixed = list(na.rm = TRUE, alpha = 0.3),
    nrow = 1, ncol = NULL
){
    value <- sample_id <- NULL
    assert_is_subset(subgroupvar, svars(object))
    dt <- sumexp_to_longdt(object, svars = c(subgroupvar))
    dt$transfo <- 'input'
    for (transfo in transformations){
        dt1 <- sumexp_to_longdt(get(transfo)(object), svars = c(subgroupvar))
        dt1$transfo <- transfo
        dt %<>% rbind(dt1)
    }
    dt$transfo %<>% factor(c('input', transformations))
    plot_data(dt, geom_density, x = value, group = sample_id,
            color = NULL, fill = !!sym(subgroupvar), ..., fixed = fixed) +
    facet_wrap(vars(transfo), scales = "free", nrow = nrow, ncol = ncol)
}


plot_transformation_violins <- function(
    object,
    subgroupvar = 'subgroup',
    transformations = c('quantnorm', 'vsn', 'zscore', 'invnorm'),
    ...,
    fixed = list(na.rm=TRUE)
){
    value <- sample_id <- NULL
    assert_is_subset(subgroupvar, svars(object))
    dt <- sumexp_to_longdt(object, svars = subgroupvar)
    dt$transfo <- 'input'
    for (transfo in transformations){
        dt1 <- sumexp_to_longdt(get(transfo)(object), svars = subgroupvar)
        dt1$transfo <- transfo
        dt %<>% rbind(dt1)
    }
    dt$transfo %<>% factor(c('input', transformations))
    plot_data(dt, geom_violin, x = sample_id, y = value, group = sample_id, 
                color = NULL, fill = !!sym(subgroupvar), ..., fixed = fixed) +
    facet_grid(rows = vars(!!sym(subgroupvar)), cols = vars(transfo), scales = "free") +
    coord_flip()
}


# file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
# object <- read_maxquant_proteingroups(file)
# plot_transformation_biplots(object, transformations = c('quantnorm', 'zscore', 'invnorm'))
plot_transformation_biplots <- function(
    object,
    subgroupvar = 'subgroup',
    transformations = c('quantnorm', 'vsn', 'zscore', 'invnorm'),
    method = c('pca', 'pls')[1], by = 'sample_id',
    dims = 1:2, color = subgroupvar, sep = FITSEP, ...,
    fixed = list(shape = 15, size = 3), nrow = 1, ncol = NULL
){
    . <- NULL
    assert_is_subset(subgroupvar, svars(object))
    assert_is_a_string(method)
    assert_is_subset(method, c('pca', 'pls'))
    assert_are_same_length(dims, 1:2)
    assert_is_numeric(dims)
    str_elem <- c(pca = 'by', pls = 'subgroupvar')
    xy <- paste0('effect', sep, get(str_elem[method]), sep, method, dims)
    x <- xy[1]; y <- xy[2]
    tmpobj <- object
    tmpobj %<>% get(method)(ndim = max(dims), verbose = FALSE)
    scoredt <- sdt(tmpobj) %>% cbind(transfo = 'input')
    mdidx <- paste0(get(str_elem[method]), sep, method)
    xvariance <- round(metadata(tmpobj)[[mdidx]][[paste0('effect', dims[1])]])
    yvariance <- round(metadata(tmpobj)[[mdidx]][[paste0('effect', dims[2])]])
    scoredt$transfo <- switch(
      method,
      pca = sprintf('input : %d + %d %%', xvariance, yvariance),
      pls = sprintf('input : %d %%'     , xvariance))
    for (transfo in transformations){
        tmpobj <- get(transfo)(object)
        tmpobj %<>% get(method)(dims = dims, verbose = FALSE)
        xvariance <- round(metadata(tmpobj)[[ mdidx ]][[ paste0('effect', dims[1]) ]])
        yvariance <- round(metadata(tmpobj)[[ mdidx ]][[ paste0('effect', dims[2]) ]])
        tmpdt <- sdt(tmpobj)
        tmpdt$transfo <- switch(
          method,
          pca = sprintf('%s : %d + %d %%', transfo, xvariance, yvariance),
          pls = sprintf('%s : %d %%',      transfo, xvariance))
        scoredt %<>% rbind(tmpdt)
    }
    scoredt$transfo %<>% factor(unique(.))
    p <- plot_data(
      scoredt, x = !!sym(x), y = !!sym(y), color = !!sym(color), ...,
                    fixed = fixed)
    p + facet_wrap(vars(transfo), nrow = nrow, ncol = ncol, scales = "free")
}


