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
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' object %<>% filter_samples(SampleGroup %in% c('t1', 't2'), verbose = TRUE)
#' rm_singleton_samples(object, svar = 'Subject_ID')
#' rm_unmatched_samples(object, subgroupvar = 'SampleGroup', block = 'Subject_ID')
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
rm_singleton_samples <- function(object, svar = 'subgroup', verbose = TRUE){
    selectedsamples <- sdt(object)[, .SD[.N>1], by = svar][['sample_id']]
    n <- length(selectedsamples)
    if (verbose & n < ncol(object)){
        message('\t\tRetain ', length(selectedsamples), '/',
                ncol(object), ' samples with replicated ', svar) }
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
#' require(magrittr)
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' object %<>% filter_medoid(by = 'subgroup', verbose=TRUE)
#' @export
filter_medoid <- function(object, by = NULL, verbose = FALSE){
    if (!requireNamespace('ICSNP', quietly = TRUE)){
        message("`BiocManager::install('ICSNP')`. Then re-run.")
        return(object) }
    if (is.null(by))  return(.filter_medoid(object, verbose=verbose))
    object %<>% split_by_svar(!!sym(by))
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
        assays(perturbs)[[ass]] %<>% sweep(1, assays(controls)[[ass]]) }
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
#'     require(magrittr)
#'     file <- download_data('atkin18.metabolon.xlsx') 
#'     object0 <- read_metabolon(file, plot=FALSE)
#'     pca(object0, plot=TRUE, color=SET)
#' 
#' # subtract_baseline: takes medoid of baseline samples if multiple
#'     object <- subtract_baseline(object0, block='SUB', subgroupvar='SET')
#'     pca(object, plot=TRUE, color=SET)
#' 
#' # subtract_pairs: optimized for many blocks
#'     object <- subtract_pairs(   object0, block='SUB', subgroupvar='SET')
#'     pca(object, plot=TRUE, color=SET)
#' 
#' # subtract differences
#'     object <- subtract_differences(object0, block='SUB', subgroupvar='SET')
#'     values(object) %<>% na_to_zero()
#'     pca(object, plot=TRUE, color=SET)
#' @export 
subtract_baseline <- function(
    object, subgroupvar, subgroupctr = slevels(object, subgroupvar)[1], 
    block = NULL, assaynames = setdiff(assayNames(object), 'weights'), 
    verbose = TRUE
){
    if (verbose){ 
        message("\t\tSubtract controls")
        message("\t\t\tcontrols  : ", subgroupvar, "=", subgroupctr," (medoid)")
        if (!is.null(block))  message("\t\t\tin block  : ", block)
        message("\t\t\tfor assays: ", paste0(assaynames, collapse = ', '))
    }
    objlist <- if (is.null(block)) list(object) else split_by_svar(
                                                        object, !!sym(block)) 
    objlist %<>% lapply(.subtract_baseline, 
                        subgroupvar = subgroupvar, subgroupctr = subgroupctr, 
                        assaynames = assaynames)
    do.call(BiocGenerics::cbind, objlist)
}


#' @rdname subtract_baseline
#' @export
subtract_pairs <- function(
    object, subgroupvar, subgroupctr = slevels(object, subgroupvar)[1], block,
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
    splitobjects <- split_by_svar(object, !!sym(subgroupvar))
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
    dt <- sumexp_to_long_dt(object, fvars=fvars0, svars=c(subgroupvar, block))
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
#' @param  assay    character vector : assays for which to perform transformation
#' @param  pseudo   number           : pseudo value to be added prior to transformation
#' @param  verbose  TRUE/FALSE       : whether to msg
#' @return Transformed sumexp
#' @examples
#' require(magrittr)
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE, impute=FALSE)
#'
#' object                       %>% plot_sample_densities()
#' invnorm(object)              %>% plot_sample_densities()
#'
#' object                       %>% plot_sample_densities()
#' quantnorm(object)            %>% plot_sample_densities()
#'
#' object                       %>% plot_sample_densities()
#' vsn(object)                  %>% plot_sample_densities()
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
    assert_is_subset(assay, assayNames(object))
    . <- NULL
    for (ass in assay){
        i <- match(ass, assayNames(object))
        if (verbose)  message('\t\tAdd pseudo value ', pseudo) ;  assays(object)[[i]] %<>% add(pseudo) 
        if (verbose)  message('\t\tLogarithmize (log2) ', ass);      assays(object)[[i]] %<>% log2()
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
    if (verbose)  message('\t\tZscore')
    values(object) %<>% scale()
    object
}


#' Center samples
#' @param object   SummarizedExperiment
#' @param selector logical vector (length = nrow(object))
#' @param fun      aggregation function (string)
#' @param verbose  TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' require(matrixStats)
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE, impute=FALSE)
#' fdata(object)$housekeeping <- FALSE
#' fdata(object)$housekeeping[order(rowVars(values(object)))[1:100]] <- TRUE
#' values(object)[, object$subgroup=='Adult'] %<>% add(5)
#' plot_sample_densities(object)
#' plot_sample_densities(center(object))
#' plot_sample_densities(center(object, housekeeping))
#' @export
center <- function(object, selector = rep(TRUE, nrow(object))==TRUE,
                    fun = 'median', verbose = TRUE
){
    selector <- enexpr(selector)
    selector <- rlang::eval_tidy(selector, data = fdata(object))
    if (verbose)  message(fun, ' center samples on ', 
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
    subgroup = sym('subgroup'),
    transformations = c('quantnorm', 'vns' , 'zscore', 'invnorm'),
    ...,
    fixed = list(na.rm=TRUE, alpha=0.3)
){
    value <- sample_id <- NULL
    subgroup <- enquo(subgroup); subgroupvar <- as_name(subgroup)
    assert_is_subset(subgroupvar, svars(object))
    dt <- sumexp_to_long_dt(object, svars = c(subgroupvar))
    dt$transfo <- 'input'
    for (transfo in transformations){
        dt1 <- sumexp_to_long_dt(
                get(transfo)(object), svars = c(subgroupvar))
        dt1$transfo <- transfo
        dt %<>% rbind(dt1)
    }
    dt$transfo %<>% factor(c('input', transformations))
    plot_data(dt, geom_density, x = value, group = sample_id,
            color = NULL, fill = !!subgroup, ..., fixed = fixed) +
    facet_wrap(~transfo, scales = "free", nrow=1)
}


plot_transformation_violins <- function(
    object,
    subgroup = sym('subgroup'),
    transformations = c('quantnorm', 'vsn', 'zscore', 'invnorm'),
    ...,
    fixed = list(na.rm=TRUE)
){
    value <- sample_id <- NULL
    subgroupvar <- as_name(subgroup); subgroup <- enquo(subgroup)
    assert_is_subset(subgroupvar, svars(object))
    dt <- sumexp_to_long_dt(object, svars = subgroupvar)
    dt$transfo <- 'input'
    for (transfo in transformations){
        dt1 <- sumexp_to_long_dt(
                get(transfo)(object), svars = subgroupvar)
        dt1$transfo <- transfo
        dt %<>% rbind(dt1)
    }
    dt$transfo %<>% factor(c('input', transformations))
    plot_data(dt, geom_violin, x = sample_id, y = value, group = sample_id, 
                color = NULL, fill = !!subgroup, ..., fixed = fixed) +
    facet_grid(rows = vars(!!subgroup), cols = vars(transfo), scales = "free") +
    coord_flip()
}


plot_transformation_biplots <- function(
    object,
    subgroup = sym('subgroup'),
    transformations = c('quantnorm', 'vsn', 'zscore', 'invnorm'),
    method = 'pca', xdim = 1, ydim = 2, color = !!enquo(subgroup), ...,
    fixed = list(shape=15, size=3)
){
    . <- NULL
    subgroup <- enquo(subgroup)
    subgroupvar <- as_name(subgroup)
    assert_is_subset(subgroupvar, svars(object))
    color <- enquo(color)
    xstr <- paste0(method, xdim)
    ystr <- paste0(method, ydim)

    object %<>% get(method)(ndim=max(xdim, ydim), verbose = FALSE)
    scoredt <- sdata(object) %>% cbind(transfo = 'input')
    var1 <- metadata(object)[[method]][[xstr]]
    var2 <- metadata(object)[[method]][[ystr]]
    scoredt$transfo <- sprintf('input : %d + %d %%', var1, var2)
    for (transfo in transformations){
        tmpobj <- get(transfo)(object) %>%
                get(method)(ndim=max(xdim, ydim), verbose=FALSE)
        var1 <- metadata(tmpobj)[[method]][[xstr]]
        var2 <- metadata(tmpobj)[[method]][[ystr]]
        tmpdt <- sdata(tmpobj)
        tmpdt$transfo <- sprintf('%s : %d + %d %%', transfo, var1, var2)
        scoredt %<>% rbind(tmpdt)
    }
    scoredt$transfo %<>% factor(unique(.))
    p <- plot_data( scoredt, x=!!sym(xstr), y=!!sym(ystr), color=!!color, ...,
                    fixed = fixed)
    p + facet_wrap(~transfo, nrow=1, scales = "free")
}


#==============================================================================
#
#                    explore_transformations
#
#==============================================================================

#' Explore transformations
#' @param object          SummarizedExperiment
#' @param subgroup        subgroup (sym)
#' @param transformations vector
#' @param method          'pca', 'pls', 'sma', or 'lda'
#' @param xdim            number (default 1)
#' @param ydim            number (default 2)
#' @param ...             passed to plot_data
#' @return grid object
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' invert <- c('EM_E', 'EM_BM', 'BM_E')
#' object <- read_proteingroups(file, invert_subgroups = invert, plot=FALSE)
#' explore_transformations(object)
#' @export
explore_transformations <- function(
    object, 
    subgroup = subgroup,
    transformations = c('quantnorm', 'vsn', 'zscore', 'invnorm'),
    method='pca', xdim=1, ydim=2, ...
){
    subgroup <- enquo(subgroup)
    subgroupvar <- as_name(subgroup)
    assert_is_subset(subgroupvar, svars(object))
    p1 <- plot_transformation_densities(
            object, subgroup=!!subgroup, transformations, ...)
    p2 <- plot_transformation_biplots(
            object, subgroup=!!subgroup, transformations, method, 
            xdim1 = xdim, ydim = ydim, ...)

    p3 <- arrangeGrob(
        p1 + theme(legend.position='none'),
        p2  + theme(legend.position='none'),
        nrow=2, right = gglegend(p2))
    grid.draw(p3)
    invisible(p3)
}

