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

# Filter for medoid sample
filter_medoid <- function(object){
    if (ncol(object)==1)  return(object)
    medoid <- which.medoid(exprs(object))
    object[, medoid]
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
#' It \cannot handle multiple baseline samples, but has instead been optimized
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
#'     exprs(object) %<>% na_to_zero()
#'     pca(object, plot=TRUE, color=SET)
#' @export 
subtract_baseline <- function(
    object, subgroupvar, subgroupctr = slevels(object, subgroupvar)[1], 
    block = NULL, assaynames = setdiff(assayNames(object), 'weights'), 
    verbose = TRUE
){
    if (verbose){ 
        cmessage("\t\tSubtract controls")
        cmessage("\t\t\tcontrols  : %s=%s (medoid)", subgroupvar, subgroupctr)
        if (!is.null(block))  cmessage("\t\t\tin block  : %s", block)
        cmessage("\t\t\tfor assays: %s", paste0(assaynames, collapse = ', '))
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
    assaynames = setdiff(assayNames(object), 'weights'), verbose = TRUE
){
# Report
    if (verbose){ 
        cmessage("\t\tSubtract pairs")
        cmessage("\t\t\tcontrols  : %s=%s (medoid)", subgroupvar, subgroupctr)
        if (!is.null(block))  cmessage("\t\t\tin block  : %s", block)
        cmessage("\t\t\tfor assays: %s", paste0(assaynames, collapse = ', '))
    }
# Ensure single ref per block
    sdata1 <- sdata(object)[, c('sample_id', subgroupvar, block)]
    sdata1 %<>% data.table()
    singlerefperblock <- sdata1[,sum(get(subgroupvar)==subgroupctr)==1,by=block]$V1
    assert_is_identical_to_true(all(singlerefperblock))
# Subtract ref
    splitobjects <- split_by_svar(object, !!sym(subgroupvar))
    refobj <- splitobjects[[subgroupctr]]
    splitobjects %<>% extract(-1)
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
    if (verbose){ 
        cmessage("\t\tSubtract differences")
        if (!is.null(block))  cmessage("\t\t\tin block  : %s", block)
        cmessage("\t\t\tfor assays: %s", paste0(assaynames, collapse = ', '))
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
    diffdt <- diffdt[, subgroups[-1], with=FALSE] - 
              diffdt[, subgroups[-n], with=FALSE]
    names(diffdt) %<>% paste0(' - ', subgroups[-n])
    newdt %<>% cbind(diffdt)
    
    newdt %<>% data.table::melt.data.table(
                id.vars =  c(fvars0, block), variable.name = subgroupvar)
    data.table::setorderv(newdt, c('feature_id', block, subgroupvar))
    newdt[, sample_id := paste0(get(block), '.', get(subgroupvar))]
    newobject <- dt2sumexp(newdt)
    assayNames(newobject)[1] <- assayNames(object)[1]
    newobject
}


#' Log2 transform
#' @param object SummarizedExperiment
#' @param verbose TRUE or FALSE
#' @return Transformed sumexp
#' @examples
#' require(magrittr)
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#'
#' plot_sample_densities(object)
#' plot_sample_densities(invnorm(object))
#'
#' plot_sample_densities(object)
#' plot_sample_densities(quantnorm(object))
#'
#' plot_sample_densities(object)
#' plot_sample_densities(zscore(object))
#'
#' plot_sample_densities(object)
#' plot_sample_densities(exp2(object))
#' plot_sample_densities(log2transform(exp2(object)))
#' @export
log2transform <- function(object, verbose = FALSE){
    if (verbose)  message('\t\tLog2 transform')
    exprs(object) %<>% log2()
    object
}

#' @rdname log2transform
#' @export
exp2 <- function(object, verbose = FALSE){
    if (verbose)  message('\t\tExp2 transform')
    exprs(object) %<>% magrittr::raise_to_power(2, .)
    object
}


#' @rdname log2transform
#' @export
zscore <- function(object, verbose = FALSE){
    if (verbose)  message('\t\tZscore')
    exprs(object) %<>% scale()
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
#' file <- download_data('billing19.proteingroups.txt')
#' select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
#' select %<>% paste0('_STD')
#' object <- read_proteingroups(file, select_subgroups = select, plot=FALSE)
#' object %<>% extract(, order(object$subgroup))
#' fdata(object)$housekeeping <- FALSE
#' fdata(object)$housekeeping[order(rowVars(exprs(object)))[1:100]] <- TRUE
#' exprs(object)[, object$subgroup=='BM00_STD'] %<>% add(5)
#' gridExtra::grid.arrange(plot_sample_densities(object),
#'                         plot_sample_densities(center(object)),
#'                         plot_sample_densities(center(object, housekeeping)))
#' @export
center <- function(object, selector = rep(TRUE, nrow(object))==TRUE,
                    fun = 'median', verbose = TRUE
){
    selector <- enexpr(selector)
    selector <- rlang::eval_tidy(selector, data = fdata(object))
    if (verbose) cmessage('%s center samples on %d features',
                        fun, nrow(object[selector, ]))
    correction_factors <- apply(exprs(object[selector, ]), 2, fun, na.rm=TRUE)
    correction_factors[is.na(correction_factors)] <- 0
    exprs(object) %<>% sweep(2, correction_factors)
    object
}


#' @rdname log2transform
#' @export
quantnorm <- function(object, verbose = FALSE){
    if (verbose)  message('\t\tQuantnorm')
    exprs(object) %<>% limma::normalizeBetweenArrays()
    object
}


#' @rdname log2transform
#' @export
invnorm <- function(object, verbose = FALSE){
    if (verbose)  message('Invnorm')
    exprs(object) %<>% apply(2, transform_to_fitting_normal)
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
        message("BiocManager::install('MASS'). Then re-run.")
        return(x)
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
    legend <- tmp$grobs[[leg]]
    return(legend)
}


plot_transformation_densities <- function(
    object,
    transformations = c('quantnorm', 'zscore', 'invnorm'),
    x = value, fill = subgroup,
    group = if (ncol(object) > 16)  subgroup  else  sample_id, ...,
    fixed = list(na.rm=TRUE, alpha=0.3)
){
    dt <- sumexp_to_long_dt(object)
    dt$transfo <- 'input'
    for (transfo in transformations){
        dt1 <- sumexp_to_long_dt(get(transfo)(object))
        dt1$transfo <- transfo
        dt %<>% rbind(dt1)
    }
    dt$transfo %<>% factor(c('input', transformations))
    plot_data(dt, geom_density, x = !!enquo(x), group = !!enquo(group),
            color = NULL, fill = !!enquo(fill), ..., fixed = fixed) +
    facet_wrap(~transfo, scales = "free", nrow=1)
}


plot_transformation_violins <- function(
    object,
    transformations = c('quantnorm', 'zscore', 'invnorm'),
    x = subgroup, y = value, fill = subgroup,
    group = if (ncol(object) > 16)  subgroup  else  sample_id, ...,
    fixed = list(na.rm=TRUE)
){
    dt <- sumexp_to_long_dt(object)
    dt$transfo <- 'input'
    for (transfo in transformations){
        dt1 <- sumexp_to_long_dt(get(transfo)(object))
        dt1$transfo <- transfo
        dt %<>% rbind(dt1)
    }
    dt$transfo %<>% factor(c('input', transformations))
    plot_data(dt, geom_violin, x = !!enquo(x), y = !!enquo(y),
            group = !!enquo(group), color = NULL, fill = !!enquo(fill),
            ..., fixed = fixed) +
    facet_grid(subgroup~transfo, scales = "free") +
    coord_flip()
}


plot_transformation_biplots <- function(
    object,
    transformations = c('quantnorm', 'zscore', 'invnorm'),
    method = 'pca', xdim = 1, ydim = 2, color = subgroup, ...,
    fixed = list(shape=15, size=3)
){
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
#' @param transformations vector
#' @param method          'pca', 'pls', 'sma', or 'lda'
#' @param xdim            number (default 1)
#' @param ydim            number (default 2)
#' @param ...             passed to plot_data
#' @return grid object
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' explore_transformations(object)
#' @export
explore_transformations <- function(
    object, transformations = c('quantnorm', 'zscore', 'invnorm'),
    method='pca', xdim=1, ydim=2, ...
){
    p1 <- plot_transformation_densities(object, transformations, ...)
    p2 <- plot_transformation_biplots(
            object, transformations, method, xdim1 = xdim, ydim = ydim, ...)

    grid.arrange(arrangeGrob(p1 + theme(legend.position='none'),
                            p2  + theme(legend.position='none'), nrow=2),
                            gglegend(p2), ncol=2, widths = c(8, 1))

}

