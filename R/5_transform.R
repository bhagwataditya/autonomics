#==============================================================================
#
#                           log2transform()
#                           zscore()
#                           quantnorm()
#                           invnorm()
#
#==============================================================================


#' Log2 transform
#' @param object SummarizedExperiment
#' @param verbose TRUE or FALSE
#' @return Transformed sumexp
#' @examples
#' require(magrittr)
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, pca=FALSE, lmfit=FALSE, plot=FALSE)
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
#' object <- read_proteingroups(
#'             file, select_subgroups = select, pca=FALSE, lmfit=FALSE, plot=FALSE)
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
#' object <- read_metabolon(file, pca=FALSE, lmfit=FALSE, plot = FALSE)
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

