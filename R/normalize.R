#==============================================================================
#
#                           .log2transform()
#                           .zscore()
#                           .quantnorm()
#                           .invnorm()
#
#==============================================================================


#' Log2 transform exprs
#' @param object SummarizedExperiment
#' @param object SummarizedExperiment
#' @param assays 'exprs', etc.
#' @param verbose TRUE or FALSE
#' @return Updated SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, log2 = FALSE, plot = FALSE)
#' object %<>% .log2transform()
#' .exp2(object)
#' .zscore(object)
#' .quantnorm(object)
#' .invnorm(object)
#' @noRd
.log2transform <- function(object, assays = 'exprs', verbose = FALSE){
    for (assay in assays){
        if (verbose)  message('\t\tLog2 ', assay)
        SummarizedExperiment::assays(object)[[assay]]  %<>% log2()
    }
    object
}

.exp2 <- function(object, assays = 'exprs', verbose = FALSE){
    for (assay in assays){
        if (verbose)  message('\t\tExp2 ', assay)
        SummarizedExperiment::assays(object)[[assay]]  %<>%
            magrittr::raise_to_power(2, .)
    }
    object
}


.zscore <- function(object, assays = 'exprs', verbose = FALSE){
    for (assay in assays){
        if (verbose)  message('\t\tZscore ', assay)
        SummarizedExperiment::assays(object)[[assay]]  %<>% scale()
    }
    object
}

.quantnorm <- function(object, assays = 'exprs', verbose = FALSE){
    for (assay in assays){
        if (verbose)  message('\t\tQuantnorm ', assay)
        SummarizedExperiment::assays(object)[[assay]] %<>%
            limma::normalizeBetweenArrays()
    }
    object
}


.invnorm <- function(object, assays = 'exprs', verbose = FALSE){
    for (assay in assays){
        if (verbose)  message('Invnorm ', assay)
        SummarizedExperiment::assays(object)[[assay]] %<>%
            apply(2, transform_to_fitting_normal)
    }
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
    . <- NULL
    x %<>% extract(!is.na(.) & !is.infinite(.))
    fitdistr(x, 'normal')[['estimate']]
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
#                       plot_transformation_densities()
#
#==============================================================================

#' Plot transformation effects
#' @param object          SummarizedExperiment
#' @param geom            function
#' @param x               svar mapped to x
#' @param group           svar mapped to group
#' @param ...             additional svar to aesthetic mappings
#' @param tranformations  character vector
#' @param fixed           list
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' plot_transformation_violins(object)
#' plot_transformation_scores(object)
#' @export
plot_transformation_violins <- function(
    object,
    transformations = c('.quantnorm', '.zscore', '.invnorm'),
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


plot_transformation_densities <- function(
    object,
    transformations = c('.quantnorm', '.zscore', '.invnorm'),
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


#' @rdname plot_transformation_violins
#' @export
plot_transformation_scores <- function(
    object,
    transformations = c('.quantnorm', '.zscore', '.invnorm'),
    method = 'pca', xdim = 1, ydim = 2, color = subgroup, ...,
    fixed = list(shape=15, size=3)
){
    color <- enquo(color)
    xstr <- paste0(method, xdim)
    ystr <- paste0(method, ydim)
    add_projection <- get(paste0('.add_', method))

    object %<>% add_projection(verbose = FALSE, ndim=max(xdim, ydim))
    scoredt <- sdata(object) %>% cbind(transfo = 'input')
    var1 <- metadata(object)[[method]][[xstr]]
    var2 <- metadata(object)[[method]][[ystr]]
    scoredt$transfo <- sprintf('input : %d + %d %%', var1, var2)
    for (transfo in transformations){
        tmpobj <- get(transfo)(object) %>% add_projection(verbose=FALSE, ndim=max(xdim, ydim))
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

gglegend<-function(p){
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

plot_transformations <- function(
    object, transformations = c('.quantnorm', '.zscore', '.invnorm'),
    method='pca', xdim=1, ydim=2, ...
){
    p1 <- plot_transformation_densities(object, transformations, ...)
    p2 <- plot_transformation_scores(
            object, transformations, method, xdim1 = xdim, ydim = ydim, ...)

    gridExtra::grid.arrange(gridExtra::arrangeGrob(
                                p1 + theme(legend.position='none'),
                                p2 + theme(legend.position='none'), nrow=2),
                            gglegend(p2), ncol=2, widths = c(8, 1))

}

#==============================================================================
#
#                       log2transform
#                       invnorm
#                       quantnorm
#                       zscore
#
#==============================================================================


#' Inverse normal transform samples
#' @param object SummarizedExperiment
#' @param plot  TRUE (default) or FALSE
#' @param ...   passed on to plot_transformations
#' @return normalized SummarizedExperimen
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE, log2=FALSE)
#'
#' object %<>% log2transform()
#' invnorm(object)
#' quantnorm(object)
#' zscore(object)
#' @export
log2transform <- function(
    object, assays = 'exprs', plot = TRUE, verbose = TRUE, ...
){
    if (plot) print(plot_transformations(object, '.log2transform', ...))
    .log2transform(object, assays=assays, verbose=verbose)
}


#' @rdname log2transform
#' @export

invnorm <- function(
    object, assays = 'exprs', plot = TRUE, verbose = TRUE, ...
){
    if (plot) print(plot_transformations(object, '.invnorm', ...))
    .invnorm(object, assays=assays, verbose=verbose)
}


#' @rdname log2transform
#' @export
quantnorm <- function(
    object, assays = 'exprs', plot = TRUE, verbose = TRUE, ...
){
    if (plot) print(plot_transformations(object, '.quantnorm', ...))
    .quantnorm(object, assays=assays, verbose=verbose)
}


#' @rdname log2transform
#' @export
zscore <- function(
    object, assays = 'exprs', plot = TRUE, verbose = TRUE, ...
){
    if (plot) print(plot_transformations(object, '.zscore', ...))
    .zscore(object, assays=assays, verbose=verbose)
}


#' @rdname log2transform
#' @export
exp2 <- function(
    object, assays = 'exprs', plot = TRUE, verbose = TRUE, ...
){
    if (plot) print(plot_transformations(object, '.log2transform', ...))
    .exp2(object, assays=assays, verbose=verbose)
}


