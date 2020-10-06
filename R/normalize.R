#==============================================================================
#
#                           .zscore()
#                           .quantnorm()
#                           .invnorm()
#
#==============================================================================

.zscore <- function(object){
    exprs(object) %<>% scale()
    object
}

.quantnorm <- function(object){
    exprs(object) %<>% limma::normalizeBetweenArrays()
    object
}


.invnorm <- function(object){
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

#' Plot transformation densities
#'
#' Plot sample densities for alternative transformations
#'
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
#' plot_transformation_densities(object)
#' @export
plot_transformation_densities <- function(
    object,
    transformations = c('.quantnorm', '.zscore', '.invnorm'),
    geom = geom_violin, x = subgroup, y = value,
    group = if (ncol(object) > 16)  subgroup  else  sample_id, ...,
    fixed = list(na.rm=TRUE)
){
    dt <- sumexp_to_long_dt(object)
    dt$transfo <- 'input'
    for (transfo in transformations){
        object <- get(transfo)(object)
        dt1 <- sumexp_to_long_dt(object)
        dt1$transfo <- transfo
        dt %<>% rbind(dt1)
    }
    dt$transfo %<>% factor(c('input', transformations))
    plot_data(dt, geom, x = !!enquo(x), y = !!enquo(y), group = !!enquo(group), color = NULL,
              fill = subgroup, ..., fixed = fixed) +
    facet_grid(subgroup~transfo, scales = "free") +
    coord_flip()
}


#==============================================================================
#
#                       invnorm
#
#==============================================================================


#' Inverse normal transform samples
#' @param object SummarizedExperiment
#' @param plot  TRUE (default) or FALSE
#' @param ...   additional aesthetic mappings
#' @return normalized SummarizedExperimen
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' invnorm(object)
#' quantnorm(object)
#' zscore(object)
#' @export
invnorm <- function(object, plot = TRUE, ...){
    if (plot) print(plot_transformation_densities(object, '.invnorm', ...))
    .invnorm(object)
}


#' @rdname invnorm
#' @export
quantnorm <- function(object, plot = TRUE, ...){
    if (plot) print(plot_transformation_densities(object, '.quantnorm', ...))
    .quantnorm(object)
}


#' @rdname invnorm
#' @export
zscore <- function(object, plot = TRUE, ...){
    if (plot) print(plot_transformation_densities(object, '.zscore', ...))
    .zscore(object)
}


