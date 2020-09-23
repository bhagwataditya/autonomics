#==============================================================================
#
#                                invnorm()
#
#==============================================================================


#' Inverse normal transform samples
#' @param object SummarizedExperiment
#' @return normalized SummarizedExperimen
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' invnorm(object)
#' @export
invnorm <- function(object){
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





