#=============================================================================
#
#                 plot_sample_densities()
#                 plot_sample_violins()
#                 plot_sample_boxplots()
#
#=============================================================================

#' Plot sample densities
#' @param object      SummarizedExperiment
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param ...         mapped aesthetics
#' @param fixed       fixed aesthetics
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' plot_sample_densities(object)
#' plot_sample_densities(object, color = subgroup, fill = NULL)
#' @export
plot_sample_densities <- function(
    object,
    fill       = subgroup,
    color      = NULL,
    ...,
    fixed = list(alpha = 0.5, na.rm = TRUE)
){
    dt <- sumexp_to_long_dt(object, svars = svars(object))
    fill <- enquo(fill)
    color <- enquo(color)
    value <- NULL
    plot_data(
        dt, geom = geom_density, x = value, fill = !!fill,
        color= !!color, ..., fixed = fixed)
}


#' Plot sample violins
#' @param object      SummarizedExperiment
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param ...         mapped aesthetics
#' @param fixed       fixed aesthetics
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' plot_sample_violins(object)
#' @export
plot_sample_violins <- function(
    object,
    fill       = subgroup,
    color      = NULL,
    ...,
    fixed = list(na.rm=TRUE)
){
    dt <- sumexp_to_long_dt(object, svars = svars(object))
    fill <- enquo(fill)
    color <- enquo(color)
    sample_id <- value <- NULL
    plot_data(
        dt, geom = geom_violin, x = sample_id, y = value, fill = !!fill,
        color= !!color, ..., fixed      = fixed)
}



#' Plot sample boxplots
#' @param object      SummarizedExperiment
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param ...         mapped aesthetics
#' @param fixed       fixed aesthetics
#' @return  ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' plot_sample_boxplots(object)
#' @export
plot_sample_boxplots <- function(
    object,
    fill       = subgroup,
    color      = NULL,
    ...,
    fixed = list(na.rm=TRUE)
){
    dt <- sumexp_to_long_dt(object, svars = svars(object))
    fill <- enquo(fill)
    color <- enquo(color)
    sample_id <- value <- NULL
    plot_data(
        dt, geom = geom_boxplot, x = sample_id, y = value, fill = !!fill,
        color = !!color, ..., fixed      = fixed)
}


