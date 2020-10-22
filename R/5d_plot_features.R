#=============================================================================
#
#                 plot_feature_boxplots()
#
#=============================================================================

#' Plot features
#' @param object      SummarizedExperiment
#' @param geom        geom_point, geom_boxplot, etc.
#' @param x           svar mapped to x
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param ...         mapped aesthetics
#' @param fixed       fixed aesthetics
#' @param theme       ggplot theme specifications
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% pca()
#' idx <- order(abs(fdata(object)$pca1), decreasing=TRUE)[1:9]
#' object %<>% extract(idx, )
#' plot_feature_boxplots(object)
#' plot_feature_profiles(object)
#' @export
plot_features <- function(
    object,
    geom,
    x          = subgroup,
    fill       = subgroup,
    color      = subgroup,
    ...,
    fixed = list(na.rm=TRUE),  #element_text(angle=90, vjust=0.5),
    theme = list(axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank())
){
    fill  <- enquo(fill)
    color <- enquo(color)
    x     <- enquo(x)
    dt <- sumexp_to_long_dt(
            object, svars = svars(object), fvars = fvars(object))
    value <- NULL
    p <- plot_data(
            dt, geom = geom, x = !!x, y = value, fill = !!fill,
            color = !!color, ..., fixed = fixed)
    p <- p + facet_wrap(~ feature_name, scales = 'free_y')
    p <- p + do.call(ggplot2::theme, theme)
    p
}

#' @rdname plot_features
#' @export
plot_feature_boxplots <- function(...){
    plot_features(geom = geom_boxplot, color = NULL, ...)
}


#' @rdname plot_features
#' @export
plot_feature_profiles <- function(...){
    plot_features(geom = geom_point, ...)
}


