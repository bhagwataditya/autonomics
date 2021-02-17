#=============================================================================
#
#                    add_color_scale
#                    add_fill_scale
#
#==============================================================================

#' Add color scale
#' @param object   SummarizedExperiment
#' @param color    symbol: svar mapped to color
#' @param show     TRUE or FALSE (default)
#' @param verbose  TRUE or FALSE (default)
#' @return default color values vector
#' @examples
#' # STEMCELL RATIOS
#'     file <- download_data('billing16.proteingroups.txt')
#'     invert_subgroups <- c('E_EM','BM_E', 'BM_EM')
#'     object <- read_proteingroups(
#'                file, invert_subgroups = invert_subgroups, plot=FALSE)
#'     p <- plot_sample_densities(object)
#'     add_color_scale(p, data=sdata(object))
#' # STEMCELL INTENSITIES
#'    file <- download_data('billing16.proteingroups.txt')
#'    object <- read_proteingroups(
#'               file, quantity = 'Intensity labeled', plot=FALSE)
#'    add_color_scale(object)
#' @noRd
add_color_scale <- function(p, color = subgroup, data){
# Assert
    assert_is_data.frame(data)
    color <- enquo(color)
# Colors
    if (!rlang::quo_is_null(color)){
        color_var <- as_name(color)
        assert_is_subset(color_var, names(data))
        values0 <- data[[color_var]]
        if (!is.numeric(values0)){
            if (is.character(values0)) values0 %<>% factor()
            levels0 <- levels(values0)
            colors0 <- make_colors(levels0, sep = guess_sep(levels0))
            p <- p + scale_color_manual(values = colors0)
        }
    }
# Return
    return(p)
}


add_fill_scale <- function(p, fill = subgroup, data){
# Assert
    assert_is_data.frame(data)
    fill <- enquo(fill)
# Colors
    if (!rlang::quo_is_null(fill)){
        fillstr <- as_name(fill)
        assert_is_subset(fillstr, names(data))
        values0 <- data[[fillstr]]
        if (!is.numeric(values0)){
            levels0 <- unique(values0)
            colors0 <- make_colors(levels0, sep = guess_sep(levels0))
            p <- p + scale_fill_manual(values = colors0)
        }
    }
# Return
    return(p)
}


make_colors <- function(
    varlevels, sep = guess_sep(varlevels), show=FALSE,
    verbose = FALSE
){
    makefun <- make_onefactor_colors
    if (!is.null(sep)){            # consistent separator
        if (length(varlevels)>2){  # 3+ samples
            n1 <- length(unique(split_extract(varlevels, 1, sep)))
            n2 <- length(unique(split_extract(varlevels, 2, sep)))
            if (n1>1 & n2>1){             # 2+ huevar levels
                makefun <- make_twofactor_colors
            }
        }
    }
    makefun(varlevels, sep = sep, show = show, verbose = verbose)
}

#' Create default ggplot colors for factor levels
#' @param varlevels  string vector
#' @param show           TRUE/FALSE
#' @param verbose        TRUE/FALSE
#' @return string vector: elements = colors, names = factor levels
#' @author John Colby
#' @references https://stackoverflow.com/questions/8197559
#' @noRd
make_onefactor_colors <- function(
    varlevels, show, verbose = TRUE, sep = NULL
){
    n <- length(varlevels)
    hues <- seq(15, 375, length = n + 1)
    colors <- hcl(h = hues, l = 65, c = 100)[seq_len(n)] %>%
                    set_names(varlevels)
    if (show) pie(rep(1, length(colors)), names(colors),
                    col = colors)
    if (verbose)  cmessage('\t\tMake default ggplot colors')
    colors
}


#' Make composite colors
#' @param varlevels string vector
#' @param sep     string
#' @param show    TRUE/FALSE: show colors in pie plot?
#' @param verbose TRUE/FALSE
#' @return named string vector (elements = colors, names = color_var levels)
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' varlevels <- subgroup_levels(object)
#' make_twofactor_colors(varlevels, show = TRUE)
#' @noRd
make_twofactor_colors <- function(
    varlevels, sep  = guess_sep(varlevels), show = FALSE, verbose = TRUE
){
    # Assert
    assert_has_no_duplicates(varlevels)
    assert_is_not_null(sep)
    if (verbose) cmessage('\t\tMake composite colors')

    # Satisfy CHECK
    subgroup <- V1 <- V2 <- color <- hue <- luminance <- NULL

    # Split into components
    #    * V1: first n-1 components => will be mapped to hue
    #    * V2: last component       => will be mapped to luminance
    # This approach works also when more than two components are present
    # It is therefore used instead of split_values()
    V1  <-  stri_split_fixed(varlevels, sep) %>%
            vapply( function(x) paste0(x[-length(x)], collapse = sep),
                    character(1))
    V2  <-  stri_split_fixed(varlevels, sep) %>%
            vapply(function(x) x[length(x)], character(1))
    V1levels <- sort(unique(V1))
    V2levels <- sort(unique(V2))
    n1 <- length(V1levels)
    n2 <- length(V2levels)
    hues <- seq(15, 375, length = n1 + 1)[seq_len(n1)] %>% set_names(V1levels)

    colors <- character(0)
    for (i in seq_along(hues)){
        colors  %<>%  c(sequential_hcl(
                                n2, h = hues[[i]], power = 1, c = c(50, 100),
                                l = c(90, 30)) %>%
                            set_names(paste0(V1levels[[i]], sep, V2levels)))
    }
    if (show) pie(rep(1, length(colors)), names(colors),
                col = colors)

    return(colors)
}

#=============================================================================
#
#     plot_data()
#
#==============================================================================

#' Plot data
#' @param data        data.frame'
#' @param geom        geom_point, etc.
#' @param color       variable mapped to color (symbol)
#' @param fill        variable mapped to fill (symbol)
#' @param ...         mapped aesthetics
#' @param fixed       fixed  aesthetics (list)
#' @param theme       list with ggplot theme specifications
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% pca()
#' data <- sdata(object)
#' plot_data(data, x = pca1, y = pca2)
#' plot_data(data, x = pca1, y = pca2, color = TIME_POINT)
#' data$TIME <- as.numeric(substr(data$TIME_POINT, 2, 3))
#' plot_data(data, x = pca1, y = pca2, color = TIME)
#' plot_data(data, x = pca1, y = pca2, color = NULL)
#'
#' fixed <- list(shape = 15, size = 3)
#' plot_data(data, x = pca1, y = pca2, fixed=fixed)
#' @author Aditya Bhagwat, Johannes Graumann
#' @export
plot_data <- function(
    data, geom = geom_point, color = subgroup, fill = !!enquo(color), ...,
    fixed = list(), theme = list()
){
    color <- enquo(color)
    fill  <- enquo(fill)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))

    p <- ggplot(data    = data,  # https://stackoverflow.com/a/55816211
                mapping = eval(expr(aes(color=!!color, fill=!!fill, !!!dots))))
    p <- p + do.call(geom, fixed)
    p <- p + theme_bw()
    p <- add_color_scale(p, !!color, data)
    p <- add_fill_scale( p, !!fill, data)
    p <- p + do.call(ggplot2::theme, theme)

    p
}


#==============================================================================
#
#                  add_highlights
#
#==============================================================================
add_highlights <- function(p, hl, geom = geom_point, fixed_color = "black") {
    hl <- enquo(hl)
    if (quo_is_null(hl)) return(p)
    hlstr <- as_name(hl)
    hl_df <- p$data[get(hlstr)==TRUE]
    args <- list(data = hl_df)
    if (identical(geom, geom_point)) {
        many_hl <- length(unique(args$data$feature_name)) > 6
        if (many_hl) args$data$feature_name <- hlstr
        args %<>% c(
            list(aes(shape = feature_name), size = rel(3), color = fixed_color))
    }
    p <- p + do.call(geom, args)
    if (identical(geom, geom_point)) p <- p +
        labs(shape = if (many_hl) NULL else hlstr) +
        guides(fill = guide_legend(override.aes = list(shape = NA)))
    p
}


#=============================================================================
#
#                plot_densities
#                    plot_sample_densities()
#                    plot_feature_densities
#
#=============================================================================


#' Plot sample/feature densities
#'
#' @param object      SummarizedExperiment
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param group       svar mapped to group
#' @param fixed       fixed aesthetics
#' @param subsetter   subsetter for showing a subset of samples/features
#' @seealso \code{\link{plot_sample_violins}},
#'          \code{\link{plot_sample_boxplots}}
#' @return  ggplot object
#' @examples
#' # Read data
#'     require(magrittr)
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#'     object %<>% extract(, order(.$subgroup))
#' # Plot distributions
#'     plot_sample_densities(object)
#'     plot_feature_densities(object)
#' @export
plot_densities <- function(object, group, fill, color = NULL,
    fixed = list(alpha = 0.5, na.rm = TRUE)
){
# Process
    assert_is_all_of(object, 'SummarizedExperiment')
    fill  <- enquo(fill)
    color <- enquo(color)
    group <- enquo(group)
    value <- NULL
    fillstr  <- if (quo_is_null(fill))  character(0) else as_name(fill)
    colorstr <- if (quo_is_null(color)) character(0) else as_name(color)
    groupstr <- if (quo_is_null(group)) character(0) else as_name(group)
# Prepare
    plotvars <- unique(c(fillstr, colorstr, groupstr))
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    assert_is_identical_to_true(is_uniquely_empty(plottedsvars, plottedfvars))
    dt <- sumexp_to_long_dt(object, svars = plottedsvars, fvars = plottedfvars)
# Plot
    plot_data(dt, geom = geom_density, x = value, fill = !!fill,
        color = !!color, group = !!group, fixed = fixed)
}

is_uniquely_empty <- function(x, y){
    is_empty <- assertive::is_empty
    ( is_empty(x) | !is_empty(y)) | (!is_empty(x) |  is_empty(y))
}

#' @rdname plot_densities
#' @export
plot_sample_densities <- function(
    object,
    fill = subgroup,
    color = NULL,
    group = sample_id,
    fixed = list(alpha=0.5, na.rm=TRUE),
    subsetter = if (ncol(object)<100){  seq_len(ncol(object))
                } else {                sample(ncol(object), 9)}
){
    plot_densities( object[, subsetter],
                    fill  = !!enquo(fill),
                    color = !!enquo(color),
                    group = !!enquo(group),
                    fixed = fixed ) +
    ggtitle("samples")
}


#' @rdname plot_densities
#' @export
plot_feature_densities <- function(
    object,
    fill = feature_id,
    color = NULL,
    group = feature_id,
    fixed = list(alpha=0.5, na.rm=TRUE),
    subsetter = if (nrow(object)<100){  seq_len(nrow(object))
                } else {                sample(nrow(object), 9)}
){
    plot_densities( object[subsetter, ],
                    fill  = !!enquo(fill),
                    color = !!enquo(color),
                    group = !!enquo(group),
                    fixed = fixed) +
    ggtitle("features")
}

#==============================================================================
#
#               plot_violins()
#                   plot_sample_violins()
#                   plot_feature_violins
#
#==============================================================================

#' Plot sample/feature violins
#'
#' @param object      SummarizedExperiment
#' @param x           svar mapped to x
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param group       svar mapped to group
#' @param facet       svar mapped to facets
#' @param highlight   fvar expressing which feature should be highlighted
#' @param fixed       fixed aesthetics
#' @return  ggplot object
#' @seealso \code{\link{plot_sample_densities}},
#'          \code{\link{plot_sample_boxplots}}
#' @examples
#' # data
#'     require(magrittr)
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#'     object %<>% extract(, order(.$subgroup))
#'     control_features <- c('biotin','phosphate')
#'     fdata(object) %<>% cbind(control=.$feature_name %in% control_features)
#' # plot
#'     plot_violins(object[1:12, ], x=feature_id, fill=feature_id)
#'     plot_feature_violins(object[1:12, ])
#'     plot_sample_violins(object[, 1:12],  highlight = control)
#'     plot_subgroup_violins(object[1:4, ])
#' @export
plot_violins <- function(object, x, fill, color = NULL, group = NULL,
    facet = NULL, highlight = NULL, fixed = list(na.rm=TRUE)
){
# Process
    assert_is_all_of(object, 'SummarizedExperiment')
    x         <- enquo(x)
    fill      <- enquo(fill)
    color     <- enquo(color)
    group     <- enquo(group)
    highlight <- enquo(highlight)
    facet     <- enquo(facet)
    sample_id <- value <- NULL
    xstr         <- as_name(x)
    fillstr      <- if (quo_is_null(fill))      character(0) else as_name(fill)
    colorstr     <- if (quo_is_null(color))     character(0) else as_name(color)
    highlightstr <- if (quo_is_null(highlight)) character(0) else as_name(
                                                                    highlight)
# Prepare
    plotvars <- unique(c('feature_name', fillstr, colorstr, highlightstr))
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    dt <- sumexp_to_long_dt(object, svars = plottedsvars, fvars = plottedfvars)
# Plot
    p <- plot_data(dt, geom = geom_violin, x = !!x, y = value,
                fill = !!fill, color= !!color, group=!!group, fixed = fixed)
    p %<>% add_highlights(!!highlight, geom = geom_point)
    p <- p + facet_wrap(vars(!!facet), scales = "free_y")
    # Finish
    breaks <- unique(dt[[xstr]])
    if (length(breaks)>50) breaks <- dt[, .SD[1], by = fillstr][[xstr]]
    p <- p + xlab(NULL) + scale_x_discrete(breaks = breaks) +
        theme(axis.text.x = element_text(angle=90, hjust=1))
# Return
    p
}

#' @rdname plot_violins
#' @export
plot_sample_violins <- function(
    object, x = sample_id, fill = subgroup, color = NULL, highlight = NULL,
    fixed=list(na.rm=TRUE)
) plot_violins(
    object, x=!!enquo(x), fill=!!enquo(fill), color=!!enquo(color),
    highlight=!!enquo(highlight), fixed = fixed
)


#' @rdname plot_violins
#' @export
plot_feature_violins <- function(
    object, x = feature_id, fill = feature_name, color = NULL, highlight = NULL,
    fixed = list(na.rm=TRUE)
) plot_violins(
    object, x=!!enquo(x), fill=!!enquo(fill), color=!!enquo(color),
    highlight=!!enquo(highlight), fixed = fixed
)

#' @rdname plot_violins
#' @export
plot_subgroup_violins <- function(
    object, x = subgroup, fill = subgroup, color = NULL, highlight = NULL,
    facet = feature_id, fixed = list(na.rm=TRUE)
) plot_violins(
    object, x = !!enquo(x), fill=!!enquo(fill), color=!!color,
    facet=!!enquo(facet), highlight=!!enquo(highlight), fixed=fixed
)


#==============================================================================
#
#                   plot_boxplots()
#                       plot_sample_boxplots()
#                       plot_feature_boxplots
#
#==============================================================================


#' Plot sample/feature boxplots
#'
#' @param object      SummarizedExperiment
#' @param x           svar mapped to x
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param facet       svar mapped to facet
#' @param highlight   fvar expressing which feature should be highlighted
#' @param fixed       fixed aesthetics
#' @return  ggplot object
#' @seealso \code{\link{plot_sample_densities}},
#'          \code{\link{plot_sample_violins}}
#' @examples
#' # data
#'     require(magrittr)
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#'     object %<>% extract(, order(.$subgroup))
#'     fdata(object) %<>% cbind(
#'                         control=.$feature_name %in% c('biotin','phosphate'))
#' # plot
#'     plot_boxplots(object[1:9,], x = feature_id, fill = feature_id)
#'     plot_boxplots(object[,1:9], x = sample_id,  fill = sample_id )
#'     plot_feature_boxplots(object[1:9, ])
#'     plot_sample_boxplots(object[, 1:12])
#'     plot_sample_boxplots(object[, 1:12], highlight = control)
#'     plot_subgroup_boxplots(object[1:2, ])
#' @export
plot_boxplots <- function(object, x, fill, color = NULL, facet = NULL,
    highlight = NULL, fixed = list(na.rm=TRUE)
){
# Assert/Process
    assert_is_all_of(object, "SummarizedExperiment")
    sample_id <- value <- NULL
    x         <- enquo(x)
    fill      <- enquo(fill)
    color     <- enquo(color)
    highlight <- enquo(highlight)
    facet     <- enquo(facet)
    xstr         <- as_name(x)
    fillstr      <- if (quo_is_null(fill))   character(0) else as_name(fill)
    colorstr     <- if (quo_is_null(color))  character(0) else as_name(color)
    highlightstr <- if (quo_is_null(highlight)) character(0) else as_name(
                                                                    highlight)
# Prepare
    plotvars <- unique(c('feature_name', xstr, fillstr, colorstr, highlightstr))
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    dt <- sumexp_to_long_dt(object, svars = plottedsvars, fvars = plottedfvars)
# Plot
    p <- plot_data(dt, geom = geom_boxplot, x = !!x, y = value,
                fill = !!fill, color = !!color, fixed = fixed)
    p %<>% add_highlights(!!highlight, geom = geom_point, fixed_color="darkred")
    p <- p + facet_wrap(vars(!!facet), scales = 'free_y')
# Finish
    breaks <- unique(dt[[xstr]])
    if (length(breaks)>50) breaks <- dt[, .SD[1], by = fillstr][[xstr]]
    p <- p + xlab(NULL) + scale_x_discrete(breaks = breaks) +
        theme(axis.text.x = element_text(angle=90, hjust=1))
# Return
    p
}


#' @rdname plot_boxplots
#' @export
plot_sample_boxplots <- function(
    object, x = sample_id, fill = subgroup, color = NULL, highlight = NULL,
    fixed = list(na.rm=TRUE)
) plot_boxplots(
    object, x = !!enquo(x), fill=!!enquo(fill), color=!!color,
    highlight=!!enquo(highlight), fixed=fixed
)


#' @rdname plot_boxplots
#' @export
plot_feature_boxplots <- function(
    object, x = feature_id, fill = feature_id, color = NULL, highlight = NULL,
    fixed = list(na.rm=TRUE)
) plot_boxplots(
    object, x = !!enquo(x), fill=!!enquo(fill), color=!!color,
    highlight=!!enquo(highlight), fixed=fixed
)

#' @rdname plot_boxplots
#' @export
plot_subgroup_boxplots <- function(
    object, x = subgroup, fill = subgroup, color = NULL, highlight = NULL,
    facet = feature_id, fixed = list(na.rm=TRUE)
) plot_boxplots(
    object, x = !!enquo(x), fill=!!enquo(fill), color=!!color,
    facet=!!enquo(facet), highlight=!!enquo(highlight), fixed=fixed
)


#' Plot samples
#'
#' Plots sample densities and scores
#' @param object    SummarizedExperiment
#' @param x         svar mapped to biplot x (sym, default pca1)
#' @param y         svar mapped to biplot y *sym, default pca2)
#' @param color     svar mapped to biplot color and density fill
#' @param nloadings n loadings added to plot
#' @param plot      whether to print plot
#' @return gtable returned by \code{\link[gridExtra]{arrangeGrob}}
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, pca=TRUE, fit='limma', plot = FALSE)
#' plot_samples(object)
#'
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, pca=TRUE, fit='limma', plot = FALSE)
#' plot_samples(object)
#' @export
plot_samples <- function(
    object, x=pca1, y=pca2, color=subgroup, nloadings=0, plot=TRUE
){
    pcaplot <-  biplot( object, x=pca1, y=pca2, color=!!enquo(color),
                        nloadings = nloadings) +
                theme(legend.position='right', legend.title=element_blank())
    samples    <- plot_sample_densities(object)  +
                    theme(legend.position='none') + xlab(NULL) + ylab(NULL)
    features   <- plot_feature_densities(object) +
                    theme(legend.position='none') + xlab(NULL) + ylab(NULL)
    detections <- plot_summarized_detections(object) +
                    theme(legend.position='none')
    volcanoes  <- plot_volcano(object, ntop = 0)
    # Cairo::CairoPNG(
    #    'contrastogram.png', width=480*7, height=480*7, dpi=400, pointsize=10)
    # plot_contrastogram(object)
    # dev.off()
    # contrastogram <- grid::rasterGrob(
    #    png::readPNG('contrastogram.png'), interpolate = TRUE)
    colors  <- gglegend(pcaplot)
    pcaplot <- pcaplot + theme(legend.position='none',
                axis.text.x = element_blank(), axis.text.y = element_blank())

    samples <- samples + theme(legend.position='none',
                                axis.text.y = element_blank())
    #features <- features   +
    #            theme(legend.position='none', axis.text.y = element_blank())
    detections <- detections + theme(legend.position='none')
    volcanoes  <- volcanoes  + theme(legend.position='none')
    pp <- grid.arrange( colors, #features, ncol=1),
                        arrangeGrob(detections, samples, pcaplot, ncol=1),
                        volcanoes,
                        ncol=3,
                        widths = c(2,4,6))
    if (plot) grid.draw(pp)
    pp
}




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
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, pca=TRUE, plot = FALSE)
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

# @rdname plot_features
# @export
#plot_feature_boxplots <- function(...){
#    plot_features(geom = geom_boxplot, color = NULL, ...)
#}


#' @rdname plot_features
#' @export
plot_feature_profiles <- function(...){
    plot_features(geom = geom_point, ...)
}

