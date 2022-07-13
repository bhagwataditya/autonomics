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
#'     add_color_scale(p, color=subgroup, data=sdata(object))
#' # STEMCELL INTENSITIES
#'    file <- download_data('billing16.proteingroups.txt')
#'    object <- read_proteingroups(
#'               file, quantity = 'Intensity labeled', plot=FALSE)
#'     p <- plot_sample_densities(object)
#'    add_color_scale(p, color=subgroup, data = sdata(object))
#' @noRd
add_color_scale <- function(p, color, data, palette = NULL){
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
            if (is.null(palette))  palette <- make_colors(levels0, sep = guess_sep(levels0))
            p <- p + scale_color_manual(values = palette, na.value = 'gray80')
        }
    }
# Return
    return(p)
}


add_fill_scale <- function(p, fill, data, palette=NULL){
# Assert
    assert_is_data.frame(data)
    fill <- enquo(fill)
# Colors
    if (!rlang::quo_is_null(fill)){
        fillstr <- as_name(fill)
        assert_is_subset(fillstr, names(data))
        values0 <- data[[fillstr]]
        if (!is.numeric(values0)){
            levels0 <- as.character(unique(values0))
            if (is.null(palette)) palette <- make_colors(levels0, sep = guess_sep(levels0))
            p <- p + scale_fill_manual(values = palette, na.value = 'gray80')
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
            n1 <- length(unique(split_extract_fixed(varlevels, sep, 1)))
            n2 <- length(unique(split_extract_fixed(varlevels, sep, 2)))
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
    if (verbose)  message('\t\tMake default ggplot colors')
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
#' varlevels <- slevels(object, 'Group')
#' make_twofactor_colors(varlevels, show = TRUE)
#' @noRd
make_twofactor_colors <- function(
    varlevels, sep  = guess_sep(varlevels), show = FALSE, verbose = TRUE
){
    # Assert
    assert_has_no_duplicates(varlevels)
    assert_is_not_null(sep)
    if (verbose)  message('\t\tMake composite colors')

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
#' @param palette     color palette (named character vector)
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
    data, geom = geom_point, color = NULL, fill = NULL, 
    ..., palette = NULL, fixed = list(), theme = list()
){
    color <- enquo(color)
    fill  <- enquo(fill)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))

    p <- ggplot(data    = data,  # https://stackoverflow.com/a/55816211
                mapping = eval(expr(aes(color=!!color, fill=!!fill, !!!dots))))
    p <- p + do.call(geom, fixed)
    p <- p + theme_bw()
    p <- add_color_scale(p, !!color, data, palette=palette)
    p <- add_fill_scale( p, !!fill,  data, palette=palette)
    p <- p + do.call(ggplot2::theme, theme)

    p
}


#==============================================================================
#
#                  add_highlights
#
#==============================================================================
add_highlights <- function(p, x, hl, geom = geom_point, fixed_color = "black") {
    feature_name <- value <- NULL
    x <- enquo(x)
    hl <- enquo(hl)
    if (quo_is_null(hl)) return(p)
    hlstr <- as_name(hl)
    hl_df <- p$data[get(hlstr)==TRUE]
    args <- list(data = hl_df)
    if (identical(geom, geom_point)) {
        many_hl <- length(unique(args$data$feature_name)) > 6
        if (many_hl) args$data$feature_name <- hlstr
        args %<>% c(
            list(aes(shape = feature_name, x = !!x, y = value), size = rel(3), color = fixed_color))
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
#' @param group       svar mapped to group
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param facet       svar mapped to facet
#' @param nrow        number of facet rows
#' @param ncol        number of facet cols
#' @param dir         'h' (horizontal) or 'v' (vertical)
#' @param labeller    e.g. ggplot2::label_value
#' @param palette     named character vector
#' @param fixed       fixed aesthetics
#' @seealso \code{\link{plot_sample_violins}},
#'          \code{\link{plot_sample_boxplots}}
#' @return  ggplot object
#' @examples
#' # Read data
#'     require(magrittr)
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#'     object %<>% extract(, order(.$Group))
#' # Plot distributions
#'     plot_sample_densities(object, fill = Group)
#'     plot_feature_densities(object)
#' @export
plot_densities <- function(
    object, group, fill, color = NULL, 
    facet = NULL, nrow = NULL, ncol = NULL, dir = 'h', labeller = label_value,
    palette = NULL, fixed = list(alpha = 0.5, na.rm = TRUE)
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
    p <- plot_data(dt, geom = geom_density, x = value, fill = !!fill,
            color = !!color, group = !!group, palette = palette, fixed = fixed)
    p <- p + facet_wrap(
            facet, nrow = nrow, ncol = ncol, dir = dir, labeller = labeller, 
            scales = "free_y")
}

is_uniquely_empty <- function(x, y){
    is_empty <- assertive::is_empty
    ( is_empty(x) | !is_empty(y)) | (!is_empty(x) |  is_empty(y))
}

#' @rdname plot_densities
#' @export
plot_sample_densities <- function(
    object,
    group   = sample_id,
    fill    = sample_id,
    color   = NULL,
    facet   = NULL, nrow = NULL, ncol = NULL, dir = 'h', labeller = label_value,
    palette = NULL,
    fixed   = list(alpha=0.5, na.rm=TRUE)
){
    if (ncol(object)>100)  object %<>% 
        extract(, round(seq(1, ncol(object), length.out = 9)))  # prevent crash
    plot_densities(
        object,
        group   = !!enquo(group),
        fill    = !!enquo(fill),
        color   = !!enquo(color),
        facet   = facet, nrow = nrow, ncol = ncol, dir = dir, labeller = labeller,
        palette = palette,
        fixed   = fixed ) +
    ggtitle("samples")
}

feature_id <- NULL
#' @rdname plot_densities
#' @export
plot_feature_densities <- function(
    object,
    group   = feature_id,
    fill    = feature_id,
    color   = NULL,
    facet   =  NULL, nrow = NULL, ncol = NULL, dir = 'h', labeller = label_value,
    palette = NULL,
    fixed   = list(alpha=0.5, na.rm=TRUE)
){
    if (nrow(object)>100)  object %<>% 
        extract(round(seq(1, ncol(object), length.out = 9)), )   # prevent crash
    plot_densities( 
        object,
        group   = !!enquo(group),
        fill    = !!enquo(fill),
        color   = !!enquo(color),
        facet   = facet,  nrow = nrow, ncol = ncol, dir = dir, labeller = labeller,
        palette = palette,
        fixed   = fixed) +
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
#' @param subgroup    subgroup svar
#' @param x           svar mapped to x
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param group       svar mapped to group
#' @param facet       svar mapped to facets
#' @param highlight   fvar expressing which feature should be highlighted
#' @param palette     named character vector with colors
#' @param fixed       fixed aesthetics
#' @return  ggplot object
#' @seealso \code{\link{plot_sample_densities}},
#'          \code{\link{plot_sample_boxplots}}
#' @examples
#' # data
#'     require(magrittr)
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#'     object %<>% extract(, order(.$Group))
#'     control_features <- c('biotin','phosphate')
#'     fdata(object) %<>% cbind(control=.$feature_name %in% control_features)
#' # plot
#'     plot_violins(object[1:12, ], x=feature_id, fill=feature_id)
#'     plot_feature_violins(object[1:12, ])
#'     plot_sample_violins(object[, 1:12],  highlight = control)
#'     plot_subgroup_violins(object[1:4, ], subgroup = Group)
#' @export
plot_violins <- function(object, x, fill, color = NULL, group = NULL,
    facet = NULL, highlight = NULL, palette = NULL, fixed = list(na.rm=TRUE)
){
# Process
    assert_is_all_of(object, 'SummarizedExperiment')
    x         <- enquo(x)
    fill      <- enquo(fill)
    color     <- enquo(color)
    group     <- enquo(group)
    highlight <- enquo(highlight)
    sample_id <- value <- NULL
    xstr         <- as_name(x)
    fillstr      <- if (quo_is_null(fill))      character(0) else as_name(fill)
    colorstr     <- if (quo_is_null(color))     character(0) else as_name(color)
    highlightstr <- if (quo_is_null(highlight)) character(0) else as_name(
                                                                    highlight)
# Prepare
    plotvars <- unique(c('feature_name', fillstr, colorstr, highlightstr))
    if (!is.null(facet))   plotvars %<>% c(vapply(facet, as_name, character(1)))
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    dt <- sumexp_to_long_dt(object, svars = plottedsvars, fvars = plottedfvars)
# Plot
    p <- plot_data(dt, geom = geom_violin, x = !!x, y = value,
                fill = !!fill, color= !!color, group=!!group, 
                palette = palette, fixed = fixed)
    p %<>% add_highlights(x=!!x, hl=!!highlight, geom = geom_point)
    p <- p + facet_wrap(facet, scales = "free")
    # Finish
    breaks <- unique(dt[[xstr]])
    if (length(breaks)>50) breaks <- dt[, .SD[1], by = fillstr][[xstr]]
    p <- p + xlab(NULL) + scale_x_discrete(breaks = breaks) +
        theme(axis.text.x = element_text(angle=90, hjust=1))
# Return
    p
}

sample_id <- NULL
#' @rdname plot_violins
#' @export
plot_sample_violins <- function(
    object, x = sample_id, fill = sample_id, color = NULL, facet = NULL, 
    highlight = NULL, fixed=list(na.rm=TRUE)
) plot_violins(
    object, x = !!enquo(x), fill = !!enquo(fill), color = !!enquo(color),
    facet = facet, highlight = !!enquo(highlight), fixed = fixed
)

feature_id <- feature_name <- NULL
#' @rdname plot_violins
#' @export
plot_feature_violins <- function(
    object, x = feature_id, fill = feature_name, color = NULL, facet = NULL, 
    highlight = NULL, fixed = list(na.rm=TRUE)
) plot_violins(
    object, x = !!enquo(x), fill = !!enquo(fill), color = !!enquo(color),
    facet = facet, highlight = !!enquo(highlight), fixed = fixed
)

subgroup <- NULL
#' @rdname plot_violins
#' @export
plot_subgroup_violins <- function(
    object, subgroup, x = !!enquo(subgroup), fill = !!enquo(subgroup), 
    color = NULL, highlight = NULL, facet = vars(feature_id), 
    fixed = list(na.rm=TRUE)
) plot_violins(
    object, x = !!enquo(x), fill = !!enquo(fill), color = !!color,
    facet = facet, highlight = !!enquo(highlight), fixed = fixed
)


#==============================================================================
#
#                   plot_boxplots()
#
#==============================================================================


#' Plot boxplots
#'
#' @param object     SummarizedExperiment
#' @param assay      string
#' @param subgroup   subgroup svar symbol
#' @param x          svar mapped to x
#' @param fill       svar mapped to fill
#' @param color      svar mapped to color
#' @param block      svar used to connect points
#' @param facet      svar mapped to facet
#' @param scales     'free', 'fixed', 'free_x', 'free_y'
#' @param nrow       number of facet rows
#' @param ncol       number of facet columns 
#' @param page       number of facet pages: \code{\link[ggforce]{facet_wrap_paginate}}
#' @param highlight  fvar expressing which feature should be highlighted
#' @param jitter     whether to add jittered data points
#' @param hlevels    xlevels for which to plot horizontal lines
#' @return  ggplot object
#' @seealso \code{\link{plot_sample_densities}},
#'          \code{\link{plot_sample_violins}}
#' @examples
#' # Read
#'     require(magrittr)
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#'     controlfeatures <- c('biotin','phosphate')
#'     fdata(object) %<>% cbind(control=.$feature_name %in% controlfeatures)
#' # Plot
#'     vars <- ggplot2::vars
#'     plot_boxplots(object[,1:9],  x = sample_id,  fill = sample_id )
#'     plot_boxplots(object[1:9,],  x = feature_id, fill = feature_id)
#'     plot_boxplots(object[1:9, ], x = SET, fill = SET, facet = vars(feature_id))
#'     plot_boxplots(object[1, ], x = SET, fill = SET, jitter=TRUE)
#'     plot_boxplots(object[1, ], x = SET, fill = SET, jitter=TRUE, block = SUB)
#'     plot_feature_boxplots(object[1:9, ])
#'     plot_sample_boxplots(object[, 1:12])
#'     plot_sample_boxplots(object[, 1:12], highlight = control)
#'     plot_subgroup_boxplots(object[1:2, ], subgroup = SET)
#'     plot_subgroup_boxplots(object[1:2, ], subgroup = SET, block = SUB)
#' @export
plot_boxplots <- function(
    object, assay = assayNames(object)[1], x, fill, color = NULL, block = NULL,
    facet = NULL, scales = 'free_y', nrow = NULL, ncol = NULL, page = 1, 
    labeller = 'label_value', highlight = NULL, 
    jitter = FALSE, palette = NULL, hlevels = NULL, ...
){
# Assert/Process
    assert_is_all_of(object, "SummarizedExperiment")
    if (nrow(object)==0)  return(ggplot2::ggplot())
    if (!is.null(facet)) assert_is_all_of(facet, 'quosures') # vars(feature_id)
    sample_id <- value <- medianvalue <- present <- direction <- NULL
    x         <- enquo(x)
    fill      <- enquo(fill)
    color     <- enquo(color)
    block     <- enquo(block)
    highlight <- enquo(highlight)
    xstr      <- as_name(x)
    fillstr  <- if (quo_is_null(fill))   character(0) else  as_name(fill)
    colorstr <- if (quo_is_null(color))  character(0) else  as_name(color)
    blockstr <- if (quo_is_null(block))  character(0) else  as_name(block)
    facetstr <- if (is.null(facet))      character(0) else  vapply(
                                                facet, as_name, character(1))
    highlightstr <- if (quo_is_null(highlight)) character(0) else as_name(
                                                highlight)
# Prepare
    plotvars <- unique(c('feature_name', xstr, fillstr, colorstr, blockstr, 
                        highlightstr, facetstr))
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    dt <- sumexp_to_long_dt(
            object, assay = assay, svars = plottedsvars, fvars = plottedfvars)
    dt[, medianvalue := median(value, na.rm=TRUE), by = c('feature_id', xstr)]
    
# Plot
    p <- ggplot(dt) + theme_bw() 
    if (nrow(dt)==0) return(p)
    p <- p + facet_wrap_paginate(facets = facet, scales = scales, 
                                 nrow = nrow, ncol = ncol, page = page, labeller = labeller)
    if (!quo_is_null(block)){
        dt[, direction := get(fillstr)[which.max(value)], by = c(facetstr, blockstr)]
        p <- p + geom_line(aes(x=!!x, y=value, color=direction, group=!!block), na.rm = TRUE)
        p <- add_color_scale(p, direction, data=dt, palette=palette) }
    outlier.shape <- if (jitter) NA else 19
    p <- p + geom_boxplot(aes(x = !!x, y = value, fill = !!fill, color = !!color), outlier.shape = outlier.shape, na.rm = TRUE)
    p <- add_color_scale(p, !!color, data=dt, palette=palette)
    p <- add_fill_scale( p, !!fill,  data=dt, palette=palette)
    p %<>% add_highlights(x=!!x, hl=!!highlight, geom = geom_point, fixed_color="darkred")
    # Add hline
    if (!is.null(hlevels)){
        mediandt <- unique(dt[, 
                              unique(c('feature_id', xstr, 'medianvalue', facetstr)), with=FALSE])
        mediandt[, present := FALSE]
        mediandt[get(xstr) %in% hlevels, present := TRUE]
        p <- p + geom_hline(data=mediandt, 
                            aes(yintercept = medianvalue, color=!!fill, alpha=present), 
                            linetype='longdash') }
    # Add jitter
    if (jitter) p <- p + geom_jitter(aes(x=!!x, y=value),
                                     position = position_jitter(width=.1, height=0), size=0.5, na.rm = TRUE)
    # Finish and Return
    breaks <- unique(dt[[xstr]])
    if (length(breaks)>50) breaks <- dt[, .SD[1], by = fillstr][[xstr]]
    p <- p + xlab(NULL) + scale_x_discrete(breaks = breaks) + 
        guides(color = 'none', alpha = 'none') +
        theme(axis.text.x = element_text(angle=90, hjust=1))
    p
}

#============================================================================
#
#                       plot_sample_boxplots
#                       plot_feature_boxplots
#                       plot_subgroup_boxplots
#                       plot_contrast_boxplots
#
#============================================================================

#' @rdname plot_boxplots
#' @export
plot_sample_boxplots <- function(
    object, assay = assayNames(object)[1], 
    x = sample_id, fill = sample_id, color = NULL, highlight = NULL,
    palette = NULL, nrow=NULL, ncol=NULL, page=1, labeller = 'label_value'
){
    x         <- enquo(x)
    fill      <- enquo(fill)
    color     <- enquo(color)
    highlight <- enquo(highlight)
    plot_boxplots(
        object, assay = assay, 
        x = !!x, fill = !!fill, color = !!color,
        highlight = !!highlight, palette = palette,  
        nrow = nrow, ncol = ncol, page = page, labeller = labeller)
}


feature_id <- NULL
#' @rdname plot_boxplots
#' @export
plot_feature_boxplots <- function(
    object, assay = assayNames(object)[1],
    x = feature_id, fill = feature_id, color = NULL, highlight = NULL,
    palette = NULL, fixed = list(na.rm=TRUE), 
    nrow = NULL, ncol = NULL, page = 1, labeller = 'label_value'
){
    x         <- enquo(x)
    fill      <- enquo(fill)
    color     <- enquo(color)
    highlight <- enquo(highlight)
    plot_boxplots(
        object, assay = assay,
        x = !!x, fill = !!fill, color = !!color,
        highlight = !!highlight, palette = palette,  
        nrow = nrow, ncol = ncol, page = page, labeller = labeller)
}


#' @rdname plot_boxplots
#' @export
plot_subgroup_boxplots <- function(
    object, assay = assayNames(object)[1],
    subgroup, x = !!enquo(subgroup), fill = !!enquo(subgroup), 
    color = NULL, block = NULL, highlight = NULL, jitter = TRUE, 
    facet = vars(feature_id), scales = 'free_y', nrow = NULL, ncol = NULL, 
    page = 1, labeller = 'label_value',
    palette = NULL, fixed = list(na.rm=TRUE), hlevels = NULL
){
    x         <- enquo(x)
    fill      <- enquo(fill)
    color     <- enquo(color)
    block     <- enquo(block)
    highlight <- enquo(highlight)
    p <- plot_boxplots(
            object, assay = assay,
            x = !!x, fill = !!fill, color = !!color,
            block = !!block, highlight = !!highlight, 
            facet = facet, scales = scales, nrow = nrow,
            ncol = ncol, page = page,  labeller = labeller, 
            jitter = jitter, palette = palette, 
            fixed = fixed, hlevels = hlevels)
    #p <- p + stat_summary(fun='mean', geom='line', aes(group=1), color='black', size=0.7, na.rm = TRUE)
    p
}


extract_coef_features <- function(
    object, 
    coef = setdiff(coefs(object), 'Intercept')[1],
    fit  = fits(object)[1], 
    fdrcutoff, 
    ntop
){
    pvals   <- p(     object, coef = coef, fit = fit)[, 1]
    fdrs    <- fdr(   object, coef = coef, fit = fit)[, 1]
    effects <- effect(object, coef = coef, fit = fit)[, 1]
    signs   <- sign(  object, coef = coef, fit = fit)[, 1]
    
    upfeatures <- names(    sort(pvals[signs>0 & fdrs < fdrcutoff]))
    dnfeatures <- names(rev(sort(pvals[signs<0 & fdrs < fdrcutoff])))
    upfeatures %<>% utils::head(ntop %/% 2)
    dnfeatures %<>% utils::tail(ntop-length(upfeatures))
    
    object %<>% extract(c(dnfeatures, upfeatures), )
}

format_coef_vars <- function(
    object, coef = setdiff(coefs(object), 'Intercept')[1], fit = fits(object)[1]
){
    effectvars <- effectvar(object, coef = coef, fit = fit)
    pvars      <- pvar(     object, coef = coef, fit = fit)
    fdrvars    <- fdrvar(   object, coef = coef, fit = fit)
    for (var in c(effectvars, pvars, fdrvars)){
        fdt(object)[[var]] %<>% formatC(format='e', digits=0)
    }
    object
    
}


#' Plot coef boxplots
#'
#' @param object       SummarizedExperiment
#' @param assay        string
#' @param subgroup     symbol: subgroup svar
#' @param block        symbol: block var
#' @param downlevels   character vector
#' @param uplevels     character vector
#' @param fit          'limma', 'lm', 'lme', 'lmer', 'wilcoxon'
#' @param coef     string
#' @param fdrcutoff    number
#' @param palette      color palette (named character vector)
#' @param title        string
#' @param ylab         NULL or string
#' @param nrow         number
#' @param ncol         number
#' @param labeller     string or function
#' @param scales       'free', 'free_x', 'free_y', 'fixed'
#' @examples 
#' require(magrittr)
#' require(grid)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' object %<>% fit_limma(subgroupvar='SampleGroup', block='Subject_ID')
#' object %<>% fit_lme(  subgroupvar='SampleGroup', block='Subject_ID')
#' object$SampleId   %<>% factor()
#' object$Subject_ID %<>% factor()
#' plot_coef_boxplots(
#'     object, subgroup = SampleGroup, block = Subject_ID, 
#'     fit='limma', coef = 'SampleGroupt1')
#' plot_coef_boxplots(
#'     object, subgroup=SampleGroup, block = Subject_ID, fit='limma', 
#'     coef = 'SampleGroupt1', fdrcutoff=0.05)
#' plot_coef_boxplots(
#'     object, subgroup = SampleGroup, block = Subject_ID, 
#'      fit = c('limma', 'lme'), contrast = 'SampleGroupt1', 
#'      fdrcutoff = 0.05)
#' @export
plot_coef_boxplots <- function(
    object, assay = assayNames(object)[1], 
    subgroup, block = NULL, fit = fits(object)[1], 
    coef = setdiff(coefs(object), 'Intercept')[1],
    facet = vars(feature_id, !!sym(paste('fdr', coef, fit, sep = FITSEP))),
    fdrcutoff = 0.05, palette = NULL, title = coef, ylab = NULL, 
    nrow = NULL, ncol = NULL, ntop = 4,
    labeller = 'label_value', scales = 'free_y'
){
# Order/Extract on p value
    . <- NULL
    subgroup  <- enquo(subgroup); subgroupvar <- if (quo_is_null(subgroup)) character(0) else as_name(subgroup)
    block     <- enquo(block);    blockvar    <- if (quo_is_null(block))    character(0) else as_name(block)
    facetvars <- vapply(facet, as_name, character(1))
    svars0 <- c(subgroupvar, blockvar)
    fdrvar    <- paste('fdr',    coef, fit, sep = FITSEP)
    pvar      <- paste('p',      coef, fit, sep = FITSEP)
    effectvar <- paste('effect', coef, fit, sep = FITSEP)
    fvars0 <- c(facetvars, fdrvar, pvar, effectvar)
    object %<>% extract_coef_features(coef = coef, fit = fit, fdrcutoff = fdrcutoff, ntop = ntop)
    object %<>% format_coef_vars(     coef = coef, fit = fit)
# Plot
    dt <- sumexp_to_long_dt(object, assay = assay, svars = svars0, fvars = fvars0)
    mediandt <- summarize_median(dt, subgroupvar)
    contrastsubgroup <- NULL
    #mediandt[, contrastsubgroup := get(subgroupvar) %in% c(uplevels, downlevels)]
    facetstr <- vapply(facet, as_name, character(1))
    p <- plot_subgroup_boxplots(
            object, subgroup = !!subgroup, block = !!block, x = !!subgroup, fill = !!subgroup, 
            facet = facet, scales = scales, nrow = nrow, 
            ncol = ncol, labeller = labeller, palette  = palette)  + 
        theme_bw() + xlab(NULL) + ggtitle(title) + ylab(ylab) + 
        # geom_hline( data = mediandt, linetype = 'longdash',
        #        aes(yintercept=!!sym('value'), alpha=contrastsubgroup, 
        #            color = !!sym(subgroupvar))) + 
        guides(alpha = 'none')
    p
# Color
    #if (is.null(palette))  palette <- make_colors(levels(object[[subgroupvar]]))
    #facetdt <- unique(object[, c(facetstr, effectvar, subgroupvar), with=FALSE])
    #facetdt <- rbind(
    #    facetdt[ as.numeric(get(fdrvar   )) < fdrcutoff & as.numeric(get(effectvar)) > 0, .SD[get(subgroupvar) %in% uplevels  ]],
    #    facetdt[ as.numeric(get(fdrvar   )) < fdrcutoff & as.numeric(get(effectvar)) < 0, .SD[get(subgroupvar) %in% downlevels]])
    #facetpalette <- facetdt[, palette[get(subgroupvar)] %>% set_names(feature_id)]
    #g <- color_facets(p, facetpalette)
# Return
    #grid::grid.draw(g)
    #invisible(g)
}

summarize_median <- function(dt, subgroupvar){
    value <- NULL
    mediandt <- copy(dt)
    mediandt[,value:=median(value, na.rm=TRUE), by=c('feature_id', subgroupvar)]
    mediandt[,sample_id := NULL]
    unique(mediandt)
}

color_facets <- function(p, palette){
    . <- NULL
    g <- ggplot_gtable(ggplot_build(p))
    facets <- which(grepl('strip-', g$layout$name))
    for (facet in facets){
        if (!is(g$grobs[[facet]], 'zeroGrob')){
            for (strip in seq_along(g$grobs[[facet]]$grobs)){
                content <- g$grobs[[facet]]$grobs[[strip]]$children[[2]]$children[[1]]$label
                if (grepl('feature_id', content)){
                    content %<>% gsub('feature_id: ', '', ., fixed=TRUE) 
                    g$grobs[[facet]]$grobs[[strip]]$children[[2]]$children[[1]]$label <- content
                    feature <- content
                }
                if (grepl('fdr', content)){
                    content %<>% gsub('^[^:]+[:] ', '', .)
                    content %<>% as.numeric()
                    if (content < 0.05){
                        g$grobs[[facet]]$grobs[[strip]]$children[[1]]$gp$fill <- palette[[feature]] } } } } }
    g
}

which.names <- function(x) names(x[x])

expand_into_vector <- function(value, templatevector){
    y <- rep(value, length(templatevector))
    names(y) <- templatevector
    y
}



#===========================================================================
#
#                plot_svar_boxplots
#
#===========================================================================


.prep_subgroup_contrasts <- function(
    object,
    svar,
    contrast,
    xvar,
    fillvar,
    blockvar,
    fit = fits(object) [1],
    nfeature = 1
){
    contr <- contrast
    sumdt <- summarize_fit(object, fit = fit)[contrast==contr]
    pvar   <- paste('p',   contrast, fit, sep = FITSEP)
    fdrvar <- paste('fdr', contrast, fit, sep = FITSEP)
    topfid <- fnames(object)[order(fdata(object)[[pvar]])[nfeature]]
    topfname <- as.character(fdata(object)[topfid, ]$feature_name)
    subject <- object[topfid, ]
    subject$xinteraction <- interaction(subject[[xvar]], subject[[svar]])
    plotdt <- sumexp_to_long_dt(subject,
                                svars = c('xinteraction', fillvar, blockvar),
                                fvars = c('feature_name', pvar, fdrvar))
    setnames(plotdt, pvar,   'p')
    setnames(plotdt, fdrvar, 'fdr')
    plotdt$facet <- contrast
    plotdt
}


prep_subgroup_contrasts <- function(
    object, svar, contrasts, xvar, fillvar, blockvar, nfeature
){
    contrastlevels <- summarize_fit(object)[rev(order(ndown+nup))]$contrast
    contrastlevels %<>% as.character()
    . <- NULL
    plotdt <- mapply(
        .prep_subgroup_contrasts,
        svar = svar,
        contrast = contrasts,
        MoreArgs = list(object   = object,
                        xvar     = xvar,
                        fillvar  = fillvar,
                        blockvar = blockvar, 
                        nfeature = nfeature),
        SIMPLIFY = FALSE)
    plotdt %<>% data.table::rbindlist()
    plotdt %<>% data.table::setorder(p)
    plotdt[, fdr := paste0('FDR ', formatC(fdr, format='e', digits=2))]
    plotdt$facet %<>% factor(contrastlevels)
    plotdt
}


plot_svar_boxplots <- function(...){
    .Deprecated('plot_subgroup_contrasts')
}



#' Plot subgroup contrasts
#' @param object    SummarizedExperiment
#' @param svar      string
#' @param contrasts character vector
#' @param xvar      string
#' @param fillvar   string
#' @param blockvar  string
#' @param nfeature  number of plotted features per contrast
#' @param nrow      number of plot rows
#' @param ncol      number of plot columns
#' @param palette   named character vector : color palette
#' @return ggplot 
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% fit_limma(subgroupvar = 'SET')
#' plot_subgroup_contrasts(object, subgroup = SET, block = SUB)
#' @export
plot_subgroup_contrasts <- function(
    object, 
    subgroup, 
    block     = NULL, 
    contrasts = coefs(object), #, svars = as_name(!!enquo(subgroup))),
    x         = !!enquo(subgroup), 
    fill      = !!enquo(subgroup),
    jitter    = FALSE,
    nfeature  = 1, 
    nrow      = NULL, 
    ncol      = NULL, 
    palette   = NULL
){
    facet <- palette <- NULL
    subgroup <- enquo(subgroup)
    x        <- enquo(x)
    fill     <- enquo(fill)
    block    <- enquo(block)

    subgroupstr <- as_name(subgroup)
    xstr     <- as_name(x)
    fillstr  <- as_name(fill)
    blockstr <- if (quo_is_null(block)) NULL else as_name(block)

    plotdt <- prep_subgroup_contrasts(
                object      = object,
                svar        = subgroupstr,
                contrasts   = contrasts,
                xvar        = xstr,
                fillvar     = fillstr,
                blockvar    = blockstr, 
                nfeature    = nfeature)
    plot_boxplots(
        plotdt,
        x       = !!x,
        fill    = !!fill,
        block   = !!block,
        facet   = vars(facet, feature_name, fdr),
        jitter  = jitter,
        scales  = 'free',
        nrow    = nrow,
        ncol    = ncol,
        palette = palette)

}


#=============================================================================
#
#                 plot_feature_points()
#
#=============================================================================

#' Plot features
#' @param object      SummarizedExperiment
#' @param subgroup    subgroup svar
#' @param block       block svar
#' @param x           svar mapped to x
#' @param color       svar mapped to color
#' @param group       svar mapped to group
#' @param facet       svar mapped to facets
#' @param nrow        number of rows
#' @param scales      'free_y' etc. 
#' @param ...         mapped aesthetics
#' @param palette     color palette (named character vector)
#' @param fixed       fixed aesthetics
#' @param theme       ggplot theme specifications
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit='limma', plot = FALSE)
#' idx <- order(fdata(object)$F.p.limma)[1:9]
#' object %<>% extract(idx, )
#' plot_sample_boxplots(  object)
#' plot_feature_boxplots( object)
#' plot_subgroup_boxplots(object, subgroup = SET)
#' plot_subgroup_points(  object, subgroup = SET)
#' plot_subgroup_points(  object, subgroup = SET, block = SUB)
#' @export
plot_subgroup_points <- function(
    object, subgroup, block = NULL, x = !!enquo(subgroup), 
    color = !!enquo(subgroup), group = !!enquo(block), 
    facet = vars(feature_id), nrow = NULL, scales = 'free_y', ...,
    palette = NULL,
    fixed = list(na.rm=TRUE),  #element_text(angle=90, vjust=0.5),
    theme = list(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
){
    x     <- enquo(x)
    color <- enquo(color)
    group <- enquo(group)
    
    dt <- sumexp_to_long_dt(
            object, svars = svars(object), fvars = fvars(object))
    value <- NULL
    p <- plot_data(dt, geom = geom_point, x = !!x, y = value, color = !!color, 
                    group = !!group, ..., palette = palette, fixed = fixed)
    if (!quo_is_null(enquo(block)))  p <- p + geom_line()
    p <- p + facet_wrap(facets = facet, scales = scales, nrow = nrow)
    p <- p + do.call(ggplot2::theme, theme)
    p
}


#=========================================================
#
#           plot_venn_heatmap
#           plot_venn
#           plot_contrast_venn
#               list2mat
#
#=========================================================


#' list to matrix
#' @param x list
#' @return matrix
#' @examples
#' x <- list(roundfruit = c('apple', 'orange'), redfruit = c('apple', 'strawberry'))
#' list2mat(x)
#' @export
list2mat <- function(x){
    uni <- unique(Reduce(union, x))
    mat <- matrix(0, nrow = length(uni), ncol = length(x),
           dimnames = list(uni, names(x)))
    for (i in seq_along(x))  mat[x[[i]], names(x)[i]] <- 1
    mat
}


#' Plot venn heatmap
#' @examples
#' x <- list(roundfruit = c('apple', 'orange'), redfruit = c('apple', 'strawberry'))
#' plot_venn_heatmap(x)
#' @export
plot_venn_heatmap <- function(x){
    if (!requireNamespace('pheatmap', quietly = TRUE)){
        stop("`BiocManager::install('pheatmap')`")
    }
    assert_is_list(x)
    x %<>% list2mat()
    pctmat <- matrix(0, nrow = ncol(x), ncol = ncol(x), dimnames = list(colnames(x), colnames(x)))
    nmat   <- matrix(0, nrow = ncol(x), ncol = ncol(x), dimnames = list(colnames(x), colnames(x)))
    for (cl1 in colnames(x)){
    for (cl2 in colnames(x)){
        set1 <- rownames(x)[x[, cl1]==1]
        set2 <- rownames(x)[x[, cl2]==1]
        nmat[  cl2, cl1] <- length(intersect(set1, set2))
        pctmat[cl2, cl1] <- length(intersect(set1, set2)) / min(length(set1), length(set2))
        pctmat[cl1, cl2] <- length(intersect(set1, set2)) / min(length(set2), length(set2))
    }
    }
    pheatmap::pheatmap(pctmat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = nmat)
}


#' Plot venn
#' @param x list
#' @examples
#' x <- list(roundfruit = c('apple', 'orange'), redfruit = c('apple', 'strawberry'))
#' plot_venn(x)
#' @export
plot_venn <- function(x){
    assert_is_list(x)
    limma::vennDiagram(list2mat(x))
}


#' Plot contrast venn
#' @param isfdr matrix(nrow, ncontrast): -1 (down), +1 (up)
#' @return nothing returned
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' object %<>% fit_wilcoxon(subgroupvar='SampleGroup', block = 'Subject_ID')
#' object %<>% fit_limma(   subgroupvar='SampleGroup', block = 'Subject_ID')
#' isfdr <- is_sig(object, contrast = 't3-t2')
#' plot_contrast_venn(isfdr)
#' @export
plot_contrast_venn <- function(issig, colors = NULL){
    assert_is_matrix(issig)
    layout(matrix(c(1,2), nrow=2))
    vennDiagram(issig, include='up',   mar = rep(0,4), show.include=TRUE, circle.col = colors)
    vennDiagram(issig, include='down', mar = rep(0,4), show.include=TRUE, circle.col = colors)
}

#' Plot binary matrix
#' @param mat 
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' mat <- sdt(object)[, .(replicate, Group)]
#' mat$present <- 1
#' mat %<>% data.table::dcast(replicate ~ Group, value.var  = 'present', fill = 0)
#' mat %<>% dt2mat()
#' plot_matrix(mat)
#' @return
#' @export
plot_matrix <- function(mat){
    nr <- nrow(mat)
    nc <- ncol(mat)
    values <- unique(c(mat)) %>% setdiff(0) %>% as.character()
    colors <- make_colors(values)
    colors %<>% unname()
    colors %<>% c('white', .)
    
    image(t(mat %>% extract(seq(nrow(.), 1), )), col = colors,  axes = FALSE)
    axis(side = 1, labels =     colnames(mat),  at = seq(0, by = 1, length.out = nc)/(nc-1), las = 1, tick = FALSE)
    axis(side = 2, labels = rev(rownames(mat)), at = seq(0, by = 1, length.out = nr)/(nr-1), las = 1, tick = FALSE)
    box()
    par(mar = c(5,5,4,2))
    abline(h = (0.5:(nr-0.5))/(nr-1), v = (0.5:(nc-0.5))/(nc-1), col = 'gray30')
}
