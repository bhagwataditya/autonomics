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
#'     object <- read_proteingroups(file, invert_subgroups = invert_subgroups, plot = FALSE)
#'     p <- plot_sample_densities(object)
#'     add_color_scale(p, color = 'subgroup', data = sdt(object))
#' # STEMCELL INTENSITIES
#'    file <- download_data('billing16.proteingroups.txt')
#'    object <- read_proteingroups(file, quantity = 'Intensity labeled', plot = FALSE)
#'     p <- plot_sample_densities(object)
#'    add_color_scale(p, color = 'subgroup', data = sdt(object))
#' @noRd
add_color_scale <- function(p, color, data, palette = NULL){
# Assert
    assert_is_data.frame(data)
    assert_is_subset(color, names(data))
# Colors
    if (!is.null(color)){
        values0 <- data[[color]]
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


add_fill_scale <- function(p, fill, data, palette = NULL){
# Assert
    assert_is_data.frame(data)
# Colors
    if (!is.null(fill)){
        assert_is_subset(fill, names(data))
        values0 <- data[[fill]]
        if (!is.numeric(values0)){
            levels0 <- as.character(unique(values0))
            if (is.null(palette)) palette <- make_colors(levels0, sep = guess_sep(levels0))
            p <- p + scale_fill_manual(values = palette, na.value = 'gray80')
        }
    }
# Return
    return(p)
}

make_sample_palette <- function(object){
    dt <- sdt(object)[, c('sample_id', 'subgroup'), with = FALSE]
    palette <- make_colors(levels(dt$subgroup))
    palette <- data.table(subgroup = names(palette), color = palette)
    dt %<>% merge(palette, by = 'subgroup')
    palette <- dt$color
    names(palette) <- dt$sample_id
    palette
}

make_subgroup_palette <- function(object){
    make_colors(subgroup_levels(object))
}
make_svar_palette <- function(object, svar){ 
    if (is.null(svar)) return(NULL)
    make_colors(slevels(object, svar))
}
make_fvar_palette <- function(object, fvar){
    if (is.null(fvar)) return(NULL)
    make_colors(flevels(object, fvar))
}
make_var_palette <- function(object, var){
    if (is.null(var)) return(NULL)
    if (var %in% svars(object)){        make_svar_palette(object, var)
    } else if (var %in% fvars(object)){ make_fvar_palette(object, var) }
}

#' Make colors
#' @param varlevels character vector
#' @param sep       string
#' @param show      TRUE or FALSE: whether to plot
#' @param verbose   TRUE or FALSE: whether to msg
#' @examples 
#' make_colors(c('A',   'B',   'C',  'D'  ), show = TRUE)
#' make_colors(c('A.1', 'B.1', 'A.2','B.2'), show = TRUE)
#' @export
make_colors <- function(
    varlevels, sep = guess_sep(varlevels), show = FALSE, verbose = FALSE
){
# Numeric colors
    if (is.numeric(varlevels)){
        colors <- brewer.pal(length(varlevels), 'YlOrRd')
        names(colors) <- varlevels
        if (show) pie(rep(1, length(colors)), names(colors), col = colors)
        return(colors)
    }
# # Twofactor colors
#     if (!is.null(sep)){            # consistent separator
#         if (length(varlevels)>2){  # 3+ samples
#             n1 <- length(unique(split_extract_fixed(varlevels, sep, 1)))
#             n2 <- length(unique(split_extract_fixed(varlevels, sep, 2)))
#             if (n1>1 & n2>1){             # 2+ huevar levels
#                 return(make_twofactor_colors(
#                     varlevels, sep = sep, show = show, verbose = verbose))
#             }
#         }
#     }
# Onefactor colors
    return(make_onefactor_colors(
        varlevels, sep = sep, show = show, verbose = verbose))
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
    varlevels, show = FALSE, verbose = TRUE, sep = NULL
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
    for (i in seq_along(hues)){  # https://stackoverflow.com/a/5738083
        basecolor  <- hcl(h = hues[[i]], c = 100, l = 50)
        newcolors <- colorRampPalette(c('white', basecolor))(n2+1)[-1]
        names(newcolors) <- paste0(V1levels[[i]], sep, V2levels)
        colors %<>% c(newcolors)
    }
    if (show) pie(rep(1, length(colors)), names(colors), col = colors)

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
                mapping = eval(expr(aes(color = !!color, 
                                        fill  = !!fill, 
                                        !!!dots))))
    p <- p + do.call(geom, fixed)
    p <- p + theme_bw()
    colorstr <- if (quo_is_null(color)) NULL else as_name(color)
    fillstr  <- if (quo_is_null(fill))  NULL else as_name(fill)
    p <- add_color_scale(p, colorstr, data, palette = palette)
    p <- add_fill_scale( p, fillstr,  data, palette = palette)
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
    if (is.null(hl)) return(p)
    hl_df <- p$data[get(hl)==TRUE]
    args <- list(data = hl_df)
    if (identical(geom, geom_point)) {
        many_hl <- length(unique(args$data$feature_name)) > 6
        if (many_hl) args$data$feature_name <- hl
        args %<>% c(list(aes(shape = feature_name, x = !!sym(x), y = value), size = rel(3), color = fixed_color))
    }
    p <- p + do.call(geom, args)
    if (identical(geom, geom_point)) p <- p +
        labs(shape = if (many_hl) NULL else hl) +
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


#' Plot sample/feature distributions
#'
#' @param object      SummarizedExperiment
#' @param group       svar (string)
#' @param fill        svar (string)
#' @param color       svar (string)
#' @param facet       svar (character vector)
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
#' # Data
#'     require(magrittr)
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#'     object %<>% extract(, order(.$subgroup))
#'     
#' # Sample distributions
#'     plot_sample_densities(object)
#'     plot_sample_violins(  object, facet = 'SET')
#'     plot_sample_boxplots( object, facet = 'SET')
#'     
#' # Feature distributions
#'     plot_feature_densities(object)
#'     plot_feature_violins(  object)
#'     plot_feature_boxplots( object)
#' @export
plot_densities <- function(
    object, group, fill, color = NULL, 
    facet = NULL, nrow = NULL, ncol = NULL, dir = 'h', scales = 'free_y', 
    labeller = label_value, 
    palette = NULL, fixed = list(alpha = 0.8, na.rm = TRUE)
){
# Assert / Process
    assert_is_all_of(object, 'SummarizedExperiment')
                            assert_is_a_string(group)
    if (!is.null(fill))     assert_is_a_string(fill)
    if (!is.null(color))    assert_is_a_string(color)
    if (!is.null(facet))    assert_is_a_string(facet)
    if (!is.null(nrow))     assert_is_a_number(nrow)
    if (!is.null(ncol))     assert_is_a_number(ncol)
    if (!is.null(palette))  assert_is_character(palette)
                            assert_is_list(fixed)
                            assert_is_subset(group, c(svars(object), fvars(object)))
    if (!is.null(fill))     assert_is_subset(fill,  c(svars(object), fvars(object))) 
    if (!is.null(color))    assert_is_subset(color, c(svars(object), fvars(object)))
    if (!is.null(facet))    assert_is_subset(facet, c(svars(object), fvars(object)))
                            assert_is_subset(dir, c('h', 'v'))
# Prepare
    plotvars <- group
    if (!is.null(fill))   plotvars %<>% c(fill)  %>% unique()
    if (!is.null(color))  plotvars %<>% c(color) %>% unique()
    if (!is.null(facet))  plotvars %<>% c(facet) %>% unique()
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    assert_is_identical_to_true(is_uniquely_empty(plottedsvars, plottedfvars))
    if (!is.null(fill))  object[[fill]] %<>% num2char()
    dt <- sumexp_to_longdt(object, svars = plottedsvars, fvars = plottedfvars)
# Plot
    groupsym <- if (is.null(group))  quo(NULL) else sym(group)
    fillsym  <- if (is.null(fill ))  quo(NULL) else sym(fill)
    colorsym <- if (is.null(color))  quo(NULL) else sym(color)
    p <- plot_data(dt, geom = geom_density, x = value, fill = !!fillsym,
            color = !!colorsym, group = !!groupsym, palette = palette, fixed = fixed)
    if (!is.null(facet))  p <- p + facet_wrap(
            facet, nrow = nrow, ncol = ncol, dir = dir, labeller = labeller, 
            scales = scales)
    p
}

is_uniquely_empty <- function(x, y){
    is_empty <- assertive::is_empty
    ( is_empty(x) | !is_empty(y)) | (!is_empty(x) |  is_empty(y))
}

#' @rdname plot_densities
#' @export
plot_sample_densities <- function(
    object,
    group    = 'sample_id',
    fill     = 'sample_id',
    color    = NULL, 
    n        = 100,
    facet    = NULL, nrow = NULL, ncol = NULL, dir = 'h', scales = 'free_y',
    labeller = label_value,
    palette  = NULL,
    fixed    = list(alpha = 0.8, na.rm = TRUE)
){
    object %<>% extract_samples_evenly(n)
    plot_densities(
        object,
        group    = group,
        fill     = fill,
        color    = color,
        facet    = facet, nrow = nrow, ncol = ncol, dir = dir, scales = scales,
        labeller = labeller, 
        palette  = palette, 
        fixed    = fixed ) +
    ggtitle("Sample Densities")
}

feature_id <- NULL
#' @rdname plot_densities
#' @export
plot_feature_densities <- function(
    object,
    group   = 'feature_id',
    fill    = 'feature_id',
    color   = NULL,
    n       = 9,
    facet   = NULL, nrow = NULL, ncol = NULL, dir = 'h', scales = 'free', 
    labeller = label_value, palette = NULL, 
    fixed = list(alpha = 0.8, na.rm = TRUE)
){
    object %<>% extract_features_evenly(n)
    plot_densities( 
        object,
        group    = group,
        fill     = fill,
        color    = color,
        facet    = facet,  nrow = nrow, ncol = ncol, dir = dir, scales = scales, 
        labeller = labeller, 
        palette  = palette, 
        fixed    = fixed) +
    ggtitle("Feature Densities")
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
#' @param assay       string
#' @param subgroup    subgroup svar
#' @param x           svar (string)
#' @param fill        svar (string)
#' @param color       svar (string)
#' @param group       svar (string)
#' @param facet       svar (character vector)
#' @param highlight   fvar expressing which feature should be highlighted (string)
#' @param palette     named color vector (character vector)
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
#'     fdata(object) %<>% cbind(control = .$feature_name %in% control_features)
#' # plot
#'     plot_violins(object[1:12, ], x = 'feature_id', fill = 'feature_id')
#'     plot_feature_violins(object[1:12, ])
#'     plot_sample_violins(object[, 1:12],  highlight = 'control')
#'     plot_subgroup_violins(object[1:4, ], subgroup = 'subgroup')
#' @export
plot_violins <- function(object, assay = assayNames(object)[1], x, fill, 
    color = NULL, group = NULL, facet = NULL, nrow = NULL, ncol = NULL, 
    dir = 'h', scales = "free", labeller = label_value, highlight = NULL, 
    palette = NULL, fixed = list(na.rm = TRUE)
){
# Process
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(assay, assayNames(object))
                              assert_is_a_string(x)
                              assert_is_a_string(fill)
    if (!is.null(color))      assert_is_a_string(color)
    if (!is.null(group))      assert_is_a_string(group)
    if (!is.null(facet))      assert_is_a_string(facet)
    if (!is.null(highlight))  assert_is_a_string(highlight)
                              assert_is_subset(x,         c(svars(object), fvars(object)))
                              assert_is_subset(fill,      c(svars(object), fvars(object)))
    if (!is.null(color))      assert_is_subset(color,     c(svars(object), fvars(object)))
    if (!is.null(group))      assert_is_subset(group,     c(svars(object), fvars(object)))
    if (!is.null(facet))      assert_is_subset(facet,     c(svars(object), fvars(object)))
    if (!is.null(highlight))  assert_is_subset(highlight, c(svars(object), fvars(object)))
    assert_is_list(fixed)
# Prepare
    plotvars <- c('feature_name')
                              plotvars %<>% c(x)         %>% unique()
                              plotvars %<>% c(fill)      %>% unique()
    if (!is.null(color))      plotvars %<>% c(color)     %>% unique()
    if (!is.null(highlight))  plotvars %<>% c(highlight) %>% unique()
    if (!is.null(facet))      plotvars %<>% c(facet)     %>% unique()
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    dt <- sumexp_to_longdt(object, assay = assay, svars = plottedsvars, fvars = plottedfvars)
    dtsum <- dt[, .(median = median(value, na.rm = TRUE), 
                       iqr =    IQR(value, na.rm = TRUE) ), by = x]
# Plot
    xsym         <- sym(x)
    fillsym      <- sym(fill)
    groupsym     <- if (is.null(group))  quo(NULL)  else  sym(group)
    colorsym     <- if (is.null(color))  quo(NULL)  else  sym(highlight)
    p <- plot_data(dt, geom = geom_violin, x = !!xsym, y = value,
                fill = !!fillsym, color = !!colorsym, group = !!groupsym, 
                palette = palette, fixed = fixed)
    #p <- p + geom_point(data = dtsum, aes(x = !!xsym, y = median))
    p <- p + geom_boxplot(width = 0.1, na.rm = TRUE)
    #p <- p + geom_errorbar(
    #    data    = dtsum, 
    #    mapping = aes(x = !!xsym, ymin = median-iqr, ymax = median+iqr, y = median), 
    #    width   = 0)
    p %<>% add_highlights(x = x, hl = highlight, geom = geom_point)
    if (!is.null(facet))  p <- p + facet_wrap(facet, nrow = nrow, ncol = ncol, 
                              dir = dir, scales = scales, labeller = labeller)
    # Finish
    breaks <- unique(dt[[x]])
    if (length(breaks)>50) breaks <- dt[, .SD[1], by = fill][[x]]
    p <- p + xlab(NULL) + scale_x_discrete(breaks = breaks) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Return
    p
}


#' @rdname plot_violins
#' @export
plot_feature_violins <- function(
    object, assay = assayNames(object)[1], x = 'feature_id', 
    fill = 'feature_id', color = NULL, n = 9,
    facet = NULL, nrow = NULL, ncol = NULL, dir = 'h', scales = 'free',
    labeller = label_value, highlight = NULL, fixed = list(na.rm = TRUE)
){
    object %<>% extract_features_evenly(n)
    plot_violins(
        object, assay = assay, x = x, fill = fill, color = color, facet = facet, 
        nrow = nrow, ncol = ncol, dir  = dir, labeller = labeller,
        highlight = highlight, fixed = fixed) + 
        ggtitle('Feature Violins')
}


#' @rdname plot_violins
#' @export
plot_sample_violins <- function(
    object, assay = assayNames(object)[1], x = 'sample_id', fill = 'sample_id', 
    color = NULL, n = 100, facet = NULL, nrow = NULL, ncol = NULL, dir = 'h', 
    scales = 'free', labeller = label_value, highlight = NULL, fixed = list(na.rm = TRUE)
){
    object %<>% extract_samples_evenly(n)
    plot_violins(
        object, assay = assay, x = x, fill = fill, color = color, facet = facet, 
        nrow = nrow, ncol = ncol, dir = dir, scales = scales, 
        labeller = labeller, highlight = highlight, fixed = fixed) + 
    ggtitle('Sample Violins')
}


extract_evenly <- function(l, p){
    round(seq(1, l, length.out = p))
}

extract_samples_evenly <- function(object, n){
    if (n < ncol(object)){
        object %<>% extract(, extract_evenly(ncol(object), n))}
    object
}

extract_features_evenly <- function(object, n){
    if (n < nrow(object)){
        object %<>% extract(extract_evenly(nrow(object), n), )
    }
    object
}


subgroup <- NULL
#' @rdname plot_violins
#' @export
plot_subgroup_violins <- function(
    object, assay = assayNames(object)[1], subgroup, x = 'subgroup', 
    fill = 'subgroup', color = NULL, highlight = NULL, facet = 'feature_id', 
    fixed = list(na.rm = TRUE)
){
    plot_violins(
       object, assay = assay(object), x = x, fill = fill, color = color,
        facet = facet, highlight = highlight, fixed = fixed) + 
    ggtitle('Subgroup violins')
}


#==============================================================================
#
#                   plot_boxplots()
#
#==============================================================================


#' Plot boxplots
#'
#' @param object        SummarizedExperiment
#' @param assay         string
#' @param subgroup      svar (string)
#' @param x             svar (string)
#' @param fill          svar (string)
#' @param color         svar (string)
#' @param block         svar used to connect points (string)
#' @param facet         svar (string)
#' @param scales       'free', 'fixed', 'free_x', 'free_y'
#' @param nrow          number of facet rows
#' @param ncol          number of facet columns 
#' @param page          number of facet pages: \code{\link[ggforce]{facet_wrap_paginate}}
#' @param highlight     fvar expressing which feature should be highlighted (string)
#' @param points        TRUE or FALSE
#' @param jitter        jitter width (number)
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
#'     fdata(object) %<>% cbind(control = .$feature_name %in% controlfeatures)
#' # Plot
#'     vars <- ggplot2::vars
#'     plot_boxplots(object[,1:9],  x = 'sample_id',  fill = 'sample_id' )
#'     plot_boxplots(object[1:9,],  x = 'feature_id', fill = 'feature_id')
#'     plot_boxplots(object[1:9, ], x = 'SET', fill = 'SET', facet = 'feature_id')
#'     plot_boxplots(object[1, ],   x = 'SET', fill = 'SET', jitter = TRUE)
#'     plot_boxplots(object[1, ],   x = 'SET', fill = 'SET', block = 'SUB')
#'     plot_feature_boxplots(object[1:9, ])
#'     plot_sample_boxplots(object[, 1:12])
#'     plot_sample_boxplots(object[, 1:12], highlight = 'control')
#'     plot_subgroup_boxplots(object[1:2, ], subgroup = 'SET')
#'     plot_subgroup_boxplots(object[1:2, ], subgroup = 'SET', block = 'SUB')
#' @export
plot_boxplots <- function(
    object, assay = assayNames(object)[1], 
    x = 'subgroup',  fill = 'subgroup', color = NULL, 
    block = NULL, facet = NULL, scales = 'free_y', nrow = NULL, ncol = NULL, 
    page = 1, labeller = 'label_value', highlight = NULL, 
    points = if (is.null(block)) FALSE else TRUE, 
    jitter = if (is.null(block)) 0.1 else 0,
    fillpalette  = make_var_palette(object, fill), 
    colorpalette = make_var_palette(object, color),
    hlevels = NULL, ...
){
# Assert/Process
    assert_is_all_of(object, "SummarizedExperiment")
    if (nrow(object)==0)  return(ggplot2::ggplot())
    if (!is.null(x))          assert_is_a_string(x)
    if (!is.null(fill))       assert_is_a_string(fill)
    if (!is.null(color))      assert_is_a_string(color)
    if (!is.null(block))      assert_is_a_string(block)
    if (!is.null(facet))      assert_is_character(facet)
    if (!is.null(nrow))       assert_is_a_number(nrow)
    if (!is.null(ncol))       assert_is_a_number(ncol)
    if (!is.null(highlight))  assert_is_a_string(highlight)
    if (!is.null(x))          assert_is_subset(x,          c(svars(object), fvars(object)))
    if (!is.null(fill))       assert_is_subset(fill,       c(svars(object), fvars(object)))
    if (!is.null(color))      assert_is_subset(color,      c(svars(object), fvars(object)))
    if (!is.null(block))      assert_is_subset(block,      c(svars(object), fvars(object)))
    if (!is.null(facet))      assert_is_subset(facet,      c(svars(object), fvars(object)))
    if (!is.null(highlight))  assert_is_subset(highlight, c(svars(object), fvars(object)))
                              assert_is_subset(scales, c('fixed', 'free', 'free_x', 'free_y'))
# Prepare
    xsym     <- sym(x)
    fillsym  <- if (is.null(fill))   quo(NULL) else  sym(fill)
    colorsym <- if (is.null(color))  quo(NULL) else  sym(color)
    blocksym <- if (is.null(block))  quo(NULL) else  sym(block)
    plotvars <- 'feature_name'
    if (!is.null(x))          plotvars %<>% c(x)         %>% unique()
    if (!is.null(fill))       plotvars %<>% c(fill)      %>% unique()
    if (!is.null(color))      plotvars %<>% c(color)     %>% unique()
    if (!is.null(block))      plotvars %<>% c(block)     %>% unique()
    if (!is.null(highlight))  plotvars %<>% c(highlight) %>% unique()
    if (!is.null(facet))      plotvars %<>% c(facet) %>% unique()
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    if (!is.null(x))   object[[x]] %<>% num2char()
    dt <- sumexp_to_longdt(object, assay = assay, svars = plottedsvars, fvars = plottedfvars)
    dt[, medianvalue := median(value, na.rm = TRUE), by = c('feature_id', x)]
    for (facetvar in facet){ 
        names(dt) %<>% stri_replace_first_fixed(facetvar, make.names(facetvar))
        facet %<>% stri_replace_first_fixed(facetvar, make.names(facetvar))
    } # otherwise facet_wrap_paginate thinks `fdr~coef~limma` is a formula
# Plot
    p <- ggplot(dt)
    if (!is.null(facet)) p <- p + facet_wrap_paginate(
        facets = facet, scales = scales, nrow = nrow, ncol = ncol, 
        page = page, labeller = labeller)
    outlier.shape <- if (points) NA else 19
    p <- p + geom_boxplot(aes(x = !!xsym, y = value, fill = !!fillsym, color = !!colorsym), 
                          outlier.shape = outlier.shape, na.rm = TRUE)
    p <- add_color_scale(p, color, data = dt, palette = colorpalette)
    p <- add_fill_scale( p, fill,  data = dt, palette = fillpalette)
    # Connect blocks
    if (!is.null(block)){
        byvar <- block
        if (!is.null(facet)) byvar %<>% c(facet)
        #dt[, direction := get(fill)[which.max(value)], by = byvar]
        p <- p + geom_line(aes(x = !!xsym, y = value, color = !!xsym, group = !!blocksym), na.rm = TRUE) # color = direction
        p <- add_color_scale(p, x, data = dt, palette = fillpalette)    # 'direction'
    }
    # Points
    if (points){
        p <- p + geom_jitter(aes(x = !!xsym, y = value),
                position = position_jitter(width = jitter, height = 0), size = 0.5, na.rm = TRUE)
    }
    p %<>% add_highlights(x = x, hl = highlight, geom = geom_point, fixed_color = "darkred")
    # Add hline
    if (!is.null(hlevels)){
        mediandt <- unique(dt[, unique(c('feature_id', x, 'medianvalue', facet)), with = FALSE])
        mediandt[, present := FALSE]
        mediandt[get(x) %in% hlevels, present := TRUE]
        p <- p + geom_hline(data = mediandt, 
                            aes(yintercept = medianvalue, color = !!fill, alpha = present), 
                            linetype = 'longdash') }
    # Finish and Return
    breaks <- unique(dt[[x]])
    if (length(breaks)>50) breaks <- dt[, .SD[1], by = fill][[x]]
    p <- p + xlab(NULL) + scale_x_discrete(breaks = breaks) + 
        guides(color = 'none', alpha = 'none') +
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p
}

#============================================================================
#
#                       plot_sample_boxplots
#                       plot_feature_boxplots
#                       plot_subgroup_boxplots
#                       plot_top_boxplots
#
#============================================================================

#' @rdname plot_boxplots
#' @export
plot_sample_boxplots <- function(
    object, assay = assayNames(object)[1], 
    x = 'sample_id', fill = 'sample_id', color = NULL, highlight = NULL,
    points = FALSE, jitter = 0.1, palette = NULL, n = 100,
    facet = NULL, scales = 'free_x', nrow = NULL, ncol = NULL, page = 1, 
    labeller = 'label_value'
){
    object %<>% extract_samples_evenly(n)
    plot_boxplots(
        object, assay = assay, 
        x = x, fill = fill, color = color,
        highlight = highlight, points = points, 
        jitter = jitter, palette = palette,  
        facet = facet, scales = scales, nrow = nrow, ncol = ncol, page = page, 
        labeller = labeller) + 
    ggtitle('Sample Boxplots')
}


#' @rdname plot_boxplots
#' @export
plot_feature_boxplots <- function(
    object, assay = assayNames(object)[1],
    x = 'feature_id', fill = 'feature_id', color = NULL, highlight = NULL,
    points = FALSE, jitter = 0.1, palette = NULL, n = 9,
    facet = NULL, scales = 'free_y', nrow = NULL, ncol = NULL, page = 1, 
    labeller = 'label_value'
){
    object %<>% extract_features_evenly(n)
    plot_boxplots(
        object, assay = assay,
        x = x, fill = fill, color = color,
        highlight = highlight, points = points, 
        jitter = jitter, palette = palette,  
        facet = facet, scales = scales, nrow = nrow, ncol = ncol, page = page, 
        labeller = labeller) + 
    ggtitle('Feature Boxplots')
}


#' @rdname plot_boxplots
#' @export
plot_subgroup_boxplots <- function(
    object, assay = assayNames(object)[1],
    subgroup = 'subgroup', x = subgroup, fill = subgroup, 
    color = NULL, block = NULL, highlight = NULL, jitter = TRUE, n = 9,
    facet = 'feature_id', scales = 'free_y', nrow = NULL, ncol = NULL, 
    page = 1, labeller = 'label_value',
    palette = make_subgroup_palette(object), fixed = list(na.rm=TRUE), hlevels = NULL
){
    plot_boxplots(
        object, assay = assay,
        x = x, fill = fill, color = color,
        block = block, highlight = highlight, 
        facet = facet, scales = scales, nrow = nrow,
        ncol = ncol, page = page,  labeller = labeller, 
        jitter = jitter, palette = palette, 
        fixed = fixed, hlevels = hlevels) + 
    ggtitle(('Subgroup Boxplots'))
    #p <- p + stat_summary(fun='mean', geom='line', aes(group=1), color='black', size=0.7, na.rm = TRUE)
}


#' Filter coefficient features
#' @param object      SummarizedXExperiment
#' @param fit         string
#' @param coef        string
#' @param effectsize  effectsize cutoff
#' @param p           p cutoff
#' @param fdr         fdr cutoff
#' @param ntop        ntop cutoff
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, fit = 'limma')
#' nrow(object)                                  # 4534 features
#' nrow(filter_top_features(object, p = 1))      # 3772   p != NA
#' nrow(filter_top_features(object))             # 1961   p < 0.05 
#' nrow(filter_top_features(object, fdr = 0.05)) # 1622 fdr < 0.05 
#' @export
filter_top_features <- function(
    object, 
    fit         = fits(object)[1], 
    coefficient = default_coefficient(object, fit = fit),
    effectsize  = 0,
    p           = 0.05,
    fdr         = 1,
    ntop        = Inf
){
# Assert
    if (length(coefficient)==0)  return(object)
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fit,  fits(object))
    assert_is_subset(coefficient, coefficients(object, fit = fit))
    assert_is_a_number(effectsize)
    assert_is_a_number(p)
    assert_is_a_number(fdr)
    assert_is_a_number(ntop)
    
# Filter
    pvar0 <- pvar(object, fit = fit, coef = coefficient)   # available
    idx <- !is.na(fdt(object)[[pvar0]])
    object %<>% extract(idx, )
    
    ntop %<>% min(nrow(object))                            # top
    idx <- order(fdt(object)[[pvar0]])[1:ntop]
    object %<>% extract(idx, )
                                                           # significant
    idx <- (abs(autonomics::effect(object, fit = fit, coef = coefficient)[, 1]) > effectsize)  & 
           (    autonomics::p(     object, fit = fit, coef = coefficient)[, 1]  < p  )         & 
           (    autonomics::fdr(   object, fit = fit, coef = coefficient)[, 1]  < fdr)
    object %<>% extract(idx, )
    
# Return
    object
}

format_coef_vars <- function(
    object, coef = setdiff(coefficients(object), 'Intercept')[1], fit = fits(object)[1]
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
#' @param jitter       TRUE or FALSE
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
#' plot_top_boxplots(
#'     object, subgroup = SampleGroup, block = Subject_ID, 
#'     fit='limma', coef = 'SampleGroupt1')
#' plot_top_boxplots(
#'     object, subgroup=SampleGroup, block = Subject_ID, fit='limma', 
#'     coef = 'SampleGroupt1', fdrcutoff=0.05)
#' plot_top_boxplots(
#'     object, subgroup = SampleGroup, block = Subject_ID, 
#'      fit = c('limma', 'lme'), contrast = 'SampleGroupt1', 
#'      fdrcutoff = 0.05)
#' @export
plot_top_boxplots <- function(
    object, assay = assayNames(object)[1], 
    subgroup = 'subgroup', block = NULL, fit = fits(object)[1], 
    coef = setdiff(coefficients(object), 'Intercept')[1],
    facet = c('feature_id', paste('fdr', coef, fit, sep = FITSEP)),
    fdrcutoff = 0.05, 
    jitter = if (is.null(block)) 0.1 else 0,
    palette = NULL, title = coef, ylab = NULL, 
    nrow = NULL, ncol = NULL, ntop = 4,
    labeller = 'label_value', scales = 'free_y'
){
# Order/Extract on p value
    . <- NULL
    subgroupsym  <- if (is.null(subgroup)) quo(NULL) else sym(subgroup)
    blocksym     <- if (is.null(block))    quo(NULL) else sym(block)
    
    svars0 <- c(subgroup, block)
    fdrvar    <- paste('fdr',    coef, fit, sep = FITSEP)
    pvar      <- paste('p',      coef, fit, sep = FITSEP)
    effectvar <- paste('effect', coef, fit, sep = FITSEP)
    
    fvars0 <- c(facet, fdrvar, pvar, effectvar)
    object %<>% filter_top_features(coef = coef, fit = fit, fdr = fdrcutoff, ntop = ntop)
    object %<>% format_coef_vars(     coef = coef, fit = fit)
# Prepare
    dt <- sumexp_to_longdt(object, assay = assay, svars = svars0, fvars = fvars0)
    mediandt <- summarize_median(dt, subgroup)
    contrastsubgroup <- NULL
    #mediandt[, contrastsubgroup := get(subgroupvar) %in% c(uplevels, downlevels)]
    facetstr <- vapply(facet, as_name, character(1))
    idx <- order(as.numeric(fdt(object)[[pvar]]))  # order on significance
    object %<>% extract(idx, )
# Plot
    p <- plot_subgroup_boxplots(
            object, subgroup = subgroup, block = block, x = subgroup, 
            fill = subgroup, facet = facet, scales = scales, nrow = nrow, 
            ncol = ncol, labeller = labeller, jitter = jitter, palette  = palette)  + 
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
    plotdt <- sumexp_to_longdt(subject,
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
    contrasts = coefficients(object), #, svars = as_name(!!enquo(subgroup))),
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
#' object <- read_metabolon(file, fit = 'limma', plot = FALSE)
#' idx <- order(fdata(object)$F.p.limma)[1:9]
#' object %<>% extract(idx, )
#' plot_sample_boxplots(  object)
#' plot_feature_boxplots( object)
#' plot_subgroup_boxplots(object, subgroup = 'SET')
#' plot_subgroup_points(  object, subgroup = 'SET')
#' plot_subgroup_points(  object, subgroup = 'SET', block = 'SUB')
#' @export
plot_subgroup_points <- function(
    object, subgroup = 'subgroup', block = NULL, x = subgroup, 
    color = subgroup, group = block, 
    facet = 'feature_id', nrow = NULL, scales = 'free_y', ...,
    palette = NULL,
    fixed = list(na.rm=TRUE),  #element_text(angle=90, vjust=0.5),
    theme = list(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
){
    dt <- sumexp_to_longdt(object, svars = svars(object), fvars = fvars(object))
    value <- NULL
    xsym     <- if (is.null(x))     quo(NULL) else sym(x)
    colorsym <- if (is.null(color)) quo(NULL) else sym(color)
    groupsym <- if (is.null(group)) quo(NULL) else sym(group)
    blocksym <- if (is.null(block)) quo(NULL) else sym(block)
    
    p <- plot_data(dt, geom = geom_point, 
                   x = !!xsym, y = value, color = !!colorsym, 
                   group = !!groupsym, ..., palette = palette, fixed = fixed)
    if (!is.null(block))  p <- p + geom_line()
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
    for (i in seq_along(x))  mat[x[[i]], i] <- 1
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
#' @return no return (base R plot)
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

#' Plot model 
#' @param object ´SummarizedExperiment
#' @return ggplot
#' @examples
#' require(magrittr)
#' file <- download_data('billing19.proteingroups.txt')
#' subgroups <- paste0(c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'), '_STD')
#' object <- read_proteingroups(file, subgroups = subgroups)
#' object$subgroup %<>% substr(1,3)
#' plot_design(object)
#' @export
plot_design <- function(object){
    designmat <- create_design(object, subgroupvar = 'subgroup', drop = TRUE)
    rownames(designmat) <- object$subgroup
    designmat %<>% unique()
    subgroups <- subgroup_levels(object)
    designmat %<>% extract(subgroups, )
    coefficients <- colnames(designmat)
    ymat <- matrix(seq_along(subgroups), nrow = ncol(designmat), ncol = 1)
    betamat <- solve(designmat) %*% ymat
    betamat[1,1] <- 1 # not strictly required, but plot is nicer if Intercept 
                      # is 1 unit long (in MASS:contr.sdif it gets much longer, 
                      # I think to maintain orthogonality of design)
    plotdt <- data.table(subgroup = subgroups, 
                         coef     = coefficients, 
                         x        = seq_along(subgroups),
                         yend     = seq_along(subgroups),
                         y        = seq_along(subgroups) - betamat[, 1])
    arrow <- arrow(length = unit(0.15, 'in'))
    
    ggplot(plotdt) + theme_bw() + 
    geom_segment(aes(x = x-0.05, xend = x+0.05, y = yend, yend = yend)) + 
    geom_text(   aes(x = x+0.1, y = yend, label = subgroup), hjust = 0) + 
    geom_segment(aes(x = x, xend = x, y = y, yend = yend), arrow = arrow) + 
    geom_label(  aes(x = x, y = y + (yend-y)/2, label = coef), parse = TRUE) +
    xlab(NULL) + ylab(NULL) + 
    theme_void()
    # theme(axis.text = element_blank(), 
    #       panel.grid.major.x = element_blank(), 
    #       panel.grid.minor.x = element_blank())
}

#' Plot top heatmap
#' @param object    SummarizedExperiment
#' @param fvar      string
#' @param lowColor  color for downregulation (Default color yellow)
#' @param highColor color for upregulation (Default color royal blue)
#' @param midColor  color for midpoint (Default color gray)
#' @param width     number: width  in inches (Default width 2 inches)
#' @param height    number: height in inches (Default height 4 inches)
#' @param fontsize  number: font size for feature name (Default size 0 for not showing feature names)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, fit = 'limma')
#' @export
plot_top_heatmap <- function(
    object,
    fit         = fits(object)[1],
    coefficient = default_coefficient(object, fit = fit),
    effectsize  = 0,
    p           = 0.5,
    fdr         = 1,
    ntop        = 9,
    fvar        = 'feature_id',
    lowColor    = "yellow",
    highColor   = "blue1",
    midColor    = "gray",
    width       = 2,
    height      = 4,
    fontsize    = 0
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fit,                 fits(object))
    assert_is_subset(coefficient, coefficients(object))
    
    object %<>% filter_top_features(
        fit        = fit,             coefficient = coefficient,
        effectsize = effectsize,      p           = p,
        fdr        = fdr,             ntop        = ntop)
  
   autonomics.import::fnames(object) <- object %>%
                                        autonomics.import::fvalues(fvar) %>%
                                        autonomics.support::uniquify()
   autonomics.support::cmessage('\t\tImpute NA(N) values')
   object %<>% autonomics.preprocess::impute(method = 'impute.QRILC')
   x <- autonomics.import::exprs(object)
   hres <- stats::hclust(stats::as.dist(1 - stats::cor(t(x), method="pearson"))) #compute correlation for matrix
   h.row.order <- hres$order # clustering on row
   row.order.df <- x[h.row.order,]
   scaled.x = t(scale(t(row.order.df),center=T)) #scaling dataframe
   df.molten.dat <- reshape2::melt(scaled.x) #resturcturing dataframe

   myPlot <- ggplot2::ggplot(data = df.molten.dat, ggplot2::aes_string(x = 'Var2', y = 'Var1', fill = 'value')) +
             ggplot2::geom_tile() +
             #ggplot2::facet_wrap(~ branch, scales = 'free') +
             ggplot2::scale_fill_gradient2(low = lowColor, mid = midColor, high = highColor) +
             ggplot2::theme(axis.text.x   = ggplot2::element_text(angle = 90, hjust = 1, size=width),
                            axis.text.y   = ggplot2::element_text(size=fontsize),
                            axis.ticks.length = grid::unit(.01, "cm"),
                            axis.ticks    = ggplot2::element_blank(),
                            axis.title.x  = ggplot2::element_blank(),
                            axis.title.y  = ggplot2::element_blank(),
                            legend.title  = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank(),
                            legend.justification = "top",
                            legend.key.width =   grid::unit(0.2, "cm"),
                            legend.key.height =  grid::unit(0.3,"cm"),
                            legend.text =        ggplot2::element_text(size=fontsize),
                            legend.box.margin =  ggplot2::margin(-10,-10,-10,-12))

   myPlot
}
