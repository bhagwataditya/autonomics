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
#'     object <- read_maxquant_proteingroups(file, invert_subgroups = invert_subgroups, plot = FALSE)
#'     p <- plot_sample_densities(object)
#'     add_color_scale(p, color = 'subgroup', data = sdt(object))
#' # STEMCELL INTENSITIES
#'    file <- download_data('billing16.proteingroups.txt')
#'    object <- read_maxquant_proteingroups(file, quantity = 'Intensity labeled', plot = FALSE)
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
    if (is.null(svar))               return(NULL)
    if (is.numeric(object[[svar]]))  return(NULL)
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
    p <- p + do.call(ggplot2::theme, {{theme}})

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
#' @param assay       string
#' @param group       svar (string)
#' @param fill        svar (string)
#' @param color       svar (string)
#' @param facet       svar (character vector)
#' @param nrow        number of facet rows
#' @param ncol        number of facet cols
#' @param dir         'h' (horizontal) or 'v' (vertical)
#' @param labeller    e.g. label_value
#' @param palette     named character vector
#' @param fixed       fixed aesthetics
#' @seealso \code{\link{plot_sample_violins}},
#'          \code{\link{boxplot_samples}}
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
#'     boxplot_samples( object, facet = 'SET')
#'     
#' # Feature distributions
#'     plot_feature_densities(object)
#'     plot_feature_violins(  object)
#'     plot_feature_boxplots( object)
#' @export
plot_densities <- function(
    object, assay = assayNames(object)[1], group, fill, color = NULL, 
    facet = NULL, nrow = NULL, ncol = NULL, dir = 'h', scales = 'free_y', 
    labeller = label_value, 
    palette = NULL, fixed = list(alpha = 0.8, na.rm = TRUE)
){
# Assert / Process
    assert_is_valid_sumexp(object)
    assert_scalar_subset(assay, assayNames(object))
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
    dt <- sumexp_to_longdt(object, assay = assay, svars = plottedsvars, fvars = plottedfvars)
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
    assay    = assayNames(object)[1],
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
        assay    = assay,
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
    assay = assayNames(object)[1],
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
        assay    = assay,
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
#' @param n           number
#' @param group       svar (string)
#' @param facet       svar (character vector)
#' @param nrow        NULL or number
#' @param ncol        NULL or number
#' @param dir         'h' or 'v' : are facets filled horizontally or vertically ?
#' @param scales      'free', 'free_x', 'free_y', or 'fixed'
#' @param labeller    label_both or label_value
#' @param highlight   fvar expressing which feature should be highlighted (string)
#' @param palette     named color vector (character vector)
#' @param fixed       fixed aesthetics
#' @return  ggplot object
#' @seealso \code{\link{plot_exprs}},
#'          \code{\link{boxplot_samples}}
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
#               extract_coef_features
#                   .extract_p_features
#                   .extract_fdr_features
#                   .extract_effectsize_features
#                       ..extract_statistic_features
#                   .extract_sign_features
#                   .extract_n_features
#
#==============================================================================


cmessage <- function(pattern, ...) message(sprintf(pattern, ...))

..extract_statistic_features <- function(
    object, coefs, statistic, comparer, threshold, 
    fit = fits(object)[1], combiner = '|', verbose  = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    if (is.null(fit))   return(object)
    if (is.null(coefs)) return(object)
    assert_scalar_subset(statistic, c('p', 'fdr', 'effect', 'effectsize'))
    assert_scalar_subset(comparer,  c('<', '>', '=='))
    assert_is_a_number(threshold)
    assert_is_subset(fit, fits(object))
    assert_is_subset(coefs, autonomics::coefs(object, fit = fit))
    assert_scalar_subset(combiner, c('|', '&'))
    assert_is_a_bool(verbose)
# Filter
    fun <- getFromNamespace(sprintf('%smat', statistic), 'autonomics')
    x <- fun(object, fit = fit, coefs = coefs)
    if (is.null(x))  return(object)
    idx <- get(comparer)(x, threshold)
    idx[is.na(idx)] <- FALSE
    fun <- function(y) Reduce(get(combiner), y)
    idx %<>% apply(1, fun)
    idx %<>% unname()
# Return
    n0 <- length(idx)
    n1 <- sum(idx, na.rm = TRUE)
    if (verbose & n1<n0){
        combiner <- paste0(' ', combiner, ' ')
        cmessage('\t\t\tRetain %d/%d features: %s(%s) %s %s', 
                n1, n0, statistic, paste0(coefs, collapse = combiner), 
                comparer, as.character(threshold))  
    }
    object[idx, ]
}

#' @rdname extract_coef_features
#' @export
.extract_p_features <- function(
    object, coefs, p = 0.05, fit = fits(object), combiner = '|',  verbose = TRUE
){
    assert_is_fraction(p)
    ..extract_statistic_features(
        object    = object,        coefs     = coefs,
        statistic = 'p',           comparer  = '<',   
        threshold = p,             fit       = fit,
        combiner  = combiner,      verbose   = verbose )
}

#' @rdname extract_coef_features
#' @export
.extract_fdr_features <- function(
    object, coefs, fdr = 0.05, fit = fits(object), combiner = '|',  verbose = TRUE
){
    assert_is_fraction(fdr)
    ..extract_statistic_features(
        object    = object,        coefs     = coefs,
        statistic = 'fdr',         comparer  = '<',
        threshold = fdr,           fit       = fit,
        combiner  = combiner,      verbose   = verbose )
}


#' @rdname extract_coef_features
#' @export
.extract_effectsize_features <- function(
    object, coefs, effectsize = 1, fit = fits(object), combiner = '|',  verbose = TRUE
){
    assert_weakly_positive_number(effectsize)
    ..extract_statistic_features(
        object    = object,        coefs     = coefs, 
        statistic = 'effectsize',  comparer  = '>',
        threshold = effectsize,    fit       = fit,
        combiner  = combiner,      verbose   = verbose )
}

#' @rdname extract_coef_features
#' @export
.extract_sign_features <- function(
    object, 
    coefs, 
    sign, 
    fit      = fits(object)[1], 
    combiner = '|',
    verbose  = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(sign, c(-1, +1))
    if (is.null(fit))    return(object)
    if (is.null(coefs))  return(object)
# Filter
    x <- autonomics::effectmat(object, fit = fit, coefs = coefs)
    idx <- unname(apply(sign(x), 1, function(y)  Reduce(get(combiner), sign(y) %in% sign) ))
# Return
    n0 <- length(idx)
    n1 <- sum(idx, na.rm = TRUE)
    if (verbose & n1<n0){
        combiner <- paste0(' ', combiner, ' ')
        cmessage('\t\t\tRetain %d/%d features: sign(%s) %%in%% c(%s)', 
            n1, n0, paste0(coefs, collapse = combiner), paste0(sign,  collapse = ','))
    }
    object[idx, ]
}

#' Order on p 
#' @param object   SummarizedExperiment
#' @param coefs    string vector: subset of `coefs(object)`
#' @param fit      string vector: subset of `fits(object)`
#' @param verbose  TRUE or FALSE
#' @examples 
#' # Read
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#' # no limma
#'     object %<>% order_on_p()  # unchanged
#'     object %<>% order_on_t()
#' # with limma
#'     object %<>% fit_limma(coefs = c('t1', 't2', 't3'))
#'     (assays(order_on_p(object))[[1]] == assays(order_on_t(object))[[1]])
#'     identical(order_on_p(object) , order_on_t(object))
#' @return SummarizedExperiment
order_on_p <- function(
    object, 
    fit = autonomics::fits(object), 
    coefs = autonomics::coefs(object, fit = fit), 
    combiner = '|',
    verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    if (is.null(fit))  return(object)
    assert_is_subset(fit,   autonomics::fits( object))
    assert_is_subset(coefs, autonomics::coefs(object, fit = fit))
    assert_scalar_subset(combiner, c('|', '&'))
    assert_is_a_bool(verbose)
# Order    
    pmat <- autonomics::pmat(object, fit = fit, coefs = coefs)
    if (is.null(pmat))  return(object)
    if (verbose)   cmessage("\t\tp-order features on: %s (%s)", 
                            paste0(fit,   collapse = ', '), 
                            paste0(coefs, collapse = ', '))
    if (combiner == '|')  idx <- order(matrixStats::rowMins(pmat))
    if (combiner == '&')  idx <- order(matrixStats::rowMaxs(pmat))
# Return
    object[idx, ]
}


#' @rdname order_on_p
#' @export
order_on_effect <- function(
    object, 
    fit = autonomics::fits(object),
    coefs = autonomics::coefs(object, fit = fit),
    combiner = '|', 
    verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    if (is.null(fit))  return(object)
    assert_is_subset(fit,   autonomics::fits( object))
    assert_is_subset(coefs, autonomics::coefs(object, fit = fit))
    assert_scalar_subset(combiner, c('|', '&'))
    assert_is_a_bool((verbose))
# Order
    effectmat <- autonomics::effectmat(object, fit = fit, coefs = coefs)
    if (verbose)   cmessage("\t\tt-order features on: %s (%s)", 
                            paste0(fit,   collapse = ', '), 
                            paste0(coefs, collapse = ', '))
    if (combiner == '|')  idx <- order(matrixStats::rowMaxs(abs(effectmat)), decreasing = TRUE)
    if (combiner == '&')  idx <- order(matrixStats::rowMins(abs(effectmat)), decreasing = TRUE)
# Return
    object[idx, ]
}


#' @rdname extract_coef_features
#' @export
.extract_n_features <- function(object, coefs, combiner = '|', n, fit = fits(object)[1], verbose = TRUE){
# Assert
    assert_is_valid_sumexp(object)
    assert_positive_number(n)
# Filter
    object %<>% order_on_effect(fit = fit, coefs = coefs, combiner = combiner, verbose = FALSE)  # pls
    object %<>% order_on_p(     fit = fit, coefs = coefs, combiner = combiner, verbose = FALSE)  # glm
    n %<>% min(nrow(object))
    idx <- c(rep(TRUE, n), rep(FALSE, nrow(object)-n))
    n0 <- length(idx)
    n1 <- sum(idx, na.rm = TRUE)
    if (verbose & n1<n0){
        combiner <- paste0(' ', combiner, ' ')
        y <- paste0(coefs, collapse = combiner)
        cmessage('\t\t\tRetain %d/%d features: p(%s) or effect(%s) in best %d', 
                 n1, n0, y, y, n)
    }
# Return
    object[idx, ]
}


#' Extract coefficient features
#' @param object      SummarizedXExperiment
#' @param fit         subset of fits(object)
#' @param coefs       subset of coefs(object)
#' @param combiner    '|' or '&': how to combine multiple fits/coefs
#' @param p           p threshold
#' @param fdr         fdr threshold
#' @param effectsize  effectsize threshold
#' @param sign        effect sign
#' @param n           number of top features (Inf means all)
#' @param verbose     TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' # Read and Fit
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     object %<>% fit_limma()
#' # Single coef
#'     object0 <- object
#'     object %<>% .extract_p_features(         coefs = 't1', p = 0.05)
#'     object %<>% .extract_fdr_features(       coefs = 't1', fdr = 0.05)
#'     object %<>% .extract_effectsize_features(coefs = 't1', effectsize = 1)
#'     object %<>% .extract_sign_features(      coefs = 't1', sign = -1)
#'     object %<>% .extract_n_features(         coefs = 't1', n = 4)
#'     object <- object0
#'     object %<>%  extract_coef_features(
#'         coefs = 't1', p = 0.05, fdr = 0.05, effectsize = 1, sign = -1, n = 4)
#' # Multiple coefs
#'     object <- object0
#'     object %<>% .extract_p_features(         coefs = c('t1', 't2'), p = 0.05)
#'     object %<>% .extract_fdr_features(       coefs = c('t1', 't2'), fdr = 0.05)
#'     object %<>% .extract_effectsize_features(coefs = c('t1', 't2'), effectsize = 1)
#'     object %<>% .extract_sign_features(      coefs = c('t1', 't2'), sign = -1)
#'     object %<>% .extract_n_features(         coefs = c('t1', 't2'), n = 4)
#'     object <- object0
#'     object %<>%  extract_coef_features(
#'         coefs = c('t1', 't2'), p = 0.05, fdr = 0.05, effectsize = 1, sign = -1, n = 4)
#' @export
extract_coef_features <- function(
    object, 
    fit         = fits(object)[1], 
    coefs       = default_coefs(object, fit = fit),
    combiner    = '|',
    p           = 1, 
    fdr         = 1, 
    effectsize  = 0, 
    sign        = c(-1,+1), 
    n           = 4,
    verbose     = TRUE
){
# Assert
# Filter
    object %<>% .extract_p_features(         coefs = coefs, p = p,                   fit = fit, combiner = combiner, verbose = verbose)
    object %<>% .extract_fdr_features(       coefs = coefs, fdr = fdr,               fit = fit, combiner = combiner, verbose = verbose)
    object %<>% .extract_effectsize_features(coefs = coefs, effectsize = effectsize, fit = fit, combiner = combiner, verbose = verbose)
    object %<>% .extract_sign_features(      coefs = coefs, sign = sign,             fit = fit, combiner = combiner, verbose = verbose)
    object %<>% .extract_n_features(         coefs = coefs, n = n,                   fit = fit, combiner = combiner, verbose = verbose)
# Return
    object
}



format_coef_vars <- function(
    object, fit = fits(object)[1], coefs = default_coefs(object, fit = fit)[1]
){
    effectvars <- effectvar(object, coefs = coefs, fit = fit)
    pvars      <- pvar(     object, coefs = coefs, fit = fit)
    fdrvars    <- fdrvar(   object, coefs = coefs, fit = fit)
    for (var in c(effectvars, pvars, fdrvars)){
        fdt(object)[[var]] %<>% formatC(format='e', digits=0)
        fdt(object)[[var]] %<>% as.character()
        fdt(object)[[var]] %<>% paste0(split_extract_fixed(var, FITSEP, 2), ' : ',  
                                       split_extract_fixed(var, FITSEP, 1), ' = ', .)
    }
    object
}

#' Add facetvars
#' @param object  SummarizedExperiment
#' @param fit     string
#' @param coefs   string vector
#' @return  SummarizedExperiment
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma')
#' fdt(object)
#' fdt(add_facetvars(object))
#' @export
add_facetvars <- function(
    object, fit = fits(object)[1], coefs = default_coefs(object, fit = fit)
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(fit, fits(object))
    assert_is_subset(coefs, autonomics::coefs(object))
# Add
    for (i in seq_along(coefs)){
               pvar <- autonomics::pvar(     object, fit = fit, coefs = coefs)
             fdrvar <- autonomics::fdrvar(   object, fit = fit, coefs = coefs)
          effectvar <- autonomics::effectvar(object, fit = fit, coefs = coefs)
           facetvar <- paste0('facet.', coefs[[i]])
        assertive.sets::assert_are_disjoint_sets(facetvar, fvars(object))
        if (!is.null(pvar))            pvalues <- fdt(object)[[     pvar]] %>% formatC(format = 'e', digits = 0) %>% as.character() 
        if (!is.null(fdrvar))        fdrvalues <- fdt(object)[[   fdrvar]] %>% formatC(format = 'e', digits = 0) %>% as.character()
        if (!is.null(effectvar))  effectvalues <- fdt(object)[[effectvar]] %>% round(3)  %>% as.character()
        fdt(object)[[facetvar]] <- 
            if (is.null(pvar)){ sprintf('%s : %s', coefs[[i]], effectvalues)
            } else {            sprintf('%s : %s (%s)', coefs[[i]], fdrvalues, pvalues) 
            }
    }
# Return
    object
}


#==============================================================================
#
#               plot_exprs_per_coef
#                   plot_exprs
#                       .plot_exprs
#
#==============================================================================

.plot_exprs <- function(
    object, assay, geom, x, fill, color, shape, size, block, linetype, 
    highlight, facet, scales, nrow, ncol, page, labeller, 
    pointsize, jitter, colorpalette, fillpalette, hlevels, 
    title, subtitle, xlab, ylab, theme
){
# Prepare
    xsym        <- sym(x)
    fillsym     <- if (is.null(fill))      quo(NULL) else  sym(fill)
    colorsym    <- if (is.null(color))     quo(NULL) else  sym(color)
    shapesym    <- if (is.null(shape))     quo(NULL) else  sym(shape)
    sizesym     <- if (is.null(size))      quo(NULL) else  sym(size)
    blocksym    <- if (is.null(block))     quo(NULL) else  sym(block)
    linetypesym <- if (is.null(linetype))  quo(NULL) else  sym(linetype)
    plotvars <- 'feature_name'
    if (!is.null(x))          plotvars %<>% c(x)         %>% unique()
    if (!is.null(fill))       plotvars %<>% c(fill)      %>% unique()
    if (!is.null(color))      plotvars %<>% c(color)     %>% unique()
    if (!is.null(shape))      plotvars %<>% c(shape)     %>% unique()
    if (!is.null(size))       plotvars %<>% c(size)      %>% unique()
    if (!is.null(block))      plotvars %<>% c(block)     %>% unique()
    if (!is.null(linetype))   plotvars %<>% c(linetype)  %>% unique()
    if (!is.null(highlight))  plotvars %<>% c(highlight) %>% unique()
    if (!is.null(facet))      plotvars %<>% c(facet)     %>% unique()
    plottedsvars <- intersect(plotvars, svars(object))
    plottedfvars <- intersect(plotvars, fvars(object))
    # if (!is.null(x))   object[[x]] %<>% num2char()
    dt <- sumexp_to_longdt(object, assay = assay, svars = plottedsvars, fvars = plottedfvars)
    dt[, medianvalue := median(value, na.rm = TRUE), by = c('feature_id', x)]
    for (facetvar in facet){ 
        names(dt) %<>% stri_replace_first_fixed(facetvar, make.names(facetvar))
        facet %<>% stri_replace_first_fixed(facetvar, make.names(facetvar))
    } # otherwise facet_wrap_paginate thinks `fdr~coef~limma` is a formula
# Initialization
    p <- ggplot(dt) + theme_bw() + xlab(xlab) + ylab(ylab) + ggtitle(title, subtitle = subtitle)
    if (!is.numeric(dt[[x]]))  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    if (!is.null(facet))   p <- p + facet_wrap_paginate(facets = facet, 
        scales = scales, nrow = nrow, ncol = ncol, page = page, labeller = labeller)
# Boxplots/Points
    if (geom == 'boxplot'){
        outlier.shape <- if (pointsize==0) NA else 19
        mapping <- aes(x = !!xsym, y = value, fill = !!fillsym)
        p <- p + geom_boxplot(mapping = mapping, outlier.shape = outlier.shape, na.rm = TRUE)
        if (pointsize > 0){
            mapping <- aes(x = !!xsym, y = value)
            position <- position_jitter(width = jitter, height = 0)
            p <- p + geom_jitter(mapping = mapping, position = position, size = pointsize, na.rm = TRUE)
        }
    } else {
        mapping <- aes(x = !!xsym, y = value, color = !!colorsym, shape = !!shapesym, size = !! sizesym)
        p <- p + geom_point(mapping = mapping, na.rm = TRUE)
    }
    p <- add_color_scale(p, color, data = dt, palette = colorpalette)
    p <- add_fill_scale( p, fill,  data = dt, palette = fillpalette)
# Lines
    if (!is.null(block)){   
        byvar <- block
        if (!is.null(facet)) byvar %<>% c(facet)
        mapping <- aes(x = !!xsym, y = value, color = !!colorsym, group = !!blocksym, linetype = !!linetypesym)
        p <- p + geom_line(mapping = mapping, na.rm = TRUE)      # color = direction
    }
# Highlights (points)
    p %<>% add_highlights(x = x, hl = highlight, geom = geom_point, fixed_color = "darkred")
# Hlines
    if (!is.null(hlevels)){
        mediandt <- unique(dt[, unique(c('feature_id', x, 'medianvalue', facet)), with = FALSE])
        mediandt[, present := FALSE]
        mediandt[get(x) %in% hlevels, present := TRUE]
        mapping <- aes(yintercept = medianvalue, color = !!fill, alpha = present)
        p <- p + geom_hline(data = mediandt, mapping = mapping, linetype = 'longdash') 
    }
# Finish
    if (!is.numeric(dt[[x]])){
        breaks <- unique(dt[[x]])
        if (length(breaks)>50)  breaks <- dt[, .SD[1], by = fill][[x]]
        p <- p + scale_x_discrete(breaks = breaks) + guides(alpha = 'none')
    }
    if (!is.null(theme))  p <- p + theme
    p
}


#' Plot exprs
#' @param object        SummarizedExperiment
#' @param dim          'samples'   (per-sample distribution across features), \cr
#'                     'features' (per-feature distribution across samples ) or 
#'                     'both'        (subgroup distribution faceted per feature)
#' @param assay         string: value in assayNames(object)
#' @param geom          'boxplot' or 'point'
#' @param x                     x svar
#' @param color             color svar: points, lines
#' @param fill               fill svar: boxplots
#' @param shape             shape svar
#' @param size               size svar
#' @param block             group svar
#' @param linetype       linetype svar
#' @param highlight     highlight svar
#' @param fit          'limma', 'lm', 'lme', 'lmer', 'wilcoxon'
#' @param coefs         subset of coefs(object) to consider in selecting top
#' @param p             fraction: p   cutoff
#' @param fdr           fraction: fdr cutoff
#' @param facet         string: fvar mapped to facet
#' @param n             number of samples (dim = 'samples') or features (dim = 'features' or 'both') to plot
#' @param nrow          number of rows in faceted plot (if dim = 'both)
#' @param ncol          number of cols in faceted plot (if dim = 'both')
#' @param scales        'free_y', 'free'x', 'fixed'
#' @param labeller      string or function
#' @param pointsize     number
#' @param jitter        jitter width (number)
#' @param fillpalette   named character vector: fill palette
#' @param colorpalette  named character vector: color palette
#' @param hlevels       xlevels for which to plot hlines
#' @param title         string
#' @param subtitle      string
#' @param xlab          string
#' @param ylab          string
#' @param theme         ggplot2::theme(...) or NULL
#' @param file          NULL or filepath
#' @param width         inches
#' @param height        inches
#' @return ggplot object
#' @seealso \code{\link{plot_sample_densities}},
#'          \code{\link{plot_sample_violins}}
#' @examples 
#' # Without limma
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     plot_exprs(object, block = 'SUB', title = 'Subgroup Boxplots')
#'     plot_exprs(object, dim = 'samples')
#'     plot_exprs(object, dim = 'features', block = 'sample_id')
#' # With limma 
#'     object %<>% fit_limma(block = 'SUB')
#'     plot_exprs(object, block = 'SUB')
#'     plot_exprs(object, block = 'SUB', coefs = c('t1', 't2', 't3'))
#'     plot_exprs_per_coef(object, x = 'SET', block = 'SUB')
#' # Points
#'     plot_exprs(object, geom = 'point', block = 'SUB')
#' # Add highlights
#'     controlfeatures <- c('biotin','phosphate')
#'     fdt(object) %<>% cbind(control = .$feature_name %in% controlfeatures)
#'     plot_exprs(object, dim = 'samples', highlight = 'control')
#' # Multiple pages
#'     plot_exprs(object, block = 'SUB', n = 4, nrow = 1, ncol = 2)
#' @export
plot_exprs <- function(
    object, 
    dim          = 'both', 
    assay        = assayNames(object)[1], 
    x            = switch(dim, both = 'subgroup', features = 'feature_id', samples = 'sample_id'),  
    geom         = default_geom(object, x = x),
    color        = x, # points/lines
    fill         = x, # boxplots
    shape        = NULL,
    size         = NULL,
    block        = NULL, 
    linetype     = NULL,
    highlight    = NULL, 
    fit          = fits(object)[1], 
    coefs        = default_coefs(object, fit = fit)[1],
    combiner     = '|',
    p            = 1,
    fdr          = 1, 
    facet        = if (dim=='both')  'feature_id' else NULL,
    n            = 4,
    ncol         = NULL,
    nrow         = NULL,
    scales       = 'free_y', 
    labeller     = 'label_value', 
    pointsize    = if (is.null(block)) 0 else 0.5, 
    jitter       = if (is.null(block)) 0.1 else 0,
    fillpalette  = make_var_palette(object, fill), 
    colorpalette = make_var_palette(object, color),
    hlevels      = NULL, 
    title        = switch(dim, both = sprintf('%s~%s', paste0(coefs, collapse = ', '), fit), 
                               features = 'Feature Boxplots', 
                               samples  =  'Sample Boxplots'), 
    subtitle     = coefs,
    xlab         = NULL, 
    ylab         = 'value', 
    theme        = ggplot2::theme(plot.title = element_text(hjust = 0.5)), 
    file         = NULL,
    width        = 7, 
    height       = 7
){
# Assert
    assert_is_valid_sumexp(object)
    if (nrow(object)==0)      return(ggplot())
    assert_scalar_subset(dim, c('features', 'samples', 'both'))
    assert_scalar_subset(assay, assayNames(object))
    assert_scalar_subset(geom, c('boxplot', 'point'))
                              assert_scalar_subset(x,         c(svars(object), fvars(object)))
    if (!is.null(color))      assert_scalar_subset(color,     c(svars(object), fvars(object)))
    if (!is.null(fill))       assert_scalar_subset(fill,      c(svars(object), fvars(object)))
    if (!is.null(shape))      assert_scalar_subset(shape,     c(svars(object), fvars(object)))
    if (!is.null(size))       assert_scalar_subset(size,      c(svars(object), fvars(object)))
    if (!is.null(block))      assert_scalar_subset(block,     c(svars(object), fvars(object)))
    if (!is.null(linetype))   assert_scalar_subset(linetype,  c(svars(object), fvars(object)))
    if (!is.null(highlight))  assert_scalar_subset(highlight, c(svars(object), fvars(object)))
    if (!is.null(facet))      assert_is_subset(facet,         c(svars(object), fvars(object)))
    if (!is.null(nrow))       assert_is_a_number(nrow)
    if (!is.null(ncol))       assert_is_a_number(ncol)
    if (!is.null(facet))      assert_is_subset(scales, c('fixed', 'free', 'free_x', 'free_y'))
# Extract
    if        (dim == 'samples' ){   n %<>% min(ncol(object));  object %<>% extract_samples_evenly(n)
    } else if (dim == 'features'){   n %<>% min(nrow(object));  object %<>% extract_features_evenly(n)
    } else if (dim == 'both'){       n %<>% min(nrow(object)); 
        if (is.null(coefs)){         object %<>% extract_features_evenly(n) 
        } else {                     object %<>% extract_coef_features(fit = fit, coefs = coefs, combiner = combiner, p = p, fdr = fdr, n = n)
                                     object %<>% add_facetvars(fit = fit, coefs = coefs)
                                     facet %<>% c(sprintf('facet.%s', coefs))
                                     #object %<>% format_coef_vars(coefs = coefs, fit = fit) 
        }
    }
# Plot
    if (is.null(ncol))  ncol <- ceiling(sqrt(n))
    if (is.null(nrow))  nrow <- ceiling(n/ncol )   # https://stackoverflow.com/a/60110740
    npages <- if (dim == 'samples' ) 1  else  ceiling(nrow(object) / nrow / ncol)
    if (!is.null(file))   pdf(file, width = width, height = height)
    for (i in seq_len(npages)){
        p <- .plot_exprs(
            object,
            assay         = assay,             geom        = geom,
            x             = x,                 fill        = fill,
            color         = color,             shape       = shape,
            size          = size,
            block         = block,             linetype    = linetype,
            highlight     = highlight,         facet       = facet,     
            scales        = scales,            nrow        = nrow,
            ncol          = ncol,              page        = i,
            labeller      = labeller,          pointsize   = pointsize,  
            jitter        = jitter, 
            colorpalette  = colorpalette,      fillpalette = fillpalette, 
            hlevels       = hlevels,           title       = title,
            subtitle      = subtitle,
            xlab          = xlab,              ylab        = ylab,
            theme         = theme
        )
        if (npages>1)  print(p)
    }
    if (!is.null(file)){ dev.off(); file  
    } else {             p }
}


#' @rdname plot_exprs
#' @export
plot_boxplots <- function(...){
    .Deprecated('plot_exprs')
    plot_exprs(...)
}

#' @rdname plot_exprs
#' @export
boxplot_subgroups <- function(...){
    .Deprecated('plot_exprs')
    plot_exprs(...)
}

#' @rdname plot_exprs
#' @export
plot_sample_boxplots <- function(...){
    .Deprecated('plot_exprs(object, dim = "samples")')
    plot_exprs(..., dim = 'samples')
}

#' @rdname plot_exprs
#' @export
boxplot_samples <- function(...){
    .Deprecated('plot_exprs(object, dim = "samples")')
    plot_exprs(..., dim = 'samples')
}

#' @rdname plot_exprs
#' @export
plot_feature_boxplots <- function(...){
    .Deprecated('plot_exprs(object, dim = "features")')
    plot_exprs(..., dim = 'features')
}

#' @rdname plot_exprs
#' @export
boxplot_features <- function(...){
    .Deprecated('plot_exprs(object, dim = "features")')
    plot_exprs(..., dim = 'features')
}

#' @rdname plot_exprs
#' @export
plot_exprs_per_coef <- function(
    object, 
    x        = 'subgroup',
    geom     = default_geom(object, x),
    block    = NULL,
    coefs    = default_coefs(object),
    orderbyp = TRUE,
    title    = x,
    subtitle = coefs,
    nrow     = 1, 
    ncol     = NULL, 
    theme    = ggplot2::theme(legend.position = 'bottom', 
                              legend.title    = element_blank(), 
                              plot.title      = element_text(hjust = 0.5), 
                              plot.subtitle   = element_text(hjust = 0.5))
){
    assert_is_valid_sumexp(object)
    if (orderbyp){
        idx <- order(vapply(coefs, function(x)  min(pmat(object, coefs = x)), numeric(1)))
        coefs %<>% extract(idx)
        if (length(x)        > 1)         x %<>% extract(idx)
        if (length(geom)     > 1)      geom %<>% extract(idx)
        if (length(title)    > 1)     title %<>% extract(idx)
        if (length(subtitle) > 1)  subtitle %<>% extract(idx)
    }
    grobs <- mapply(
        plot_exprs, 
        x        = x, 
        geom     = geom,
        coefs    = coefs, 
        title    = title,
        MoreArgs = list(object = object, block = block, n = 1, theme = theme), 
        SIMPLIFY = FALSE)
    gridExtra::grid.arrange(grobs = grobs, nrow = nrow)
}

#' Default geom
#' @param object SummarizedExperiment
#' @param x      svar
#' @return character vector
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' svars(object)
#' default_geom(object, x = c('AGE', 'T2D'))
#' default_geom(object, x = 'AGE')
#' @export
default_geom <- function(object, x){
    sdt0 <- sdt(object)[, x, with = FALSE]
    y <- vapply(sdt0, class, character(1))
    y %<>% unname()
    y <- c(numeric = 'point', factor = 'boxplot', character = 'boxplot')[y]
    names(y) <- x
    y
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
    p <- p + do.call(ggplot2::theme, {{theme}})
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
#' object <- read_maxquant_proteingroups(file, subgroups = subgroups)
#' object$subgroup %<>% substr(1,3)
#' plot_design(object)
#' @export
plot_design <- function(object){
    designmat <- create_design(object, subgroupvar = 'subgroup', drop = TRUE)
    rownames(designmat) <- object$subgroup
    designmat %<>% unique()
    subgroups <- subgroup_levels(object)
    designmat %<>% extract(subgroups, )
    coefs <- colnames(designmat)
    ymat <- matrix(seq_along(subgroups), nrow = ncol(designmat), ncol = 1)
    betamat <- solve(designmat) %*% ymat
    betamat[1,1] <- 1 # not strictly required, but plot is nicer if Intercept 
                      # is 1 unit long (in MASS:contr.sdif it gets much longer, 
                      # I think to maintain orthogonality of design)
    plotdt <- data.table(subgroup = subgroups, 
                         coef     = coefs, 
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
#' @param object           SummarizedExperiment
#' @param assay            string: one of assayNames(object)
#' @param fit             'limma', 'lm', 'lme(r)', 'wilcoxon'
#' @param coef             string: one of coefs(object)
#' @param effectsize       number: effectsize filter
#' @param p                number: p    filter
#' @param fdr              number: fdr  filter
#' @param n                number: n filter
#' @param cluster_samples  TRUE or FALSE
#' @param flabel           string: feature label
#' @param group            sample groupvar
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file, fit = 'limma')
#' plot_top_heatmap(object)
#' plot_top_heatmap(object, cluster_samples = FALSE)
#' @export
plot_top_heatmap <- function(
    object,
    assay            = assayNames(object)[1],
    fit              = fits(object)[1],
    coef             = default_coefs(object, fit = fit)[1],
    effectsize       = 0,
    p                = 1,
    fdr              = 0.05,
    n                = 100,
    cluster_features = is.null(coef),
    cluster_samples  = FALSE,
    flabel           = intersect(c('gene', 'feature_id'), fvars(object))[1], 
    group            = 'subgroup'
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(assay,         assayNames(object))
    assert_is_a_string(fit)
    assert_is_subset(  fit, fits(object))
    if (!is.null(coef))  assert_is_a_string(coef)
    if (!is.null(coef))  assert_is_subset(  coef, coefs(object))
    assert_is_a_number(effectsize)
    assert_is_a_number(p)
    assert_is_a_number(fdr)
    assert_is_a_number(n)
    assert_is_a_string(flabel)
    assert_is_subset(  flabel, fvars(object))
    assert_is_subset(group,  svars(object))
# Filter: significant features
    object0 <- object
    if (is.null(coef)){   object %<>% extract_features_evenly(n)
    } else {              object %<>% extract_coef_features(
                            fit = fit, coefs = coef, effectsize = effectsize, p = p, fdr = fdr, n = n) }
# Zscore
    assays(object)[[assay]] %<>% t() %>% scale(center = TRUE, scale = TRUE) %>% t()
    assays(object)[[assay]] %<>% na_to_zero()
# Order features                                # in an edge case one of the groups had no obs
    idx <- rowSds(assays(object)[[assay]]) > 0  # still limma::lmFit produced a p value - limma bug ?
    object %<>% extract(idx, )                  # leads to a 0 variance error in the next line
    if (cluster_features){
        idx <- hclust(as.dist(1-cor(t(assays(object)[[assay]]))))$order
        object  %<>% extract(  idx , )                          # order features
    }
    if (!is.null(coef)){
        idx <- effectmat(object, fit = fit, coef = coef)[, 1] < 0; down <- object[idx, ]  # split down/up
        idx <- effectmat(object, fit = fit, coef = coef)[, 1] > 0;   up <- object[idx, ]
        object <- rbind(rev(down), rev(up))
    }
# Add pvalues
    if (!is.null(coef)){
        pvar   <- paste('p',   coef, fit, sep = FITSEP)
        fdrvar <- paste('fdr', coef, fit, sep = FITSEP)
        pvalues   <- fdt(object)[[  pvar]] %>% formatC(format = 'e', digits = 0) %>% as.character() 
        fdrvalues <- fdt(object)[[fdrvar]] %>% formatC(format = 'e', digits = 0) %>% as.character()
        fdt(object)[[flabel]] %<>% paste0('  ', pvalues, '  ', fdrvalues)
    }
    fdt(object)[[flabel]] %<>% factor(unique(.))                # fix order
# Order samples
    if (cluster_samples){
        idx <- matrixStats::colSds(assays(object)[[assay]]) > 0 # this ad-hoc dropping of samples is undesirable
        object %<>% extract(, idx)
        idx <- hclust(as.dist(1-cor(assays(object)[[assay]])))$order
        object %<>% extract(, idx)                              # order samples
    }
    object %<>% split_samples(group)                            # split by group
    object %<>% Reduce(cbind, .)                                # cbind
    sdt(object)$sample_id %<>% factor(unique(.))                # fix order
# Prepare
    dt <- sumexp_to_longdt(object, assay = assay, fvars = flabel)
    setnames(dt, 'value', 'z-score')
    vlines <- 0.5 + c(0, cumsum(table(object[[group]])))
    if (!is.null(coef)){
        hlines <- 0.5 + c(0, sum(effectmat(object, fit = fit, coef = coef)[, 1] < 0), nrow(object))
    }
# Plot
    p <- ggplot(data = dt, aes(x = sample_id, y = !!sym(flabel), fill = `z-score`)) +
        geom_tile() +
        theme_minimal() + xlab(NULL) + ylab(NULL) + 
        scale_x_discrete(position = 'top') + 
        theme(axis.text.x = element_text(angle = 90, hjust = 0)) + 
        scale_fill_gradient2(low = '#ff5050', high = '#009933', na.value = 'white') + 
        geom_vline(xintercept = vlines)
    if (!is.null(coef)){
        p <- p + geom_hline(yintercept = hlines)
    }
    p
}
