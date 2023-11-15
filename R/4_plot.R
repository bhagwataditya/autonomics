#=============================================================================
#
#                    add_color_scale
#                    add_fill_scale
#
#==============================================================================

#' Add color scale
#' @param object   SummarizedExperiment
#' @param color    string: svar mapped to color
#' @param show     TRUE or FALSE (default)
#' @param verbose  TRUE or FALSE (default)
#' @return default color values vector
#' @examples
#' # STEMCELL RATIOS
#'     file <- download_data('billing16.proteingroups.txt')
#'     invert_subgroups <- c('E_EM','BM_E', 'BM_EM')
#'     object <- read_maxquant_proteingroups(file, invert = invert_subgroups)
#'     p <- plot_sample_densities(object)
#'     add_color_scale(p, color = 'subgroup', data = sdt(object))
#' # STEMCELL INTENSITIES
#'    file <- download_data('billing16.proteingroups.txt')
#'    object <- read_maxquant_proteingroups(file, quantity = 'labeledintensity')
#'    p <- plot_sample_densities(object)
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
    return(make_onefactor_colors(varlevels, show = show, verbose = verbose))
}

#' Create default ggplot colors for factor levels
#' @param varlevels  string vector
#' @param h          start hue
#' @param l          luminance
#' @param show       TRUE/FALSE
#' @param verbose    TRUE/FALSE
#' @return string vector: elements = colors, names = factor levels
#' @author John Colby
#' @references https://stackoverflow.com/questions/8197559
#' @noRd
make_onefactor_colors <- function(
    varlevels, h = 15, l = 65, show = FALSE, verbose = TRUE
){
    n <- length(varlevels)
    hues <- seq(h, h + 360, length = n + 1)
    colors <- hcl(h = hues, l = l, c = 100)[seq_len(n)] %>%
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
#' varlevels <- slevels(object, 'subgroup')
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
        # OLD IMPLEMENTATION
        #     colors  %<>%  c(sequential_hcl(
        #                        n2, h = hues[[i]], power = 1, c = c(50, 100),
        #                        l = c(90, 30)) %>%
        #                    set_names(paste0(V1levels[[i]], sep, V2levels)))
        basecolor  <- hcl(h = hues[[i]], c = 100, l = 50)
        newcolors <- grDevices::colorRampPalette(c('white', basecolor))(n2+1)[-1]
        names(newcolors) <- paste0(V1levels[[i]], sep, V2levels)
        colors %<>% c(newcolors)
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
#' @param linetype    variable mapped to linetype (symbol)
#' @param ...         mapped aesthetics
#' @param palette     color palette (named character vector)
#' @param fixed       fixed  aesthetics (list)
#' @param theme       list with ggplot theme specifications
#' @return ggplot object
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% pca()
#' data <- sdt(object)
#' plot_data(data, x = `effect~sample_id~pca1`, y = `effect~sample_id~pca2`)
#' plot_data(data, x = `effect~sample_id~pca1`, y = `effect~sample_id~pca2`, color = TIME_POINT)
#' data$TIME <- as.numeric(substr(data$TIME_POINT, 2, 3))
#' plot_data(data, x = `effect~sample_id~pca1`, y = `effect~sample_id~pca2`, color = TIME)
#' plot_data(data, x = `effect~sample_id~pca1`, y = `effect~sample_id~pca2`, color = NULL)
#' fixed <- list(shape = 15, size = 3)
#' plot_data(data, x = `effect~sample_id~pca1`, y = `effect~sample_id~pca2`, fixed = fixed)
#' @author Aditya Bhagwat, Johannes Graumann
#' @export
plot_data <- function(
    data, geom = geom_point, color = NULL, fill = NULL, linetype = NULL, ..., palette = NULL, 
    fixed = list(), theme = list()
){
    color <- enquo(color)
    fill  <- enquo(fill)
    linetype <- enquo(linetype)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))

    p <- ggplot(data    = data,  # https://stackoverflow.com/a/55816211
                mapping = eval(expr(aes(color=!!color, fill=!!fill, linetype = !!linetype, !!!dots))))
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
#' @param linetype    svar (string)
#' @param facet       svar (character vector)
#' @param n           number
#' @param nrow        number of facet rows
#' @param ncol        number of facet cols
#' @param dir         'h' (horizontal) or 'v' (vertical)
#' @param scales      'free', 'fixed', 'free_y'
#' @param labeller    e.g. label_value
#' @param palette     named character vector
#' @param fixed       fixed aesthetics
#' @seealso \code{\link{plot_sample_violins}},
#'          \code{\link{plot_sample_boxplots}}
#' @return  ggplot object
#' @examples
#' # Data
#'     file <- download_data('atkin.metabolon.xlsx')
#'     object <- read_metabolon(file, plot = FALSE)
#'     object %<>% extract(, order(.$subgroup))
#'     
#' # Sample distributions
#'     plot_sample_densities(object)
#'     plot_sample_violins(  object, facet = 'Time')
#'     plot_sample_boxplots(object)
#'     plot_exprs(object)
#'     plot_exprs(object, dim = 'samples', x = 'subgroup', facet = 'Time')
#'     
#' # Feature distributions
#'     plot_feature_densities(object)
#'     plot_feature_violins(  object)
#'     plot_feature_boxplots( object)
#' @export
plot_densities <- function(
    object, assay = assayNames(object)[1], group, fill, color = NULL, linetype = NULL,
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
    value <- NULL
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
    groupsym    <- if (is.null(group))    quo(NULL) else sym(group)
    fillsym     <- if (is.null(fill ))    quo(NULL) else sym(fill)
    colorsym    <- if (is.null(color))    quo(NULL) else sym(color)
    linetypesym <- if (is.null(linetype)) quo(NULL) else sym(linetype)
    p <- plot_data(dt, geom = geom_density, x = value, fill = !!fillsym,
            color = !!colorsym, linetype = !!linetypesym, group = !!groupsym, palette = palette, fixed = fixed)
    if (!is.null(facet))  p <- p + facet_wrap(
            facet, nrow = nrow, ncol = ncol, dir = dir, labeller = labeller, 
            scales = scales)
    p
}

is_uniquely_empty <- function(x, y){
    ( is_empty(x) | !is_empty(y)) | (!is_empty(x) |  is_empty(y))
}

#' @rdname plot_densities
#' @export
plot_sample_densities <- function(
    object,
    assay    = assayNames(object)[1],
    group    = 'sample_id',
    fill     = if ('subgroup' %in% svars(object)) 'subgroup' else  'sample_id',
    color    = NULL, 
    linetype = NULL,
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
        linetype = linetype,
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
    assay    = assayNames(object)[1],
    fill     = 'feature_id',
    group    = fill,
    color    = NULL,
    linetype = NULL,
    n        = 9,
    facet    = NULL, nrow = NULL, ncol = NULL, dir = 'h', scales = 'free', 
    labeller = label_value, palette = NULL, 
    fixed    = list(alpha = 0.8, na.rm = TRUE)
){
    object %<>% extract_features_evenly(n)
    plot_densities( 
        object,
        assay    = assay,
        group    = group,
        fill     = fill,
        color    = color,
        linetype = linetype,
        facet    = facet,  nrow = nrow, ncol = ncol, dir = dir, scales = scales, 
        labeller = labeller, 
        palette  = palette, 
        fixed    = fixed ) +
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
#'          \code{\link{plot_densities}}
#' @examples
#' # data
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file)
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
    value <- NULL
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
    colorsym     <- if (is.null(color))  quo(NULL)  else  sym(color)
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
    object,
    assay = assayNames(object)[1],
    x     = 'sample_id', 
    fill  = if ('subgroup' %in% svars(object)) 'subgroup'  else  'sample_id', 
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
        object %<>% extract(rowSums(!is.na(values(object))) > 2, )
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
       object, assay = assay, x = x, fill = fill, color = color,
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
#' @param fit      string vector: subset of `fits(object)`
#' @param coefs    string vector: subset of `coefs(object)`
#' @param combiner '|' or '&'
#' @param verbose  TRUE or FALSE
#' @examples 
#' # Read
#'     file <- download_data('atkin.metabolon.xlsx')
#'     object <- read_metabolon(file)
#' # no limma
#'     object %<>% order_on_p()  # unchanged
#' # with limma
#'     object %<>% fit_limma(coefs = c('t1', 't2', 't3'))
#'     object %<>% order_on_p()
#' @return SummarizedExperiment
#' @export
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
#'     file <- download_data('atkin.metabolon.xlsx')
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
#' file <- download_data('atkin.metabolon.xlsx')
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
               pvar <- autonomics::pvar(     object, fit = fit, coefs = coefs[i])
             fdrvar <- autonomics::fdrvar(   object, fit = fit, coefs = coefs[i])
          effectvar <- autonomics::effectvar(object, fit = fit, coefs = coefs[i])
           facetvar <- paste0('facet.', coefs[[i]])
        assert_are_disjoint_sets(facetvar, fvars(object))
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


# Prepare
    xsym        <- sym(x)
    fillsym     <- if (is.null(fill))      quo(NULL) else  sym(fill)
    colorsym    <- if (is.null(color))     quo(NULL) else  sym(color)
    shapesym    <- if (is.null(shape))     quo(NULL) else  sym(shape)
    sizesym     <- if (is.null(size))      quo(NULL) else  sym(size)
    alphasym    <- if (is.null(alpha))     quo(NULL) else  sym(alpha)
    blocksym    <- if (is.null(block))     quo(NULL) else  sym(block)
    linetypesym <- if (is.null(linetype))  quo(NULL) else  sym(linetype)
    plotvars <- 'feature_name'
    if (!is.null(x))          plotvars %<>% c(x)         %>% unique()
    if (!is.null(fill))       plotvars %<>% c(fill)      %>% unique()
    if (!is.null(color))      plotvars %<>% c(color)     %>% unique()
    if (!is.null(shape))      plotvars %<>% c(shape)     %>% unique()
    if (!is.null(size))       plotvars %<>% c(size)      %>% unique()
    if (!is.null(alpha))      plotvars %<>% c(alpha)     %>% unique()
    if (!is.null(block))      plotvars %<>% c(block)     %>% unique()
    if (!is.null(linetype))   plotvars %<>% c(linetype)  %>% unique()
    if (!is.null(highlight))  plotvars %<>% c(highlight) %>% unique()
    if (!is.null(facet))      plotvars %<>% c(facet)     %>% unique()
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
    object, x = sample_id, fill = sample_id, color = NULL, highlight = NULL,
    fixed = list(na.rm=TRUE)
) plot_boxplots(
    object, x = !!enquo(x), fill=!!enquo(fill), color=!!color,
    highlight=!!enquo(highlight), fixed=fixed
)


feature_id <- NULL
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
    object, subgroup, x = !!enquo(subgroup), fill = !!enquo(subgroup), 
    color = NULL, highlight = NULL, facet = feature_id, fixed = list(na.rm=TRUE)
) plot_boxplots(
    object, x = !!enquo(x), fill = !!enquo(fill), color = !!color,
    facet = !!enquo(facet), highlight = !!enquo(highlight), fixed = fixed
)



#=============================================================================
#
#                 plot_feature_boxplots()
#
#=============================================================================

#' Plot features
#' @param object      SummarizedExperiment
#' @param geom        geom_point, geom_boxplot, etc.
#' @param subgroup    subgroup svar
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
#' plot_subgroup_boxplots(object, subgroup=Group)
#' plot_feature_profiles( object, subgroup=Group)
#' @export
plot_features <- function(
    object,
    geom,
    subgroup, 
    x     = !!enquo(subgroup),
    fill  = !!enquo(subgroup),
    color = !!enquo(subgroup),
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

