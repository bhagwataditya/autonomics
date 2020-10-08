#=============================================================================
#
#                    default_colorscale
#
#==============================================================================

#' Default color values
#' @param object   SummarizedExperiment
#' @param color    symbol: svar mapped to color
#' @param show     TRUE or FALSE (default)
#' @param verbose  TRUE or FALSE (default)
#' @return default color values vector
#' @examples
#' # STEMCELL RATIOS
#'     file <- download_data('stemcells.proteinGroups.txt')
#'     invert_subgroups <- c('E_EM','BM_E', 'BM_EM')
#'     object <- read_proteingroups(
#'                 file, invert_subgroups = invert_subgroups, plot = FALSE)
#'     default_colorscale(object, show = TRUE)
#'
#' # STEMCELL INTENSITIES
#'    file <- download_data('stemcells.proteinGroups.txt')
#'    object <- read_proteingroups(
#'                 file, quantity = 'Intensity labeled', plot = FALSE)
#'    default_colorscale(object, show = TRUE)
#'
#' # GLUTAMINASE
#'    file <- download_data('glutaminase.metabolon.xlsx')
#'    object <- read_metabolon(file, plot = FALSE)
#'    default_colorscale(object, show = TRUE)
#' @noRd
default_colorscale <- function(
    object, color = subgroup, show = FALSE, verbose = FALSE
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    color_var <- as_string(ensym(color))
    assert_is_subset(color_var, svars(object))

# Extract
    color_var_levels <- slevels(object, color_var)
    possible_separators  <- if (contains_ratios(object)){ c('.',' ')
    } else {                      c('.',' ', '_') }
    sep <- guess_sep(object, color_var,
                possible_separators = possible_separators, verbose=FALSE)
    make_colorscale(color_var_levels, sep, show, verbose)
}


make_colorscale <- function(
    color_var_levels, sep = guess_sep(color_var_levels), show=FALSE,
    verbose = FALSE
){
    if (is.null(color_var_levels)){
        return(make_gg_colorscale('default', show = show, verbose = verbose))

    } else if (is.null(sep)){
        return(make_gg_colorscale(
                color_var_levels, show = show, verbose = verbose))

    } else {
        if (verbose) cmessage('\t\tMake composite colors')
        return(make_composite_colorscale(
                color_var_levels, sep = sep, show = show, verbose = verbose))
    }
}



#' Create default ggplot colors for factor levels
#' @param factor_levels  string vector
#' @param show           TRUE/FALSE
#' @param verbose        TRUE/FALSE
#' @return string vector: elements = colors, names = factor levels
#' @author John Colby
#' @references https://stackoverflow.com/questions/8197559
#' @noRd
make_gg_colorscale <- function(factor_levels, show, verbose = TRUE) {
    n <- length(factor_levels)
    hues <- seq(15, 375, length = n + 1)
    color_levels <- hcl(h = hues, l = 65, c = 100)[seq_len(n)] %>%
                    set_names(factor_levels)
    if (show) pie(rep(1, length(color_levels)), names(color_levels),
                    col = color_levels)
    if (verbose)  cmessage('\t\tMake default ggplot colors')
    color_levels
}


#' Make composite colors
#' @param svalues string vector
#' @param sep     string
#' @param show    TRUE/FALSE: show colors in pie plot?
#' @param verbose TRUE/FALSE
#' @return named string vector (elements = colors, names = color_var levels)
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' svalues <- subgroup_levels(object)
#' make_composite_colorscale(svalues, show = TRUE)
#' @noRd
make_composite_colorscale <- function(
    svalues, sep  = guess_sep(svalues), show = FALSE, verbose = TRUE
){
    # Assert
    assert_is_not_null(sep)

    # Satisfy CHECK
    subgroup <- V1 <- V2 <- color <- hue <- luminance <- NULL

    # Split into components
    #    * V1: first n-1 components => will be mapped to hue
    #    * V2: last component       => will be mapped to luminance
    # This approach works also when more than two components are present
    # It is therefore used instead of autonomics.import::split_values()
    V1  <-  stri_split_fixed(svalues, sep) %>%
            vapply( function(x) paste0(x[-length(x)], collapse = sep),
                    character(1))
    V2  <-  stri_split_fixed(svalues, sep) %>%
            vapply(function(x) x[length(x)], character(1))
    V1levels <- sort(unique(V1))
    V2levels <- sort(unique(V2))
    n1 <- length(V1levels)
    n2 <- length(V2levels)
    hues <- seq(15, 375, length = n1 + 1)[seq_len(n1)] %>% set_names(V1levels)

    color_levels <- character(0)
    for (i in seq_along(hues)){
        color_levels  %<>%  c(sequential_hcl(
                                n2, h = hues[[i]], power = 1, c = c(50, 100),
                                l = c(90, 30)) %>%
                            set_names(paste0(V1levels[[i]], sep, V2levels)))
    }
    if (show) pie(rep(1, length(color_levels)), names(color_levels),
                col = color_levels)

    return(color_levels)
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
#' @param colorscale  vector(names = svarlevels, values = colordefs)
#' @param fillscale   vector(names = svarlevels, values = colordefs)
#' @param ...         mapped aesthetics
#' @param fixed       fixed  aesthetics (list)
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% add_pca()
#' data <- sdata(object)
#' plot_data(data, x = pca1, y = pca2)
#' plot_data(data, x = pca1, y = pca2, color = TIME_POINT)
#' plot_data(data, x = pca1, y = pca2, color = NULL)
#'
#' fixed <- list(shape = 15, size = 3)
#' plot_data(data, x = pca1, y = pca2, fixed=fixed)
#'
#' guides <- list(color=FALSE, fill=FALSE)
#' plot_data(data, x = pca1, y = pca2, fixed=fixed, guides = guides)
#' @author Aditya Bhagwat, Johannes Graumann
#' @export
plot_data <- function(
    data,
    geom = geom_point,
    color      = subgroup,
    fill       = !!enquo(color),
    colorscale = make_colorscale(unique(eval_tidy(enquo(color), data))),
    fillscale  = make_colorscale(unique(eval_tidy(enquo(fill), data))),
    ...,
    fixed = list(),
    theme = list()
){
    color <- enquo(color)
    fill  <- enquo(fill)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))

    p <- ggplot(data    = data,  # https://stackoverflow.com/a/55816211
                mapping = eval(expr(aes(color=!!color, fill=!!fill, !!!dots))))
    p <- p + do.call(geom, fixed)
    p <- p + theme_bw()
    if (!rlang::quo_is_null(color)) p <- p + scale_color_manual(values = colorscale)
    if (!rlang::quo_is_null(fill )) p <- p + scale_fill_manual( values = fillscale)
    p <- p + do.call(ggplot2::theme, theme)

    p
}


#=============================================================================
#
#                       plot_sample_scores()
#                       plot_pca()
#                       plot_pls()
#                       plot_lda()
#                       plot_sma()
#
#==============================================================================


#' @rdname plot_sample_scores
#' @export
plot_pca <- function(
    object, xdim = 1, ydim = 2, color = subgroup,
    colorscale = default_colorscale(object, !!enquo(color)),
    feature_label = feature_name, ...,  fixed = list(shape=15, size=3),
    nloadings = 1
){
    plot_sample_scores(
        object, method = 'pca', xdim = xdim, ydim = ydim,
        color = !!enquo(color), colorscale = colorscale,
        feature_label = !!enquo(feature_label), ..., fixed = fixed,
        nloadings=nloadings)
}


#' @rdname plot_sample_scores
#' @export
plot_pls <- function(
    object, xdim = 1, ydim = 2, color = subgroup,
    colorscale = default_colorscale(object, !!enquo(color)),
    feature_label = feature_name, ...,  fixed = list(shape=15, size=3),
    nloadings = 1
){
    plot_sample_scores(
        object, method = 'pls', xdim = xdim, ydim = ydim,
        color = !!enquo(color), colorscale = colorscale,
        feature_label = !!enquo(feature_label), ..., fixed = fixed,
        nloadings=nloadings)
}


#' @rdname plot_sample_scores
#' @export
plot_lda <- function(
    object, xdim = 1, ydim = 2, color = subgroup,
    colorscale = default_colorscale(object, !!enquo(color)),
    feature_label = feature_name, ...,  fixed = list(shape=15, size=3),
    nloadings = 1
){
    plot_sample_scores(
        object, method = 'lda', xdim = xdim, ydim = ydim,
        color = !!enquo(color), colorscale = colorscale,
        feature_label = !!enquo(feature_label), ..., fixed = fixed,
        nloadings=nloadings)
}


#' @rdname plot_sample_scores
#' @export
plot_sma <- function(
    object, xdim = 1, ydim = 2, color = subgroup,
    colorscale = default_colorscale(object, !!enquo(color)),
    feature_label = feature_name, ...,  fixed = list(shape=15, size=3),
    nloadings = 1
){
    plot_sample_scores(
        object, method = 'sma', xdim = xdim, ydim = ydim,
        color = !!enquo(color), colorscale = colorscale,
        feature_label = !!enquo(feature_label), ..., fixed = fixed,
        nloadings=nloadings)
}


add_colorscale <- function(p, color, colorscale){
    if (!rlang::quo_is_null(enquo(color))){
        p <- p + scale_color_manual(values = colorscale)
    }
    p
}


#' Plot sample scores
#' @param object         SummarizedExperiment
#' @param method         'pca', 'pls', 'lda', or 'sma'
#' @param xdim           number (default 1)
#' @param ydim           number (default 2)
#' @param color          svar mapped to color (symbol)
#' @param colorscale     vector(names = svarlevels, values = colordefs)
#' @param ...            additional svars mapped to aesthetics
#' @param feature_label  fvar mapped to (loadings) label
#' @param fixed          fixed plot aesthetics
#' @param nloadings      number of loadings per half-axis to plot
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' plot_pca(object)
#' plot_pca(object, xdim=3, ydim=4)
#' plot_pca(object, nloadings = 0)
#' plot_pca(object, color = TIME_POINT)
#' plot_pca(object, color = TIME_POINT, xdim=3, ydim=4)
#' plot_pca(object, color = NULL)
#' @export
plot_sample_scores <- function(
    object, method='pca', xdim = 1, ydim = 2,
    color = subgroup, colorscale = default_colorscale(object, !!enquo(color)),
    feature_label = feature_name,
    ...,
    fixed = list(shape=15, size=3), nloadings = 1
){
    object %<>% add_projection(method, ndim=max(xdim, ydim), verbose = TRUE)
    color <- enquo(color)
    xstr <- paste0(method, xdim)
    ystr <- paste0(method, ydim)
    x     <- sym(xstr)
    y     <- sym(ystr)
    feature_label <- enquo(feature_label)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))

    xlab  <- paste0(xstr, ' : ', metadata(object)[[method]][[xstr]],'% ')
    ylab  <- paste0(ystr, ' : ', metadata(object)[[method]][[ystr]],'% ')

    p <- ggplot() + theme_bw() + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    p %<>% add_loadings(object, !!x, !!y, label = !!feature_label, nloadings = nloadings)
    p %<>% add_scores(object, !!x, !!y, color = !!color, !!!dots, fixed = fixed)
    p %<>% add_colorscale(!!color, colorscale)

    p
}

add_scores <- function(
    p, object, x = pca1, y = pca2, color = subgroup, ...,
    fixed = list(shape=15, size=3)
){
    x     <- enquo(x)
    y     <- enquo(y)
    color <- enquo(color)

    p + layer(  geom = 'point',
                mapping = aes(x = !!x, y = !!y, color = !!color, ...),
                stat    = "identity",
                data    = sdata(object),
                params  = fixed,
                position= 'identity')

}

utils::globalVariables('feature_name')

add_loadings <- function(
    p, object, x = pca1, y = pca2, label = feature_name, nloadings = 1
){
# Process args
    if (nloadings==0) return(p)
    x     <- enquo(x)
    y     <- enquo(y)
    label <- enquo(label)
    xstr <- rlang::as_name(x)
    ystr <- rlang::as_name(y)
# Loadings
    xloadings <- fdata(object)[[xstr]]
    yloadings <- fdata(object)[[ystr]]
    idx <- unique(c(headtail(order(xloadings, na.last=NA), nloadings),
                    headtail(order(yloadings, na.last=NA), nloadings)))
# Scale loadings to scoreplot
    xscores <- sdata(object)[[xstr]]
    yscores <- sdata(object)[[ystr]]
    maxscore <- min(abs(min(c(xscores, yscores, na.rm=TRUE))),
                    abs(max(c(xscores, yscores, na.rm=TRUE))), na.rm=TRUE)
    scorefactor <- maxscore/max(abs(c(xloadings, yloadings)),  na.rm=TRUE)

    plotdt <- fdata(object)
    plotdt[[xstr]] %<>% multiply_by(scorefactor)
    plotdt[[ystr]] %<>% multiply_by(scorefactor)
    plotdt %<>% extract(idx, )

# Plot
    feature_name <- NULL
    if (!'feature_name' %in% names(plotdt)){
        setnames(plotdt, 'feature_id', 'feature_name')}
    p + layer(  geom     = 'segment',
                mapping  = aes(x=0, y=0, xend=!!x, yend=!!y),
                stat     = "identity",
                data     = plotdt,
                params   = list(alpha = 0.1, size=1, na.rm = TRUE),#params   = list(alpha = 0.05, size=3),
                position = "identity") +
        layer(  geom     = "text",
                mapping  = aes(x = !!x, y = !!y, label = !!label),
                stat     = "identity",
                data     = plotdt,
                params   = list(alpha = 0.5, na.rm = TRUE),
                position ='identity')
}



headtail <- function(x, n){
    c(x[seq(1, n)], x[seq(length(x)+1-n, length(x))])
}


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
#' @param colorscale  named vector (names = varlevels, values = colors)
#' @param fillscale   named vector (names = varlevels, values = colors)
#' @param ...         mapped aesthetics
#' @param fixed       fixed aesthetics
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' plot_sample_densities(object)
#' plot_sample_densities(object, color = NULL)
#' @export
plot_sample_densities <- function(
    object,
    fill       = subgroup,
    color      = NULL,
    colorscale = default_colorscale(object, color = !!enquo(color)), # default_colorscale on NULL ?
    fillscale  = default_colorscale(object, color = !!enquo(fill)),
    ...,
    fixed = list(alpha = 0.5, na.rm = TRUE)
){
    dt <- sumexp_to_long_dt(object, svars = svars(object))
    fill <- enquo(fill)
    color <- enquo(color)
    value <- NULL
    plot_data(  dt,
                geom       = geom_density,
                x          = value,
                fill       = !!fill,
                color      = !!color,
                colorscale = colorscale,
                fillscale  = fillscale,
                ...,
                fixed      = fixed)
}


#' Plot sample violins
#' @param object      SummarizedExperiment
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param colorscale  named vector (names = varlevels, values = colors)
#' @param fillscale   named vector (names = varlevels, values = colors)
#' @param ...         mapped aesthetics
#' @param fixed       fixed aesthetics
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
    colorscale = default_colorscale(object, color = !!enquo(color)), # default_colorscale on NULL ?
    fillscale  = default_colorscale(object, color = !!enquo(fill)),
    ...,
    fixed = list(na.rm=TRUE)
){
    dt <- sumexp_to_long_dt(object, svars = svars(object))
    fill <- enquo(fill)
    color <- enquo(color)
    sample_id <- value <- NULL
    plot_data(  dt,
                geom       = geom_violin,
                x          = sample_id,
                y          = value,
                fill       = !!fill,
                color      = !!color,
                colorscale = colorscale,
                fillscale  = fillscale,
                ...,
                fixed      = fixed)
}



#' Plot sample boxplots
#' @param object      SummarizedExperiment
#' @param fill        svar mapped to fill
#' @param color       svar mapped to color
#' @param colorscale  named vector (names = varlevels, values = colors)
#' @param fillscale   named vector (names = varlevels, values = colors)
#' @param ...         mapped aesthetics
#' @param fixed       fixed aesthetics
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
    colorscale = default_colorscale(object, color = !!enquo(color)), # default_colorscale on NULL ?
    fillscale  = default_colorscale(object, color = !!enquo(fill)),
    ...,
    fixed = list(na.rm=TRUE)
){
    dt <- sumexp_to_long_dt(object, svars = svars(object))
    fill <- enquo(fill)
    color <- enquo(color)
    sample_id <- value <- NULL
    plot_data(  dt,
                geom       = geom_boxplot,
                x          = sample_id,
                y          = value,
                fill       = !!fill,
                color      = !!color,
                colorscale = colorscale,
                fillscale  = fillscale,
                ...,
                fixed      = fixed)
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
#' @param colorscale  named vector (names = varlevels, values = colors)
#' @param fillscale   named vector (names = varlevels, values = colors)
#' @param ...         mapped aesthetics
#' @param fixed       fixed aesthetics
#' @param theme       ggplot theme specifications
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% add_pca()
#' object %<>% extract(order(abs(fdata(.)$pca1, decreasing = TRUE)[1:9]), )
#' plot_feature_boxplots(object)
#' plot_feature_profiles(object)
#' @export
plot_features <- function(
    object,
    geom,
    x          = subgroup,
    fill       = subgroup,
    color      = subgroup,
    colorscale = default_colorscale(object, color = !!enquo(color)), # default_colorscale on NULL ?
    fillscale  = default_colorscale(object, color = !!enquo(fill)),
    ...,
    fixed = list(na.rm=TRUE),
    theme = list(axis.text.x  = element_blank(),#element_text(angle=90, vjust=0.5),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
){
    fill  <- enquo(fill)
    color <- enquo(color)
    x     <- enquo(x)
    dt <- sumexp_to_long_dt(object, svars = svars(object), fvars = fvars(object))
    value <- NULL
    p <- plot_data(  dt,
                geom       = geom,
                x          = !!x,
                y          = value,
                fill       = !!fill,
                color      = !!color,
                colorscale = colorscale,
                fillscale  = fillscale,
                ...,
                fixed      = fixed)
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


#==============================================================================
#
#                        plot_detects_per_subgroup
#
#==============================================================================

#' Plot detects per subgroup
#'
#' Plot number of detects, partial detects, and nondetects per subgroup
#'
#' @param object SummarizedExperiment
#' @param group  svar (symbol)
#' @param colorscale vector(names = svarlevels, values = colordefs)
#' @return ggplot object
#' @examples
#' # STEMCELLS
#'     file <- download_data('stemcells.proteinGroups.txt')
#'     invert_subgroups <- c('E_EM', 'E_BM', 'EM_BM')
#'     object <- read_proteingroups(file, invert_subgroups = invert_subgroups,
#'                                   impute = FALSE, plot = FALSE)
#'     plot_detects_per_subgroup(object)
#' @export
plot_detects_per_subgroup <- function(
    object, group = subgroup,
    colorscale = default_colorscale(object, !!ensym(group))
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    group_var <- as_string(ensym(group))
    assert_is_subset(group_var, svars(object))
    assert_is_subset(slevels(object, group_var), names(colorscale))
    variable <- subgroup <- value <- . <- NULL
# Prepare datatable
    split_objects  <- split_by_svar(object, group_var)
    imps     <- vapply(split_objects, count_imputes,    numeric(1))
    nons     <- vapply(split_objects, count_nondetects, numeric(1))
    partials <- vapply(split_objects, count_pardetects, numeric(1)) - nons
    fulls    <- nrow(object) - partials - imps - nons
    plot_dt  <- data.table(subgroup = factor(names(split_objects),
                                            slevels(object, group_var)))
    if (any(fulls   !=0))            plot_dt %<>% extract(, fulls   := fulls)
    if (any(partials!=0))            plot_dt %<>% extract(, partials:= partials)
    if (any(nons!=0) | any(imps!=0)) plot_dt %<>% extract(, nons := nons + imps)
    plot_dt %<>% data.table::melt(id.vars = 'subgroup')
# Set order of variable levels
    variable_levels <- c('nons', 'partials', 'fulls') %>%
                        extract(.%in% plot_dt$variable)
    plot_dt$variable %<>% factor(variable_levels)
# Plot
    title <- 'fulldetects  |  partialdetects'
    if ( all(imps==0) & !all(nons==0))  title %<>% paste0('  |  nondetects')
    if (!all(imps==0) &  all(nons==0))  title %<>% paste0('  |  imputes')
    plot_dt$subgroup %<>% factor(rev(levels(.)))
    ggplot(plot_dt, aes(x = subgroup, y = value, fill = subgroup,
                        group = variable)) +
    ggtitle(title) + theme_bw() +
    geom_col(color = 'black', position = position_stack()) +
    scale_fill_manual(values = colorscale) +
    geom_text(aes(label=value, x = subgroup),
                    position = position_stack(vjust=0.5),
                    size = rel(3)) +
    theme(#axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust=1),
                panel.grid.major.x = element_blank(),
                panel.border       = element_blank(),
                plot.title         = element_text(hjust = 0.5)) +
    xlab(NULL) +
    ylab(NULL) +
    coord_flip()
}

count_imputes <- function(object){
    sum(rowAlls(is_imputed(object)))
}

count_nondetects <- function(object){
    sum(rowAlls(is.na(exprs(object))))
}

count_pardetects <- function(object){
    sum(rowAnys(is.na(exprs(object))))
}

