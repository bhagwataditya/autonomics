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
#'     file <- download_data('stemcells.proteinGroups.txt')
#'     invert_subgroups <- c('E_EM','BM_E', 'BM_EM')
#'     object <- read_proteingroups(
#'                 file, invert_subgroups = invert_subgroups, plot = FALSE)
#'     p <- plot_sample_densities(object)
#'     add_color_scale(p, data=sdata(object))
#'
#' # STEMCELL INTENSITIES
#'    file <- download_data('stemcells.proteinGroups.txt')
#'    object <- read_proteingroups(
#'                 file, quantity = 'Intensity labeled', plot = FALSE)
#'    add_color_scale(object)
#'
#' # GLUTAMINASE
#'    require(magrittr)
#'    file <- download_data('glutaminase.metabolon.xlsx')
#'    object <- read_metabolon(file, plot = FALSE)
#'    object %<>% pca()
#'    plot_data(sdata(object), x=pca1, y=pca2, color=TIME_POINT)
#' @noRd
add_color_scale <- function(p, color = subgroup, data){
# Assert
    assert_is_data.frame(data)
    color <- enquo(color)
# Colors
    if (!rlang::quo_is_null(color)){
        color_var <- as_string(ensym(color))
        assert_is_subset(color_var, names(data))
        values0 <- data[[color_var]]
        if (!is.numeric(values0)){
            levels0 <- unique(values0)
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
        fill_var <- as_string(ensym(fill))
        assert_is_subset(fill_var, names(data))
        values0 <- data[[fill_var]]
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
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' varlevels <- subgroup_levels(object)
#' make_twofactor_colors(varlevels, show = TRUE)
#' @noRd
make_twofactor_colors <- function(
    varlevels, sep  = guess_sep(varlevels), show = FALSE, verbose = TRUE
){
    # Assert
    assertive::assert_has_no_duplicates(varlevels)
    assert_is_not_null(sep)
    if (verbose) cmessage('\t\tMake composite colors')

    # Satisfy CHECK
    subgroup <- V1 <- V2 <- color <- hue <- luminance <- NULL

    # Split into components
    #    * V1: first n-1 components => will be mapped to hue
    #    * V2: last component       => will be mapped to luminance
    # This approach works also when more than two components are present
    # It is therefore used instead of autonomics.import::split_values()
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
#' file <- download_data('glutaminase.metabolon.xlsx')
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
#' @return ggplot object
#' @examples
#' # STEMCELLS
#'     file <- download_data('stemcells.proteinGroups.txt')
#'     invert_subgroups <- c('E_EM', 'E_BM', 'EM_BM')
#'     object <- read_proteingroups(file, invert_subgroups = invert_subgroups,
#'                                   impute = FALSE, plot = FALSE)
#'     plot_detects_per_subgroup(object)
#' @export
plot_detects_per_subgroup <- function(object, group = subgroup){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    group_var <- as_string(ensym(group))
    assert_is_subset(group_var, svars(object))
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

    plot_data(plot_dt, geom_col,
              x = subgroup, y = value, fill = subgroup, group = variable,
              fixed = list(position = position_stack(), color="black")) +
    #ggplot(plot_dt, aes(x = subgroup, y = value, fill = subgroup,
    #                    group = variable)) +
    ggtitle(title) + theme_bw() +
    #geom_col(color = 'black', position = position_stack()) +
    #scale_fill_manual(values = colorscale) +
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

