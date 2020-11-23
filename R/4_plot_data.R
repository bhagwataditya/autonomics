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
#'                 file, invert_subgroups = invert_subgroups, plot = FALSE)
#'     p <- plot_sample_densities(object)
#'     add_color_scale(p, data=sdata(object))
#'
#' # STEMCELL INTENSITIES
#'    file <- download_data('billing16.proteingroups.txt')
#'    object <- read_proteingroups(
#'                 file, quantity = 'Intensity labeled', plot = FALSE)
#'    add_color_scale(object)
#'
#' # GLUTAMINASE
#'    require(magrittr)
#'    file <- download_data('halama18.metabolon.xlsx')
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
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
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



