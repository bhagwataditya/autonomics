#=============================================================================
#
#                    default_color_values
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
#'     default_color_values(object, show = TRUE)
#'
#' # STEMCELL INTENSITIES
#'    file <- download_data('stemcells.proteinGroups.txt')
#'    object <- read_proteingroups(
#'                 file, quantity = 'Intensity labeled', plot = FALSE)
#'    default_color_values(object, show = TRUE)
#'
#' # GLUTAMINASE
#'    file <- download_data('glutaminase.metabolon.xlsx')
#'    object <- read_metabolon(file, plot = FALSE)
#'    default_color_values(object, show = TRUE)
#' @noRd
default_color_values <- function(
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
    make_colors(color_var_levels, sep, show, verbose)
}


make_colors <- function(
    color_var_levels, sep = guess_sep(color_var_levels), show=FALSE,
    verbose = FALSE
){
    if (is.null(color_var_levels)){
        return(make_gg_colors('default'))

    } else if (is.null(sep)){
        return(make_gg_colors(color_var_levels, show = show, verbose = verbose))

    } else {
        if (verbose) cmessage('\t\tMake composite colors')
        return(make_composite_colors(
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
make_gg_colors <- function(factor_levels, show, verbose = TRUE) {
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
#' make_composite_colors(svalues, show = TRUE)
#' @noRd
make_composite_colors <- function(
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
#' @param data          data.frame'
#' @param geom          geom_point, etc.
#' @param color         variable mapped to color (symbol)
#' @param color_values  vector(names = svarlevels, values = colordefs)
#' @param ...           mapped aesthetics
#' @param fixed         fixed  aesthetics (list)
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% add_pca(plot = FALSE)
#' data <- sdata(object)
#' plot_data(data, x = pca1, y = pca2)
#' plot_data(data, x = pca1, y = pca2, fixed = list(size=3))
#' plot_data(data, x = pca1, y = pca2, color = TIME_POINT, fixed=list(size=3))
#' @author Aditya Bhagwat, Johannes Graumann
#' @export
plot_data <- function(
    data, geom = geom_point, color = subgroup,
    color_values = make_colors(eval_tidy(color, data)), ..., fixed = list()
){
    color <- enquo(color)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))
    # https://stackoverflow.com/a/55816211
    p <- ggplot(data = data, mapping = eval(expr(aes(color=!!color, !!!dots))))
    p <- p + do.call(geom, fixed)
    p <- p + theme_bw()
    if (!is.null(color_values)){
        p <- p + scale_color_manual(values = color_values) }
    p
}

#=============================================================================
#
#     plot_sample_scores()
#
#==============================================================================


#' Plot sample scores
#' @param object  SummarizedExperiment
#' @param method  string: 'pca', 'pls', 'lda', 'sma'
#' @param xdim    number (default 1): x axis dimension
#' @param ydim    number (default 2): y axis dimension
#' @param color   svar mapped to color (symbol)
#' @param color_values vector(names = svarlevels, values = colordefs)
#' @param ...     additional svars mapped to aesthetics
#' @param fixed   fixed plot aesthetics
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% add_pca(plot = FALSE, ndim = 4)  # Principal Component Analysis
#' plot_sample_scores(object, 'pca')
#' plot_sample_scores(object, 'pca', color = TIME_POINT)
#' plot_sample_scores(object, 'pca', color = TIME_POINT, xdim=3, ydim=4)
#' @export
plot_sample_scores <- function(
    object, method, xdim = 1, ydim = 2, color = subgroup,
    color_values = default_color_values(object, !!ensym(color)), ...,
    fixed = list(shape=15, size=3)
){
    x <- paste0(method, xdim)
    y <- paste0(method, ydim)
    xlab  <- paste0(x, ' : ', metadata(object)[[method]][[x]], '% ')
    ylab  <- paste0(y, ' : ', metadata(object)[[method]][[y]], '% ')

    p <- plot_data( sdata(object), x = !!sym(x), y = !!sym(y),
                    color = !!ensym(color), color_values = color_values, ...,
                    fixed = fixed)
    p <- p + ggplot2::xlab(xlab)
    p <- p + ggplot2::ylab(ylab)
    p
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
#' @param color_values vector(names = svarlevels, values = colordefs)
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
    color_values = default_color_values(object, !!ensym(group))
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    group_var <- as_string(ensym(group))
    assert_is_subset(group_var, svars(object))
    assert_is_subset(slevels(object, group_var), names(color_values))
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
    scale_fill_manual(values = color_values) +
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

