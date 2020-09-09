#' default color_var
#' @param object SummarizedExperiment
#' @return default value of color_var
#' @examples
#' # STEMCELLS
#'     file <- download_data('stemcells.proteinGroups.txt')
#'     object <- read_proteingroups(file,
#'                   invert_subgroups = c('E_EM','BM_E', 'BM_EM'), plot = FALSE)
#'     default_color_var(object)
#'
#' # GLUAMINASE
#'     file <- download_data('glutaminase.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     default_color_var(object)
#' @noRd
default_color_var <- function(object){
   if (     'block'    %in% svars(object))  'block'
   else if ('subgroup' %in% svars(object))  'subgroup'
   else                                     NULL
}


#' Default color values
#' @param object     SummarizedExperiment
#' @param color_var  string: svar mapped to color
#' @param show       logical
#' @return default color values vector
#' @examples
#' # STEMCELL RATIOS
#'     file <- download_data('stemcells.proteinGroups.txt')
#'     invert <- c('E_EM','BM_E', 'BM_EM')
#'     object <- read_proteingroups(
#'                  file, invert_subgroups = invert_subgroups, plot = FALSE)
#'     default_color_values(object, show = TRUE)
#'
#' # STEMCELL INTENSITIES
#'    file <- download_data('stemcells.proteinGroups.txt')
#'    object <- read_proteingroups(
#'                 file, quantity = 'Intensity labeled', plot = FALSE)
#'    default_color_values(object, show = TRUE)
#'
#' # GLUTAMINASE
#'     file <- download_data('glutaminase.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     default_color_values(object, show = TRUE)
#' @noRd
default_color_values <- function(
    object, color_var = default_color_var(object), show=FALSE, verbose=FALSE
){

# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(color_var, svars(object))

# Extract
    color_var_levels <- slevels(object, color_var)
    possible_separators  <- if (contains_ratios(object)){ c('.',' ')
                            } else {                      c('.',' ', '_') }
    sep <- guess_sep(object, color_var,
                    possible_separators = possible_separators, verbose=FALSE)

# Create color values
    color_values <- if (is.null(color_var_levels)){
                        make_gg_colors('default')
                    } else if (is.null(sep)){
                        make_gg_colors(color_var_levels, show = show,
                            verbose = verbose)
                    } else {
                       make_composite_colors(color_var_levels, sep = sep,
                            show = show, verbose = verbose)
                    }
# Return
    if (verbose) cmessage('\t\tMake composite colors')
    color_values
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
    color_levels <- hcl(h = hues, l = 65, c = 100)[1:n] %>%
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
    hues <- seq(15, 375, length = n1 + 1)[1:n1] %>% set_names(V1levels)

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


