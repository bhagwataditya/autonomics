#==============================================================================
#
#          plot_volcano
#              make_volcano_dt
#
#==============================================================================

default_coefs <- function(featuredt, fit = fits(featuredt)){
    if (length(fit)==0) return(NULL)    # none
    y <- autonomics::coefs(featuredt, fit = fit)   # intercept
    if (length(y)==1)   return(y)
    y %<>% setdiff('Intercept')                 # intercept + others
    y
}

#' Bin continuous variable
#' @param object numeric or SummarizedExperiment
#' @param fvar   string or NULL
#' @param probs  numeric
#' @param ... (S3 dispatch)
#' @return  factor vector
#' @examples 
#' # Numeric vector
#'     object <- rnorm(10, 5, 1)
#'     bin(object)
#' # SummarizedExperiment
#'     file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'     fdt(object <- read_maxquant_proteingroups(file))
#'     fdt(bin(object, 'pepcounts'))
#' @export
bin <- function(object, ...) UseMethod('bin')


#' @rdname bin
#' @export
bin.logical <- function(object, ...)  object

#' @rdname bin
#' @export
bin.character <- function(object, ...) object

#' @rdname bin
#' @export
bin.factor <- function(object, ...) object

#' @rdname bin
#' @export
bin.numeric <- function(object, probs = c(0, 0.33, 0.66, 1), ...){
    breaks <- quantile(object, probs = probs)
    breaks[1] %<>% subtract(1e-7)      # avoid smallest number from falling outside of bin
    object %<>% cut(breaks)                 # explicit breaks avoid negative bin
    levels(object) %<>% substr(2, nchar(.)) #    https://stackoverflow.com/questions/47189232
    levels(object) %<>% split_extract_fixed(',', 1)
    levels(object) %<>% paste0('>', .)
    object
}

#' @rdname bin
#' @export
bin.SummarizedExperiment <- function(
    object, fvar, probs = c(0, 0.33, 0.66, 1), ...
){
    if (is.null(fvar))  return(object)
    fdt(object)[[fvar]] %<>% bin()
    object
}


#' Add assay means
#' @param object SummarizedExperiment or NULL
#' @param assay  string
#' @param bin    TRUE or FALSE
#' @return SummarizedExperiment
#' @examples 
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' fdt(object) %<>% extract(, 1:2)
#' fdt(object)
#' object %<>% add_assay_means(SummarizedExperiment::assayNames(.))
#' fdt(object)
#' @export
add_assay_means <- function(
    object, 
    assay = assayNames(object)[1], 
    bin = TRUE
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(assay, assayNames(object))
# Add
    if (is.null(assay))  return(object)
    for (.assay in assay){
        if (is.numeric(assays(object)[[.assay]])){
            fdt(object)[[.assay]] <- rowMeans(assays(object)[[.assay]], na.rm = TRUE)
        }
    }
# Return
    object
}



#' Add adjusted pvalues
#' 
#' @param featuredt  fdt(object)
#' @param method    'fdr', 'bonferroni', ... (see `p.adjust.methods`)
#' @param fit       'limma', 'lm', 'lme', 'lmer'
#' @param coefs      coefficient (string)
#' @param verbose    TRUE or FALSE
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' fdt(object) %<>% extract(, 1:2)
#' object %<>% fit_limma(coef = 'Adult')
#' object %<>% extract(order(fdt(.)$`p~Adult~limma`), )
#'  fdt(object)
#' (fdt(object) %<>% add_adjusted_pvalues('fdr'))
#' (fdt(object) %<>% add_adjusted_pvalues('fdr'))      # smart enough not to add second column
#' (fdt(object) %>% add_adjusted_pvalues('bonferroni'))
#' @return SummarizedExperiment
#' @export
add_adjusted_pvalues <- function(
    featuredt, 
       method,
          fit = fits(featuredt),
        coefs = default_coefs(featuredt, fit = fit),
      verbose = TRUE
){
# Assert
    assert_is_data.table(featuredt)
    assert_is_subset(method, stats::p.adjust.methods)
# Reset
    fdrvars <- paste(method, coefs, fit, sep = FITSEP)
    fdrvars %<>% intersect(names(featuredt))
    featuredt[ , (fdrvars) := NULL ]
# Compute
    sep <- guess_fitsep(featuredt)
    adjdt <- pdt(featuredt, fit = fit, coef = coefs)
    adjdt <- adjdt[, lapply(.SD, p.adjust, method = method), .SDcols = names(adjdt)[-1] ]
    names(adjdt) %<>% stri_replace_first_regex(sprintf('^p%s', sep), sprintf('%s%s', method, sep))
    featuredt %<>% cbind(adjdt)
    featuredt
}


    
#' Create volcano datatable
#' @param object  SummarizedExperiment
#' @param fit    'limma', 'lme', 'lm', 'wilcoxon'
#' @param coefs   character vector: coefs for which to plot volcanoes
#' @param shape   fvar or NULL
#' @param size    fvar or NULL
#' @param alpha   fvar or NULL
#' @param label   fvar or NULL
#' @return data.table
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file, impute = TRUE, fit = 'limma')
#' make_volcano_dt(object, fit = 'limma', coefs = 'Adult')
#' @export
make_volcano_dt <- function(
    object,
       fit = fits(fdt(object))[1], 
     coefs = default_coefs(fdt(object), fit = fit)[1],
     shape = 'imputed', 
     size  = NULL, 
     alpha = NULL,
     label = 'feature_id'
){
# Assert    
    assert_is_all_of(object, "SummarizedExperiment")
    sep <- guess_fitsep(fdt(object))
    assert_any_are_matching_regex(fvars(object), paste0('^p', sep))
    assert_is_subset(fit, fits(fdt(object)))
    assert_is_subset(coefs, autonomics::coefs(fdt(object), fit = fit))
    if (!is.null(shape)){ assert_is_subset(shape, fvars(object)); object %<>% bin(shape) }
    if (!is.null(size) ){ assert_is_subset(size,  fvars(object)); object %<>% bin(size)  }
    if (!is.null(alpha)){ assert_is_subset(alpha, fvars(object)); object %<>% bin(alpha) }
    if (!is.null(label))  assert_is_subset(label, fvars(object))
    fdt(object) %<>% add_adjusted_pvalues('bonferroni', fit = fit, coefs = coefs)
# Prepare
    bon <- direction <- effect <- fdr <- mlp <- p <- significance <- NULL
    idvars <- 'feature_id'
    if ('control' %in% fvars(object))  idvars %<>% c('control')
    if (!is.null(label))  idvars %<>% union(label)
    if (!is.null(shape))  idvars %<>% union(shape)
    if (!is.null(size))   idvars %<>% union(size)
    if (!is.null(alpha))  idvars %<>% union(alpha)
    valuevars  <-  effectvar(fdt(object), coef = coefs, fit = fit)  # elminate similar function pvars etc.
    valuevars %<>%  c(  pvar(fdt(object), coef = coefs, fit = fit))

    dt <- fdt(object)[, c(idvars, valuevars), with = FALSE]
    dt %<>% melt.data.table(id.vars = idvars)
    dt %<>% tidyr::separate_wider_delim(cols = 'variable', names = c('quantity', 'coef', 'fit'), delim = sep)
    idvars %<>% c('coef', 'fit')
    #dt %<>% dcast.data.table(feature_id+feature_name+coef+fit ~ quantity, value.var = 'value')
    dt %<>% tidyr::pivot_wider(id_cols = tidyselect::all_of(idvars), names_from = 'quantity', values_from = 'value')
    dt %<>% data.table()
    dt$coef %<>% factor(coefs)
    dt[, fdr := p.adjust(p, method = 'fdr'),        by = c('fit', 'coef')]
    dt[, bon := p.adjust(p, method = 'bonferroni'), by = c('fit', 'coef')]
    dt[, mlp := -log10(p)]
    dt %<>% extract(!is.na(effect) & !is.na(p))
    dt[, direction := 'unchanged']
    dt[p>0.05, direction := 'unchanged']
    dt[p<=0.05 & effect>0, direction := 'up']
    dt[p<=0.05 & effect<0, direction := 'down']
    dt[]
}


#' Plot volcano
#' @param object         SummarizedExperiment
#' @param fit           'limma', 'lme', 'lm', 'wilcoxon'
#' @param coefs          character vector
#' @param facet          character vector
#' @param shape          fvar  (string)
#' @param size           fvar  (string)
#' @param alpha          fvar  (string)
#' @param label          fvar  (string)
#' @param max.overlaps   number: passed to ggrepel
#' @param features       feature ids (character vector): features to encircle 
#' @param nrow           number: no of rows in plot
#' @param p              number: p cutoff for labeling
#' @param fdr            number: fdr cutoff for labeling
#' @param xndown         x position of ndown labels
#' @param xnup           x position of nup labels
#' @param title          string or NULL
#' @return ggplot object
#' @examples
#' # Regular Usage
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object <- read_metabolon(file)
#'     object %<>% fit_limma()
#'     object %<>% fit_lm()
#'     plot_volcano(object, coefs = 't3', fit = 'limma')                   # single contrast
#'     plot_volcano(object, coefs = c('t2', 't3'), fit = 'limma')          # multiple contrasts
#'     plot_volcano(object, coefs = c('t2', 't3'), fit = c('limma', 'lm')) # multiple contrasts and methods
#' 
#' # When nothing passes FDR
#'     fdt(object) %<>% add_adjusted_pvalues('fdr', fit = 'limma',coefs = 't3')
#'     object %<>% extract( fdrvec(fdt(object), fit = 'limma', coef = 't3') > 0.05, )
#'     plot_volcano(object, coefs = 't3', fit = 'limma')
#' 
#' # Additional mappings
#'     file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'     object <- read_maxquant_proteingroups(file, impute = TRUE)
#'     object %<>% fit_limma()
#'     plot_volcano(object)
#'     plot_volcano(object, label = 'gene')
#'     plot_volcano(object, label = 'gene', size = 'log2maxlfq')
#'     plot_volcano(object, label = 'gene', size = 'log2maxlfq', alpha = 'pepcounts')
#'     plot_volcano(object, label = 'gene', features = c('hmbsb'))
#' @export
plot_volcano <- function(
          object,
             fit = fits(fdt(object))[1], 
           coefs = default_coefs(fdt(object), fit = fit)[1],
           facet = if (is_scalar(fit)) 'coef' else c('fit', 'coef'),
           shape = if ('imputed' %in% fvars(object)) 'imputed' else NULL, 
            size = NULL,
           alpha = NULL,
           label = 'feature_id', #if ('gene' %in% fvars(object)) 'gene' else 'feature_id', 
    max.overlaps = 10,
        features = NULL,
            nrow = length(fit),
               p = 0.05, 
             fdr = 0.05,
          xndown = NULL,
            xnup = NULL,
           title = NULL
){
# Assert
    assert_is_valid_sumexp(object)
    if (nrow(object)==0){  cmessage('nrow(object)=0 : return NULL');   return(NULL) }
    assert_is_a_number(nrow)
    bon <- effect <- direction <- mlp <- ndown <- nup <- significance <- NULL
    singlefeature <- yintercept <- NULL
    facet %<>% lapply(sym)
    facet <- vars(!!!facet)
# Volcano 
    plotdt <- make_volcano_dt(object, fit = fit, coefs = coefs, 
                  label = label, shape = shape, size = size, alpha = alpha)
    g <- ggplot(plotdt) + facet_wrap(facet, nrow = nrow)
    g <- g + theme_bw() + theme(panel.grid = element_blank())
    g <- g + xlab('log2(FC)') + ylab('-log10(p)') + ggtitle(title)
    shapesym <- if (is.null(shape))  quo(NULL) else sym(shape)
    sizesym  <- if (is.null(size))   quo(NULL) else sym(size)
    alphasym <- if (is.null(alpha))  quo(NULL) else sym(alpha)
    colorvalues <- c(down = '#ff5050', unchanged = 'grey', up = '#009933')
    g <- g + geom_point(data = plotdt, 
                        mapping = aes(x = effect, y = mlp, color = direction, 
                                      shape = !!shapesym, alpha = !!alphasym, 
                                      size = !!sizesym), 
                        na.rm = TRUE) + 
             scale_color_manual(values = colorvalues)
    if (!is.null(size))   g <- g + scale_size_manual( values = 1:3)
    if (!is.null(alpha))  g <- g + scale_alpha_manual(values = c(0.3, 0.5, 1))
# Significance lines
    if (is.null(xndown))  xndown <- min(plotdt$effect)
    if (is.null(xnup  ))  xnup   <- max(plotdt$effect)
    dy <- 0.03*(max(plotdt$mlp) - min(plotdt$mlp))
    minn <- function(x) if (length(x)==0) return(Inf) else min(x, na.rm = TRUE) 
        # Avoid warning 'no non-missing arguments to min; returning Inf'
    summarydt <- plotdt[ ,
        .( 
            significance = c('p = 0.05',   'fdr = 0.05',             'bon = 0.05'), 
            yintercept   = c( -log10(0.05), -log10(fdr2p(c(0.05, fdr))[1]),   -log10(0.05/(length(feature_id)+1))), # +1 because dummy feature is added
            ndown        = c( sum(p <0.05 & effect < 0),  sum(fdr <0.05 & effect < 0),  sum(bon <0.05 & effect < 0)),
            nup          = c( sum(p <0.05 & effect > 0),  sum(fdr <0.05 & effect > 0),  sum(bon <0.05 & effect > 0))
        ),
        by = c('fit', 'coef')]
    g <- g + geom_hline(data = summarydt, mapping = aes(yintercept = yintercept, linetype = significance), color = 'gray30')
    g <- g + ggplot2::scale_linetype_manual(values = c(`p = 0.05` = 3, `fdr = 0.05` = 2, `bon = 0.05` = 1))
    do_geom_label <- function(mapping, label.size){  
                        geom_label(data       = summarydt, 
                                   mapping    = mapping, 
                                   label.size = label.size, 
                                   hjust      = 0.5, 
                                   fill       = alpha(c('white'), 0.6), 
                                   color      = 'gray30')  }
  # g <- g + do_geom_label(mapping = aes(x = xsignificance, y = yintercept,      label = significance), label.size = 0.5)
    g <- g + do_geom_label(mapping = aes(x = xndown,        y = yintercept + dy, label = ndown       ), label.size = NA)
    g <- g + do_geom_label(mapping = aes(x = xnup,          y = yintercept + dy, label = nup         ), label.size = NA)
    g <- g + guides(color = 'none')
# Labels
    if (!is.null(label)){
        idx <- plotdt$fdr < fdr & plotdt$p < p
        labeldt <- plotdt[idx]
        if (!is.null(features)){
            seldt <- copy(plotdt)
            seldt[, singlefeature := feature_id]
            seldt %<>% uncollapse(singlefeature, sep = ';')
            seldt %<>% extract(singlefeature %in% features)
            seldt[, singlefeature := NULL]
            seldt %<>% unique()
        }
        g <- g + geom_label_repel(     data = labeldt, 
                                    mapping = aes(x = effect, y = mlp, label = !!sym(label), color = direction), 
                                    # color = 'black', 
                                 label.size = NA, 
                                       fill = alpha(c('white'), 0.6),
                                      na.rm = TRUE, 
                                show.legend = FALSE,
                               max.overlaps = max.overlaps )
        if (!is.null(features)){
            g <- g + geom_label_repel( data = seldt, 
                                    mapping = aes(x = effect, y = mlp, label = !!sym(label)), 
                                      color = 'black', 
                                 label.size = NA, 
                                       fill = alpha(c('white'), 0.6),
                                      na.rm = TRUE,
                                show.legend = FALSE,
                               max.overlaps = max.overlaps )
            
            g <- g + geom_point(       data = seldt, 
                                    mapping = aes(x = effect, y = mlp), 
                                      shape = 1,
                                       size = 4,
                                      color = 'black')
        }
    }
    g
}


#' Detect fixed patterns in collapsed strings
#' @param x         vector with collapsed strings
#' @param patterns  vector with fixed patterns (strings)
#' @param sep       collapse separator (string) or NULL (if uncollapsed)
#' @return boolean vector
#' @examples 
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' x <- fdt(object)$uniprot
#' patterns <- c('A0A0R4IKT8', 'Q7T3G6')
#' table(stri_detect_fixed_in_collapsed(x = x, patterns = patterns, sep = ';'))
#' @export
stri_detect_fixed_in_collapsed <- function(x, patterns, sep){
# Assert
    assert_is_character(x)
    assert_is_character(patterns)
    assert_is_a_string(sep)
    detect <- x0 <- NULL
# Detect
    dt <- data.table( x0 = x, x = x)
    if (!is.null(sep))  dt %<>% uncollapse(x, sep = sep)
    dt[, detect := x %in% patterns]
    dt %<>% extract(, .(detect = any(detect)), by = x0)
    dt <- dt[x, on = 'x0']
# Return
    dt$detect
}


#' Map fvalues
#' @param object   SummarizedExperiment
#' @param fvalues  uncollapsed string vector
#' @param from     string (fvar)
#' @param to       string (svar)
#' @param sep      collapse separator
#' @return string vector
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' fdt(object)
#' map_fvalues(object, c('Q6DHL5', 'Q6PFS7'), from = 'uniprot', to = 'feature_id', sep = ';')
#' @export
map_fvalues <- function(
    object, 
    fvalues, 
    from    = 'uniprot', 
    to      = 'feature_id', 
    sep     = ';'){
# Assert 
    assert_is_valid_sumexp(object)
    assert_is_character(fvalues)
    assert_scalar_subset(from, fvars(object))
    assert_scalar_subset(to,   fvars(object))
# Map
    x <- fdt(object)[[from]]
    idx <- stri_detect_fixed_in_collapsed(x, fvalues, sep )
    y <- fdt(object)[[to ]][idx]
# Return
    y
}


#' fdr to p
#' @param fdr fdr values
#' @examples
#' # Read/Fit
#'    file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'    object <- read_metabolon(file)
#'    object %<>% fit_limma()
#'    pcol <- pvar(fdt(object), fit = 'limma', coef = 't3')
#'    object %<>% extract(order(fdt(.)[[pcol]]), )
#'    object %<>% extract(1:10, )
#'    fdt(object) %<>% extract(, 1)
#'    object %<>% fit_limma(coefs = 't3')
#' # fdr2p
#'    fdt(object)[[pcol]]
#'    fdt(object)[[pcol]] %>% p.adjust(method = 'fdr')
#'    fdt(object)[[pcol]] %>% p.adjust(method = 'fdr') %>% fdr2p()
#' @export
fdr2p <- function(fdr){
    idx <- order(fdr)
    p <- fdr
    p[idx] <- fdr[idx]*seq_along(fdr[idx])/length(fdr[idx])
    p
}

