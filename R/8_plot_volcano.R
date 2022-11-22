#==============================================================================
#
#          plot_volcano
#              make_volcano_dt
#                  top_down/top_up
#                      nmax/nmin
#
#==============================================================================

invwhich <- function(indices, totlength) is.element(seq_len(totlength), indices)

#' Return nth max (min) value in vector
#'
#' Orders a vector and returns n'th ordered value.
#' When vector length is smaller than n, returns last value.
#'
#' @param x numeric vector
#' @param n integer
#' @return value
#' @examples
#' nmax(c(1,2,3,4,5), 2)
#' nmin(c(1,2,3,4,5), 2)
#' @noRd
nmax <- function(x, n){
    . <- NULL
    sort(x, decreasing = TRUE) %>% extract(min(length(.), n))
}
nmin <- function(x, n){
    . <- NULL
    sort(x) %>% extract(min(length(.), n))
}

top_down <- function(effect, fdr, mlp, ntop){
    fdr_ok   <- fdr  < 0.05
    coef_ok  <- effect < -1
    coef_top <- if (any(fdr_ok)) {  effect < nmin(effect[fdr_ok], ntop + 1)
                } else {           rep(FALSE, length(effect))            }
    mlp_top  <- if (any(coef_ok)) { mlp  > nmax(mlp[coef_ok], ntop + 1)
                } else {           rep(FALSE, length(effect))            }
    fdr_ok & coef_ok & (coef_top | mlp_top)
}

#' @examples
#' file <- download_data("billing16.proteingroups.txt")
#' invert_subgroups <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'           file, invert_subgroups=invert_subgroups, fit='limma', plot=FALSE)
#' effect <-      limma(object)[,1,'effect']
#' fdr    <-      limma(object)[,1,'fdr']
#' mlp    <- -log(limma(object)[,1,'p'])
#' ntop   <- 3
#' table(top_up(effect, fdr, mlp, ntop))
#' table(top_down(effect, fdr, mlp, ntop))
#' @noRd
top_up <- function(effect, fdr, mlp, ntop){
    fdr_ok   <- fdr  < 0.05
    coef_ok  <- effect >  0 # currently no filter
    coef_top <- if (any(fdr_ok)) {  effect > nmax(effect[fdr_ok], ntop + 1)
                } else {          rep(FALSE, length(effect)) }
    mlp_top  <- if (any(coef_ok)) { mlp > nmax(mlp[coef_ok], ntop + 1)
                } else {           rep(FALSE, length(effect)) }
    fdr_ok & coef_ok & (coef_top | mlp_top)
}


melt_contrastdefs <- function(contrastdefmat){
    facetrow <- NULL
    contrastdefdt <- data.table(contrastdefmat, facetrow = "")
    if (!is.null(rownames(contrastdefmat))) contrastdefdt[,
                                        facetrow := rownames(contrastdefmat)]
    melt.data.table(
        contrastdefdt,
        id.vars       = 'facetrow',
        variable.name = 'facetcol',
        value.name    = 'contrast', by = 'contrast')
}


default_coef <- function(object, fit = fits(object)){
    if (length(fit)==0) return(character(0))    # none
    y <- autonomics::coefs(object, fit = fit)   # intercept
    if (length(y)==1)  return(y)
    y %<>% setdiff('Intercept')                 # intercept + others
    y
}

#' Bin continuous variable
#' @param object numeric or SummarizedExperiment
#' @param fvar   string or NULL
#' @param probs  numeric
#' @return  factor vector
#' @examples 
#' # Numeric vector
#'     object <- rnorm(10, 5, 1)
#'     bin(object)
#' # SummarizedExperiment
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     fdt(object <- read_proteingroups(file))
#'     fdt(object %<>% add_assay_means('pepcounts'))
#'     fdt(bin(object, 'pepcounts'))
#' @export
bin <- function(object, ...) UseMethod('bin')


#' @rdname bin
#' @export
bin.logical <- function(object)  object

#' @rdname bin
#' @export
bin.numeric <- function(object, probs = c(0, 0.33, 0.66, 1)){
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
bin.SummarizedExperiment <- function(object, fvar, probs = c(0, 0.33, 0.66, 1)){
    if (is.null(fvar))  return(object)
    fdt(object)[[fvar]] %<>% bin()
    object
}


#' Add assay means
#' @param object SummarizedExperiment or NULL
#' @param assay string
#' @return SummarizedExperiment
#' @examples 
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file)
#' fdt(object)
#' fdt(add_assay_means(object))
#' @export
add_assay_means <- function(
    object, 
    assay = assayNames(object)[1], 
    bin = TRUE
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(assay, assayNames(object))
    assertive.sets::are_disjoint_sets(assay, fvars(object))
# Add
    if (is.null(assay))  return(object)
    for (.assay in assay){
        if (.assay %in% fvars(object))  fdt(object)[[.assay]] <- NULL
        fdt(object)$placeholder <- rowMeans(assays(object)[[assay]], na.rm = TRUE)
        fvars(object) %<>% stri_replace_first_fixed('placeholder', paste0('mean.', assay))
    }
# Return
    object
}

 
#' Add adjusted pvalues
#' 
#' @param object  SummarizedExperiment
#' @param method 'fdr', 'bonferroni', ... (see `p.adjust.methods`)
#' @param fit    'limma', 'lm', 'lme', 'lmer'
#' @param coef   coefficient (string)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file)
#' object %<>% fit_limma()
#' object %<>% extract(order(fdt(.)$`p~Adult~limma`), )
#' fdt(object)
#' fdt(object %>% add_adjusted_pvalues('bonferroni'))
#' @return SummarizedExperiment
#' @export
add_adjusted_pvalues <- function(
    object, method, fit = fits(object)[1], coef = default_coef(object, fit)
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(fit, fits(object))
    assert_is_subset(coef, coefs(object, fit))
    assert_is_subset(method, stats::p.adjust.methods)
# Compute
    for (.fit  in fit){
    for (.coef in coef){
    for (.method in method){
        pdt <- fdt(object)[ , pvar(object, coef = .coef, fit = .fit) , with = FALSE]
        pdt[, names(pdt) := lapply(.SD, p.adjust, method = .method), .SDcols = names(pdt)]
        names(pdt) %<>% stri_replace_first_fixed('p~', paste0(method, '~'))
    }}}
# Add
    fdt(object) %<>% cbind(pdt)
    object
}

    
#' Create volcano datatable
#' @param object  SummarizedExperiment
#' @param fit     'limma', 'lme', 'lm', 'wilcoxon'
#' @param coef   character vector: coefs for which to plot volcanoes
#' @return data.table
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = TRUE, fit = 'limma', plot = FALSE)
#' make_volcano_dt(object, fit = 'limma', coef = 'Adult')
#' @export
make_volcano_dt <- function(
    object, fit = fits(object)[1], coef = default_coef(object, fit = fit),
    shape = 'imputed', size = NULL, alpha = NULL,
    label = 'feature_id'
){
# Assert    
    assert_is_all_of(object, "SummarizedExperiment")
    assertive::assert_any_are_matching_regex(fvars(object), paste0('^p', FITSEP))
    assert_is_subset(fit, fits(object))
    assert_is_subset(coef, autonomics::coefs(object, fit))
    if (!is.null(shape)){ assert_is_subset(shape, fvars(object)); object %<>% bin(shape) }
    if (!is.null(size) ){ assert_is_subset(size,  fvars(object)); object %<>% bin(size)  }
    if (!is.null(alpha)){ assert_is_subset(alpha, fvars(object)); object %<>% bin(alpha) }
    if (!is.null(label))  assert_is_subset(label, fvars(object))
    object %<>% add_adjusted_pvalues('bonferroni', fit, coef)
# Prepare
    effect <- p <- mlp <- significance <- fdr <- NULL
    idvars <- 'feature_id'
    if ('control' %in% fvars(object))  idvars %<>% c('control')
    if (!is.null(label))  idvars %<>% union(label)
    if (!is.null(shape))  idvars %<>% union(shape)
    if (!is.null(size))   idvars %<>% union(size)
    if (!is.null(alpha))  idvars %<>% union(alpha)
    valuevars  <-  effectvar(object, coef = coef, fit = fit)  # elminate similar function pvars etc.
    valuevars %<>%  c(  pvar(object, coef = coef, fit = fit))

    dt <- fdt(object)[, c(idvars, valuevars), with = FALSE]
    dt %<>% melt.data.table(id.vars = idvars)
    dt %<>% tidyr::separate(.data$variable, into = c('quantity', 'coef', 'fit'), sep = FITSEP)
    idvars %<>% c('coef', 'fit')
    #dt %<>% dcast.data.table(feature_id+feature_name+coef+fit ~ quantity, value.var = 'value')
    dt %<>% tidyr::pivot_wider(id_cols = all_of(idvars), names_from = 'quantity', values_from = 'value')
    dt %<>% data.table()
    dt$coef %<>% factor(coef)
    dt[, fdr := p.adjust(p, method = 'fdr'),        by = c('fit', 'coef')]
    dt[, bon := p.adjust(p, method = 'bonferroni'), by = c('fit', 'coef')]
    dt[, mlp := -log10(p)]
    dt %<>% extract(!is.na(effect) & !is.na(p))
    dt[, direction := 'unchanged']
    dt[p>0.05, direction := 'unchanged']
    dt[p<=0.05 & effect>0, direction := 'up']
    dt[p<=0.05 & effect<0, direction := 'down']
}

    

#' Plot volcano
#' @param object    SummarizedExperiment
#' @param fit      'limma', 'lme', 'lm', 'wilcoxon'
#' @param coef      character vector
#' @param facet     svars mapped to facet
#' @param label     label fvar (string)
#' @param shape     shape fvar (string)
#' @param size      size  fvar (string)
#' @param alpha     alpha fvar (string)
#' @param max.overlaps  number: passed to ggrepel
#' @param features  character vector: features to plot 
#' @param nrow      number: no of rows in plot
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = TRUE)
#' object %<>% fit_limma()
#' plot_volcano(object)
#' plot_volcano(object, label = 'genesymbol')
#' plot_volcano(object, label = 'genesymbol', size = 'mean.log2.LFQ.intensity')
#' plot_volcano(object, label = 'genesymbol', size = 'mean.log2.LFQ.intensity', alpha = 'mean.pepcounts')
#' plot_volcano(object, label = 'genesymbol', features = c('hmbsb'))
#' object %<>% fit_lm()
#'
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, impute = TRUE, fit = 'limma', plot = FALSE)
#' plot_volcano(object, coef = c('t1', 't2', 't3'))
#' object %<>% fit_lm(subgroupvar = 'SET')
#' plot_volcano(object, coef = c('t1', 't2', 't3'), nrow=2)
#' plot_volcano(object, coef = c('t1', 't2', 't3'), fit='lm')
#' @export
plot_volcano <- function(
    object,
    fit   = fits(object), 
    coef  = default_coef(object, fit[1]),
    facet = if (is_scalar(fit)) 'coef' else c('fit', 'coef'),
    shape = if ('imputed' %in% fvars(object)) 'imputed' else NULL, 
    size  = NULL,
    alpha = NULL,
    label = 'feature_id', 
    max.overlaps = 10,
    features  = NULL,
    nrow  = length(fit)
){
# Assert
    assert_is_a_number(nrow)
    effect <- mlp <- NULL
    facet %<>% lapply(sym)
    facet <- vars(!!!facet)
# Volcano 
    plotdt <- make_volcano_dt(object, fit = fit, coef = coef, 
                  label = label, shape = shape, size = size, alpha = alpha)
    p <- ggplot(plotdt) + facet_wrap(facet, nrow = nrow)
    p <- p + theme_bw() + theme(panel.grid = element_blank())
    p <- p + xlab('log2(FC)') + ylab('-log10(p)') + ggtitle('volcano')
    shapesym <- if (is.null(shape))  quo(NULL) else sym(shape)
    sizesym  <- if (is.null(size))   quo(NULL) else sym(size)
    alphasym <- if (is.null(alpha))  quo(NULL) else sym(alpha)
    colorvalues <- c(down = '#ff5050', unchanged = 'grey', up = '#009933')
    p <- p + geom_point(data = plotdt, 
                        mapping = aes(x = effect, y = mlp, color = direction, 
                                      shape = !!shapesym, alpha = !!alphasym, 
                                      size = !!sizesym), 
                        na.rm = TRUE) + 
             scale_color_manual(values = colorvalues)
    if (!is.null(size))   p <- p + scale_size_manual( values = 1:3)
    if (!is.null(alpha))  p <- p + scale_alpha_manual(values = c(0.3, 0.5, 1))
# Significance lines
    xmin <- min(plotdt$effect)
    xmax <- max(plotdt$effect)
    dy <- 0.03*(max(plotdt$mlp) - min(plotdt$mlp))
    minn <- function(x) if (length(x)==0) return(Inf) else min(x, na.rm = TRUE) 
        # Avoid warning 'no non-missing arguments to min; returning Inf'
    summarydt <- plotdt[, 
        .( 
            label      = c('p = 0.05',                 'fdr = 0.05',                 'bon = 0.05'),
            yintercept = c( minn(mlp[p<0.05]),          minn(mlp[fdr<0.05]),          minn(mlp[bon<0.05])),
            down       = c( sum(p <0.05 & effect < 0),  sum(fdr <0.05 & effect < 0),  sum(bon <0.05 & effect < 0)),
            up         = c( sum(p <0.05 & effect > 0),  sum(fdr <0.05 & effect > 0),  sum(bon <0.05 & effect > 0))
        ),
        by = c('fit', 'coef')
    ]
    summarydt %<>% extract(!is.infinite(yintercept))
    p <- p + geom_hline(data = summarydt, mapping = aes(yintercept = yintercept), color = 'gray30', linetype = 'longdash')
    do_geom_label <- function(mapping, label.size){  
                        geom_label(data       = summarydt, 
                                   mapping    = mapping, 
                                   label.size = label.size, 
                                   hjust      = 0.5, 
                                   fill       = alpha(c('white'), 0.6), 
                                   color      = 'gray30')  }
    p <- p + do_geom_label(mapping = aes(x = 0,    y = yintercept,      label = label), label.size = 0.5)
    p <- p + do_geom_label(mapping = aes(x = xmin, y = yintercept + dy, label = down ), label.size = NA)
    p <- p + do_geom_label(mapping = aes(x = xmax, y = yintercept + dy, label = up   ), label.size = NA)
    p <- p + guides(color = 'none')
# Labels
    if (!is.null(label)){
        labeldt <- plotdt[fdr<0.05]
        if (!is.null(features)){
            seldt <- copy(plotdt)
            seldt[, singlefeature := feature_id]
            seldt %<>% separate_rows(singlefeature) %>% data.table()
            seldt %<>% extract(singlefeature %in% features)
            seldt[, singlefeature := NULL]
            seldt %<>% unique()
        }
        p <- p + ggrepel::geom_label_repel(data = labeldt, 
                    aes(x = effect, y = mlp, label = !!sym(label), color = direction), #color = 'black', 
                    label.size = NA, fill = alpha(c('white'), 0.6),
                    na.rm = TRUE, show.legend = FALSE, max.overlaps = max.overlaps)
        if (!is.null(features)){
            p <- p + ggrepel::geom_label_repel(data = seldt, 
                        aes(x = effect, y = mlp, label = !!sym(label)), color = 'black', 
                        label.size = NA, fill = alpha(c('white'), 0.6),
                        na.rm = TRUE, show.legend = FALSE, max.overlaps = max.overlaps)
            p <- p + geom_point(data = seldt, aes(x = effect, y = mlp), shape = 1, size = 4, color = 'black')
        }
    }
    p
}
