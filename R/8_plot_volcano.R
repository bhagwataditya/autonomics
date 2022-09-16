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


#' Create volcano datatable
#' @param object  SummarizedExperiment
#' @param fit     'limma', 'lme', 'lm', 'wilcoxon'
#' @param coef   character vector: coefs for which to plot volcanoes
#' @param ntop    no of top features to be annotated
#' @return data.table
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = TRUE, fit='limma', plot=FALSE)
#' make_volcano_dt(object, fit = 'limma', coef = 'subgroupAdult')
#' @export
make_volcano_dt <- function(
    object, fit = fits(object)[1], coef = autonomics::coefs(object)[1],
    label = 'feature_id', ntop = 3
){
    effect <- p <- mlp <- topdown <- topup <- significance <- fdr <- NULL
    id.vars <- c('feature_id', label, 'imputed', 'control')
    id.vars %<>% intersect(fvars(object))
    value.vars  <-  effectvar(object, coef = coef, fit = fit)  # elminate similar function pvars etc.
    value.vars %<>%  c(  pvar(object, coef = coef, fit = fit))
    value.vars %<>%  c(fdrvar(object, coef = coef, fit = fit))
    value.vars %<>%  c(bonvar(object, coef = coef, fit = fit))
    
    dt <- fdt(object)[, c(id.vars, value.vars), with = FALSE]
    dt %<>% melt.data.table(id.vars = id.vars)
    dt %<>% tidyr::separate(
                .data$variable, into = c('quantity', 'coef', 'fit'), sep = FITSEP)
    id.vars %<>% c('coef', 'fit')
    #dt %<>% dcast.data.table(feature_id+feature_name+coef+fit ~ quantity, value.var = 'value')
    dt %<>% tidyr::pivot_wider(
            id_cols = id.vars, names_from = 'quantity', values_from = 'value')
    dt %<>% data.table()
    #dt <- extract_fit_dt(object, fit)
    #dt %<>% merge(melt_contrastdefs(contrastdefmat), by = 'coef')
    dt %<>% extract(!is.na(effect) & !is.na(p))
    dt[, mlp  := -log10(p)]

    # Prepare volcano datatable
    # Note: Using effect <= 0 (rather than effect <0) is required.
    # Otherwise the (very few) features with effect=0 will have no effect for
    # 'significance'
    by <- intersect(c('coef', 'imputed', 'fit'), names(dt))
    #dt[,topdown := top_down(effect, fdr, mlp, ntop), by = by]
    #dt[,topup   := top_up(  effect, fdr, mlp, ntop), by = by]
    dt[effect <= 0,              significance := 'down']
    dt[effect <= 0 & p   < 0.05, significance := 'down: p < 0.05']
    dt[effect <= 0 & fdr < 0.05, significance := 'down: fdr < 0.05']
    dt[effect <= 0 & bon < 0.05, significance := 'down: bon < 0.05']
    dt[effect >  0,              significance := 'up']
    dt[effect >  0 & p   < 0.05, significance := 'up: p < 0.05']
    dt[effect >  0 & fdr < 0.05, significance := 'up: fdr < 0.05']
    dt[effect >  0 & bon < 0.05, significance := 'up: bon < 0.05']
    # dt[topdown == TRUE,        significance := 'down: top']
    # dt[topup == TRUE,        significance :=   'up: top']
    dt$significance %<>% factor(c(
        'down: bon < 0.05', 'down: fdr < 0.05', 'down: p < 0.05', 'down',
        'up', 'up: p < 0.05', 'up: fdr < 0.05', 'up: bon < 0.05'))
    dt[]
}

#' Plot volcano
#' @param object    SummarizedExperiment
#' @param fit      'limma', 'lme', 'lm', 'wilcoxon'
#' @param coef     character vector
#' @param label     fvar for labeling top features (string)
#' @param features  character vector: features to plot 
#' @param ntop      number: n top features to be annotated
#' @param nrow      number: no of rows in plot
#' @param intercept TRUE/FALSE: plot also intercept?
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = TRUE, fit = 'limma', plot = FALSE)
#' plot_volcano(object)
#' plot_volcano(object, label = 'genesymbol')
#' plot_volcano(object, label = 'genesymbol', features = c('F1QDE4', 'Q503D2'))
#' object %<>% fit_lm()
#' plot_volcano(object)
#'
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, impute=TRUE, fit='limma', plot=FALSE)
#' plot_volcano(object, coef = c('t1', 't2', 't3'))
#' object %<>% fit_lm(subgroupvar = 'SET')
#' plot_volcano(object, coef = c('t1', 't2', 't3'), nrow=2)
#' plot_volcano(object, coef = c('t1', 't2', 't3'), fit='lm')
#' @export
plot_volcano <- function(object,
    fit       = fits(object), 
    coef     = autonomics::coefs(object, fit[1]),
    label     = 'feature_id', 
    max.overlaps = 10,
    features  = NULL,
    nrow      = length(fit),
    intercept = identical("Intercept", coef)
){
# Assert/Process
    assert_is_all_of(object, "SummarizedExperiment")
    assert_is_subset(fit, fits(object))
    assert_is_subset(coef, autonomics::coefs(object, fit))
    if (!is.null(label)){
        assert_is_a_string(label)
        assert_is_subset(label, fvars(object))
    }
    assert_is_a_number(ntop)
    assert_is_a_number(nrow)
    assert_is_a_bool(intercept)
    if (!intercept) coef %<>% setdiff('Intercept')
    topup <- topdown <- effect <- mlp <- facetrow <- facetcol <- NULL
# Bonferroni
    pdt <- fdt(object)[ , pvar(object, coef = coef, fit = fit) , with = FALSE]
    pdt[, names(pdt) := lapply(.SD, p.adjust, method = 'bonferroni'), .SDcols = names(pdt)]
    names(pdt) %<>% stri_replace_first_fixed('p~', 'bon~')
    fdt(object) %<>% cbind(pdt)
# Prepare
    plotdt <- make_volcano_dt(object, fit = fit, coef = coef, ntop = ntop, label = label)
    bondt <- plotdt[bon <= 0.05]
    if (!is.null(features)){
        seldt <- copy(plotdt)
        seldt[, singlefeature := feature_id]
        seldt %<>% separate_rows(singlefeature) %>% data.table()
        seldt %<>% extract(singlefeature %in% features)
        seldt[, singlefeature := NULL]
        seldt %<>% unique()
    }
    colorvalues <- c(hcl(h =   0, l = c(30,  60, 80, 100), c = 100), # 20 70 100
                     hcl(h = 120, l = c(100, 80, 60,  30), c = 100)) # 100 70 20
    names(colorvalues) <- levels(plotdt$significance)
# Plot
    imputed <- NULL # fallback when plotdt misses "imputed"
    significance <- NULL
    p <- ggplot(plotdt) + 
         theme_bw() + 
         facet_wrap(fit~coef, nrow = nrow) +
         geom_point(aes(x = effect,y = mlp, color = significance, shape = imputed), na.rm = TRUE) + 
         scale_color_manual(values = colorvalues, name = NULL) + 
         xlab('log2(FC)') +
         ylab('-log10(p)') +
         ggtitle('volcano')#+
        
    if (!is.null(label)){
        p <- p + ggrepel::geom_label_repel(data = bondt, 
                    aes(x = effect, y = mlp, label = !!sym(label), color = significance), #color = 'black', 
                    label.size = NA, fill = alpha(c('white'), 1),
                    na.rm = TRUE, show.legend = FALSE, max.overlaps = max.overlaps)
        if (!is.null(features)){
            p <- p + ggrepel::geom_label_repel(data = seldt, 
                        aes(x = effect, y = mlp, label = !!sym(label)), color = 'black', 
                        label.size = NA, fill = alpha(c('white'), 1),
                        na.rm = TRUE, show.legend = FALSE, max.overlaps = max.overlaps, )
            p <- p + geom_point(data = seldt, aes(x = effect, y = mlp), shape = 1, size = 4, color = 'black')
        }
    }
    p
        #guides(color = 'none')
}
