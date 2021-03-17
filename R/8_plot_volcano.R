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
    coef_top <- if (any(fdr_ok)){  effect < nmin(effect[fdr_ok], ntop+1)
                } else {           rep(FALSE, length(effect))            }
    mlp_top  <- if (any(coef_ok)){ mlp  > nmax(mlp[coef_ok], ntop+1)
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
    coef_top <- if(any(fdr_ok)){  effect > nmax(effect[fdr_ok], ntop+1)
                } else {          rep(FALSE, length(effect)) }
    mlp_top  <- if (any(coef_ok)){ mlp > nmax(mlp[coef_ok], ntop+1)
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
#' @param object          SummarizedExperiment
#' @param fit            'limma', 'lme', 'lm', 'wilcoxon'
#' @param contrastdefmat  contrastdef matrix
#' @param ntop            no of top features to be annotated
#' @return data.table
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, fit='limma', plot=FALSE)
#' make_volcano_dt(object, fit = 'limma')
#' @export
make_volcano_dt <- function(
    object, fit, contrastdefmat = contrastdefs(object)[[1]], ntop = 3
){
    effect <- p <- mlp <- topdown <- topup <- significance <- fdr <- NULL
    dt <- extract_fit_dt(object, fit)
    dt %<>% merge(melt_contrastdefs(contrastdefmat), by = 'contrast')
    dt %<>% extract(!is.na(effect) & !is.na(p))
    dt[, mlp  := -log10(p)]

    # Prepare volcano datatable
    # Note: Using effect <= 0 (rather than effect <0) is required.
    # Otherwise the (very few) features with effect=0 will have no effect for
    # 'significance'
    by <- intersect(c('contrast', 'imputed'), names(dt))
    dt[,topdown := top_down(effect, fdr, mlp, ntop), by=by]
    dt[,topup   := top_up(  effect, fdr, mlp, ntop), by=by]
    dt[effect<=0,            significance := 'down']
    dt[effect> 0,            significance :=   'up']
    dt[effect<=0 & fdr<0.05, significance := 'down: fdr<0.05']
    dt[effect> 0 & fdr<0.05, significance :=   'up: fdr<0.05']
    dt[topdown==TRUE,        significance := 'down: top']
    dt[topup  ==TRUE,        significance :=   'up: top']
    dt[,significance := factor(significance, c(
        'down: top','down: fdr<0.05','down','up','up: fdr<0.05','up: top'))]
    dt[]
}

#' Plot volcano
#' @param object         SummarizedExperiment
#' @param fit           'limma', 'lme', 'lm', 'wilcoxon'
#' @param contrastdefs   contrastdef vector / matrix / list
#' @param label          fvar for labeling top features
#' @param ntop           number: n top features to be annotated
#' @return ggplot object
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, fit='limma', plot=FALSE)
#' plot_volcano(object)
#' @export
plot_volcano <- function(object, 
    fit = intersect(names(metadata(object)), TESTS)[1], 
    contrastdefs = autonomics::contrastdefs(object)[[1]], 
    label = feature_name, ntop = 1
){
# Assert/Process
    assert_is_all_of(object, "SummarizedExperiment")
    assert_is_a_string(fit)
    assert_is_subset(fit, TESTS)
    assert_is_subset(fit, names(metadata(object)))
    if (is.null(contrastdefs)) contrastdefs <- colnames(metadata(object)[[fit]])
    if (is.character(contrastdefs)) contrastdefs %<>% contrvec2mat()
    if (is.list(contrastdefs))      contrastdefs %<>% extract2(1)
    assert_is_matrix(contrastdefs)
    topup <- topdown <- effect <- mlp <- facetrow <- facetcol <- NULL
    label <- enquo(label)
# Prepare
    plotdt <- make_volcano_dt(object, fit, contrastdefs, ntop = ntop)
    txtdt  <- copy(plotdt)[topup==TRUE | topdown==TRUE]
    colorvalues <-c(hcl(h=  0, l=c(20, 70, 100), c=100), # 20 70 100
                    hcl(h=120, l=c(100, 70, 20), c=100)) # 100 70 20
    names(colorvalues) <- levels(plotdt$significance)
# Plot
    imputed <- NULL # fallback when plotdt misses "imputed"
    significance <- NULL
    p <- ggplot(plotdt) + facet_grid(
        rows = if (all(stri_isempty(plotdt$facetrow))) NULL else vars(facetrow),
        cols = if (all(stri_isempty(plotdt$facetcol))) NULL else vars(facetcol),
        scales = 'fixed') +
    geom_point(aes(x=effect,y=mlp,color=significance,shape=imputed),na.rm=TRUE)
    if (!quo_is_null(label)) p <- p + geom_text_repel(
                        data = txtdt,
                        aes(x=effect, y=mlp, label=!!label, color=significance),
                        #hjust = 'outward',
                        na.rm = TRUE,
                        show.legend = FALSE)#,
                        #direction = 'x'
    p + theme_bw() +
        scale_color_manual(values = colorvalues, name = NULL) +
        xlab('log2(FC)') +
        ylab('-log10(p)') +
        ggtitle(fit)#+
        #guides(color=FALSE)
}
