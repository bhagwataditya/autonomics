
#' Analyze
#' @param object     SummarizedExperiment
#' @param pca        TRUE / FALSE: perform pca ?
#' @param pls        TRUE / FALSE: perform pls ?
#' @param fit        linmod engine: 'limma', 'lm', 'lme(r)', 'lmer', 'wilcoxon'
#' @param formula    model formula
#' @param drop       TRUE / FALSE : drop varname in designmat ?
#' @param codingfun  factor coding function
#' \itemize{
#'     \item contr.treatment:          intercept = y0,     coefi = yi - y0
#'     \item contr.treatment.explicit: intercept = y0,     coefi = yi - y0
#'     \item code_control:             intercept = ymean,  coefi = yi - y0
#'     \item contr.diff:               intercept = y0,     coefi = yi - y(i-1)
#'     \item code_diff:                intercept = ymean,  coefi = yi - y(i-1)
#'     \item code_diff_forward:        intercept = ymean,  coefi = yi - y(i+)
#'     \item code_deviation:           intercept = ymean,  coefi = yi - ymean (drop last)
#'     \item code_deviation_first:     intercept = ymean,  coefi = yi - ymean (drop first)
#'     \item code_helmert:             intercept = ymean,  coefi = yi - mean(y0:(yi-1))
#'     \item code_helmert_forward:     intercept = ymean,  coefi = yi - mean(y(i+1):yp)
#' }
#' @param contrasts  model coefficient contrasts of interest: string vector or NULL
#' @param coefs      model coefficients          of interest: string vector or NULL
#' @param block      model blockvar
#' @param weightvar  NULL or name of weight matrix in assays(object)
#' @param plot       TRUE / FALSE
#' @param label      fvar
#' @param palette    NULL or colorvector
#' @param verbose    TRUE / FALSE: message?
#' @return SummarizedExperiment
#' @examples 
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% analyze()
#' @export
analyze <- function(
    object,
    pca          = TRUE,
    pls          = TRUE,
    fit          = 'limma',
    formula      = default_formula(object),
    drop         = varlevels_dont_clash(object, all.vars(formula)),
    codingfun    = contr.treatment, 
    contrasts    = NULL,
    coefs        = colnames(create_design(object, formula = formula, drop = drop)),
    block        = NULL,
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL,
    plot         = pca & !is.null(fit),
    label        = 'feature_id',
    palette      = NULL,
    verbose      = TRUE
){
    # Analyze
    if (is.null(palette))       palette <- make_svar_palette(object, svar = all.vars(formula)[1])
    if (pca)  object %<>% pca(verbose = verbose, plot = FALSE)
    if (pls)  object %<>% pls(by = all.vars(formula)[1], verbose = FALSE)
    for (curfit in fit){
        fitfun <- get(paste0('fit_', curfit))
        if (is.null(formula)) formula <- default_formula(object)
        if (is.null(coefs))   coefs <- colnames(create_design(object, formula = formula, drop = drop, verbose = FALSE))
        object %<>% fitfun(
            formula      = formula,       drop         = drop,
            codingfun    = codingfun,     contrasts    = contrasts,
            coefs        = coefs,         block        = block,
            weightvar    = weightvar,     verbose      = verbose,
            plot         = FALSE) 
    }
    # Plot/Return
    if (plot)  plot_summary(object, fit = fit, formula = formula, block = block, label = label, palette = palette)
    object
}


#' Plot summary
#' @param object   SummarizedExperiment
#' @param fit      linmod engine : 'limma', 'lm', 'lme', 'lmer' or 'wilcoxon'
#' @param formula  model formula
#' @param block    NULL or svar
#' @param label    fvar
#' @param palette  NULL or colorvector
#' @examples 
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% pca()
#' object %<>% pls(by = 'subgroup')
#' object %<>% fit_limma()
#' plot_summary(object, block = 'Subject')
#' @export
plot_summary <- function(
    object, 
    fit = 'limma', 
    formula = default_formula(object), 
    block = NULL, 
    label = 'feature_id', 
    palette = make_svar_palette(object, svar = svar)
){
# Assert
    assert_is_valid_sumexp(object)
    cmessage('\t\tPlot summary')
    if (!all(c('pca1', 'pca2') %in% fits(object))){ 
        cmessage('\t\t\tFirst `pca(object)`, then `plot_summary(object)`'); return(NULL) }
    if (!all(c('pls1', 'pls2') %in% fits(object))){
        cmessage('\t\t\tFirst `pls(object)`, then `plot_summary(object)`'); return(NULL) }
    assert_is_subset(fit, fits(object))
    assert_any_are_matching_regex(fvars(object), fit)
    assert_is_formula(formula)
# Plot
    svar <- all.vars(formula)[1]
    detections <- plot_subgroup_nas(object, by = svar,
                    palette  = palette) + ggtitle('Detections') + xlab(NULL) + 
                    theme(plot.title = element_text(hjust = 0.5))
    pcaplot <- biplot(object, method = 'pca', color = svar, colorpalette = palette) + 
               guides(color = 'none', linetype = 'none') + # scale_x_continuous(expand = c(0.2, 0.2)) + 
               theme(axis.text.x  = element_blank(), axis.text.y  = element_blank(), 
                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
                     plot.title   = element_text(hjust = 0.5))
    plsplot <- biplot(object, method = 'pls', by = svar, color = svar, colorpalette = palette) + 
               guides(color = 'none', linetype = 'none') + # scale_x_continuous(expand = c(0.2, 0.2)) +
               theme(axis.text.x  = element_blank(), axis.text.y  = element_blank(), 
                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
                     plot.title   = element_text(hjust = 0.5))
    samples  <- plot_top_samples(object, svar = svar, palette = palette)
    features <- plot_top_features(object, fit = fit)
    pca1exprs <- plot_exprs(object, fit = 'pca1',  block = block, n = 1, subtitle = 'X1', title = NULL, verbose = FALSE, color = svar, fill = svar, x = svar)
    pca2exprs <- plot_exprs(object, fit = 'pca2',  block = block, n = 1, subtitle = 'X2', title = NULL, verbose = FALSE, color = svar, fill = svar, x = svar)
    pls1exprs <- plot_exprs(object, fit = 'pls1',  block = block, n = 1, subtitle = 'X1', title = NULL, verbose = FALSE, color = svar, fill = svar, x = svar)
    pls2exprs <- plot_exprs(object, fit = 'pls2',  block = block, n = 1, subtitle = 'X2', title = NULL, verbose = FALSE, color = svar, fill = svar, x = svar)
    coefs <- default_coefs(object, fit = fit)
    glm1exprs <- plot_exprs(object, fit = fit, coefs = coefs[1], block = block, n = 1, nrow = 1, subtitle = coefs[1], title = NULL, verbose = FALSE, color = svar, fill = svar, x = svar)
    pca1exprs <- pca1exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    pca2exprs <- pca2exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    pls1exprs <- pls1exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    pls2exprs <- pls2exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    glm1exprs <- glm1exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    glm1volcano <- plot_top_volcano(object, fit = fit, coef = coefs[1], label = label) + guides(shape = 'none', linetype = 'none')
# Arrange
    layout_matrix <- matrix(c(1, 2, 3, 
                              4, 5, 6, 
                              7, 8, 6, 
                              9,10,11), nrow = 4, byrow = TRUE)
    grid.arrange(detections, samples, features, 
                 pcaplot,   plsplot,   glm1volcano, 
                 pca1exprs, pls1exprs,#glm1volcano, 
                 pca2exprs, pls2exprs, glm1exprs,
                 layout_matrix = layout_matrix)
}


plot_top_samples <- function(object, svar, palette = make_svar_palette(object, svar)){
    plot_sample_densities( object, n = 4, fill = svar, palette = palette) + 
        ggtitle('Samples') + theme_void() + coord_flip() + 
        theme(plot.title       = element_text(hjust = 0.5), 
              legend.position  = 'left', 
              legend.direction = 'vertical', 
              legend.title     = element_blank())
}

plot_top_features <- function(object, fit, label){
    idx1 <- order(abs(loadings(object, method = 'pca', dim = 1)), decreasing = TRUE)[1]
    idx2 <- order(abs(loadings(object, method = 'pls', dim = 1)), decreasing = TRUE)[1]
    pvr <- pvar(object, fit = fit, coef = default_coefs(object, fit = fit))[1]
    idx3 <- order(abs(fdt(object)[[pvr]]))[1]
    idx  <- unique(c(idx1,idx2,idx3))

    plot_feature_densities(object[idx, ], n = length(idx)) + 
        xlab(feature_id) + ggtitle('Features') + theme_void() + 
        theme(plot.title       = element_text(hjust = 0.5), 
              legend.position  = 'bottom', 
              legend.direction = 'vertical', 
              legend.title     = element_blank())
}

top_coef <- function(object, fit){
    ndown <- NULL
    summarydt <- summarize_fit(object, fit = fit)
    summarydt %<>% extract(which.max(ndown))  # nup favours intercept !
    summarydt$contrast
}

plot_top_volcano <- function(object, fit, coef, label){
    coefficient <- NULL
    summarydt <- summarize_fit(object, fit = fit)
    summarydt %<>% extract(coefficient == coef)  # nup favours intercept !
    
    plot_volcano(object, fit = fit, coefs = coef, label = label) + 
    guides(color = 'none') + 
    ggtitle(sprintf('%s~%s', fit, coef)) + 
    #ggtitle(sprintf('%d down | %d up', summarydt$ndown, summarydt$nup)) +
    theme(plot.title = element_text(hjust = 0.5))
}