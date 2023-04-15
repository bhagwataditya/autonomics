
#' Analyze
#' @param object       SummarizedExperiment
#' @param pca          TRUE / FALSE: perform pca ?
#' @param pls          TRUE / FALSE: perform pls ?
#' @param fit          glm engine: 'limma', 'lm', 'lme(r)', 'lmer', 'wilcoxon'
#' @param formula      model formula
#' @param drop         TRUE / FALSE : drop varname in designmat ?
#' @param coding       factor coding system : 'treatment', 'reference', difference', 
#'                                            'grandref', 'granddiff', 'sum', 'helmert'
#' @param contrasts    model coefficient contrasts of interest: string vector or NULL
#' @param coefs        model coefficients          of interest: string vector or NULL
#' @param block        model blockvar
#' @param weightvar    NULL or name of weight matrix in assays(object)
#' @param plot         TRUE / FALSE
#' @param palette      NULL or colorvector
#' @param verbose      TRUE / FALSE: message?
#' @return SummarizedExperiment
#' @examples 
#' require(magrittr)
#' file <- download_data('atkin18.metabolon.xlsx')
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
    coding       = 'treatment', 
    contrasts    = NULL,
    coefs        = colnames(create_design(object, formula = formula, drop = drop)),
    block        = NULL,
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL,
    plot         = pca & !is.null(fit),
    palette      = NULL,
    verbose      = TRUE
){
    # Analyze
    if (is.null(palette))       palette <- make_subgroup_palette(object)
    if (pca)  object %<>% pca(verbose = verbose, plot = FALSE)
    if (pls)  object %<>% pls(by = all.vars(formula)[1], verbose = FALSE)
    for (curfit in fit){
        fitfun <- get(paste0('fit_', curfit))
        if (is.null(formula)) formula <- default_formula(object)
        if (is.null(coefs))   coefs <- colnames(create_design(object, formula = formula, drop = drop))
        object %<>% fitfun(
            formula      = formula,       drop         = drop,
            coding       = coding,        contrasts    = contrasts,
            coefs        = coefs,         block        = block,
            weightvar    = weightvar,     verbose      = verbose,
            plot         = FALSE) 
    }
    # Plot/Return
    if (plot)  plot_summary(object, fit = fit, block = block, palette = palette)
    object
}


#' Plot summary
#' @param object       SummarizedExperiment
#' @param fit          glm engine : 'limma', 'lm', 'lme', 'lmer' or 'wilcoxon'
#' @param block        NULL or svar
#' @param palette      NULL or colorvector
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% pca()
#' object %<>% pls(by = 'subgroup')
#' object %<>% fit_limma()
#' plot_summary(object, block = 'SUB')
#' @export
plot_summary <- function(
    object, fit = 'limma', block = NULL, palette = make_subgroup_palette(object)
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(c('effect~sample_id~pca1', 'effect~sample_id~pca2'), svars(object))
    assert_is_subset(c('effect~sample_id~pca1', 'effect~sample_id~pca2'), fvars(object))
    assert_is_subset(c('effect~subgroup~pls1',  'effect~subgroup~pls2'), svars(object))
    assert_is_subset(c('effect~subgroup~pls1',  'effect~subgroup~pls2'), fvars(object))
    assertive::assert_any_are_matching_regex(fvars(object), 'limma')
# Plot
    detections <- plot_subgroup_nas(object, 
                    palette  = palette) + ggtitle('Detections') + xlab(NULL) + 
                    theme(plot.title = element_text(hjust = 0.5))
    pcaplot <- biplot(object, method = 'pca', color = 'subgroup', colorpalette = palette, nfeatures = 1) + 
               guides(color = 'none', linetype = 'none') + # scale_x_continuous(expand = c(0.2, 0.2)) + 
               theme(axis.text.x  = element_blank(), axis.text.y  = element_blank(), 
                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
                     plot.title   = element_text(hjust = 0.5))
    plsplot <- biplot(object, method = 'pls', color = 'subgroup', colorpalette = palette) + 
               guides(color = 'none', linetype = 'none') + # scale_x_continuous(expand = c(0.2, 0.2)) +
               theme(axis.text.x  = element_blank(), axis.text.y  = element_blank(), 
                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
                     plot.title   = element_text(hjust = 0.5))
    samples  <- plot_top_samples(object, palette = palette)
    features <- plot_top_features(object, fit = fit)
    pca1exprs <- plot_exprs(object, fit = 'pca1',  block = block, n = 1, subtitle = 'X1', title = NULL)
    pca2exprs <- plot_exprs(object, fit = 'pca2',  block = block, n = 1, subtitle = 'X2', title = NULL)
    pls1exprs <- plot_exprs(object, fit = 'pls1',  block = block, n = 1, subtitle = 'X1', title = NULL)
    pls2exprs <- plot_exprs(object, fit = 'pls2',  block = block, n = 1, subtitle = 'X2', title = NULL)
    coefs <- default_coefs(object, fit = fit)
    glm1exprs <- plot_exprs(object, fit = fit, coefs = coefs[1], block = block, n = 1, nrow = 1, subtitle = coefs[1], title = NULL)
    pca1exprs <- pca1exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    pca2exprs <- pca2exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    pls1exprs <- pls1exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    pls2exprs <- pls2exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    glm1exprs <- glm1exprs + guides(color = 'none', fill = 'none') + ylab(NULL)
    glm1volcano <- plot_top_volcano(object, fit = fit, coef = coefs[1]) + guides(shape = 'none', linetype = 'none')
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


plot_top_samples <- function(object, palette = make_subgroup_palette(object)){
    plot_sample_densities( object, n = 4, fill = 'subgroup', palette = palette) + 
        ggtitle('Samples') + theme_void() + coord_flip() + 
        theme(plot.title       = element_text(hjust = 0.5), 
              legend.position  = 'left', 
              legend.direction = 'vertical', 
              legend.title     = element_blank())
}

plot_top_features <- function(object, fit){
    idx1 <- order(abs(fdt(object)$`effect~sample_id~pca1`), decreasing = TRUE)[1]
    #idx2 <- order(abs(fdt(object)$`effect~sample_id~pca2`), decreasing = TRUE)[1]
    idx3 <- order(abs(fdt(object)$`effect~subgroup~pls1` ), decreasing = TRUE)[1]
    #idx4 <- order(abs(fdt(object)$`effect~subgroup~pls2` ), decreasing = TRUE)[1]
    pvr <- pvar(object, fit = fit, coef = default_coefs(object, fit = fit))[1]
    idx5 <- order(abs(fdt(object)[[pvr]]))[1]
    idx  <- unique(c(idx1,idx3,idx5))
    #idx  <- unique(c(idx1,idx2,idx3,idx4,idx5))
    
    plot_feature_densities(
        object[idx, ], n = length(idx)) + 
        xlab(feature_id) + ggtitle('Features') + theme_void() + 
        theme(plot.title       = element_text(hjust = 0.5), 
              legend.position  = 'bottom', 
              legend.direction = 'vertical', 
              legend.title     = element_blank())
}

top_coef <- function(object, fit){
    summarydt <- summarize_fit(fdt(object), fit = fit)
    summarydt %<>% extract(which.max(ndown))  # nup favours intercept !
    summarydt$contrast
}

plot_top_volcano <- function(object, fit, coef){
    summarydt <- summarize_fit(fdt(object), fit = fit)
    summarydt %<>% extract(coefficient == coef)  # nup favours intercept !
    
    plot_volcano(object, fit = fit, coef = coef, ) + 
    guides(color = 'none') + 
    ggtitle(sprintf('%s~%s', fit, coef)) + 
    #ggtitle(sprintf('%d down | %d up', summarydt$ndown, summarydt$nup)) +
    theme(plot.title = element_text(hjust = 0.5))
}