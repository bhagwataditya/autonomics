
#' Analyze
#' @param object       SummarizedExperiment
#' @param pca          whether to perform pca
#' @param fit          NULL, 'limma', 'lm', 'lme', 'lmer', or 'wilcoxon'
#' @param subgroupvar  subgroup svar
#' @param formula      model formula
#' @param block        block svar
#' @param weightvar    NULL or name of weight matrix in assays(object)
#' @param coefs        NULL or character vector: model coefficients to test
#' @param contrasts NULL or character vector: coefficient contrasts to test
#' @param verbose      whether to msg
#' @param plot         whether to plot
#' @param feature_id   string: which feature to visualize
#' @param sample_id    string: which sample to visualize
#' @param palette      color vector: values = colors, names = slevels
#' @return SummarizedExperiment
#' @examples 
#' require(magrittr)
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' object %<>% analyze(pca = TRUE, subgroupvar = 'Group', fit = 'limma')
#' @export
analyze <- function(
    object,
    pca          = TRUE,
    fit          = 'limma',
    subgroupvar  = default_subgroupvar(object),
    contrasts    = NULL,
    formula      = default_formula(object, subgroupvar, contrasts),
    drop         = varlevels_dont_clash(object, all.vars(formula)),
    block        = NULL,
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL,
    coefs        = colnames(create_design(object, formula = formula, drop = drop)),
    verbose      = TRUE,
    plot         = pca & !is.null(fit),
    feature_id   = NULL,
    sample_id    = NULL,
    palette      = NULL
){
    # Analyze
    if (is.null(subgroupvar))  subgroupvar <- default_subgroupvar(object)
    if (is.null(palette))      palette <- make_subgroup_palette(object)
    subgroup <- if (is.null(subgroupvar))  quo(NULL) else sym(subgroupvar)
    if (pca)  object %<>% pca(verbose = verbose, plot = FALSE)
    for (curfit in fit){
        fitfun <- get(paste0('fit_', curfit))
        if (is.null(formula)) formula <- default_formula(object, subgroupvar, contrasts)
        if (is.null(coefs)) coefs <- colnames(create_design(object, formula = formula, drop = drop))
        object %<>% fitfun( subgroupvar  = subgroupvar,
                            formula      = formula,
                            coefs        = coefs,
                            contrasts    = contrasts,
                            block        = block,
                            weightvar    = weightvar,
                            verbose      = verbose,
                            plot         = FALSE) 
    }
    # Plot/Return
    if (plot)  plot_summary(object, fit = fit, feature_id = feature_id, 
                            sample_id = sample_id, palette = palette)
    object
}


#' Plot summary
#' @param object       SummarizedExperiment
#' @param fit          string
#' @param coef         string
#' @param feature_id   NULL or string
#' @param sample_id    NULL or string
#' @param palette      NULL or color vector
#' @examples 
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot = FALSE)
#' plot_summary(object, fit = 'limma')
#' @export
plot_summary <- function(
    object, 
    fit         = fits(object)[1], 
    coef        = c(setdiff(coefs(object), 'Intercept'), 'Intercept')[1],
    feature_id  = NULL, 
    sample_id   = NULL, 
    palette     = make_subgroup_palette(object)
){
# Initialize
    if (is.null(sample_id)){  # most avg sample
        sample_id  <- snames(object)[which.min(abs(object$pca1))] 
    }
    if (is.null(feature_id)){
        #idx <- which.max(abs(fdata(object)$pca1))
        pvar <- paste('p', coef, fit, sep = FITSEP)
        idx <- which.min(fdt(object)[[pvar]])
        feature_id <- fnames(object)[idx]
    }
# Create plots
    detections <- plot_summarized_detections(
        object, 
        palette  = palette) + ggtitle('Detections') + xlab('subgroup') + 
        theme(plot.title = element_text(hjust = 0.5))

    pcaplot <- biplot(
        object, color = 'subgroup', x = 'pca1', y = 'pca2', palette = palette) + 
        guides(color = 'none') + ggtitle(NULL) +
        theme(axis.text.x  = element_blank(), 
              axis.text.y  = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank()) +
        scale_x_continuous(expand = c(0.2, 0.2)) + 
        ggtitle('Pca') +
        theme(plot.title = element_text(hjust = 0.5))

    sample  <- plot_top_sample(object[, sample_id],  palette = palette)
    feature <- plot_top_feature(object[feature_id,], palette = palette)
    volcano <- plot_top_volcano(object, fit = fit, coef = coef)
    # Layout
    layout_matrix <- matrix(c(1,2,3,4,5,3), nrow = 2, byrow = TRUE)
    grid.arrange(pcaplot, detections, volcano, 
                 sample, feature, layout_matrix = layout_matrix)
}


plot_top_sample <- function(object, palette = make_subgroup_palette(object)){
    plot_sample_densities(
        object, 
        fill = 'subgroup', facet = 'sample_id', palette = palette,
        fixed = list(alpha=1, na.rm = TRUE), labeller = label_value) + 
        guides(fill = 'none') + ggtitle('Sample') +
        xlab(assayNames(object)[1]) + ylab(NULL) + 
        theme(axis.text.y  = element_blank(), 
              axis.ticks.y = element_blank(), 
              plot.title = element_text(hjust = 0.5))
}

plot_top_feature <- function(
    object, palette = make_subgroup_palette(object)
){
    plot_subgroup_boxplots(
        object, 
        subgroup = subgroup, 
        fill     = subgroup, 
        facet    = vars(feature_id),
        palette  = palette, 
        labeller = label_value) + 
    ggtitle('Feature') + guides(fill = 'none') + 
    ylab(assayNames(object)[1]) +
    theme(plot.title = element_text(hjust = 0.5))
}

top_coef <- function(object, fit){
    summarydt <- summarize_fit(object, fit = fit)
    summarydt %<>% extract(which.max(ndown))  # nup favours intercept !
    summarydt$contrast
}

plot_top_volcano <- function(object, fit, coef){
    summarydt <- summarize_fit(object, fit = fit)
    summarydt %<>% extract(contrast==coef)  # nup favours intercept !
    
    plot_volcano(object, coefs = coef) + 
    guides(color = 'none') + 
    ggtitle('Volcano') + 
    #ggtitle(sprintf('%d down | %d up', summarydt$ndown, summarydt$nup)) +
    theme(plot.title = element_text(hjust = 0.5))
}