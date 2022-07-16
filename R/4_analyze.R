
#' Analyze
#' @param object       SummarizedExperiment
#' @param pca          whether to perform pca
#' @param fit          NULL, 'limma', 'lm', 'lme', 'lmer', or 'wilcoxon'
#' @param subgroupvar  subgroup svar
#' @param formula      model formula
#' @param block        block svar
#' @param weightvar    NULL or name of weight matrix in assays(object)
#' @param coefs        NULL or character vector: model coefficients to test
#' @param contrastdefs NULL or character vector: coefficient contrasts to test
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
#' object %<>% analyze(pca=TRUE, subgroupvar = 'Group', fit='limma')
#' @export
analyze <- function(
        object,
        pca          = TRUE,
        fit          = 'limma',
        subgroupvar  = default_subgroupvar(object),
        formula      = default_formula(object, subgroupvar, fit),
        block        = NULL,
        weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL,
        coefs        = colnames(create_design(object, formula = formula)),
        contrastdefs = NULL,
        verbose      = TRUE,
        plot         = pca & !is.null(fit),
        feature_id   = NULL,
        sample_id    = NULL,
        palette      = make_palette(object)
){
    # Analyze
    if (is.null(subgroupvar))  subgroupvar <- default_subgroupvar(object)
    subgroup <- if (is.null(subgroupvar))  quo(NULL) else sym(subgroupvar)
    if (pca)  object %<>% pca(verbose = verbose, plot = FALSE)
    for (curfit in fit){
        fitfun <- get(paste0('fit_', curfit))
        if (is.null(formula)) formula <- default_formula(object,subgroupvar,fit)
        if (is.null(coefs)) coefs <- colnames(create_design(object, formula = formula))
        object %<>% fitfun( subgroupvar  = subgroupvar,
                            formula      = formula,
                            coefs        = coefs,
                            contrastdefs = contrastdefs,
                            block        = block,
                            weightvar    = weightvar,
                            verbose      = verbose,
                            plot         = FALSE) 
    }
    # Plot/Return
    if (plot)  plot_summary(
        object, subgroupvar = subgroupvar, fit = fit, 
        feature_id = feature_id, sample_id = sample_id, palette = palette)
    object
}


#' Plot summary
#' @param object       SummarizedExperiment
#' @param subgroupvar  string
#' @param fit          string
#' @param coef         string
#' @param feature_id   NULL or string
#' @param sample_id    NULL or string
#' @param palette      NULL or color vector
#' @examples 
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot = FALSE)
#' plot_summary(object, subgroupvar = 'subgroup', fit = 'limma')
#' @export
plot_summary <- function(
    object, 
    subgroupvar = 'subgroup',
    fit         = fits(object)[1], 
    coef        = c(setdiff(coefs(object), 'Intercept'), 'Intercept')[1],
    feature_id  = NULL, 
    sample_id   = NULL, 
    palette     = make_palette(object)
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
    detections   <- plot_top_detections(object, subgroupvar = subgroupvar, palette = palette)
    sampledistr  <- plot_top_sample(object[, sample_id], 
                                subgroupvar = subgroupvar, palette = palette) + 
                    ggtitle(assayNames(object)[1]) + 
                    theme(plot.title = element_text(hjust = 0.5))
    featuredistr <- plot_top_feature(object[feature_id,])
    boxplot <- plot_top_boxplot(object[feature_id,], subgroupvar = subgroupvar, 
                    palette = palette)
    volcano <- plot_top_volcano(object, fit = fit, coef = coef)
    # Layout
    layout_matrix <- matrix(c(1,2,5,3,4,5), nrow = 2, byrow = 2)
    grid.draw(grid.arrange(
        detections, sampledistr, featuredistr,  boxplot, volcano, 
        layout_matrix = layout_matrix))
}

plot_top_detections <- function(
    object, subgroupvar, palette = make_palette(object)
){
    detections <- plot_summarized_detections(
                    object, 
                    subgroup = !!sym(subgroupvar), 
                    fill     = !!sym(subgroupvar), 
                    palette  = palette)
        # numeric colorscale cant be overridden
    pcaplot <-
        biplot(object, color = !!sym(subgroupvar), x = pca1, y = pca2, palette = palette) + 
        guides(color = 'none') + ggtitle(NULL) +
        theme(axis.text.x  = element_blank(), 
              axis.text.y  = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank()) +
        scale_x_continuous(expand = c(0.2, 0.2))
    
    detections + 
        annotation_custom(ggplotGrob(pcaplot), 
                          xmin = 0.2*max(detections$data$xmax), 
                          xmax = 0.8*max(detections$data$xmax), 
                          ymin = 0.1*max(detections$data$ymax), 
                          ymax = 0.7*max(detections$data$ymax))
}

plot_top_sample <- function(object, subgroupvar, palette = make_palette(object)){
    plot_sample_densities(
        object, 
        fill = !!sym(subgroupvar), facet = vars(sample_id), palette = palette,
        fixed = list(alpha=1)) + 
        guides(fill = 'none') + ggtitle(NULL) +
        xlab(NULL) + ylab(NULL) + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

plot_top_feature <- function(object){
    plot_feature_densities(
        object,
        facet = vars(feature_id), fixed = list(fill = 'grey80')) + 
        ggtitle(NULL) + xlab(NULL) + ylab(NULL) + 
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
        coord_flip()
}

plot_top_boxplot <- function(
    object, subgroupvar, palette = make_palette(object)
){
    plot_subgroup_boxplots(
        object, 
        subgroup = !!sym(subgroupvar), 
        fill     = !!sym(subgroupvar), 
        facet    = vars(feature_id),
        palette  = palette) + 
        ggtitle(NULL) + guides(fill = 'none') + ylab(NULL)
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
    ggtitle(sprintf('%d down | %d up', summarydt$ndown, summarydt$nup)) +
    theme(plot.title = element_text(hjust = 0.5))
}