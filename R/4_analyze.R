
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
#' @return SummarizedExperiment
#' @examples 
#' require(magrittr)
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' object %<>% analyze(pca=TRUE, subgroupvar = 'Group', fit='limma')
#' @export
analyze <- function(
        object,
        pca = TRUE,
        fit = 'limma',
        subgroupvar = default_subgroupvar(object),
        formula = default_formula(object, subgroupvar, fit),
        block = NULL,
        weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL,
        coefs    = colnames(create_design(object, formula = formula)),
        contrastdefs = contrast_coefs(object, formula),
        verbose = TRUE,
        plot    = pca & !is.null(fit)
){
    # Analyze
    if (is.null(subgroupvar))  subgroupvar <- default_subgroupvar(object)
    subgroup <- if (is.null(subgroupvar))  quo(NULL) else sym(subgroupvar)
    if (pca)  object %<>% pca(verbose = verbose, plot = FALSE)
    for (curfit in fit){
        fitfun <- get(paste0('fit_', curfit))
        if (is.null(formula)) formula <- default_formula(object,subgroupvar,fit)
        if (is.null(coefs)) coefs <- colnames(create_design(
            object, formula = formula))
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
    if (plot)  plot_summary(object, subgroupvar, fit)
    object
}


#' Plot summary
#' @param object SummarizedExperiment
#' @examples 
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot = FALSE)
#' plot_summary(object, subgroupvar = 'subgroup', fit = 'limma')
#' @export
plot_summary <- function(object, subgroupvar, fit){
    detections   <- plot_top_detections(object, subgroupvar)
    sampledistr  <- plot_top_sample( object, subgroupvar)
    featuredistr <- plot_top_feature(object)
    boxplot      <- plot_top_boxplot(object, subgroupvar)
    volcano      <- plot_top_volcano(object, fit)
    layout_matrix <- matrix(c(1,2,5,3,4,5), nrow = 2, byrow = 2)
    grid.draw(grid.arrange(detections,   sampledistr, 
                           featuredistr,  boxplot, volcano, 
                           layout_matrix = layout_matrix))
}

plot_top_detections <- function(object, subgroupvar){
    detections <- plot_summarized_detections(
        object, subgroup = !!sym(subgroupvar), fill = !!sym(subgroupvar))
    pcaplot <- biplot(object, color = !!sym(subgroupvar), x = pca1, y = pca2) + 
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
                          ymin = 0.2*max(detections$data$ymax), 
                          ymax = 0.8*max(detections$data$ymax))
}

plot_top_sample <- function(object, subgroupvar){
    palette <- make_colors(slevels(object, subgroupvar))
    idx <- which.min(abs(object$pca1))   # most average sample
    
    plot_sample_densities(
        object[, idx], 
        fill = !!sym(subgroupvar), facet = vars(sample_id), palette = palette,
        fixed = list(alpha=1)) + 
        guides(fill = 'none') + ggtitle(NULL) + coord_flip() + 
        xlab(NULL) + ylab(NULL) + 
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plot_top_feature <- function(object){
    idx <- which.max(abs(fdata(object)$pca1))
    
    plot_feature_densities(
        object[idx, ],
        facet = vars(feature_id), fixed = list(fill = 'grey80')) + 
        ggtitle(NULL) + xlab(NULL) + ylab(NULL) + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

plot_top_boxplot <- function(object, subgroupvar){
    idx <- which.max(abs(fdata(object)$pca1))
    palette <- make_colors(slevels(object, subgroupvar))
    plot_subgroup_boxplots(
        object[idx, ], 
        subgroup = subgroup, fill = subgroup, facet = vars(feature_id),
        palette = palette) + 
        ggtitle(NULL) + guides(fill = 'none') + ylab(NULL)
}

plot_top_volcano <- function(object, fit){
    selectedcoef <- summarize_fit(object, fit = fit)
    selectedcoef %<>% extract(which.max(ndown))  # nup favours intercept !
    selectedcoef %<>% extract2('contrast')
    selectedcoef %<>% as.character()
    plot_volcano(object, coefs = selectedcoef) + 
        guides(color = 'none') + ggtitle(NULL)
}