#=============================================================================
#
#             plot_contrastogram
#                 compute_connections
#
#=============================================================================


 true_names <- function(x) names(x)[ x]
false_names <- function(x) names(x)[!x]

fit_limma_contrastogram <- function(object, subgroupvar, design){
    colcontrasts <- contrast_subgroup_cols(object, subgroupvar)
    rowcontrasts <- contrast_subgroup_rows(object, subgroupvar)
    contrasts <-  c( c(t(colcontrasts)), c(t(rowcontrasts)))
    object %<>% fit_limma(design = design, contrasts = contrasts)
    object
}

compute_connections <- function(
    object, subgroupvar, design,
    colors = make_colors(slevels(object, subgroupvar), guess_sep(object))
){
# subgroup matrix, difference contrasts, limma
    fdrvalues <- fdr(object)
    effects <- effect(object)
    colnames(fdrvalues) %<>% split_extract_fixed(FITSEP, 1)
    colnames(effects) %<>% split_extract_fixed(FITSEP, 1)
    nsignif <- apply(fdrvalues < 0.05, 2, sum, na.rm=TRUE)
                #colSums( fdrvalues < 0.05, na.rm=TRUE)  # BREAKS ON SINGLE CONTR!
    nup     <- apply(fdrvalues < 0.05 & effects>0, 2, sum, na.rm=TRUE)
    ndown   <- apply(fdrvalues < 0.05 & effects<0, 2, sum, na.rm=TRUE)
# Create diagram
    sep <- guess_sep(object)
    subgroupmatrix <- subgroup_matrix(object, subgroupvar = subgroupvar)
    subgrouplevels <- c(t(subgroupmatrix))
    arrowsizes <- arrowcolors <- matrix(0,
        nrow = length(subgrouplevels), ncol = length(subgrouplevels),
        dimnames = list(subgrouplevels, subgrouplevels))
    arrowlabels <- matrix("0", nrow = nrow(arrowsizes), ncol = ncol(arrowsizes),
                        dimnames = dimnames(arrowsizes))
# Add contrast numbers
    contrastmat  <- makeContrasts(contrasts = coefs(object), levels = design)
    for (contrastname in colnames(contrastmat)){
        contrastvector <- contrastmat[, contrastname]
        to   <- true_names(contrastvector>0)
        from <- if (any(contrastvector<0)) true_names(contrastvector<0) else to
        ns <- nsignif[[contrastname]]
        nu <- nup[[contrastname]]
        nd <- ndown[[contrastname]]
        arrowsizes[ to, from] <- nu#ns
        arrowsizes[ from, to] <- nd#ns
        arrowcolors[to, from] <- colors[[to]]
        arrowcolors[from, to] <- colors[[from]]
        arrowlabels[to, from] <- if (nu>0) nu else 0
                            #paste0(nu,  " %up% phantom(.)") else "phantom(.)"
        arrowlabels[from, to] <- if (nd>0) nd else 0
                            #paste0(nd," %down% phantom(.)") else "phantom(.)"
    }
# Return
    #arrowlabels[arrowcolors==0] <- "0"
    list(arrowsizes = arrowsizes,
        arrowcolors = arrowcolors,
        arrowlabels = arrowlabels)
}


#' Plot contrastogram
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup svar
#' @param formula      formula
#' @param colors       named color vector (names = subgroups)
#' @param curve        arrow curvature
#' @return list returned by \code{\link[diagram]{plotmat}}
#' @examples
#' if (requireNamespace('diagram', quietly = TRUE)){
#'    file <- download_data('halama18.metabolon.xlsx')
#'    object <- read_metabolon(file, fit='limma', plot=FALSE)
#'    plot_contrastogram(object, subgroupvar = 'Group')
#' }
#' @export
plot_contrastogram <- function(
    object, 
    subgroupvar,
    formula = as.formula(paste0('~ 0 +', subgroupvar)),
    colors = make_colors(slevels(object, subgroupvar), guess_sep(object)),
    curve  = 0.1
){
# Initialize
    V2 <- N <- NULL
    if (!requireNamespace('diagram', quietly = TRUE)){
        stop("BiocManager::install('diagram'). Then re-run.") }
# Fit limma
    formula <- as.formula(paste0('~ 0 + ', subgroupvar))
    design <- create_design(object, formula = formula)
    colnames(design) %<>% stri_replace_first_regex(subgroupvar, '')
    object %<>% fit_limma_contrastogram(
        subgroupvar = subgroupvar, design = design)
# Compute connections
    contrastogram_matrices <- compute_connections(
        object, design = design, subgroupvar = subgroupvar, colors = colors)
    arrowsizes  <- contrastogram_matrices$arrowsizes
    arrowcolors <- contrastogram_matrices$arrowcolors
    arrowlabels <- contrastogram_matrices$arrowlabels
    widths <- scales::rescale(arrowsizes, c(0.01,30))
# Plot
    dt <- split_subgroup_levels(object, subgroupvar)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    if (all(nperrow==1)) nperrow %<>% length()
    # basedir <- file.path(tools::R_user_dir('autonomics', 'cache'), 'contrastogram')
    # dir.create(basedir)
    # pdf(file.path(basedir, 'directed_contrastogram.pdf'), 
    # width = 9, height = 9)
    arrowlabels %<>% as.data.frame()
    diagram::plotmat(A          = arrowlabels,
                    pos         = nperrow,
                    curve       = curve,
                    name        = rownames(arrowsizes),
                    relsize     = 1,
                    box.size    = 0.05,
                    box.col     = colors[rownames(arrowsizes)],
                    box.type    = 'square',
                    box.prop    = 0.8,
                    arr.lwd     = widths,
                    shadow.size = 0, # sqrt(arrowsizes)
                    arr.lcol    = arrowcolors,
                    arr.col     = arrowcolors,
                    arr.type    = 'triangle')
    #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
# Return
    object # limma!
}


