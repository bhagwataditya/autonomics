#=============================================================================
#
#             plot_contrastogram
#                 compute_connections
#
#=============================================================================


true_names <- function(x) names(x[x])

compute_connections <- function(
    object, colors = make_colors(subgroup_levels(object), guess_sep(object))
){
# subgroup matrix, difference contrasts, limma
    pvalues <- limma(object)[, , 'p',      drop=FALSE]
    effects <- limma(object)[, , 'effect', drop=FALSE]
    nsignif <- apply(pvalues < 0.05, 2, sum, na.rm=TRUE)
                #colSums( pvalues < 0.05, na.rm=TRUE)  # BREAKS ON SINGLE CONTR!
    nup     <- apply(pvalues < 0.05 & effects>0, 2, sum, na.rm=TRUE)
                #colSums((pvalues < 0.05) & (effects > 0), na.rm=TRUE)
    ndown   <- apply(pvalues < 0.05 & effects<0, 2, sum, na.rm=TRUE)
                #colSums((pvalues < 0.05) & (effects < 0), na.rm=TRUE)
# Create diagram
    sep <- guess_sep(object)
    subgroupmatrix <- subgroup_matrix(object)
    subgrouplevels <- c(t(subgroupmatrix))
    arrowsizes <- arrowcolors <- matrix(0,
        nrow = length(subgrouplevels), ncol = length(subgrouplevels),
        dimnames = list(subgrouplevels, subgrouplevels))
    arrowlabels <- matrix("0", nrow = nrow(arrowsizes), ncol = ncol(arrowsizes),
                        dimnames = dimnames(arrowsizes))
# Add contrast numbers
    design <- create_design(
        object, formula=attr(limma(object), 'formula'), verbose = FALSE)
    colcontrasts <- contrastdefs(object)[[1]]
    rowcontrasts <- contrastdefs(object)[[2]]
    contrastdefs <-  c( c(t(colcontrasts)), c(t(rowcontrasts)))
    contrastmat  <- makeContrasts(contrasts = contrastdefs, levels = design)
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
#' @param object SummarizedExperiment
#' @param colors named color vector (names = subgroups)
#' @param curve  arrow curvature
#' @return list returned by \code{\link[diagram]{plotmat}}
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, fit='limma', plot=FALSE)
#' plot_contrastogram(object)
#' @export
plot_contrastogram <- function(
    object,
    colors = make_colors(subgroup_levels(object), guess_sep(object)),
    curve  = 0.1
){
# Initialize
    V2 <- N <- NULL
# Prepare
    contrastogram_matrices <- compute_connections(object, colors = colors)
    arrowsizes  <- contrastogram_matrices$arrowsizes
    arrowcolors <- contrastogram_matrices$arrowcolors
    arrowlabels <- contrastogram_matrices$arrowlabels
    widths <- scales::rescale(arrowsizes, c(0.01,30))
# Plot
    dt <- split_subgroup_levels(object)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    if (all(nperrow==1)) nperrow %<>% length()
    #dir.create('~/autonomicscache/contrastogram')
    #pdf('~/autonomicscache/contrastogram/directed_contrastogram.pdf',
    #width = 9, height = 9)
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
}


