#==============================================================================
#
#                           getters / setters
#
#==============================================================================


#' @title Get/set limma results
#' @description Get/Set limma results
#' @param object SummarizedExperiment
#' @param value list
#' @return limma results (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'            file, invert_subgroups=inv, fit='limma', plot=FALSE)
#' dim(limma(object))
#' dim(limma(object[1:5, ]))
#' @export
setGeneric("limma", function(object)   standardGeneric("limma") )


#' @rdname limma
setMethod("limma", signature("SummarizedExperiment"),
function(object){
    limma_array <- metadata(object)$limma
    if (is.null(limma_array)) NULL else limma_array[
                                            fnames(object), , , drop=FALSE] })

#' @rdname limma
#' @export
setGeneric("limma<-", function(object, value)  standardGeneric("limma<-") )


#' @rdname limma
setReplaceMethod("limma", signature("SummarizedExperiment", "array"),
function(object, value){
    metadata(object)$limma <- value
    object  })


#' @rdname limma
setReplaceMethod("limma", signature("SummarizedExperiment", "NULL"),
function(object, value) object)


#' @title Get/set design matrix
#' @param object SummarizedExperiment
#' @param value list
#' @return design (get) or SummarizedExperiment (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'               file, invert_subgroups=inv, fit='limma', plot=FALSE)
#' design(object)
#' @export
setGeneric("design", function(object)   standardGeneric("design") )


#' @rdname design
#' @export
setGeneric("design<-",
function(object, value)  standardGeneric("design<-") )


#' @rdname design
setMethod("design", signature("SummarizedExperiment"),
function(object) metadata(object)$design)


#' @rdname design
setReplaceMethod("design", signature("SummarizedExperiment", "matrix"),
function(object, value){
    metadata(object)$design <- value
    object  })


#' @title Get/set contrastdefs
#' @param object SummarizedExperiment
#' @param value list
#' @return contrastdefs (get) or SummarizedExperiment (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'             file, invert_subgroups=inv, fit='limma', plot=FALSE)
#' contrastdefs(object)
#' @export
setGeneric("contrastdefs", function(object)   standardGeneric("contrastdefs") )


#' @rdname contrastdefs
setMethod("contrastdefs", signature("SummarizedExperiment"),
function(object) metadata(object)$contrastdefs)


#' @rdname contrastdefs
#' @export
setGeneric("contrastdefs<-",
function(object, value)  standardGeneric("contrastdefs<-") )


#' @rdname contrastdefs
setReplaceMethod("contrastdefs", signature("SummarizedExperiment", "list"),
function(object, value){
    metadata(object)$contrastdefs <- value
    object  })

#==============================================================================
#
#                       create_design
#                           single_subgroup
#                           are_factor
#                           singlelevel
#                           multilevel
#                               nlevels
#
#==============================================================================

single_subgroup <- function(object){
    assert_is_subset('subgroup', svars(object))
    length(unique(object$subgroup))==1
}


are_factor <- function(df) vapply(df, is.factor, logical(1))

nlevels <- function(object, svar){
    if (!svar %in% svars(object))  return(0)
    length(unique(object[[svar]]))
}

singlelevel <- function(object, svar)   nlevels(object, svar) ==1
multilevel  <- function(object, svar)   nlevels(object, svar) > 1

#' 
default_formula <- function(object, subgroupvar, fit){
    formula <- if (is.null(subgroupvar))        '~1'
    else if (!subgroupvar %in% svars(object))   '~1'
    else if (singlelevel(object, subgroupvar))  '~1'
    else if (fit=='limma')                      sprintf('~0 + %s', subgroupvar)
    else                                        sprintf('~ %s', subgroupvar)
    formula %<>% as.formula()
    formula
}

#' Create design
#'
#'  Create design matrix  for statistical analysis
#' @param object  SummarizedExperiment
#' @param subgroupvar subgroup svar
#' @param formula formula with svars
#' @param verbose whether to message
#' @return design matrix
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' create_design(object)
#'
#' object$subgroup <- 'billing16'
#' create_design(object)
#'
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' create_design(object)
#' create_design(object, ~ 0 + subgroup + Sex + T2D + age + bmi)
#' object$subgroup <- 'atkin18'
#' create_design(object)
#' @export
create_design <- function(
    object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar), verbose = TRUE
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(all.vars(formula), svars(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    for (var in setdiff(all.vars(formula), subgroupvar)){
        if (is.character(object[[var]])) object[[var]] %<>% factor() }
# Create design matrix
    if (verbose)  cmessage('\t\tDesign: %s', deparse(formula))
    myDesign <- model.matrix(formula,  data = sdata(object))
# Rename columns
    colnames(myDesign) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    factors <- svars(object)[are_factor(sdata(object))] # Rename intercept
    for (var in factors) colnames(myDesign) %<>% gsub(var, '', ., fixed = TRUE)
        # Fails for e.g. T2D = YES/NO: a meaningless column "YES" is created
        # For other cases it works wonderfully, so I keep it for now.
        # If it gives too many issues, roll back to doing the dropping only
        # for "subgroup" levels:
        #colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
# Validify names
    colnames(myDesign) %<>% gsub(':', '..', ., fixed = TRUE)
    colnames(myDesign) %<>% make.names()
# Return
    return(myDesign)
}


#=============================================================================
#
#               contrast_subgroups
#                   contrast_subgroup_cols
#                   contrast_subgroup_rows
#
#==============================================================================


#' Row/Col contrasts
#' @param object SummarizedExperiment
#' @return  matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object)
#' contrast_subgroup_cols(object)
#' contrast_subgroup_rows(object)
#' @export
contrast_subgroup_cols <- function(object){
    subgroupmat <- subgroup_matrix(object)
    if (is_scalar(subgroupmat))  return(subgroupmat)
    if (ncol(subgroupmat)==1) return(matrix(, ncol=0, nrow=nrow(subgroupmat)))
    colcontrasts <- matrix(  sprintf('%s - %s',
                                    subgroupmat[, -1],
                                    subgroupmat[, -ncol(subgroupmat)]),
                            nrow = nrow(subgroupmat),
                            ncol = ncol(subgroupmat)-1)
    rownames(colcontrasts) <- rownames(subgroupmat)
    colnames(colcontrasts) <- sprintf('%s - %s',
                            colnames(subgroupmat)[-1],
                            colnames(subgroupmat)[-ncol(subgroupmat)])
    colcontrasts
}


#' @rdname contrast_subgroup_cols
#' @export
contrast_subgroup_rows <- function(object){
    subgroupmat <- subgroup_matrix(object)
    if (nrow(subgroupmat)==1) return(matrix(, nrow=0, ncol=ncol(subgroupmat)))
    rowcontrasts <- matrix(  sprintf('%s - %s',
                                    subgroupmat[-nrow(subgroupmat), ],
                                    subgroupmat[-1, ]),
                            nrow = nrow(subgroupmat)-1,
                            ncol = ncol(subgroupmat))
    colnames(rowcontrasts) <- colnames(subgroupmat)
    rownames(rowcontrasts) <- sprintf('%s - %s',
                            rownames(subgroupmat)[-nrow(subgroupmat)],
                            rownames(subgroupmat)[-1])
    rowcontrasts
}


contrast_subgroups <- function(object, design){
    if (ncol(design)==1){
        list(matrix(colnames(design), nrow=1, ncol=1), 
            matrix(nrow=0, ncol=0))
    } else {
        list(contrast_subgroup_cols(object),
            contrast_subgroup_rows(object)) }
}



#==============================================================================
#
#                            fit_limma
#
#==============================================================================

contrvec2mat  <- function(contrastdefs)  matrix(
                    contrastdefs, nrow=1, dimnames=list("", contrastdefs))

contrmat2list <- function(contrastdefs)  list(colcontrasts = contrastdefs)

vectorize_contrastdefs <- function(contrastdefs){
    unname(unlist(lapply(contrastdefs, function(x) na.exclude(c(t(x))))))
}

#' Fit model and test for differential expression
#'
#' Limma results can be easily accessed with limma(object).
#' @param object       SummarizedExperiment
#' @param contrastdefs contrastdef vector / matrix / list
#' \itemize{
#' \item{c("t1-t0", "t2-t1", "t3-t2")}
#' \item{matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE)}
#' \item{list(matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE), \cr
#'      matrix(c("KD.t0-WT.t0", "KD.t1-WT.t1", "KD.t2-WT.t2", "KD.t3-WT.t3"),\cr
#'      nrow=1, byrow=TRUE))}}
#' @param formula   designmat formula
#' @param block     block svar (or NULL)
#' @param verbose   whether to msg
#' @param plot      whether to plot
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' fit_wilcoxon(object)
#' fit_lm(object)
#' fit_limma(object)
#' fit_limma(object, block = 'Subject_ID')
#' fit_lme(  object, block = 'Subject_ID')
#' fit_lmer( object, block = 'Subject_ID')
#'
#' file <- download_data('billing19.proteingroups.txt')
#' select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
#' select %<>% paste0('_STD')
#' object <- read_proteingroups(file, select_subgroups = select, plot = FALSE)
#' object %<>% impute_systematic_nondetects(plot=FALSE)
#' object %<>% fit_limma()
#' object %<>% fit_lm()
#'
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE, voom=TRUE)
#' object %<>% fit_limma(plot=FALSE)
#' weights(object) <- NULL; object %<>% fit_limma(plot=FALSE)
#' object %<>% fit_lm(plot=FALSE)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% impute_systematic_nondetects(plot=FALSE)
#' object %<>% fit_limma(plot=FALSE)
#' object %<>% fit_lm(plot=FALSE)
#' @export
fit_limma <- function(object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar, 'limma'), 
    contrastdefs = NULL, block = NULL, verbose = TRUE, plot=FALSE
){
# Set design/contrasts
    assert_is_all_of(object, 'SummarizedExperiment')
    if (verbose)  cmessage('\t\tlimma: lmFit(%s%s%s)',
        Reduce(paste, deparse(formula)),
        if(is.null(block))            '' else paste0(', block = object$ `', block, '`'),
        if(is.null(weights(object)))  '' else paste0(', weights = weights(object)'))
    design <- create_design(object, formula=formula, verbose = FALSE)
    design(object)    <- design
    if (is.null(contrastdefs)) contrastdefs <- contrast_subgroups(object,design)
    if (is.character(contrastdefs)) contrastdefs %<>% contrvec2mat()
    if (is.matrix(contrastdefs))    contrastdefs %<>% contrmat2list()
    contrastdefs(object) <- contrastdefs
# Block
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$blockcor)){
            if (verbose)  cmessage('\t\t\t\tdupcor `%s`', blockvar)
            metadata(object)$blockcor <- duplicateCorrelation(
                exprs(object), design=design, block=block
            )$consensus.correlation }}
# Lmfit
    fit <- suppressWarnings(lmFit(object = exprs(object[, rownames(design)]),
        design = design, block = block, correlation = metadata(object)$blockcor,
        weights = weights(object)))
# Contrast
    object %<>% .limmacontrast(fit)
    if (plot)  print(plot_volcano(object, fit='limma')) 
                    # plot_contrastogram(object)
    if (verbose) cmessage_df('\t\t\t%s', summarize_fit(object, 'limma'))
    return(object)
}


.limmacontrast <- function(object, fit){
    # compute contrasts
    contrastmat <- makeContrasts(
        contrasts = vectorize_contrastdefs(contrastdefs(object)),
        levels    = design(object))
    fit %<>% contrasts.fit(contrasts = contrastmat)
    limma_quantities <- if (all(fit$df.residual==0)){ c('effect', 'rank')
    } else { c('effect','rank','t','se','p','fdr','bonf') }
    limma(object) <- array( dim=c(nrow(fit),ncol(fit),length(limma_quantities)),
                            dimnames = list(feature  = rownames(fit),
                                            contrast = colnames(fit),
                                            quantity = limma_quantities))
    limma(object)[,,'effect'] <- fit$coefficients
    limma(object)[,,'rank'  ] <- apply(-abs(fit$coefficients), 2, rank)
    # perform moderated t test
    if (!all(fit$df.residual==0)){
        fit %<>% eBayes()
        pp <- fit$p.value
        limma(object)[,,'t' ] <- fit$t
        limma(object)[,,'se'] <- sqrt(fit$s2.post) * fit$stdev.unscaled
        limma(object)[,,'p' ] <- pp
        limma(object)[,,'rank'] <- apply(pp, 2, rank)
        limma(object)[,,'fdr' ] <- apply(pp, 2, p.adjust, 'fdr')
        limma(object)[,,'bonf'] <- apply(pp, 2, p.adjust, 'bonferroni')
        fdata(object)$F.limma   <- fit$F
        fdata(object)$F.p.limma <- fit$F.p
    }
    return(object)
}
