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

#' @rdname default_formula
#' @export
default_subgroupvar <- function(object){
    if ('subgroup' %in% svars(object))  'subgroup' else NULL
}

#' Create default formula
#' @param object SummarizedExperiment
#' @param subgroupvar string
#' @param fit 'limma', 'lm', 'lme', 'lmer'
#' @return formula
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <-.read_metabolon(file)
#' default_subgroupvar(object)
#' default_formula(object, fit = 'limma')
#' default_formula(object, fit = 'lm')
#' @export 
default_formula <- function(
    object, subgroupvar = default_subgroupvar(object), fit
){
    formula <- if (is.null(subgroupvar))        '~1'
    else if (!subgroupvar %in% svars(object))   '~1'
    else if (singlelevel(object, subgroupvar))  '~1'
    else if (fit %in% c('limma', 'wilcoxon'))   sprintf('~0 + %s', subgroupvar)
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
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' unique(create_design(object))
#'
#' object$subgroup <- 'billing19'
#' unique(create_design(object))
#'
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' unique(create_design(object))
#' create_design(object, formula= ~ 0 + SampleGroup + Sex + T2D + age + bmi)
#' object$subgroup <- 'atkin18'
#' unique(create_design(object))
#' @export
create_design <- function(
    object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar, fit = 'limma'), 
    verbose = TRUE
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(all.vars(formula), svars(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    for (var in all.vars(formula)){
        if (is.character(object[[var]])) object[[var]] %<>% factor() }
# Create design matrix
    if (verbose)   message('\t\tDesign: ', formula2str(formula))
    myDesign <- model.matrix(formula,  data = sdata(object))
    colnames(myDesign) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    is_factor_var <- function(x, object) is.factor(object[[x]])
    for (predictor in all.vars(formula)){
        if (is.factor(object[[predictor]]))  colnames(myDesign) %<>% 
                    stri_replace_first_fixed(predictor, '') %>% make.names() }
        # Fails for e.g. T2D = YES/NO: a meaningless column "YES" is created
        # For other cases it works wonderfully, so I keep it for now.
        # If it gives too many issues, roll back to doing the dropping only
        # for "subgroup" levels:
        #colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
# Return
    return(myDesign)
}


#=============================================================================
#
#               contrast_coefs
#                   contrast_subgroup_cols
#                   contrast_subgroup_rows
#
#==============================================================================


#' Row/Col contrasts
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup svar
#' @return  matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object, subgroupvar = 'Group')
#' contrast_subgroup_cols(object, subgroupvar = 'Group')
#' contrast_subgroup_rows(object, subgroupvar = 'Group')
#' @export
contrast_subgroup_cols <- function(object, subgroupvar){
    subgroupmat <- subgroup_matrix(object, subgroupvar)
    if (is_scalar(subgroupmat))  return(subgroupmat)
    if (ncol(subgroupmat)==1) return(matrix(, ncol=0, nrow=nrow(subgroupmat)))
    colcontrasts <- matrix(  sprintf('%s-%s',    # no space: as lm(contr.sdiff)
                                    subgroupmat[, -1],
                                    subgroupmat[, -ncol(subgroupmat)]),
                            nrow = nrow(subgroupmat),
                            ncol = ncol(subgroupmat)-1)
    rownames(colcontrasts) <- rownames(subgroupmat)
    colnames(colcontrasts) <- sprintf('%s-%s',   # no space: as lm(contr.sdiff)
                            colnames(subgroupmat)[-1],
                            colnames(subgroupmat)[-ncol(subgroupmat)])
    colcontrasts
}


#' @rdname contrast_subgroup_cols
#' @export
contrast_subgroup_rows <- function(object, subgroupvar){
    subgroupmat <- subgroup_matrix(object, subgroupvar)
    if (nrow(subgroupmat)==1) return(matrix(, nrow=0, ncol=ncol(subgroupmat)))
    rowcontrasts <- matrix(  sprintf('%s-%s',  # no space: as lm(contr.sdiff)
                                    subgroupmat[-nrow(subgroupmat), ],
                                    subgroupmat[-1, ]),
                            nrow = nrow(subgroupmat)-1,
                            ncol = ncol(subgroupmat))
    colnames(rowcontrasts) <- colnames(subgroupmat)
    rownames(rowcontrasts) <- sprintf('%s-%s', # no space: as lm(contr.sdiff)
                            rownames(subgroupmat)[-nrow(subgroupmat)],
                            rownames(subgroupmat)[-1])
    rowcontrasts
}


contrast_coefs <- function(object, formula){
    subgroupvar <- all.vars(formula)[1]
    design <- create_design(object, formula = formula, verbose = FALSE)
    if (ncol(design)==1){
        list(matrix(colnames(design), nrow=1, ncol=1), 
            matrix(nrow=0, ncol=0))
    } else if (all(design[, 1]==1)){
        list(colnames(design)[-1][seq_len(nlevels(object, subgroupvar)-1)], 
            matrix(nrow=0, ncol=0))
    } else {
        list(contrast_subgroup_cols(object, subgroupvar),
            contrast_subgroup_rows( object, subgroupvar)) }
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


.limmacontrast <- function(object, fit, formula){
    # compute contrasts
    design <- create_design(object, formula=formula, verbose = FALSE)
    contrastmat <- makeContrasts(
        contrasts = vectorize_contrastdefs(contrastdefs(object)),
        levels    = design)
    fit %<>% contrasts.fit(contrasts = contrastmat)
    limma_quantities <- if (all(fit$df.residual==0)){ c('effect', 'rank')
    } else { c('effect','rank','t','se','p','fdr','bonf') }
    limma(object) <- array( dim=c(nrow(fit),ncol(fit),length(limma_quantities)),
                            dimnames = list(feature  = rownames(fit),
                                            contrast = colnames(fit),
                                            quantity = limma_quantities))
    limma(object)[,,'effect'] <- fit$coefficients
    limma(object)[,,'rank'  ] <- apply(-abs(fit$coefficients), 2, rank)
    #names(dimnames(limma(object)))[2] <- formula2str(formula)
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


#' Fit model and test for differential expression
#'
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup variable
#' @param formula      modeling formula
#' @param contrastdefs contrastdef vector / matrix / list
#' \itemize{
#' \item{c("t1-t0", "t2-t1", "t3-t2")}
#' \item{matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE)}
#' \item{list(matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE), \cr
#'      matrix(c("KD.t0-WT.t0", "KD.t1-WT.t1", "KD.t2-WT.t2", "KD.t3-WT.t3"),\cr
#'      nrow=1, byrow=TRUE))}}
#' @param block     block svar (or NULL)
#' @param weightvar NULL or name of weight matrix in assays(object)
#' @param verbose   whether to msg
#' @param plot      whether to plot
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' object %<>% fit_limma(subgroupvar = 'SampleGroup')
#' object %<>% fit_lm(   subgroupvar = 'SampleGroup')
#' plot_venn(is_sig(object, contrast='t3-t2'))
#' 
#' S4Vectors::metadata(object)$limma <- S4Vectors::metadata(object)$lm <- NULL
#' object %<>% fit_limma(   subgroupvar = 'SampleGroup', block = 'Subject_ID')
#' object %<>% fit_wilcoxon(subgroupvar = 'SampleGroup', block = 'Subject_ID')
#' # object %<>% fit_lme(   subgroupvar = 'SampleGroup', block = 'Subject_ID')
#' # object %<>% fit_lmer(  subgroupvar = 'SampleGroup', block = 'Subject_ID')
#' plot_venn(is_sig(object, contrast='t3-t2'))
#' @export
fit_limma <- function(object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL,
    formula = default_formula(object, subgroupvar, 'limma'), 
    contrastdefs = contrast_coefs(object, formula), 
    block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    verbose = TRUE, plot=FALSE
){
# Design/contrasts
    assert_is_all_of(object, 'SummarizedExperiment')
    if (verbose)  message('\t\tlmFit(', formula2str(formula),
        if(is.null(block))     '' else paste0(' | ',block),
        if(is.null(weightvar)) '' else paste0(', weights = assays(object)$', 
                                            weightvar), ')')
    design <- create_design(object, formula=formula, verbose = FALSE)
    if (is.character(contrastdefs)) contrastdefs %<>% contrvec2mat()
    if (is.matrix(contrastdefs))    contrastdefs %<>% contrmat2list()
    contrastdefs(object) <- contrastdefs
# Block
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$dupcor)){
            if (verbose)  message('\t\t\t\tdupcor `', blockvar, '`')
            metadata(object)$dupcor <- duplicateCorrelation(
                values(object), design=design, block=block
            )$consensus.correlation }}
# Exprs/Weights
    exprmat <-  assays(object)[[1]][, rownames(design)]
    weightmat <- if (is.null(weightvar)){ NULL 
            } else {assert_is_a_string(weightvar)
                    assert_is_subset(weightvar, assayNames(object))
                    assays(object)[[weightvar]][, rownames(design)] }
# Fit
    fit <- suppressWarnings(lmFit(
                object = exprmat, design = design, 
                block = block, correlation = metadata(object)$dupcor,
                weights = weightmat))
# Contrast
    object %<>% .limmacontrast(fit, formula)
    if (plot)  print(plot_volcano(object, fit='limma')) 
    if (verbose)  message_df('\t\t\t%s', summarize_fit(object, 'limma'))
    return(object)
}


#' formula to string
#' @param formula formula
#' @return string
#' @examples 
#' formula2str(~0+subgroup)
#' @export
formula2str <- function(formula)  Reduce(paste, deparse(formula))


