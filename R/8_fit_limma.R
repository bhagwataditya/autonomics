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
    #else if (fit %in% c('limma', 'wilcoxon'))   sprintf('~0 + %s', subgroupvar)
    else if (fit %in% c('wilcoxon'))            sprintf('~0 + %s', subgroupvar)
    else                                        sprintf('~ %s', subgroupvar)
    formula %<>% as.formula()
    formula
}

character2factor <- function(x)  if (is.character(x)) factor(x) else x

create_design <- function(object, ...)  UseMethod('create_design')


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
create_design <- function(object, ...) UseMethod('create_design')

#' @rdname create_design
#' @export
create_design.SummarizedExperiment <- function(
    object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar, fit = 'limma'),
    verbose = FALSE, ...
){
    sdt <- data.table(sdata(object))
    create_design.data.table(
        sdt, subgroupvar=subgroupvar, formula=formula, verbose=verbose)
}

#' @rdname create_design
#' @export
create_design.data.table <- function(
    object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL,
    formula = default_formula(object, subgroupvar, fit = 'limma'),
    verbose = FALSE, ...
){
# Assert
    assert_is_subset(all.vars(formula), names(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    for (var in all.vars(formula)){
        if (is.character(object[[var]])) object[[var]] %<>% factor() }
# Create design matrix
    if (verbose)   message('\t\tDesign: ', formula2str(formula))
    myDesign <- model.matrix(formula, data=object)
    colnames(myDesign) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    is_factor_var <- function(x, object) is.factor(object[[x]])
    for (predictor in all.vars(formula)){
        if (is.factor(object[[predictor]]))  colnames(myDesign) %<>% 
                    stri_replace_first_fixed(predictor, '') }
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


# contrast_coefs <- function(object, formula){
#     subgroupvar <- all.vars(formula)[1]
#     design <- create_design(object, formula = formula)
#     if (ncol(design)==1){
#         list(matrix(colnames(design), nrow=1, ncol=1), 
#             matrix(nrow=0, ncol=0))
#     } else if (all(design[, 1]==1)){
#         list(colnames(design)[-1][seq_len(nlevels(object, subgroupvar)-1)], 
#             matrix(nrow=0, ncol=0))
#     } else {
#         list(contrast_subgroup_cols(object, subgroupvar),
#             contrast_subgroup_rows( object, subgroupvar)) }
# }

contrast_coefs <- function(object, formula){
    subgroupvar <- all.vars(formula)[1]
    design <- create_design(object, formula = formula)
    if (ncol(design)==1){  
        colnames(design)
    } else if (all(design[, 1]==1)){ 
        colnames(design)[-1][seq_len(nlevels(object, subgroupvar)-1)]
    } else {
        c(contrast_subgroup_cols(object, subgroupvar), 
          contrast_subgroup_rows( object, subgroupvar))
    }
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

add_fdr <- function(fitres){
    . <- NULL
    fdr <- fitres %>% extract(, stri_startswith_fixed(
                                    names(.), paste0('p', FITSEP)), with=FALSE)
    fdr[] %<>% lapply(p.adjust, method='fdr')
    names(fdr) %<>% stri_replace_first_fixed(
                paste0('p', FITSEP), paste0('fdr', FITSEP))
    fitres %<>% cbind(fdr)
    fitres
}

reset_fitres <- function(object, fit){
    . <- NULL
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_a_string(fit)
    fitcols <- fvars(object) %>% extract(stri_detect_fixed(., fit))
    fdata(object)[fitcols] <- NULL
    metadata(object)[[fit]] <- NULL
    object
}

#' Fit results separator
#' @examples
#' FITSEP
#' @export
FITSEP <- '~'

# object: SumExp
# fitres: data.table(p.contr1, p.contr2, effect.contr1, effect.contr2)
# stat:  'p', 'effect', 'fdr', 't'
# fit:   'limma', 'wilcoxon'
merge_fitres <- function(object, fitres, fit, statistic=NULL){
    . <- NULL
    fitresdt <- data.table::copy(fitres)   # dont change in original
    firstcols <- intersect(c('feature_id', 'Intercept'), names(fitresdt))
    fitresdt %<>% extract(,c(firstcols, 
                            sort(setdiff(names(.), firstcols))), with=FALSE)
    if (!is.null(statistic)) names(fitresdt)[-1] %<>% paste0(statistic,FITSEP,.)
    names(fitresdt)[-1] %<>% paste0(FITSEP, fit)
    object %<>% merge_fdata(fitresdt)
    object
}

mat2fdt <- function(mat)  mat2dt(mat, 'feature_id')

.limmacontrast <- function(object, fit, design, coefs, contrastdefs){
    object %<>% reset_fitres('limma')
    if (is.null(contrastdefs)){ 
                fit %<>% contrasts.fit(coefficients = coefs) 
    } else {    fit %<>% contrasts.fit(contrasts    = makeContrasts(
                            contrasts = contrastdefs, levels = design)) }
    limma_quantities <- if (all(fit$df.residual==0)){ 'effect'
                        } else { c('effect','fdr', 'p', 't') }
    limma(object) <- array( dim=c(nrow(fit),ncol(fit),length(limma_quantities)),
                            dimnames = list(feature  = rownames(fit),
                                            contrast = colnames(fit),
                                            quantity = limma_quantities))
    limma(object)[,,'effect'] <- fit$coefficients
    object %<>% merge_fitres(
                    mat2fdt(fit$coefficients), statistic='effect', fit='limma')
    # perform moderated t test
    if (!all(fit$df.residual==0)){
        fit %<>% eBayes()
        tvalues <- fit$t
        pvalues <- fit$p.value
        fdrvals <- apply(pvalues, 2, p.adjust, 'fdr')
        limma(object)[,,'p'  ] <- pvalues
        limma(object)[,,'fdr'] <- fdrvals
        limma(object)[,,'t'  ] <- tvalues
        # limma(object)[,,'se'  ] <- sqrt(fit$s2.post) * fit$stdev.unscaled
        # limma(object)[,,'rank'] <- apply(fit$p.value, 2, rank)
        #object %<>% merge_contrastmat(fit$t,       't')
        object %<>% merge_fitres(mat2fdt(fdrvals), statistic='fdr', fit='limma')
        object %<>% merge_fitres(mat2fdt(pvalues), statistic='p',   fit='limma')
        object %<>% merge_fitres(mat2fdt(tvalues), statistic='t',   fit='limma')
        #fdata(object)$F.limma   <- fit$F
        #fdata(object)$F.p.limma <- fit$F.p
    }
    return(object)
}


#' Fit model and test for differential expression
#'
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup variable
#' @param formula      modeling formula
#' @param coefs        NULL or character vector: model coefficients to test
#' @param contrastdefs NULL or character vector: coefficient contrasts to test
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
#' # classical: lm & limma
#'     require(magrittr)
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot=FALSE, impute=TRUE)
#'     object$SET %<>% factor(); object$SUB %<>% factor()
#'     object %<>% fit_limma(subgroupvar = 'SET')
#'     object %<>% fit_lm(subgroupvar = 'SET')
#'     plot_venn(testmat(object, coef = 't3'))
#'     
#' # blocked: limma, lme, lmer
#'     object %<>% fit_limma(subgroupvar = 'SET', block = 'SUB')
#'     object %<>% fit_lme(subgroupvar = 'SET', block = 'SUB')
#'     # object %<>% fit_lme(subgroupvar = 'SET', block = 'SUB') # slow
#'     plot_venn(testmat(object, coef = 't3', fit = c('limma', 'lme')))
#'     
#' # flexible: limma contrastdefs
#'     object %<>% fit_limma(formula=~SET,   coef='t3',            block='SUB')
#'     object %<>% fit_limma(formula=~0+SET, contrastdefs='t3-t0', block='SUB')
#'         # flexible, but only approximate
#'         # stat.ethz.ch/pipermail/bioconductor/2014-February/057682.html
#' 
#' # alternative coding: 
#'     stats::contrasts(object$SET) <- MASS::contr.sdif(levels(object$SET))
#'     object %<>% fit_limma(subgroupvar = 'SET', block = 'SUB') # backward difs
#'     stats::contrasts(object$SET) <-stats::contr.treatment(levels(object$SET))
#'     object %<>% fit_limma(subgroupvar = 'SET', block = 'SUB') # baseline difs
#'     
#' # non-parametric: wilcoxon
#'     object %<>% fit_limma(subgroupvar = 'SET', block = 'SUB')
#'     object %<>% fit_wilcoxon(subgroupvar='SET', block='SUB', 
#'                    contrastdefs=c('t1-t0', 't2-t0', 't3-t0'))
#'     plot_venn(testmat(object,coef=c('t3','t3-t0'),fit=c('limma','wilcoxon')))
#' @export
fit_limma <- function(object, 
    subgroupvar  = if ('subgroup' %in% svars(object)) 'subgroup' else NULL,
    formula      = default_formula(object, subgroupvar, 'limma'), 
    coefs        = colnames(create_design(object, formula=formula)), 
    contrastdefs = NULL,
    block        = NULL, 
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    verbose = TRUE, plot=FALSE
){
# Design/contrasts
    assert_is_all_of(object, 'SummarizedExperiment')
    if (verbose)  message('\t\tlmFit(', formula2str(formula),
        if(is.null(block))     '' else paste0(' | ',block),
        if(is.null(weightvar)) '' else paste0(', weights = assays(object)$', 
                                            weightvar), ')')
    design <- create_design(object, formula=formula)
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
    exprmat <-  values(object)[, rownames(design)]
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
    object %<>% .limmacontrast(fit, design, coefs = coefs, 
                                contrastdefs = contrastdefs)
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


