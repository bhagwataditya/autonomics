
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
#' @return formula
#' @examples 
#' # Abundances
#'     file <- download_data('atkin.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     default_formula(object)
# # Ratios
#'     file <- download_data('billing16.proteingroups.txt')
#'     object <- read_maxquant_proteingroups(file)
#'     default_formula(object)
#' @export 
default_formula <- function(object){
    if ('subgroup' %in% svars(object)){
        if (singlelevel(object, 'subgroup')){    return(~1)
        } else if (contains_ratios(object)){      return(~0+subgroup)
        } else {                                  return(~subgroup)
        }
    } else {                                      return(~1)
    }
}

character2factor <- function(x)  if (is.character(x)) factor(x) else x


#' Create design
#'
#'  Create design matrix  for statistical analysis
#' @param object       SummarizedExperiment or data.frame
#' @param formula      formula with svars
#' @param drop         whether to drop predictor names
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
#' @param verbose      whether to message
#' @param ...          required to s3ify
#' @return design matrix
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' unique(create_design(object))
#' unique(create_design(object, ~ Time))
#' unique(create_design(object, ~ Time, codingfun = contr.treatment.explicit))
#' unique(create_design(object, ~ Time, codingfun = contr.diff))
#' unique(create_design(object, ~ Time + Diabetes))
#' unique(create_design(object, ~ Time / Diabetes))
#' unique(create_design(object, ~ Time * Diabetes))
#' @export
create_design <- function(object, ...) UseMethod('create_design')


#' @rdname create_design
#' @export
create_design.SummarizedExperiment <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)), 
    codingfun = contr.treatment,
    verbose   = TRUE, 
    ...
){
    create_design.data.table(sdt(object), 
                            formula   = formula,
                            codingfun = codingfun,
                            drop      = drop,
                            verbose   = verbose)
}

#' @rdname create_design
#' @export
create_design.data.table <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)), 
    codingfun = contr.treatment,
    verbose   = TRUE, 
    ...
){
# Assert
    assert_is_subset(all.vars(formula), names(object))
    . <- NULL
# Contrast Code Factors
    object %<>% code(codingfun = codingfun, vars = all.vars(formula), verbose = verbose)
# Create design matrix
    #if (verbose)   message('\t\tDesign: ', formula2str(formula))
    object %<>% data.frame(row.names = .$sample_id)
    myDesign <- model.matrix(formula, data = object)
    colnames(myDesign) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    is_factor_var <- function(x, object) is.factor(object[[x]])
    if (drop){
        for (predictor in all.vars(formula)){
            if (is.factor(object[[predictor]]))  colnames(myDesign) %<>% 
                        stri_replace_first_fixed(predictor, '') }
            # Fails for e.g. Diabetes = YES/NO: a meaningless column "YES" is created
            # For other cases it works wonderfully, so I keep it for now.
            # If it gives too many issues, roll back to doing the dropping only
            # for "subgroup" levels:
            #colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
    }
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


