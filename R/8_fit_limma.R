
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


#' Contrast Code Factor
#' 
#' Contrast Code Factor for General Linear Model
#'
#' @param object  factor vector
#' @param vars svars
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
#' @param verbose TRUE or FALSE
#' @param n character vector
#' @param ... used for s3 dispatch
#' @return (explicitly coded) factor vector
#' @details
#' A General Linear Model contains:                                                                   \cr
#'   * An Intercept Coefficient: expressing some form of sample average                               \cr
#'   * For each numeric variable: a slope coefficient                                                 \cr
#'   * For each k-leveled factor: (k-1) Contrast Coefficients.                                        \cr
#'        The interpretation of (intercept and contrast) coefficients depends on the contrast coding function used.
#'        Several contrast coding functions are available in 'stats' and 'codingMatrices'
#'        But their (function and coefficient) namings are a bit confusing and unsystematic.
#'        Instead, the functions below offer an intuitive interface (to the otherwise powerful stats/codingMatrices packages).
#'        The names of these functions reflect the contrast coding used (treatment, backward, sum, or helmert contrasts).
#'        They also reflect the intercept interpretation (either first factor's first level or grand mean).
#'        They all produce intuitive coefficient names (e.g. 't1-t0' rather than just 't1').
#'        They all have unit scaling (a coefficient of 1 means a backward of 1).
#' @examples
#' # Coding functions
#'     x <- factor(paste0('t', 0:3))
#'     xlevels <- levels(x)
#'     contr.treatment(         xlevels)
#'     contr.treatment.explicit(xlevels)
#'     contr.diff(              xlevels)
#'     code_control(            xlevels)
#'     code_diff(               xlevels)
#'     code_diff_forward(       xlevels)
#'     code_deviation(          xlevels)
#'     code_deviation_first(    xlevels)
#'     code_helmert(            xlevels)
#'     code_helmert_forward(    xlevels)
#' 
#' # Code
#'     x %<>% code(contr.treatment)
#'     x %<>% code(contr.treatment.explicit)
#'     x %<>% code(contr.diff)
#'     x %<>% code(code_control)
#'     x %<>% code(code_diff)
#'     x %<>% code(code_diff_forward)
#'     x %<>% code(code_deviation)
#'     x %<>% code(code_deviation_first)
#'     x %<>% code(code_helmert)
#'     x %<>% code(code_helmert_forward)
#'
#' # Model
#'     file <- download_data('atkin.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     object %<>% fit_limma(codingfun = contr.treatment) # default
#'     object %<>% fit_limma(codingfun = contr.treatment.explicit)
#'     object %<>% fit_limma(codingfun = contr.diff)
#'     object %<>% fit_limma(codingfun = code_control)
#'     object %<>% fit_limma(codingfun = code_diff)
#'     object %<>% fit_limma(codingfun = code_diff_forward)
#'     object %<>% fit_limma(codingfun = code_deviation)
#'     object %<>% fit_limma(codingfun = code_deviation_first)
#'     object %<>% fit_limma(codingfun = code_helmert)
#'     object %<>% fit_limma(codingfun = code_helmert_forward)
#' @export
code <- function(object, ...)  UseMethod('code')

#' @rdname code
#' @export
code.factor <- function(object, codingfun, verbose = TRUE, ...){
    if (is.null(codingfun))  return(object)
    assert_is_function(codingfun)
    k <- length(levels(object))
    contrasts(object) <- codingfun(levels(object))
    if (verbose){
        contrastmat <- codingMatrices::mean_contrasts(contrasts(object))
        colnames(contrastmat) <- levels(object)
        rownames(contrastmat)[1] <- 'Intercept'
        names(dimnames(contrastmat)) <- c('coefficient', 'level')
        message_df('\t\t\t\t%s', contrastmat)
    }
    object
}


#' @rdname code
#' @export
code.data.table <- function(object, codingfun, vars = names(object), verbose = TRUE, ...){
    if (verbose)  cmessage('\t\t code factors')
    for (var in vars){
        if (is.character(object[[var]]))  object[[var]] %<>% factor()
        if (is.logical(  object[[var]]))  object[[var]] %<>% factor()
    }
    if (is.null(codingfun)) return(object)
    for (var in vars){
        if (is.factor(object[[var]])){
            if (verbose)  cmessage('\t\t\t%s', var)
            object[[var]] %<>% code.factor(codingfun, verbose = verbose)
        }
    }
    object
}

#' @rdname code
#' @export
contr.treatment.explicit <- function(n){
    y <- contr.treatment(n)
    colnames(y) %<>% paste0('-', n[1])
    y
}


#' @rdname code
#' @export
code_control <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    codingMatrices::code_control(n, abbreviate = FALSE)
}


#' @rdname code
#' @export
contr.diff <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    codingMatrices::contr.diff(n, abbreviate = FALSE)
}

#' @rdname code
#' @export
code_diff <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    codingMatrices::code_diff(n, abbreviate = FALSE)
}

#' @rdname code
#' @export
code_diff_forward <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    codingMatrices::code_diff_forward(n, abbreviate = FALSE)
}

#' @rdname code
#' @export
code_deviation <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    k <- length(n)
    contrastnames <- paste0(n, collapse = '+')
    contrastnames <- paste0('(', contrastnames, ')')
    contrastnames <- paste0(contrastnames, '/', length(n))
    contrastnames <- paste0(n[-k], '-', contrastnames) 
    y <- codingMatrices::code_deviation(n)
    colnames(y) <- contrastnames
    y
}

#' @rdname code
#' @export
code_deviation_first <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    k <- length(n)
    contrastnames <- paste0(n, collapse = '+')
    contrastnames <- paste0('(', contrastnames, ')')
    contrastnames <- paste0(contrastnames, '/', length(n))
    contrastnames <- paste0(n[-1], '-', contrastnames) 
    y <- codingMatrices::code_deviation_first(n)
    colnames(y) <- contrastnames
    y
}

#' @rdname code
#' @export
code_helmert <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    y <- codingMatrices::code_helmert(n) # properly scaled version of stats::contr.helmert
    for (i in seq(2, ncol(y)+1)){
        curlevel <- n[i]
        prevlevels <- n[seq(1,i-1)]
        helmertmean <- paste0(prevlevels, collapse = '+')
        if (i>2)  helmertmean <- paste0('(', helmertmean, ')/', i-1)
        colnames(y)[i-1] <- paste0(curlevel, '-', helmertmean)
    }
    y
}

#' @rdname code
#' @export
code_helmert_forward <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    y <- codingMatrices::code_helmert_forward(n) # properly scaled version of stats::contr.helmert
    k <- length(n)
    for (i in seq(1, k-1)){
        curlevel <- n[i]
        nextlevels <- n[seq(i+1,k)]
        fwdmean <- nextlevels
        if (length(nextlevels)>1){
            fwdmean %<>% paste0(collapse = '+')
            fwdmean %<>% paste0('(', ., ')')
            fwdmean %<>% paste0('/', length(nextlevels))
        }
        colnames(y)[i] <- sprintf('%s-%s',  curlevel, fwdmean)
    }
    y
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
#' object <- read_metabolon(file)
#' subgroup_matrix(object, subgroupvar = 'subgroup')
#' contrast_subgroup_cols(object, subgroupvar = 'subgroup')
#' contrast_subgroup_rows(object, subgroupvar = 'subgroup')
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
contrvec2mat  <- function(contrasts)  matrix(
                    contrasts, nrow=1, dimnames=list("", contrasts))

contrmat2list <- function(contrasts)  list(colcontrasts = contrasts)

vectorize_contrasts <- function(contrasts){
    unname(unlist(lapply(contrasts, function(x) na.exclude(c(t(x))))))
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

#' Reset fit
#' @param object  SummarizedExperiment
#' @param fit     character vector
#' @param coefs   character vector
#' @param verbose TRUE or FALSE
#' @examples 
#' file <- download_data('atkin.metabolon.xlsx')
#' (object <- read_metabolon(file))
#' object %<>% reset_fit()
#' object %<>% fit_limma() %>% reset_fit()
#' object %<>% fit_limma() %>% fit_lm() %>% reset_fit()
#' object %<>% fit_limma() %>% fit_lm() %>% reset_fit('limma')
#' @export
reset_fit <- function(
    object, 
    fit     = fits(object), 
    coefs   = autonomics::coefs(object, fit = fit), 
    verbose = TRUE
){
# Assert
    . <- NULL
    assert_is_valid_sumexp(object)
    if (is.null(fits(object)))  return(object)
    assert_is_a_bool(verbose)
# Reset
    vars <- c('effect', 'p', 'fdr', 't')
    varpat  <- paste0(vars,  collapse = '|')
    coefpat <- paste0(coefs, collapse = '|')
    fitpat  <- paste0(fit,   collapse = '|')
        
    pattern <- sprintf('^(%s)%s(%s)%s(%s)$', varpat, FITSEP, coefpat, FITSEP, fitpat)
    cols <- grep(pattern, fvars(object), value = TRUE)
    if (length(cols)>0){
        if (verbose)  cmessage('\t\tRm from fdt: %s', pattern)
        for (col in cols)  fdt(object)[[col]] <- NULL
    }
# Return
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
merge_fit <- function(object, fitres, fit, statistic = NULL){
    . <- NULL
    fitresdt <- data.table::copy(fitres)   # dont change in original
    firstcols <- intersect(c('feature_id', 'Intercept'), names(fitresdt))
    fitresdt %<>% extract(,c(firstcols, sort(setdiff(names(.), firstcols))), with = FALSE)
    if (!is.null(statistic)) names(fitresdt)[-1] %<>% paste0(statistic,FITSEP,.)
    names(fitresdt)[-1] %<>% paste0(FITSEP, fit)
    object %<>% merge_fdt(fitresdt)
    object
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


