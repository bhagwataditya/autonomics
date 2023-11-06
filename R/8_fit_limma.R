
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

mat2fdt <- function(mat)  mat2dt(mat, 'feature_id')


#' Fit model and test for differential expression
#'
#' @param object    SummarizedExperiment
#' @param formula   modeling formula
#' @param drop      TRUE or FALSE
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
#' @param design    design matrix
#' @param contrasts NULL or character vector: coefficient contrasts to test
#' @param coefs     NULL or character vector: model coefs to test
#' @param block     block svar (or NULL)
#' @param weightvar NULL or name of weight matrix in assays(object)
#' @param statvars  character vector: subset of c('effect', 'p', 'fdr', 't')
#' @param sep       string: pvar separator  ("~" in "p~t2~limma")
#' @param suffix    string: pvar suffix ("limma" in "p~t2~limma")
#' @param verbose   whether to msg
#' @param plot      whether to plot
#' @return Updated SummarizedExperiment
#' @examples
#' # Default
#'     file <- download_data('atkin.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     object %<>% fit(~ subgroup)
#' # Standard
#'     object %<>% fit_lm(        ~ subgroup)                 #     statistics default
#'     object %<>% fit_limma(     ~ subgroup)                 # bioinformatics default
#' # Blocked
#'     object %<>% fit_limma(     ~ subgroup, block = 'Subject')  #        simple random effects
#'     object %<>% fit_lme(       ~ subgroup, block = 'Subject')  #      powerful random effects
#'     object %<>% fit_lmer(      ~ subgroup, block = 'Subject')  # more powerful random effects
#' # Intuitive : alternative coding
#'     object %<>% fit_lme(       ~ subgroup, block = 'Subject', codingfun = contr.treatment.explicit)
#'     object %<>% fit_lmer(      ~ subgroup, block = 'Subject', codingfun = contr.treatment.explicit)
#'     object %<>% fit_limma(     ~ subgroup, block = 'Subject', codingfun = contr.treatment.explicit)
#' # Flexible : limma contrasts
#'     object %<>% fit_limma( ~ 0 + subgroup, block = 'Subject', contrasts = c('t1-t0'))
#'         # flexible, but only approximate
#'         # stat.ethz.ch/pipermail/bioconductor/2014-February/057682.html
#' # Non-parametric: wilcoxon
#'     object %<>% fit_wilcoxon( ~ subgroup)                # unpaired
#'     object %<>% fit_wilcoxon( ~ subgroup, block = 'Subject') # paired
#'     plot_contrast_venn(is_sig(object, contrast = 't2', fit = c('lm', 'limma')))
#'    #plot_contrast_venn(is_sig(object, contrast = 't3', fit = c('limma', 'lme')))
#' @export
fit <- function(
    object, 
    formula   = default_formula(object),
    engine    = 'limma', 
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment, 
    design    = create_design(object, formula = formula, drop = drop, codingfun = codingfun),
    contrasts = NULL,
    coefs     = if (is.null(contrasts))  colnames(design)     else NULL,
    block     = NULL,
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL,
    statvars  = c('effect', 'p', 'fdr'),
    sep       = FITSEP,
    suffix    = paste0(sep, 'limma'),
    verbose   = TRUE, 
    plot      = FALSE
){
    assert_scalar_subset(engine, c('limma', 'lme', 'lmer', 'wilcoxon', 'lm'))
    fitfun <- paste0('fit_', engine)
    get(fitfun)(object, formula = formula, drop = drop, codingfun = codingfun, 
                design = design, contrasts = contrasts, coefs = coefs, block = block, 
                weightvar = weightvar, statvars = statvars, sep = sep, suffix = suffix, 
                verbose = verbose, plot = plot)
}



#' @rdname fit
#' @export
fit_limma <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment,
    design    = create_design(object, formula = formula, drop = drop, codingfun = codingfun),
    contrasts = NULL,
    coefs     = if (is.null(contrasts))  colnames(design)     else NULL,
    block     = NULL,
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL,
    statvars  = c('effect', 'p', 'fdr'),
    sep       = FITSEP,
    suffix    = paste0(sep, 'limma'),
    verbose   = TRUE, 
    plot      = FALSE
){
    object %<>% reset_fit(fit = 'limma', coefs = coefs)
    limmadt <- .fit_limma(
        object       = object,        formula      = formula,
        drop         = drop,          codingfun    = codingfun,
        design       = design,        contrasts    = contrasts, 
        coefs        = coefs,         block        = block,
        weightvar    = weightvar,     statvars     = statvars,
        sep          = sep,           suffix       = suffix,
        verbose      = verbose)
    object %<>% merge_fdt(limmadt)
        #fdata(object)$F.limma   <- limmares$F
        #fdata(object)$F.p.limma <- limmares$F.p
    if (plot)  print(plot_volcano(object, fit = 'limma')) 
    object
}




#' Are varlevels unique
#' 
#' @param object SummarizedExperiment or data.table
#' @param vars character vector
#' @param ... required for s3 dispatch
#' @return TRUE or FALSE
#' @examples 
#' require(data.table)
#' object1 <- data.table(expand.grid(genome = c('WT', 'MUT'), treat = c('control', 'drug')))
#' object2 <- data.table(expand.grid(mutant = c('YES', 'NO'), treated = c('YES', 'NO')))
#' varlevels_dont_clash(object1)
#' varlevels_dont_clash(object2)
#' @export
varlevels_dont_clash <- function(object, ...)  UseMethod('varlevels_dont_clash')

#' @rdname varlevels_dont_clash
#' @export
varlevels_dont_clash.data.table <- function(
    object, vars = names(object), ...
){
    object                         %>% 
    extract(, vars, with = FALSE)  %>%
    lapply(factor)                 %>% 
    lapply(levels)                 %>% 
    unlist()                       %>% 
    duplicated()                   %>% 
    any()                          %>%
    magrittr::not()
}

#' @rdname varlevels_dont_clash
#' @export
varlevels_dont_clash.SummarizedExperiment <- function(
    object, vars = svars(object), ...
){
    varlevels_dont_clash.data.table(sdt(object), vars)
}


#' @rdname fit
#' @export
.fit_limma <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    codingfun    = 'treatment',
    design    = create_design(object, formula = formula, drop = drop, codingfun = codingfun),
    contrasts = NULL,
    coefs     = if (is.null(contrasts))  colnames(design) else NULL,
    block     = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    statvars  = c('effect', 'p', 'fdr'),
    sep       = FITSEP,
    suffix    = paste0(sep, 'limma'),
    verbose   = TRUE, 
    plot      = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_formula(formula)
    assert_is_subset(all.vars(formula), svars(object))
    assert_is_a_bool(drop)
    assert_is_matrix(design)
    assert_is_subset(coefs, colnames(design))
    if (!is.null(block))      assert_is_subset(block, svars(object))
    if (!is.null(weightvar))  assert_scalar_subset(weightvar, assayNames(object))
    assert_is_subset(statvars, c('effect', 'p', 'fdr', 't'))
# Design/contrasts/block/weights
    . <- NULL
    if (verbose)  message('\t\tlmFit(', formula2str(formula),
        if(is.null(block))     '' else paste0(' | ',block),
        if(is.null(weightvar)) '' else paste0(', weights = assays(object)$', 
                                            weightvar), ')')
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$dupcor)){
            if (verbose)  message('\t\t\t\tdupcor `', blockvar, '`')
            metadata(object)$dupcor <- duplicateCorrelation(
                values(object), design = design, block = block
            )$consensus.correlation }}
    design %<>% extract(intersect(snames(object), rownames(.)), , drop = FALSE) # required in mae
    exprmat <-  values(object)[, rownames(design)]
    weightmat <- if (is.null(weightvar)){ NULL 
            } else {assert_is_a_string(weightvar)
                    assert_is_subset(weightvar, assayNames(object))
                    assays(object)[[weightvar]][, rownames(design)] }
# Fit
    limmafit <- suppressWarnings(lmFit(
                    object = exprmat, design = design, 
                    block = block, correlation = metadata(object)$dupcor,
                    weights = weightmat))
# Effect
    if (is.null(contrasts)){  
                limmafit %<>% contrasts.fit(coefficients = coefs) 
    } else {    limmafit %<>% contrasts.fit(contrasts = makeContrasts(
                                    contrasts = contrasts, levels = design)) }
    limmadt <- data.table(feature_id = rownames(limmafit))
    if ('effect' %in% statvars){
        dt0 <- data.table(limmafit$coefficients)
        names(dt0) %<>% paste0('effect', sep, ., suffix)
        limmadt %<>% cbind(data.table(dt0)) 
    }
# p/t/fdr
    if (!all(limmafit$df.residual==0)){
        limmafit %<>% eBayes()
        if ('p' %in% statvars){ 
            dt0 <- data.table(limmafit$p.value)
            names(dt0) %<>% paste0('p', sep, ., suffix)
            limmadt %<>% cbind(dt0)  
        } 
        if ('t' %in% statvars){ 
            dt0 <- data.table(limmafit$t)
            names(dt0) %<>% paste0('t', sep, ., suffix)
            limmadt %<>% cbind(dt0)  
        } 
        if ('fdr' %in% statvars){
            dt0 <- data.table(apply(limmafit$p.value, 2, p.adjust, 'fdr') )
            names(dt0) %<>% paste0('fdr', sep, .,suffix)
            limmadt %<>% cbind(dt0)  
        } 
    }
# Return
    if (verbose)  message_df('\t\t\t%s',  summarize_fit(limmadt, fit = 'limma'))
    limmadt

}

old_summarize_fit <- function(object, fit = fits(object)){
    . <- NULL
    downdt <- colSums(downmat(object, fit = fit)) %>% data.table(coef = names(.), ndown = .)
    downdt %<>% tidyr::separate(
                    col = .data$coef, into = c('contrast', 'fit'), sep = FITSEP)
    
    updt <- colSums(upmat(object)) %>% data.table(coef = names(.), nup = .)
    updt %<>% tidyr::separate(
                    col = .data$coef, into = c('contrast', 'fit'), sep = FITSEP)
    
    sumdt <- merge(downdt, updt, by = c('fit', 'contrast'))
    sumdt$contrast %<>% factor()
    sumdt$contrast %<>% pull_level('Intercept')
    setorderv(sumdt, c('fit', 'contrast'))
    sumdt %<>% extract(fit, on = 'fit')
    sumdt %<>% extract(coefs(object), on = 'contrast')
    sumdt
}

pull_level <- function(x, lev){
    assert_is_factor(x)
    if (lev %in% levels(x))  x %<>% 
        factor(levels = c(lev, setdiff(levels(x), lev)))
    x
}

#' Summarize fit
#' @param fdt  fdt(object)
#' @param fit  'limma', 'lme', 'lm', 'lme', 'wilcoxon' or NULL
#' @param coefs string vector
#' @return data.table(contrast, nup, ndown)
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% fit_limma()
#' object %<>% fit_lm()
#' summarize_fit(fdt(object), coefs = c('t1', 't2', 't3'))
#' @export
summarize_fit <- function(fdt, fit = NULL, coefs = NULL){
    assert_is_data.table(fdt)
    statistic <- coefficient <- fit <- variable <- NULL
    cols <- names(fdt) %>% extract(stri_detect_fixed(., FITSEP))
    fdt %<>% extract(, c('feature_id', cols), with = FALSE)
    
    longdt <- fdt %>% melt.data.table(id.vars = 'feature_id')
    longdt[, statistic    := split_extract_fixed(variable, FITSEP, 1) %>% factor(unique(.))]
    longdt[,  coefficient := split_extract_fixed(variable, FITSEP, 2) %>% factor(unique(.))]
    longdt[,       fit    := split_extract_fixed(variable, FITSEP, 3) %>% factor(unique(.))]
    longdt[, variable := NULL]
    
    
    sumdt <- dcast.data.table(longdt, feature_id + coefficient + fit ~ statistic, value.var = 'value')
    sumdt <- sumdt[, .(
        downfdr = sum(effect < 0  & fdr < 0.05, na.rm = TRUE), 
        upfdr   = sum(effect > 0  & fdr < 0.05, na.rm = TRUE),
        downp   = sum(effect < 0  &   p < 0.05, na.rm = TRUE), 
        upp     = sum(effect > 0  &   p < 0.05, na.rm = TRUE)), by = c('coefficient', 'fit') ]
    if (!is.null(fit)){
        idx <- sumdt$fit %in% fit
        sumdt %<>% extract(idx)
    }
    if (!is.null(coefs)){
        sumdt <- sumdt[coefficient %in% coefs]
    }
    sumdt
}

#' Plot fit summary
#' @param sumdt data.table
#' @param nrow number
#' @param ncol number
#' @param order TRUE or FALSE
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% fit_lm()
#' object %<>% fit_limma(block = 'Subject')
#' sumdt <- summarize_fit(fdt(object), coefs = c('t1', 't2', 't3'))
#' plot_fit_summary(sumdt)
#' @export
plot_fit_summary <- function(sumdt, nrow = NULL, ncol = NULL, order = FALSE){
    coefficient <- downfdr <- downp <- fit <- upfdr <- upp <- NULL
    if (order){
        sumdt <- sumdt[order(downfdr+upfdr, downp+upp)]
        sumdt[, coefficient := factor(coefficient, unique(coefficient))]
    }
    ggplot(sumdt) + facet_wrap(vars(fit), nrow = nrow, ncol = ncol) + 
    geom_col(aes(y = coefficient, x = -downp),   fill = 'firebrick',   alpha = 0.3) +
    geom_col(aes(y = coefficient, x =    upp),   fill = 'forestgreen', alpha = 0.3) + 
    geom_col(aes(y = coefficient, x = -downfdr), fill = 'firebrick',   alpha = 1) +
    geom_col(aes(y = coefficient, x =    upfdr), fill = 'forestgreen', alpha = 1) + 
    geom_text(data = sumdt[  downp>0], aes(y = coefficient, x = -max(downp), label = paste0(downp, ' | ', downfdr) ), hjust = +1) + 
    geom_text(data = sumdt[    upp>0], aes(y = coefficient, x =    max(upp), label = paste0(upfdr, ' | ', upp) ), hjust = 0) + 
    xlab('count') + 
    ylab(NULL) + 
    xlim(c(-max(sumdt$downp)-100, max(sumdt$upp)+100)) + 
    #scale_x_continuous(n.breaks = 20) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


setna <- function(dt, value){
   for (j in seq_len(ncol(dt)))  set(dt, which(is.na(dt[[j]])), j, value)
    dt
}


#' formula to string
#' @param formula formula
#' @return string
#' @examples 
#' formula2str(~0+subgroup)
#' @export
formula2str <- function(formula)  Reduce(paste, deparse(formula))


