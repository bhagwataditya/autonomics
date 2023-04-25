
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
#'     file <- download_data('atkin18.metabolon.xlsx')
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
#' @param coding    contrast coding function
#' @param verbose      whether to message
#' @param ...          required to s3ify
#' @return design matrix
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' unique(create_design(object))
#' unique(create_design(object, ~ subgroup))
#' unique(create_design(object, ~ subgroup, coding = 'reference' ))
#' unique(create_design(object, ~ subgroup, coding = 'difference'))
#' unique(create_design(object, ~ subgroup + T2D))
#' unique(create_design(object, ~ subgroup / T2D))
#' unique(create_design(object, ~ subgroup * T2D))
#' @export
create_design <- function(object, ...) UseMethod('create_design')


#' @rdname create_design
#' @export
create_design.SummarizedExperiment <- function(
    object, 
    formula = default_formula(object),
    drop    = varlevels_dont_clash(object, all.vars(formula)), 
    coding  = 'treatment',
    verbose = TRUE, 
    ...
){
    create_design.data.table(sdt(object), 
                            formula = formula,
                            coding  = coding,
                            drop    = drop,
                            verbose = verbose)
}

#' @rdname create_design
#' @export
create_design.data.table <- function(
    object, 
    formula = default_formula(object),
    drop    = varlevels_dont_clash(object, all.vars(formula)), 
    coding  = 'treatment',
    verbose = TRUE, 
    ...
){
# Assert
    assert_is_subset(all.vars(formula), names(object))
    . <- NULL
# Contrast Code Factors
    object %<>% code(coding = coding, vars = all.vars(formula), verbose = verbose)
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
            # Fails for e.g. T2D = YES/NO: a meaningless column "YES" is created
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
#' @param coding  coding system \cr
#' @param vars    character vector
#'    'treatment'  :  coef = level - firstlevel, intercept = firstlevel, implicit coefnames (t2).   \cr
#'    'reference'  :  coef = level - firstlevel, intercept = firstlevel, explicit coefnames (t2-t0) \cr
#'    'difference' :  coef = level -  prevlevel, intercept = firstlevel, epxlicit coefnames (t2-t1) \cr
#'    'grandref'   :  coef = level - firstlevel, intercept = grandmean   \cr
#'    'granddiff'  :  coef = level -  prevlevel, intercept = grandmean   \cr
#'    'sum'        :  coef = level -  grandmean, intercept = grandmean   \cr
#'    'helmert'    :  coef = level - prevlevels, intercept = grandmean   \cr
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
#'        The names of these functions reflect the contrast coding used (treatment, difference, sum, or helmert contrasts).
#'        They also reflect the intercept interpretation (either first factor's first level or grand mean).
#'        They all produce intuitive coefficient names (e.g. 't1-t0' rather than just 't1').
#'        They all have unit scaling (a coefficient of 1 means a difference of 1).
#' @examples
#' # Coding functions
#'     x <- factor(paste0('t', 0:3))
#'     code_reference( levels(x))
#'     code_grandref(  levels(x))
#'     code_difference(levels(x))
#'     code_granddiff( levels(x))
#'     code_sum(       levels(x))
#'     code_helmert(   levels(x))
#' 
#' # Code
#'     require(magrittr)
#'     x %<>% code('treatment')
#'     x %<>% code('reference')
#'     x %<>% code('grandref')
#'     x %<>% code('difference')
#'     x %<>% code('granddiff')
#'     x %<>% code('sum')
#'     x %<>% code('helmert')
#'
#' # Model
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     object %<>% fit_limma(coding = 'treatment') # default
#'     object %<>% fit_limma(coding = 'reference')
#'     object %<>% fit_limma(coding = 'grandref')
#'     object %<>% fit_limma(coding = 'difference')
#'     object %<>% fit_limma(coding = 'granddiff')
#'     object %<>% fit_limma(coding = 'sum')
#'     object %<>% fit_limma(coding = 'helmert')
#' @export
code <- function(object, ...)  UseMethod('code')

#' @rdname code
#' @export
code.factor <- function(object, coding, verbose = TRUE, ...){
    
    assert_scalar_subset(coding, c('treatment', 'reference', 'difference', 
                                   'grandref', 'granddiff', 'sum', 'helmert'))
    codingfun <- switch(coding, treatment  = contr.treatment, 
                                reference  =  code_reference, 
                                difference =  code_difference, 
                                grandref   =  code_grandref, 
                                granddiff  =  code_granddiff, 
                                sum        =  code_sum,
                                helmert    =  code_helmert)
    
    
    if (is.null(coding))  return(object)
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
code.data.table <- function(object, coding, vars = names(object), verbose = TRUE, ...){
    if (verbose)  cmessage('\t\t%s code factors', stringr::str_to_title(coding))
    for (var in vars){
        if (is.character(object[[var]]))  object[[var]] %<>% factor()
        if (is.logical(  object[[var]]))  object[[var]] %<>% factor()
    }
    if (is.null(coding)) return(object)
    for (var in vars){
        if (is.factor(object[[var]])){
            if (verbose)  cmessage('\t\t\t%s', var)
            object[[var]] %<>% code.factor(coding, verbose = verbose)
        }
    }
    object
}

#' @rdname code
#' @export
code_reference <- function(n){
    y <- contr.treatment(n)
    colnames(y) %<>% paste0('-', n[1])
    y
}

#' @rdname code
#' @export
code_grandref <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    y <- codingMatrices::code_control(n)
    colnames(y) <- paste0(n[-1], '-', n[1])
    y
}

#' @rdname code
#' @export
code_difference <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    k <- length(n)
    y <- codingMatrices::contr.diff(n)
    colnames(y) <- paste0(n[-1], '-', n[-k])
    y
}

#' @rdname code
#' @export
code_granddiff <- function(n){
    if (!requireNamespace('codingMatrices', quietly = TRUE)){
        message("install.packages('codingMatrices'). Then re-run.")
        return(n) 
    }
    k <- length(n)
    y <- codingMatrices::code_diff(n)
    colnames(y) <- paste0(n[-1], '-', n[-k])
    y
}

#' @rdname code
#' @export
code_sum <- function(n){
    y <- contr.sum(n)
    k <- length(n)
    colnames(y) <- paste0(n[-k],  '-grand')
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
    colnames(y) <- paste0(n[-1],  '-helmert')
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
#' file <- download_data('atkin18.metabolon.xlsx')
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
#' @param object       SummarizedExperiment
#' @param formula      modeling formula
#' @param drop         TRUE or FALSE
#' @param coding       factor coding system: 'treatment', 'reference', 'difference', 
#'                                           'grandref',  'granddiff', 'sum', 'helmert'
#' @param design       design matrix
#' @param contrasts    NULL or character vector: coefficient contrasts to test
#' @param coefs        NULL or character vector: model coefs to test
#' @param block        block svar (or NULL)
#' @param weightvar    NULL or name of weight matrix in assays(object)
#' @param statvars     character vector: subset of c('effect', 'p', 'fdr', 't')
#' @param sep          string: pvar separator  ("~" in "p~t2~limma")
#' @param suffix       string: pvar suffix ("limma" in "p~t2~limma")
#' @param verbose      whether to msg
#' @param plot         whether to plot
#' @return Updated SummarizedExperiment
#' @examples
#' # Read
#'     require(magrittr)
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     
#' # Standard
#'     object %<>% fit_lm(        ~ subgroup)                 #     statistics default
#'     object %<>% fit_limma(     ~ subgroup)                 # bioinformatics default
#' 
#' # Blocked
#'     object %<>% fit_limma(     ~ subgroup, block = 'SUB')  #        simple random effects
#'     object %<>% fit_lme(       ~ subgroup, block = 'SUB')  #      powerful random effects
#'     object %<>% fit_lmer(      ~ subgroup, block = 'SUB')  # more powerful random effects
#'     
#' # Intuitive : alternative coding
#'     object %<>% fit_lme(       ~ subgroup, block = 'SUB', coding = 'reference')
#'     object %<>% fit_lmer(      ~ subgroup, block = 'SUB', coding = 'reference')
#'     object %<>% fit_limma(     ~ subgroup, block = 'SUB', coding = 'reference')
#'     
#' # Flexible : limma contrasts
#'     object %<>% fit_limma( ~ 0 + subgroup, block = 'SUB', contrasts = c('t1-t0'))
#'         # flexible, but only approximate
#'         # stat.ethz.ch/pipermail/bioconductor/2014-February/057682.html
#'     
#' # Non-parametric: wilcoxon
#'     object %<>% fit_wilcoxon( ~ subgroup)                # unpaired
#'     object %<>% fit_wilcoxon( ~ subgroup, block = 'SUB') # paired
#'     plot_contrast_venn(is_sig(object, contrast = 't2', fit = c('lm', 'limma')))
#'    #plot_contrast_venn(is_sig(object, contrast = 't3', fit = c('limma', 'lme')))
#'     
#' @export
fit_limma <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    coding    = 'treatment',
    design    = create_design(object, formula = formula, drop = drop, coding = coding),
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
        drop         = drop,          coding       = coding,
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
#' require(magrittr)
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


#' @rdname fit_limma
#' @export
.fit_limma <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    coding    = 'treatment',
    design    = create_design(object, formula = formula, drop = drop, coding = coding),
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
#' file <- download_data('atkin18.metabolon.xlsx')
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
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% fit_lm()
#' object %<>% fit_limma(block = 'SUB')
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


