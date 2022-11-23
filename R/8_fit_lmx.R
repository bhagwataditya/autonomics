#==============================================================================
#
#                       fit_lmx: lm / lme / lmer
#
#==============================================================================

.cbindstats <- function(
    fitres, fname = 'F value', f.pname = 'Pr(>F)', pname = 'Pr(>|t|)', 
    tname = 't value', effectname = 'Estimate', sename = 'Std. Error'
){
    . <- NULL
    fval <- stats::anova(fitres)['subgroup', fname]
    f.p <- stats::anova(fitres)['subgroup', f.pname]
    fitres %<>% summary()
    fitres %<>% stats::coefficients()
    fitres
    
    pvalues <- as.data.table(as.list(fitres[, pname,      drop = FALSE]))
    tvalues <- as.data.table(as.list(fitres[, tname,      drop = FALSE]))
    effects <- as.data.table(as.list(fitres[, effectname, drop = FALSE]))
    stderrs <- as.data.table(as.list(fitres[, sename,     drop = FALSE]))
    names(pvalues) %<>% paste0('p.', .)
    names(tvalues) %<>% paste0('t.', .)
    names(effects) %<>% paste0('effect.', .)
    names(stderrs) %<>% paste0('se.', .)
    cbind(pvalues, tvalues, effects, stderrs, F = fval, F.p = f.p )
}

.lm <- function(sd, formula, block, weights){
    formula <- as.formula(formula)  # stackoverflow.com/questions/51142338
    environment(formula) <- environment()
    fitres <- stats::coefficients(summary(lm(
                    formula   = formula, 
                    data      = sd,
                    weights   = weights,
                    na.action = stats::na.omit)))
    fitres %<>% extract(, c('Estimate', 'Pr(>|t|)', 't value'), drop=FALSE)
            #      F = stats::anova(fitres)[, 'F value'], 
            #       F.p = stats::anova(fitres)[, 'Pr(>F)'])
    colnames(fitres) %<>% stri_replace_first_fixed('Estimate', 'effect')
    # colnames(fitres) %<>% stri_replace_first_fixed('Std. Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t value',  't')
    colnames(fitres) %<>% stri_replace_first_fixed('Pr(>|t|)', 'p')
    fitmat <- matrix(fitres, nrow=1)
    colnames(fitmat) <- paste(rep(colnames(fitres), each=nrow(fitres)), 
                        rep(rownames(fitres), times =ncol(fitres)), sep=FITSEP)
    #fitmat %<>% cbind(F=0, F.p=1)
    data.table(fitmat)
}

.lme <- function(sd, formula, block, weights){
    fitres <- stats::coefficients(summary(nlme::lme(
                    fixed     = formula, 
                    random    = block, 
                    data      = sd,
                    na.action = stats::na.omit)))
    fitres %<>% extract(, c('Value', 'p-value', 't-value'), drop = FALSE)
    colnames(fitres) %<>% stri_replace_first_fixed('Value', 'effect')
    #colnames(fitres) %<>% stri_replace_first_fixed('Std.Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t-value', 't')
    colnames(fitres) %<>% stri_replace_first_fixed('p-value', 'p')
    #fitres %<>% extract(, c('effect', 'se', 't', 'p'), drop = FALSE) # drop DF
    fitmat <- matrix(fitres, nrow=1)
    colnames(fitmat) <- paste(rep(colnames(fitres), each=nrow(fitres)), 
                        rep(rownames(fitres), times = ncol(fitres)), sep=FITSEP)
    #fitmat %<>% cbind(F=0, F.p=1)
    data.table(fitmat)
}

.lmer <- function(sd, formula, block, weights){
    fitres <- lme4::lmer(
        formula   = formula,
        data      = sd,
        weights   = weights,
        na.action = stats::na.omit, 
        control   = lme4::lmerControl(
            check.conv.grad     = lme4::.makeCC(action = 'ignore', tol=2e-3),
            check.conv.singular = lme4::.makeCC(action = "ignore", tol=1e-4), #
            check.conv.hess     = lme4::.makeCC(action = 'ignore', tol=1e-6), 
        )) # https://stackoverflow.com/a/55367171
            
    fitres %<>% lmerTest::as_lmerModLmerTest()
    fitres %<>% summary() %>% stats::coefficients()
    fitres %<>% extract(, c('Estimate', 'Pr(>|t|)', 't value'), drop=FALSE)
    colnames(fitres) %<>% stri_replace_first_fixed('Estimate', 'effect')
    #colnames(fitres) %<>% stri_replace_first_fixed('Std. Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t value',  't')
    colnames(fitres) %<>% stri_replace_first_fixed('Pr(>|t|)', 'p')
    #fitres %<>% extract(, c('effect', 'se', 't', 'p'), drop = FALSE) # drop DF
    fitmat <- matrix(fitres, nrow=1)
    colnames(fitmat) <- paste(rep(colnames(fitres), each=nrow(fitres)), 
                        rep(rownames(fitres), times = ncol(fitres)), sep=FITSEP)
    #fitmat %<>% cbind(F=0, F.p=1)
    data.table(fitmat)
}

.extractstat <- function(fitres, quantity){
    idx <- stri_startswith_fixed(names(fitres), paste0(quantity, FITSEP))
    mat <- as.matrix(fitres[, idx, with=FALSE])
    rownames(mat) <- fitres$feature_id
    colnames(mat) %<>% stri_replace_first_fixed(paste0(quantity, FITSEP), '')
    mat
}

#' Statistical models supported in autonomics
#' @examples
#' TESTS
#' @export
TESTS <- c('limma','lm','lme','lmer', 'wilcoxon')

addlhs  <- function(formula)  as.formula(paste0('value ', formula2str(formula)))
droplhs <- function(formula)  as.formula(stri_replace_first_regex(
    formula2str(formula), '^value[ ]*', ''))

fit_lmx <- function(
    object, 
    fit, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula     = default_formula(object, subgroupvar, contrasts = NULL), 
    drop        = varlevels_dont_clash(object, all.vars(formula)),
    block       = NULL, 
    weightvar   = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    coefficients       = NULL, 
    verbose     = TRUE, 
    plot        = FALSE
){
# Initialize
    assert_is_a_string(fit);  assert_is_subset(fit, TESTS)
    formula %<>% addlhs()
    allx <- c(setdiff(all.vars(formula), 'value'), all.vars(block))
    assnames <- assayNames(object)[1]
    if (!is.null(weightvar)){   assert_is_character(weightvar)
                                assert_is_subset(weightvar, assayNames(object)) 
                                assnames %<>% c(weightvar)
        message('\t\t\tweights = assays(object)$', weightvar) }
    dt <- sumexp_to_longdt(object, svars = allx, assay = assnames)
    fixedx <- setdiff(allx, all.vars(block))
    #for (x in fixedx){          dt[[x]] %<>% factor()
    #                            n <- length(levels(dt[[x]]))
    #    if (n>1) stats::contrasts(dt[[x]]) <- MASS::contr.sdif(levels(dt[[x]]))}
# fit
    fitmethod <- get(paste0('.', fit))
    if (is.null(weightvar)){ weightvar <- 'weights'; weights <- NULL }
    fitres <- dt[, fitmethod(.SD, formula = formula, block = block, 
            weights = get(weightvar)), by = 'feature_id' ]
    names(fitres) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    if (drop)  for (var in all.vars(formula)){               # drop varnames
        names(fitres) %<>% stri_replace_first_fixed(var, '')  }
    fitres %<>% add_fdr()
    object %<>% reset_fitres(fit)
    object %<>% merge_fitres(fitres, fit = fit)
# extract
    quantities <- c('effect', 'fdr', 'p', 't')
    fitres <- mapply(.extractstat, quantity=quantities, 
                            MoreArgs = list(fitres=fitres), SIMPLIFY = FALSE)
    formula %<>% droplhs() %<>% formula2str()
    if (!is.null(weights))  formula %<>% paste0(', weights = assays(object)$', weightvar)
    if (verbose)  message_df('\t\t\t%s', summarize_fit(object, fit))
    if (is.null(coefficients)) coefficients <- coefficients(object, fit = fit)
    if (length(coefficients) > 1) coefficients %<>% setdiff('Intercept')
    if (plot)  print(plot_volcano(object, fit = fit, coefficients = coefficients))
    object 
}


#' @rdname fit_limma
#' @export
fit_lm <- function(
    object,
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    contrasts   = NULL,  # for interface equivalence only
    formula     = default_formula(object, subgroupvar, contrasts = NULL), 
    drop         = varlevels_dont_clash(object, all.vars(formula)),
    block       = NULL, 
    weightvar   = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    coefficients       = NULL, 
    verbose     = TRUE, 
    plot        = FALSE
){
# Assert
    for (by in all.vars(formula))  assert_is_identical_to_false(
                                has_consistent_nondetects(object, by))
        # cannot be generified into fit_lmx because for lme/lmer block is 
        # integrated in formula and also gets checked for (unnecessarily!)
# Fit    
    if (verbose)  message('\t\tlm(', formula2str(formula), ')')
    fit_lmx(
        object,                     fit          = 'lm', 
        subgroupvar = subgroupvar,  formula      = formula, 
        drop        = drop,         block        = block, 
        weightvar   = weightvar,    coefficients = coefficients, 
        verbose     = verbose,      plot         = plot)
}


#' @rdname fit_limma
#' @export
fit_lme <- function(
    object, 
    subgroupvar  = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    contrasts    = NULL,  # for interface equivalence only
    formula      = default_formula(object, subgroupvar, contrasts = NULL), 
    drop         = varlevels_dont_clash(object, all.vars(formula)),
    block        = NULL, 
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    coefficients = NULL, 
    verbose      = TRUE, 
    plot         = FALSE
){
# Assert
    . <- NULL
    if (!requireNamespace('nlme', quietly = TRUE)){
        message("BiocManager::install('nlme'). Then re-run.")
        return(object) 
    }
    for (by in all.vars(formula))  assert_is_identical_to_false(
                                has_consistent_nondetects(object, by))
        # cannot be generified into fit_lmx because for lme/lmer block is 
        # integrated in formula and also gets checked for (unnecessarily!)
    assert_is_not_null(block)
# Prepare
    if (is_a_string(block)){ 
        if (is_subset(block, svars(object)))  block %<>% sprintf('~1|%s', .)
        block %<>% as.formula() }
# Fit
    if (verbose)  message('\t\tlme(', formula2str(formula), ', ', 
                        'random = ',  formula2str(block),')')
    fit_lmx(
        object,                     fit          = 'lme', 
        subgroupvar = subgroupvar,  formula      = formula, 
        drop        = drop,         block        = block, 
        weightvar   = weightvar,    coefficients = coefficients, 
        verbose     = verbose,      plot         = plot)
}


#' @rdname fit_limma
#' @export
fit_lmer <- function(
    object, 
    subgroupvar  = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    contrasts    = NULL,  # for function equivalence only
    formula      = default_formula(object, subgroupvar, contrasts = NULL), 
    drop          = varlevels_dont_clash(object, all.vars(formula)),
    block        = NULL, 
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    coefficients = NULL, 
    verbose      = TRUE, 
    plot         = FALSE
){
# Assert
    . <- NULL
    if (!requireNamespace('lme4', quietly = TRUE)){
        message("`BiocManager::install('lme4')`. Then re-run.")
        return(object) }
    if (!requireNamespace('lmerTest', quietly = TRUE)){
        message("`BiocManager::install('lmerTest')`. Then re-run.")
        return(object) }
    for (by in all.vars(formula))  assert_is_identical_to_false(
                                has_consistent_nondetects(object, by))
        # cannot be generified into fit_lmx because for lme/lmer block is 
        # integrated in formula and also gets checked for (unnecessarily!)
# Prepare
    if (is_formula(block))  block %<>% formula2str()
    if (is_a_string(block)){ 
        if (is_subset(block, svars(object)))  block %<>% sprintf('1|%s', .)
        formula %<>% formula2str()
        formula %<>% sprintf('%s + (%s)', ., block)
        formula %<>% as.formula()
    }
# Fit
    if (verbose)  message('\t\tlmer(', formula2str(formula), ')')
    metadata(object)$lmer <- NULL
    fit_lmx(
        object,                     fit          = 'lmer', 
        subgroupvar = subgroupvar,  formula      = formula, 
        drop        = drop,         block        = NULL, 
        weightvar   = weightvar,    coefficients = coefficients, 
        verbose     = verbose,      plot         = plot)
}
