#==============================================================================
#
#                       fit_lmx: lm / lme / lmer
#
#==============================================================================

.cbindstats <- function(
    fitres, fname = 'F value', f.pname = 'Pr(>F)', pname = 'Pr(>|t|)', 
    tname = 't value', effectname = 'Estimate', sename = 'Std. Error'
){
    fval <- stats::anova(fitres)['subgroup', fname]
    f.p <- stats::anova(fitres)['subgroup', f.pname]
    fitres %<>% summary()
    fitres %<>% stats::coefficients()
    fitres
    
    pvalues <- as.data.table(as.list(fitres[, pname,      drop=FALSE]))
    tvalues <- as.data.table(as.list(fitres[, tname,      drop=FALSE]))
    effects <- as.data.table(as.list(fitres[, effectname, drop=FALSE]))
    stderrs <- as.data.table(as.list(fitres[, sename,     drop=FALSE]))
    names(pvalues) %<>% paste0('p.', .)
    names(tvalues) %<>% paste0('t.', .)
    names(effects) %<>% paste0('effect.', .)
    names(stderrs) %<>% paste0('se.', .)
    cbind(pvalues, tvalues, effects, stderrs, F = fval, F.p = f.p )
}

.lm <- function(sd, formula, block, weights){
    fitres <- stats::coefficients(summary(lm(
                    formula   = formula, 
                    data      = sd,
                    weights   = weights,
                    na.action = stats::na.omit)))
             #      F = stats::anova(fitres)[, 'F value'], 
            #       F.p = stats::anova(fitres)[, 'Pr(>F)'])
    colnames(fitres) %<>% stri_replace_first_fixed('Estimate', 'effect')
    colnames(fitres) %<>% stri_replace_first_fixed('Std. Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t value', 't')
    colnames(fitres) %<>% stri_replace_first_fixed('Pr(>|t|)', 'p')
    fitmat <- matrix(fitres, nrow=1)
    colnames(fitmat) <- paste(rep(colnames(fitres), each=nrow(fitres)), 
                        rep(rownames(fitres), times =ncol(fitres)), sep = '.')
    fitmat %<>% cbind(F=0, F.p=1)
    data.table(fitmat)
}

.lme <- function(sd, formula, block, weights){
    fitres <- stats::coefficients(summary(lme(
                    fixed     = formula, 
                    random    = block, 
                    data      = sd,
                    na.action = stats::na.omit)))
    colnames(fitres) %<>% stri_replace_first_fixed('Value', 'effect')
    colnames(fitres) %<>% stri_replace_first_fixed('Std.Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t-value', 't')
    colnames(fitres) %<>% stri_replace_first_fixed('p-value', 'p')
    fitres %<>% extract(, c('effect', 'se', 't', 'p')) # drop DF column
    fitmat <- matrix(fitres, nrow=1)
    colnames(fitmat) <- paste(rep(colnames(fitres), each=nrow(fitres)), 
                        rep(rownames(fitres), times = ncol(fitres)), sep = '.')
    fitmat %<>% cbind(F=0, F.p=1)
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
    colnames(fitres) %<>% stri_replace_first_fixed('Estimate', 'effect')
    colnames(fitres) %<>% stri_replace_first_fixed('Std. Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t value', 't')
    colnames(fitres) %<>% stri_replace_first_fixed('Pr(>|t|)', 'p')
    fitres %<>% extract(, c('effect', 'se', 't', 'p')) # drop DF column
    fitmat <- matrix(fitres, nrow=1)
    colnames(fitmat) <- paste(rep(colnames(fitres), each=nrow(fitres)), 
                        rep(rownames(fitres), times = ncol(fitres)), sep = '.')
    fitmat %<>% cbind(F=0, F.p=1)
    data.table(fitmat)
}

.extractstat <- function(fitres, quantity){
    idx <- stri_startswith_fixed(names(fitres), paste0(quantity, '.'))
    mat <- as.matrix(fitres[, idx, with=FALSE])
    rownames(mat) <- fitres$feature_id
    colnames(mat) %<>% stri_replace_first_fixed(paste0(quantity, '.'), '')
    mat
}

#' Statistical models supported in autonomics
#' @export
TESTS <- c('limma','lm','lme','lmer', 'wilcoxon')


fit_lmx <- function(object, fit, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar, fit), block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    contrastdefs = NULL, verbose = TRUE, plot =  FALSE
){
# Initialize
    assert_is_a_string(fit);  assert_is_subset(fit, TESTS)
    formula %<>% formula2str()
    formula %<>% paste0('value ', .)
    formula %<>% as.formula()
    allx <- c(setdiff(all.vars(formula), 'value'), all.vars(block))
    assnames <- assayNames(object)[1]
    if (!is.null(weightvar)){   
        assert_is_character(weightvar)
        assert_is_subset(weightvar, assayNames(object)) 
        assnames %<>% c(weightvar)
        cmessage('\t\t\tweights = assays(object)$%s', weightvar) }
    dt <- sumexp_to_long_dt(object, svars = allx, assay = assnames)
    fixedx <- setdiff(allx, all.vars(block))
    for (x in fixedx){  
        dt[[x]] %<>% factor()
        n <- length(levels(dt[[x]]))
        if (n>1) stats::contrasts(dt[[x]]) <- MASS::contr.sdif(levels(dt[[x]]))}
# fit
    fitmethod <- get(paste0('.', fit))
    if (is.null(weightvar)){ weightvar <- 'weights'; weights <- NULL }
    fitres <- dt[, fitmethod(.SD, formula=formula, block=block, 
             weights = get(weightvar)), by = 'feature_id' ]
    fdt <- fitres[, c('feature_id', 'F', 'F.p'), with=FALSE]
    names(fdt) %<>% stri_replace_first_fixed('F.p', paste0('F.p.', fit))
    names(fdt) %<>% stri_replace_first_fixed('F',   paste0('F.',   fit))
    object %<>% merge_fdata(fdt)
# extract
    for (x in setdiff(all.vars(formula), 'value'))  names(fitres) %<>% 
        stri_replace_first_fixed(x, '')
    quantities <- c('p', 't', 'effect', 'se')
    fitres <- mapply(.extractstat, quantity=quantities, 
                            MoreArgs = list(fitres=fitres), SIMPLIFY = FALSE)
    fitres$fdr  <- apply(fitres$p, 2, p.adjust, 'fdr')
    fitres$bonf <- apply(fitres$p, 2, p.adjust, 'bonf')
    fitarray <- do.call(abind::abind, c(fitres, along=3))
    colnames(fitarray) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    metadata(object)$fit <- fitarray
    names(metadata(object)) %<>% stri_replace_first_fixed('fit', fit)
    if (verbose) cmessage_df('\t\t\t%s', summarize_fit(object, fit))
    if (is.null(contrastdefs)) contrastdefs <-colnames(metadata(object)[[fit]])
    if (length(contrastdefs) > 1) contrastdefs %<>% setdiff('(Intercept)')
    if (plot)  print(plot_volcano(object, fit=fit, contrastdefs = contrastdefs)) 
    object
}


#' @rdname fit_limma
#' @export
fit_lm <- function(
    object,
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar, fit='lm'), block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    contrastdefs = NULL, verbose = TRUE, plot =  FALSE
){
    if (verbose)  cmessage('\t\tlm(%s)', formula2str(formula))
    fit_lmx(object, fit = 'lm', subgroupvar = subgroupvar, 
            formula = formula, block = block, weightvar = weightvar, 
            contrastdefs = contrastdefs, verbose = verbose, plot = plot)
}


#' @rdname fit_limma
#' @export
fit_lme <- function(
    object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar, fit='lme'), block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    contrastdefs = NULL, verbose = TRUE, plot =  FALSE
){
    assert_is_not_null(block)
    if (is_a_string(block)){ 
        if (is_subset(block, svars(object)))  block %<>% sprintf('~1|%s', .)
        block %<>% as.formula() }
    if (verbose)  cmessage('\t\tlme(%s, random = %s)', 
                        formula2str(formula), formula2str(block))
    fit_lmx(object, fit = 'lme', subgroupvar = subgroupvar, 
            formula = formula, block = block, weightvar = weightvar, 
            contrastdefs = contrastdefs, verbose = verbose, plot = plot)
}


#' @rdname fit_limma
#' @export
fit_lmer <- function(
    object, 
    fit, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = default_formula(object, subgroupvar, fit='lmer'), block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    contrastdefs = NULL, verbose = TRUE, plot =  FALSE
){
    if (is_formula(block))  block %<>% formula2str()
    if (is_a_string(block)){ 
        if (is_subset(block, svars(object)))  block %<>% sprintf('1|%s', .)
        formula <- formula2str()
        formula %<>% sprintf('%s + (%s)', ., block)
        formula %<>% as.formula()
    }
    if (verbose)  cmessage('\t\tlmer(%s)', formula2str(formula))
    fit_lmx(object, fit = 'lmer', subgroupvar = subgroupvar, 
            formula = formula, block = NULL, weightvar = weightvar, 
            contrastdefs = contrastdefs, verbose = verbose, plot = plot)
}
