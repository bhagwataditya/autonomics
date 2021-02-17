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
    pvalues <- data.table::as.data.table(as.list(fitres[, pname]))
    tvalues <- data.table::as.data.table(as.list(fitres[, tname]))
    effects <- data.table::as.data.table(as.list(fitres[, effectname]))
    stderrs <- data.table::as.data.table(as.list(fitres[, sename]))
    names(pvalues) %<>% paste0('p.', .)
    names(tvalues) %<>% paste0('t.', .)
    names(effects) %<>% paste0('effect.', .)
    names(stderrs) %<>% paste0('se.', .)
    cbind(pvalues, tvalues, effects, stderrs, F = fval, F.p = f.p )
}

.fit_lm <- function(sd, formula, block=NULL){
    fitres <- lm(  formula    = formula, 
                    data      = sd,
                    na.action = stats::na.omit)
    .cbindstats(fitres)
}

.fit_lme <- function(sd, formula, block){
    fitres <- lme(  fixed     = formula, 
                    random    = block, 
                    data      = sd,
                    na.action = stats::na.omit)
    .cbindstats( fitres, fname='F-value', f.pname='p-value', pname='p-value', 
                tname = 't-value', effectname = 'Value', sename = 'Std.Error')
}

.fit_lmer <- function(sd, formula, block = NULL){
    fitres <- lme4::lmer(
        formula   = formula,
        data      = sd,
        na.action = stats::na.omit, 
        control   = lme4::lmerControl(check.conv.singular = lme4::.makeCC(
                                            action = "ignore",  tol = 1e-4)))
                    # https://stackoverflow.com/a/55367171
    fitres %<>% lmerTest::as_lmerModLmerTest()
    .cbindstats(fitres)
}

.extractstat <- function(lmeres, quantity){
    idx <- stri_startswith_fixed(names(lmeres), paste0(quantity, '.'))
    lmemat <- as.matrix(lmeres[, idx, with=FALSE])
    rownames(lmemat) <- lmeres$feature_id
    colnames(lmemat) %<>% stri_replace_first_fixed(paste0(quantity, '.'), '')
    lmemat
}

#' Statistical models supported in autonomics
#' @export
TESTS <- c('limma','lm','lme','lmer', 'wilcoxon')


fit_lmx <- function(
    object, 
    fit, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                'value ~ %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    block = NULL, contrastdefs = NULL, verbose = TRUE, plot =  TRUE
){
# Assert
    assert_is_a_string(fit);  assert_is_subset(fit, TESTS)
# Prepare dt
    allx <- c(setdiff(all.vars(formula), 'value'), all.vars(block))
    dt <- sumexp_to_long_dt(object, svars = allx)
    fixedx <- setdiff(allx, all.vars(block))
    for (x in fixedx){  
        dt[[x]] %<>% factor()
        stats::contrasts(dt[[x]]) <- MASS::contr.sdif(levels(dt[[x]])) }
# lmefit
    fitmethod <- get(paste0('.fit_', fit))
    lmeres <- dt[, fitmethod(.SD, formula=formula, block=block),by='feature_id']
    fdt <- lmeres[, c('feature_id', 'F', 'F.p'), with=FALSE]
    setnames(fdt, c('F', 'F.p'), c('F.lme', 'F.p.lme'))
    object %<>% merge_fdata(fdt)
# lmeextract
    for (x in setdiff(all.vars(formula), 'value'))  names(lmeres) %<>% 
        stri_replace_first_fixed(x, '')
    quantities <- c('p', 't', 'effect', 'se')
    lmeres <- mapply(.extractstat, quantity=quantities, 
                            MoreArgs = list(lmeres=lmeres), SIMPLIFY = FALSE)
    lmeres$fdr  <- apply(lmeres$p, 2, p.adjust, 'fdr')
    lmeres$bonf <- apply(lmeres$p, 2, p.adjust, 'bonf')
    metadata(object)$lmx <- do.call(abind::abind, c(lmeres, along=3))
    metadata(object)$lmx %<>% extract(rownames(object),,)
    names(metadata(object)) %<>% stri_replace_first_fixed('lmx', fit)
    
    if (verbose) cmessage_df('\t\t\t%s', summarize_fit(object, fit))
    if (is.null(contrastdefs)) contrastdefs <-colnames(metadata(object)[[fit]])
    if (length(contrastdefs) > 1) contrastdefs %<>% setdiff('(Intercept)')
    if (plot)  print(plot_volcano(
        object, fit=fit, contrastdefs = contrastdefs)) 
    object
}


#' @rdname fit_limma
#' @export
fit_lm <- function(
    object,
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                'value ~ %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    block = NULL, contrastdefs = NULL, verbose = TRUE, plot =  TRUE
){
    if (verbose)  cmessage('\t\tlm(%s)', Reduce(paste, deparse(formula)))
    fit_lmx(object, fit = 'lm', subgroupvar = subgroupvar, 
            formula = formula, block = block, contrastdefs = contrastdefs, 
            verbose = verbose, plot = plot)
}


#' @rdname fit_limma
#' @export
fit_lme <- function(
    object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                'value ~ %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    block = NULL, contrastdefs = NULL, verbose = TRUE, plot =  TRUE
){
    assert_is_not_null(block)
    if (is_a_string(block)){ 
        if (is_subset(block, svars(object)))  block %<>% sprintf('~1|%s', .)
        block %<>% as.formula() }
    if (verbose)  cmessage('\t\tlme(%s, random = %s)', 
                        Reduce(paste, deparse(formula)), 
                        Reduce(paste, deparse(block)))
    fit_lmx(object, fit = 'lme', subgroupvar = subgroupvar, 
            formula = formula, block = block, contrastdefs = contrastdefs, 
            verbose = verbose, plot = plot)
}


#' @rdname fit_limma
#' @export
fit_lmer <- function(
    object, 
    fit, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                'value ~ %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    block = NULL, contrastdefs = NULL, verbose = TRUE, plot =  TRUE
){
    if (is_formula(block))  block <- Reduce(paste, deparse(block))
    if (is_a_string(block)){ 
        if (is_subset(block, svars(object)))  block %<>% sprintf('1|%s', .)
        formula <- Reduce(paste, deparse(formula))
        formula %<>% sprintf('%s + (%s)', ., block)
        formula %<>% as.formula()
    }
    if (verbose)  cmessage('\t\tlmer(%s)', Reduce(paste, deparse(formula)))
    fit_lmx(object, fit = 'lmer', subgroupvar = subgroupvar, 
            formula = formula, block = NULL, contrastdefs = contrastdefs, 
            verbose = verbose, plot = plot)
}
