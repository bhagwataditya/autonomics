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
    fitres <- nlme::lme(
                    fixed     = formula, 
                    random    = block, 
                    data      = sd,
                    na.action = stats::na.omit)
    suppressWarnings(fitres %<>% nlme:::summary.lme())  # only 2 replicates in a group -> df = 0 -> p = NaN -> warning
    fitres %<>% stats::coefficients()
    fitres %<>% extract(, c('Value', 'p-value', 't-value'), drop = FALSE)
    colnames(fitres) %<>% stri_replace_first_fixed('Value', 'effect')
    #colnames(fitres) %<>% stri_replace_first_fixed('Std.Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t-value', 't')
    colnames(fitres) %<>% stri_replace_first_fixed('p-value', 'p')
    #fitres %<>% extract(, c('effect', 'se', 't', 'p'), drop = FALSE) # drop DF
    fitmat <- matrix(fitres, nrow = 1)
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

#' Switch blockvar format
#' @param formula  formula
#' @param block    block: charactervector or formula
#' @param verbose  TRUE or FALSE
#' @examples
#' # lme
#'    block_formula_lme( block += c('SUB', 'SEX'), formula = ~ GROUP + SUBGROUP)
#'    block_formula_lme( block = ~1|SUB + 1|SEX,  formula = ~ GROUP + SUBGROUP)
#'    block_vars(~1|SUB + 1|SEX)
#' # lmer
#'    block_formula_lmer(block = c('SUB', 'SEX'), formula = ~ GROUP + SUBGROUP)
#'    block_formula_lmer(formula = ~GROUP + SUBGROUP + (1|SUB) + (1|SEX))
#'    block_vars(formula = ~GROUP + SUBGROUP + (1|SUB) + (1|SEX))
#' @export
block_formula_lme  <- function(block, formula, verbose = TRUE){
    if (is_formula(block))  return(block) # already in required format
    formula %<>% formula2str()
    block %<>% paste0('1|', .)
    block %<>% paste0(collapse = ' + ')
    block %<>% paste0('~ ', .)
    if (verbose)  cmessage('\t\tlme(%s, random = %s)',  formula, block)
    block %<>% as.formula()
    block
}

#' @rdname block_formula_lme
#' @export
block_formula_lmer <- function(formula, block, verbose = TRUE){
    if (stri_detect_fixed(formula2str(formula), '|'))  return(formula)
    formula %<>% formula2str()
    block %<>% paste0('(1|', ., ')')
    block %<>% paste0(collapse = ' + ')
    formula %<>% paste0(' + ', block)
    if (verbose)  cmessage('\t\tlmer(%s)', formula)
    formula %<>% as.formula()
    formula
}

#' @rdname block_formula_lme
#' @export
block_vars <- function(formula){
    formula %<>% formula2str()
    formula %<>% stri_replace_first_regex('^[~][ ]*', '')
    formula %<>% stri_split_regex('[ ]*[+][ ]*')
    formula %<>% extract2(1)
    formula %<>% extract(stri_detect_fixed(., '|'))
    formula %<>% stri_replace_first_regex('[ ]*[(][ ]*', '') # lmer
    formula %<>% stri_replace_first_regex('[ ]*[)][ ]*', '') # lmer
    formula %<>% split_extract_regex('[ ]*[|][ ]*', 2)
    formula
}




#' Extract full block features
#' 
#' Extract features with two blocks (or more) spanning all within-factor levels
#' @param object    SummarizedExperiment
#' @param formula   formula
#' @param block     factorvector (numeric object) or svar (sumexp object)
#' @param n         number of full factor blocks
#' @param verbose   TRUE or FALSE
#' @examples
#' # Read
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     cols <- c('sample_id', 'subgroup', 'SUB', 'SET', 'AGE', 'SEX', 'T2D')
#'     sdt(object) %<>% extract(, cols, with = FALSE)
#'     sdt(object)
#'     object$SUPERSET <- ''
#'     object$SUPERSET[object$SET == 't0'] <- 'I'
#'     object$SUPERSET[object$SET == 't1'] <- 'I'
#'     object$SUPERSET[object$SET == 't2'] <- 'II'
#'     object$SUPERSET[object$SET == 't3'] <- 'II'
#'     object$SUBSET <- ''
#'     object$SUBSET[  object$SET == 't0'] <- 'a'
#'     object$SUBSET[  object$SET == 't1'] <- 'a'
#'     object$SUBSET[  object$SET == 't2'] <- 'b'
#'     object$SUBSET[  object$SET == 't3'] <- 'b'
#' # Extract
#'     extract_full_block_features(object, formula = ~ SET,                     block = 'SUB')
#'     extract_full_block_features(object, formula = ~ SUPERSET + SUBSET,       block = 'SUB')
#'     extract_full_block_features(object, formula = ~ SUPERSET + SUBSET,       block = 'SUB')
#'     extract_full_block_features(object, formula = ~ SUPERSET + SUBSET,       block = c('SUB', 'SEX'))
#'     extract_full_block_features(object, formula = ~ SUPERSET + SUBSET + AGE, block = c('SUB', 'SEX'))
#'     extract_full_block_features(object, formula = ~ SUPERSET + SUBSET + T2D, block = c('SUB', 'SEX'))
#' @export
extract_full_block_features <- function(
    object, formula, block, n = 2, verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(all.vars(formula), svars(object))
    if (is_formula(block))    block %<>% block_vars()
    assert_is_subset(block, svars(object))
    assert_is_a_number(n)
    assert_is_a_bool(verbose)
# Identify within-block factors
    sdt0 <- sdt(object)
    sdt0 %<>% extract(, c(all.vars(formula), block), with = FALSE)                    # design/block
    sdt0 %<>% extract(, false_names(vapply(., is.numeric, logical(1))), with = FALSE) # factors
    for (curblock in block){
        # Identify in-block factors
            inblockdt <- sdt0[, lapply(.SD, function(x) paste0(unique(x), collapse = ';')), by = curblock]
            inblockdt %<>% extract(, -1, with = FALSE)
            inblockdt %<>% unique()
            idx <- vapply(inblockdt, function(x)  any(stri_detect_fixed(x, ';')), logical(1) )
            inblockdt %<>% extract(, idx, with = FALSE)
        # Identify features: at least two blocks span across every in-block factor
            for (inblockfactor in names(inblockdt)){
                dt <- sumexp_to_longdt(object, svars = c(curblock, inblockfactors), fvars = 'feature_id')
                dt %<>% extract(!is.na(value))
                dt %<>% extract(, c('feature_id', curblock, inblockfactors), with = FALSE)
                
                inblockfactorgroups <- inblockdt[[inblockfactor]] %>% extract(stri_detect_fixed(., ';'))
                for (inblocklevels in inblockfactorgroups){
                    inblocklevels %<>% stri_split_fixed(';') %>% extract2(1)
                    nlevel <- length(inblocklevels)
                    idx <- dt[[inblockfactor]] %in% inblocklevels
                    subdt <- dt[idx]
                    subdt <- subdt[, .(alllevels = length(unique(get(inblockfac))) == nlevel), by = c(curblock, 'feature_id')]
                    subdt <- subdt[, .(completeblocks = sum(alllevels)) , by = 'feature_id']
                    subdt <- subdt[completeblocks >= n]
                    if (verbose & nrow(subdt) < nrow(object)){
                        cmessage('\t\tRetain %d/%d features: %d or more %ss span %ss: %s', 
                                 nrow(subdt), nrow(object), n, curblock, inblockfactor, paste0(inblocklevels, collapse = ', '))
                        idx <- fdt(object)$feature_id %in% subdt$feature_id
                        object %<>% extract(idx, )
                    }
                }
            }
    }
# Identify completeblock features
# Extract and Return
    object
}



#' Fit lm, lme, or lmer
#' @param object       SummarizedExpriment
#' @param fit         'lm', 'lme', or 'lmer'
#' @param subgroupvar  svar
#' @param formula      formula
#' @param drop         TRUE or FALSE
#' @param block        NULL or svar
#' @param weightvar    NULL or svar
#' @param coefs        NULL or stringvector
#' @param verbose      TRUE or FALSE
#' @param plot         TRUE or FALSE
#' @return SummarizedExperiment
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% extract_full_block_features(formula = ~subgroup, block = 'SUB')
#' fit_lm(   object, formula = ~subgroup)
#' fit_limma(object, formula = ~subgroup)
#' fit_limma(object, formula = ~subgroup, block = 'SUB')
#' fit_lme(  object, formula = ~subgroup, block = 'SUB')
#' fit_lmer( object, formula = ~subgroup, block = 'SUB')
#' fit_lme(  object, formula = ~subgroup, block = ~1|SUB)
#' fit_lmer( object, formula = ~subgroup + (1|SUB))
#' @export
fit_lmx <- function(
    object, 
    fit, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula     = default_formula(object, subgroupvar, contrasts = NULL), 
    drop        = varlevels_dont_clash(object, all.vars(formula)),
    block       = NULL, 
    weightvar   = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    coefs       = NULL, 
    verbose     = TRUE, 
    plot        = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(fit, c('lm', 'lme', 'lmer'))
    assertive.types::assert_is_formula(formula)
    assert_is_subset(all.vars(formula), svars(object))
    assert_is_a_bool(drop)
    if (!is.null(weightvar)){
        assert_is_character(weightvar)
        assert_is_subset(weightvar, assayNames(object)) 
        assnames %<>% c(weightvar)
        message('\t\t\tweights = assays(object)$', weightvar) 
    }
# Filter / Customize
    obj <- object
    if (fit %in% c('lme', 'lmer')){  
        obj %<>% extract_full_block_features(formula = formula, block = block)  }
    if (       fit == 'lme'){  block   %<>% block_formula_lme(formula = formula, verbose = verbose)
    } else if (fit == 'lmer'){ formula %<>% block_formula_lmer( block   = block,   verbose = verbose)  }
# Fit
    dt <- sumexp_to_longdt(obj, svars = c(all.vars(formula), all.vars(block)), assay = assayNames(object)[1])
    fitmethod <- get(paste0('.', fit))
    if (is.null(weightvar)){ weightvar <- 'weights'; weights <- NULL }
    formula %<>% addlhs()
    fitres <- dt[, fitmethod(.SD, formula = formula, block = block, weights = get(weightvar)), by = 'feature_id' ]
    names(fitres) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    if (drop)  for (var in all.vars(formula))   names(fitres) %<>% stri_replace_first_fixed(var, '')
    fitres %<>% add_fdr()
    object %<>% reset_fitres(fit)
    object %<>% merge_fitres(fitres, fit = fit)
# extract
    quantities <- c('effect', 'fdr', 'p', 't')
    fitres <- mapply(.extractstat, quantity=quantities, 
                            MoreArgs = list(fitres=fitres), SIMPLIFY = FALSE)
    formula %<>% droplhs() %<>% formula2str()
    if (!is.null(weights))  formula %<>% paste0(', weights = assays(object)$', weightvar)
    if (verbose)  message_df('\t\t\t%s', summarize_fit(fdt(object), fit = fit))
    if (is.null(coefs)) coefs <- autonomics::coefs(object, fit = fit)
    if (length(coefs) > 1) coefs %<>% setdiff('Intercept')
    if (plot)  print(plot_volcano(object, fit = fit, coefs = coefs))
    object 
}


#' @rdname fit_lmx
#' @export
fit_lm <- function(
    object,
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    contrasts   = NULL,  # for interface equivalence only
    formula     = default_formula(object, subgroupvar, contrasts = NULL), 
    drop         = varlevels_dont_clash(object, all.vars(formula)),
    block       = NULL, 
    weightvar   = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    coefs       = NULL, 
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
        weightvar   = weightvar,    coefs        = coefs, 
        verbose     = verbose,      plot         = plot)
}


#' @rdname fit_lmx
#' @export
fit_lme <- function(
    object, 
    subgroupvar  = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    contrasts    = NULL,  # for interface equivalence only
    formula      = default_formula(object, subgroupvar, contrasts = NULL), 
    drop         = varlevels_dont_clash(object, all.vars(formula)),
    block        = NULL, 
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    coefs        = NULL, 
    verbose      = TRUE, 
    plot         = FALSE
){
# Assert
    . <- NULL
    if (!requireNamespace('nlme', quietly = TRUE)){
        message("BiocManager::install('nlme'). Then re-run.")
        return(object)   }
# Fit
    fit_lmx(
        object,                     fit          = 'lme', 
        subgroupvar = subgroupvar,  formula      = formula, 
        drop        = drop,         block        = block, 
        weightvar   = weightvar,    coefs        = coefs, 
        verbose     = verbose,      plot         = plot)
}


#' @rdname fit_lmx
#' @export
fit_lmer <- function(
    object, 
    subgroupvar  = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    contrasts    = NULL,  # for function equivalence only
    formula      = default_formula(object, subgroupvar, contrasts = NULL), 
    drop          = varlevels_dont_clash(object, all.vars(formula)),
    block        = NULL, 
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    coefs        = NULL, 
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
# Fit
    metadata(object)$lmer <- NULL
    fit_lmx(
        object,                     fit          = 'lmer', 
        subgroupvar = subgroupvar,  formula      = formula, 
        drop        = drop,         block        = block, 
        weightvar   = weightvar,    coefs        = coefs, 
        verbose     = verbose,      plot         = plot)
}
