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

.lm <- function(sd, formula, block, weights, optim = NULL){
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

.lme <- function(sd, formula, block, weights, opt = 'optim'){
    ctrl <- nlme::lmeControl(opt = opt)  # https://stats.stackexchange.com/a/40664
    fitres <- nlme::lme(
                    fixed     = formula, 
                    random    = block, 
                    data      = sd,
                    na.action = stats::na.omit, 
                    control   = ctrl)
    suppressWarnings(fitres %<>% summary())  # only 2 replicates in a group -> df = 0 -> p = NaN -> warning
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

.lmer <- function(sd, formula, block, weights, optim = NULL){
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
#'    block_formula_lme( block = c('subject', 'batch'),           formula = ~ subgroup)
#'    block_formula_lme( block = list(subject = ~1, batch = ~1),  formula = ~ subgroup)
#'    block_formula_lme( block = ~1|SUB,                          formula = ~ subgroup)
#' # lmer
#'    block_formula_lmer( block = c('subject', 'batch'),          formula = ~ subgroup)
#'    block_formula_lmer(formula = ~ subgroup + (1|subject) + (1|batch))
#' @export
block_formula_lme  <- function(block, formula, verbose = TRUE){
    if (is.list(block))     return(block)
    if (is_formula(block))  return(block)
    formula %<>% formula2str()
    y <- rep('~1', length(block))
    names(y) <- block
    y %<>% lapply(as.formula)
    if (verbose)  cmessage('\t\tlme(%s, random = %s)',  formula, capture.output(dput(y)))
    y
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




#' Extract connected features
#' 
#' Extract features for which connected factor levels are measured in at least two blocks
#' @param object           SummarizedExperiment
#' @param formula          formula
#' @param blockvars        factorvector (numeric object) or svar (sumexp object)
#' @param nconnectedblocks  minimum number of connected blocks
#' @param verbose          TRUE or FALSE
#' @examples
#' # Read
#'     file <- download_data('atkin.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     cols <- c('sample_id', 'subgroup', 'Subject', 'Time', 'Age', 'Sex', 'Diabetes')
#'     sdt(object) %<>% extract(, cols, with = FALSE)
#'     sdt(object)
#'     object$SUPERSET <- ''
#'     object$SUPERSET[object$Time == 't0'] <- 'I'
#'     object$SUPERSET[object$Time == 't1'] <- 'I'
#'     object$SUPERSET[object$Time == 't2'] <- 'II'
#'     object$SUPERSET[object$Time == 't3'] <- 'II'
#'     object$SUBSET <- ''
#'     object$SUBSET[  object$Time == 't0'] <- 'a'
#'     object$SUBSET[  object$Time == 't1'] <- 'a'
#'     object$SUBSET[  object$Time == 't2'] <- 'b'
#'     object$SUBSET[  object$Time == 't3'] <- 'b'
#' # Extract
#'     extract_connected_features(
#'         object, formula = ~ Time,                         blockvars =   'Subject')
#'     extract_connected_features(
#'         object, formula = ~ SUPERSET + SUBSET,            blockvars =   'Subject')
#'     extract_connected_features(
#'         object, formula = ~ SUPERSET + SUBSET,            blockvars = c('Subject', 'Sex'))
#'     extract_connected_features(
#'         object, formula = ~ SUPERSET + SUBSET + Age,      blockvars = c('Subject', 'Sex'))
#'     extract_connected_features(
#'         object, formula = ~ SUPERSET + SUBSET + Diabetes, blockvars = c('Subject', 'Sex'))
#' @export
extract_connected_features <- function(
    object, formula, blockvars, nconnectedblocks = 2, verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(all.vars(formula), svars(object))
    if (is_formula(blockvars))    blockvars %<>% block_vars()
    assert_is_subset(blockvars, svars(object))
    assert_is_a_number(nconnectedblocks)
    assert_is_a_bool(verbose)
    connected <- connectedblocks <- value <- NULL
# Identify within-blockvar factors
    fixedvars <- all.vars(formula)
    sdt0 <- sdt(object)[, c(fixedvars, blockvars), with = FALSE]                      # design/block vars
    sdt0 %<>% extract(, false_names(vapply(., is.numeric, logical(1))), with = FALSE) # factors
    fixedvars %<>% intersect(names(sdt0))                                             # keep only factor fixedvars
    for (blockvar in blockvars){
        # Within blockvar: identify factors with (sets of) connected levels
            connectedfactorsdt <- sdt0[, lapply(.SD, function(x) paste0(unique(x), collapse = ';')), by = blockvar, .SDcols = fixedvars]
            connectedfactorsdt %<>% extract(, -1, with = FALSE)
            connectedfactorsdt %<>% unique()
            idx <- vapply(connectedfactorsdt, function(x)  any(stri_detect_fixed(x, ';')), logical(1) )
            connectedfactorsdt %<>% extract(, idx, with = FALSE)
            connectedfactorsdt %<>% unique()
        # Identify features: at least two blocks span across every within-blockvar factor
            for (connectedfactor in names(connectedfactorsdt)){
                dt <- sumexp_to_longdt(object, svars = c(blockvar, names(connectedfactorsdt)), fvars = 'feature_id')
                dt %<>% extract(!is.na(value))
                dt %<>% extract(, c('feature_id', blockvar, names(connectedfactorsdt)), with = FALSE)
                connectedsets <- connectedfactorsdt[[connectedfactor]] 
                connectedsets %<>% extract(stri_detect_fixed(., ';'))
                for (connectedlevels in connectedsets){
                    connectedlevels %<>% stri_split_fixed(';') %>% extract2(1)
                    nlevel <- length(connectedlevels)
                    idx <- dt[[connectedfactor]] %in% connectedlevels
                    connectedlevelsdt <- dt[idx]
                    connectedlevelsdt %<>% extract(, .(connected = length(unique(get(connectedfactor))) == nlevel), by = c(blockvar, 'feature_id'))
                    connectedlevelsdt %<>% extract(, .(connectedblocks = sum(connected)) , by = 'feature_id')
                    connectedlevelsdt %<>% extract(connectedblocks >= nconnectedblocks)
                    n0 <- nrow(object); n1 <- nrow(connectedlevelsdt)
                    if (verbose & n1<n0){
                        cmessage('\t\t\tRetain %d/%d features: %d or more %s span %ss: %s', 
                                 n1, n0, nconnectedblocks, blockvar, connectedfactor, paste0(connectedlevels, collapse = ', '))
                        idx <- fdt(object)$feature_id %in% connectedlevelsdt$feature_id
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
#' @param formula      formula
#' @param drop         TRUE or FALSE
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
#' @param codingfun    coding function
#' @param coefs        NULL or stringvector
#' @param contrasts    unused. only to allow generic get(fitfun)(contrasts)
#' @param block        NULL or svar
#' @param opt          optimizer used in fit_lme: 'optim' (more robust) or 'nlminb'
#' @param weightvar    NULL or svar
#' @param statvars     character vector: subset of c('effect', 'p', 'fdr', 't')
#' @param verbose      TRUE or FALSE
#' @param plot         TRUE or FALSE
#' @return SummarizedExperiment
#' @examples 
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' fit_lm(    object, formula = ~subgroup)
#' fit_limma( object, formula = ~subgroup)
#' fit_limma( object, formula = ~subgroup, block = 'Subject')
#' fit_lme(   object, formula = ~subgroup, block = 'Subject')
#' fit_lmer(  object, formula = ~subgroup, block = 'Subject')
#' # fit_lme( object, formula = ~subgroup, block = ~1|Subject) # needs fine-tuning
#' # fit_lmer(object, formula = ~subgroup + (1|Subject))       # needs fine-tuning
#' @export
fit_lmx <- function(
    object, 
    fit, 
    formula   = default_formula(object), 
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment,
    coefs     = colnames(create_design(
                object, formula = formula, drop = drop, codingfun = codingfun, verbose = FALSE)), 
    block     = NULL, 
    opt       = 'optim',
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    statvars  = c('effect', 'p', 'fdr'),
    verbose   = TRUE, 
    plot      = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(fit, c('lm', 'lme', 'lmer'))
    assert_valid_formula(formula, object)
    assert_is_a_bool(drop)
    if (!is.null(weightvar)){
        assert_is_character(weightvar)
        assert_is_subset(weightvar, assayNames(object)) 
        message('\t\t\tweights = assays(object)$', weightvar) 
    }
    N <- value <- V1 <- NULL
# Filter / Customize
    obj <- object
    if (verbose)  cmessage('\t\tFilter features')
    obj %<>% keep_replicated_features(formula, verbose = verbose)
    obj %<>% keep_connected_blocks(   block,   verbose = verbose) # keep samples from fully connected blocks
    obj %<>% keep_connected_features( block,   verbose = verbose) # keep features with 2+ connected blocks
# Fit
    if (       fit == 'lme'){  block   %<>% block_formula_lme(formula = formula, verbose = verbose); blockvars <- names(block)
    } else if (fit == 'lmer'){ formula %<>% block_formula_lmer( block = block,   verbose = verbose); blockvars <- character(0)
    } else {                                                                                         blockvars <- block }
    fitmethod <- get(paste0('.', fit))
    if (is.null(weightvar)){ weightvar <- 'weights'; weights <- NULL }
    assays <- assayNames(object) %>% intersect(c(.[1], 'weights'))
    dt <- sumexp_to_longdt(obj, svars = c(all.vars(formula), blockvars), assay = assays)
    lhsformula <- addlhs(formula)
    fitres <- dt[, fitmethod(.SD, formula = lhsformula, block = block, weights = get(weightvar), opt = opt), by = 'feature_id' ]
    names(fitres) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    if (drop)  for (var in all.vars(formula))   names(fitres) %<>% stri_replace_first_fixed(var, '')
    pattern1 <- sprintf('^(feature_id|%s)',  paste0(statvars, collapse = '|'))   # select statvars
    pattern2 <- sprintf( '(feature_id|%s)$', paste0(coefs,    collapse = '|'))   # select coefs
    fitres <- fitres[, .SD, .SDcols = patterns(pattern1) ]
    fitres <- fitres[, .SD, .SDcols = patterns(pattern2) ]
    fitres %<>% add_fdr()
# Merge back
    object %<>% reset_fit(fit)
    object %<>% merge_fit(fitres, fit = fit)
    formula %<>% droplhs() %<>% formula2str()
    
    if (!is.null(weights))  formula %<>% paste0(', weights = assays(object)$', weightvar)
    if (verbose)  message_df('\t\t\t%s', summarize_fit(fdt(object), fit = fit, coefs = coefs))
    if (length(coefs) > 1)  coefs %<>% setdiff('Intercept')
    if (plot)  print(plot_volcano(object, fit = fit, coefs = coefs))
    object 
}


#' @rdname fit_lmx
#' @export
fit_lm <- function(
    object,
    formula   = default_formula(object), 
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment,
    block     = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    statvars  = c('effect', 'p', 'fdr'),
    coefs     = NULL, 
    contrasts = NULL,
    verbose   = TRUE, 
    plot      = FALSE
){
    
    sdt(object) %<>% code(codingfun = codingfun, vars = all.vars(formula), verbose = verbose)
    if (verbose)  message('\t\tlm(', formula2str(formula), ')')
    fit_lmx(
        object,                      fit          = 'lm', 
        formula      = formula,      drop         = drop,
        codingfun    = codingfun,    block        = block,
        weightvar    = weightvar,    statvars     = statvars,
        coefs        = coefs,        verbose      = verbose,
        plot         = plot)
}


#' @rdname fit_lmx
#' @export
fit_lme <- function(
    object, 
    formula   = default_formula(object), 
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment,
    block     = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    opt       = 'optim',
    statvars  = c('effect', 'p', 'fdr'),
    coefs     = NULL, 
    contrasts = NULL,
    verbose   = TRUE, 
    plot      = FALSE
){
# Assert
    . <- NULL
    if (!requireNamespace('nlme', quietly = TRUE)){
        message("BiocManager::install('nlme'). Then re-run.")
        return(object)   }
# Fit
    sdt(object) %<>% code(codingfun = codingfun, vars = all.vars(formula), verbose = verbose)
    fit_lmx(
        object,                      fit          = 'lme', 
        formula      = formula,      drop         = drop,
        codingfun    = codingfun,    block        = block, 
        weightvar    = weightvar,    opt          = opt,
        coefs        = coefs, 
        verbose      = verbose,      plot         = plot)
}


#' @rdname fit_lmx
#' @export
fit_lmer <- function(
    object, 
    formula   = default_formula(object), 
    drop      = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment,
    block     = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    statvars  = c('effect', 'p', 'fdr'),
    coefs     = NULL, 
    contrasts = NULL,
    verbose   = TRUE, 
    plot      = FALSE
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
    sdt(object) %<>% code(codingfun = codingfun, vars = all.vars(formula), verbose = verbose)
    fit_lmx(
        object,                    fit        = 'lmer', 
        formula    = formula,      drop       = drop,
        codingfun  = codingfun,    block      = block, 
        weightvar  = weightvar,    coefs      = coefs, 
        verbose    = verbose,      plot       = plot)
}
